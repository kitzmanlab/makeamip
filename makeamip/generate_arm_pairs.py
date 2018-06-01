#! /usr/bin/env python

import sys

import os
import os.path
import argparse
from collections import defaultdict, OrderedDict

import numpy as np
import numpy.random as rand
import re

import pandas as pd

import pysam
# import pybedtools as pbt

from functools import partial

from makeamip.capcommon import *
#     

from makeamip.capcommon import *

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--genome', dest='genome', default='/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa',
                      help='path to indexed genome file')
    
    opts.add_argument('--genomeDict', dest='genomeDict', default='/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa.dict', 
                      help='path to genome contig file dictionary, eg from Picard CreateSequenceDictionary')


    opts.add_argument('--inArmStore', dest='inArmStore')

    opts.add_argument('--outPairStore', dest='outPairStore')

    opts.add_argument('--armBaseName', default='', dest='armBaseName')
    opts.add_argument('--pairsBaseName', default='', dest='pairsBaseName')
    
    opts.add_argument('--hitTblSufifx', default='', dest='hitTblSufifx')
    opts.add_argument('--snpTblSuffix', default='', dest='snpTblSuffix')

    opts.add_argument('--generateDsMIPs', action='store_true', default=False, dest='generateDsMIPs')

    opts.add_argument('--gapFillRange', type=intRange, dest='gapFillRange')

    opts.add_argument('--tmDiffRange', type=floatRange, dest='tmDiffRange')

    opts.add_argument('--enforceGapFillOnTopOnly', action='store_true', default=False, dest='enforceGapFillOnTopOnly')

    opts.add_argument('--scoringClass',default='ScoringModelLearned',dest='scoringClass')
    opts.add_argument('--generateMipSeqFxn',default='makeseq_dsmip_eari',dest='generateMipSeqFxn')
    # machine learning scoring parameters
    opts.add_argument('--learningModel', default = 'GradientBoosting', dest = 'clfName', 
        help = 'name a classifier, GradientBoosting or AdaBoost, for machine learning algorithm')
    opts.add_argument('--trainingMode', default = 'preTrained', dest = 'mode', 
        help = 'choose to use a pre-trained model or train model on the fly using new raw data tablee')
    opts.add_argument('--trainSetPath', default = '/nfs/kitzman2/diaorch/trainsets/', dest = 'train_path', 
        help = 'set the path to the raw table and pklFiles folders')
    opts.add_argument('--setPrefix', default = 'trainset', dest = 'set_prefix', 
        help = 'file name prefix for raw table tsv files or pkl files')

    o = opts.parse_args()
    
    # from plugin_script import score_pair_model_learning
    scorer = eval(o.scoringClass)( o.genome, o.genomeDict, AFthresh=0.10, **vars(o) )
    generateMipSeqFxn =  eval(o.generateMipSeqFxn)(scorer) #getattr(scorer, o.generateMipSeqFxn)

    candsStore = pd.HDFStore( o.inArmStore, 'r' )
    outPairStore = pd.HDFStore( o.outPairStore, complib='zlib', complevel=9 )
    
    def newPairsTbl():
        tbl = OrderedDict()
        for col in ['targetidx','extarmidx','ligarmidx','chrom','start',
                    'end','strand_fill','mip_seq_5to3','score']: tbl[col]=[]
        return tbl

    tblnameExtArm,tblnameLigArm='%sextArm'%o.armBaseName,'%sligArm'%o.armBaseName
    tblExtArms = candsStore[ tblnameExtArm ]
    tblLigArms = candsStore[ tblnameLigArm ]

    tblnameExtHits,tblnameLigHits='%s_%s'%(tblnameExtArm,o.hitTblSufifx),'%s_%s'%(tblnameLigArm,o.hitTblSufifx)
    tblnameExtSnps,tblnameLigSnps='%s_%s'%(tblnameExtArm,o.snpTblSuffix),'%s_%s'%(tblnameLigArm,o.snpTblSuffix)

    tblExtHits,tblLigHits=candsStore[ tblnameExtHits ],candsStore[ tblnameLigHits ]
    tblExtSnps,tblLigSnps=candsStore[ tblnameExtSnps ],candsStore[ tblnameLigSnps ]

    tblPairsPlus = newPairsTbl()
    if o.generateDsMIPs:
        tblPairsMinus = newPairsTbl()

        tblnameLigArmMinus = candsStore.get_storer( '%sextArm'%o.armBaseName ).attrs['mippipe__dsmip_paired_table']
        tblnameExtArmMinus = candsStore.get_storer( '%sligArm'%o.armBaseName ).attrs['mippipe__dsmip_paired_table']

        tblLigArmMinus = candsStore[ tblnameLigArmMinus ]
        tblExtArmMinus = candsStore[ tblnameExtArmMinus ]

        tblnameExtMinusHits,tblnameLigMinusHits='%s_%s'%(tblnameExtArmMinus,o.hitTblSufifx),'%s_%s'%(tblnameLigArmMinus,o.hitTblSufifx)
        tblnameExtMinusSnps,tblnameLigMinusSnps='%s_%s'%(tblnameExtArmMinus,o.snpTblSuffix),'%s_%s'%(tblnameLigArmMinus,o.snpTblSuffix)

        tblExtMinusHits,tblLigMinusHits=candsStore[ tblnameExtMinusHits ],candsStore[ tblnameLigMinusHits ]
        tblExtMinusSnps,tblLigMinusSnps=candsStore[ tblnameExtMinusSnps ],candsStore[ tblnameLigMinusSnps ]

    tblnameTargets = candsStore.get_storer( '%sextArm'%o.armBaseName  ).attrs['mippipe__target_table']
    tblTargets = candsStore[ tblnameTargets ]

    nPairs=0

    nPairsSuppressed,nPairsPenalized=0,0

    for itgt in tblTargets.index:

        iextCurTgt = (tblExtArms.end + o.gapFillRange[1] - 1 >= tblTargets.start[itgt]) & \
                     (tblExtArms.end <  tblTargets.end[itgt])

        iligCurTgt = (tblLigArms.start - o.gapFillRange[1] + 1 <= tblTargets.end[itgt]) & \
                     (tblLigArms.start >= tblTargets.start[itgt])

        tblExtCurTgt = tblExtArms[ iextCurTgt ]
        tblLigCurTgt = tblLigArms[ iligCurTgt ]

        if o.generateDsMIPs:
            tblExtMinusCurTgt = tblExtArmMinus[ iligCurTgt ]
            tblLigMinusCurTgt = tblLigArmMinus[ iextCurTgt ]


        # tblExtCurTgt = tblExtArms[ tblExtArms.targetidx==itgt ]
        # tblLigCurTgt = tblLigArms[ tblLigArms.targetidx==itgt ]

        # if o.generateDsMIPs:
        #     tblExtMinusCurTgt = tblExtArmMinus[ tblExtArmMinus.targetidx==itgt ]
        #     tblLigMinusCurTgt = tblLigArmMinus[ tblLigArmMinus.targetidx==itgt ]


        if tblExtCurTgt.shape[0]==0:
            print 'warning no ext arms for target %d'%itgt
            continue

        if tblLigCurTgt.shape[0]==0:
            print 'warning no lig arms for target %d'%itgt
            continue

        print 'target %d - %d ext arms, %d lig arms'%( itgt, tblExtCurTgt.shape[0], tblLigCurTgt.shape[0] )
        
        for iextArm in tblExtCurTgt.index:

            armExt = tblExtCurTgt.ix[iextArm]
            tblLigCurElig = tblLigCurTgt[ (tblLigCurTgt.start >= armExt.end + o.gapFillRange[0] - 1)&
                                          (tblLigCurTgt.start <= armExt.end + o.gapFillRange[1] - 1)&
                                          ((armExt.arm_tm - tblLigCurTgt.arm_tm)>=o.tmDiffRange[0])&
                                          ((armExt.arm_tm - tblLigCurTgt.arm_tm)<=o.tmDiffRange[1]) ]

            hitExt = tblExtHits.ix[iextArm]
            snpExt = tblExtSnps.ix[iextArm]

            for iligArm in tblLigCurElig.index:

                armLig = tblLigCurElig.ix[ iligArm ]

                hitLig = tblLigHits.ix[iligArm]
                snpLig = tblLigSnps.ix[iligArm]

                # is gap fill within allowable lengths?
                #  no- continue
                # is the Tm diff within allowable range?
                #  no- continue

                # are we generating dsmip primers?  
                # if so, pull the paired arms
                # is that gap fill within allowable lengths, and do we care?
                #  no-continue

                if o.generateDsMIPs:
                    
                    armExtMinus,armLigMinus = tblExtMinusCurTgt.ix[iligArm],tblLigMinusCurTgt.ix[iextArm]

                    if not o.enforceGapFillOnTopOnly:
                        fillLenBtm = armLigMinus.start-armExtMinus.end-1      # = armLig.start-1-(armExt.end+1)+1 
                        tmdiffBtm = armExtMinus.arm_tm - armLigMinus.arm_tm

                        if fillLenBtm < o.gapFillRange[0] or fillLenBtm > o.gapFillRange[1]:
                            continue
                        elif tmdiffBtm < o.tmDiffRange[0] or tmdiffBtm > o.tmDiffRange[1]:
                            continue

                    hitExtMinus,snpExtMinus = tblExtMinusHits.ix[iligArm], tblExtMinusSnps.ix[iligArm]
                    hitLigMinus,snpLigMinus = tblLigMinusHits.ix[iextArm], tblLigMinusSnps.ix[iextArm]

                nPairs += 1
                if (nPairs%1000) == 0:
                    print '%d [%d/%d]..'%(nPairs,iextArm,tblExtCurTgt.shape[0])

                # if all else is good, fill in 
                corzExt = ( armExt.start, armExt.end )
                corzLig = ( armLig.start, armLig.end )

                if o.generateDsMIPs:
                    corzExtMinus = ( armExtMinus.start, armExtMinus.end )
                    corzLigMinus = ( armLigMinus.start, armLigMinus.end )
                else:
                    corzExtMinus = corzLigMinus = None

                # score pair
                if o.generateDsMIPs:
                    larmRec=[ armExt,armLig,armExtMinus,armLigMinus ]
                    lhitRec=[ hitExt,hitLig,hitExtMinus,hitLigMinus  ]
                    lsnpRec=[ snpExt,snpLig,snpExtMinus,snpLigMinus ]
                else:
                    larmRec=[ armExt,armLig,None,None ]
                    lhitRec=[ hitExt,hitLig,None,None  ]
                    lsnpRec=[ snpExt,snpLig,None,None ]

                pairScore = scorer.score_pair( larmRec, lhitRec, lsnpRec )

                suppressPair, scorePen, mipSeq = generateMipSeqFxn( larmRec, lhitRec, lsnpRec )

                if suppressPair:
                    nPairsSuppressed+=1
                    continue

                pairScore += scorePen

                if scorePen>0:
                    nPairsPenalized+=1

                # return (not suppress, penalty, seqMipCrop)

                tblPairsPlus[ 'targetidx' ].append(itgt)
                tblPairsPlus[ 'extarmidx' ].append(iextArm)
                tblPairsPlus[ 'ligarmidx' ].append(iligArm)
                tblPairsPlus[ 'chrom' ].append(tblTargets.chrom[itgt])
                tblPairsPlus[ 'start' ].append( armExt.end+1 )
                tblPairsPlus[ 'end' ].append( armLig.start-1 )
                tblPairsPlus[ 'strand_fill' ].append( 1 )


                tblPairsPlus[ 'mip_seq_5to3' ].append( mipSeq )
                tblPairsPlus[ 'score' ].append( pairScore )

                if o.generateDsMIPs:
                    tblPairsMinus[ 'targetidx' ].append(itgt)
                    tblPairsMinus[ 'extarmidx' ].append(iligArm)
                    tblPairsMinus[ 'ligarmidx' ].append(iextArm)
                    tblPairsMinus[ 'chrom' ].append(tblTargets.chrom[itgt])
                    tblPairsMinus[ 'start' ].append( armLigMinus.end+1 )
                    tblPairsMinus[ 'end' ].append( armExtMinus.start-1 )
                    tblPairsMinus[ 'strand_fill' ].append( -1 )

                    tblPairsMinus[ 'mip_seq_5to3' ].append( mipSeq )
                    tblPairsMinus[ 'score' ].append( pairScore )


    print '%d pairs suppressed'%nPairsSuppressed
    print '%d pairs penalized'%nPairsPenalized

    tblPairsPlus = pd.DataFrame( tblPairsPlus )
    outPairStore[ '%spairsPlus'%o.pairsBaseName ] = tblPairsPlus

    outPairStore.get_storer( '%spairsPlus'%o.pairsBaseName  ).attrs['mippipe__target_table'] = tblnameTargets
    outPairStore.get_storer( '%spairsPlus'%o.pairsBaseName  ).attrs['mippipe__ext_and_lig_arm_tables'] = (tblnameExtArm,tblnameLigArm)
    outPairStore.get_storer( '%spairsPlus'%o.pairsBaseName  ).attrs['mippipe__dsmip_paired_pair_table'] = '%spairsMinus'%o.pairsBaseName

    print 'made table %spairsPlus with %d entries'%(o.pairsBaseName,tblPairsPlus.shape[0])


    if o.generateDsMIPs:
        tblPairsMinus = pd.DataFrame( tblPairsMinus )
        outPairStore[ '%spairsMinus'%o.pairsBaseName ] = tblPairsMinus

        print 'made table %spairsMinus with %d entries'%(o.pairsBaseName,tblPairsMinus.shape[0])

        outPairStore.get_storer( '%spairsMinus'%o.pairsBaseName ).attrs['mippipe__target_table'] = tblnameTargets
        outPairStore.get_storer( '%spairsMinus'%o.pairsBaseName  ).attrs['mippipe__ext_and_lig_arm_tables'] = (tblnameExtArmMinus,tblnameLigArmMinus)    
        outPairStore.get_storer( '%spairsMinus'%o.pairsBaseName  ).attrs['mippipe__dsmip_paired_pair_table'] = '%spairsPlus'%o.pairsBaseName

