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

    opts.add_argument('--genome', dest='genome',
                      help='path to indexed genome file')
    
    opts.add_argument('--genomeDict', dest='genomeDict',
                      help='path to genome contig file dictionary, eg from Picard CreateSequenceDictionary')


    opts.add_argument('--inArmStore', dest='inArmStore')

    opts.add_argument('--inPairStore', dest='inPairStore')

    opts.add_argument('--outPairStore', dest='outPairStore')

    opts.add_argument('--armBaseName', default='', dest='armBaseName')
    opts.add_argument('--pairBasename', default='', dest='pairBasename')
    
    opts.add_argument('--hitTblSuffix', default='', dest='hitTblSuffix')
    opts.add_argument('--snpTblSuffix', default='', dest='snpTblSuffix')

    opts.add_argument('--scoringClass',default='ScoringModelLearned',dest='scoringClass')
    opts.add_argument('--generateMipSeqFxn',default='makeseq_dsmip_eari',dest='generateMipSeqFxn')
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



    # # scoreFxn = eval(o.scoreFxn)
    # generateMipSeqFxn = eval(o.generateMipSeqFxn)(None)

    sfGenome = pysam.Fastafile( o.genome )    
    mCtgLen = parseContigLengths( o.genomeDict )

    candsStore = pd.HDFStore( o.inArmStore, 'r' )

    inPairStore = pd.HDFStore( o.inPairStore, 'r' )

    outPairStore = pd.HDFStore( o.outPairStore, 'w', complib='zlib', complevel=9 )
    
    chromLast,seqChromLast=None,None

    tblnameExtArm,tblnameLigArm='%sextArm'%o.armBaseName,'%sligArm'%o.armBaseName
    tblExtArms = candsStore[ tblnameExtArm ]
    tblLigArms = candsStore[ tblnameLigArm ]

    tblnameExtHits,tblnameLigHits='%s_%s'%(tblnameExtArm,o.hitTblSuffix),'%s_%s'%(tblnameLigArm,o.hitTblSuffix)
    tblnameExtSnps,tblnameLigSnps='%s_%s'%(tblnameExtArm,o.snpTblSuffix),'%s_%s'%(tblnameLigArm,o.snpTblSuffix)

    tblExtHits,tblLigHits=candsStore[ tblnameExtHits ],candsStore[ tblnameLigHits ]
    tblExtSnps,tblLigSnps=candsStore[ tblnameExtSnps ],candsStore[ tblnameLigSnps ]

    tblPairsPlus = inPairStore[ '%spairsPlus'%o.pairBasename ].copy()
    tblPairsMinus = inPairStore[ '%spairsMinus'%o.pairBasename ].copy()

    # tblnameLigArmMinus = candsStore.get_storer( '%sextArm'%o.armBaseName ).attrs['mippipe__dsmip_paired_table']
    # tblnameExtArmMinus = candsStore.get_storer( '%sligArm'%o.armBaseName ).attrs['mippipe__dsmip_paired_table']

    tblnameExtArmMinus = '%sextArmMinus'%o.armBaseName
    tblnameLigArmMinus = '%sligArmMinus'%o.armBaseName 

    tblLigArmMinus = candsStore[ tblnameLigArmMinus ]
    tblExtArmMinus = candsStore[ tblnameExtArmMinus ]

    tblnameExtMinusHits,tblnameLigMinusHits='%s_%s'%(tblnameExtArmMinus,o.hitTblSuffix),'%s_%s'%(tblnameLigArmMinus,o.hitTblSuffix)
    tblnameExtMinusSnps,tblnameLigMinusSnps='%s_%s'%(tblnameExtArmMinus,o.snpTblSuffix),'%s_%s'%(tblnameLigArmMinus,o.snpTblSuffix)

    tblExtMinusHits,tblLigMinusHits=candsStore[ tblnameExtMinusHits ],candsStore[ tblnameLigMinusHits ]
    tblExtMinusSnps,tblLigMinusSnps=candsStore[ tblnameExtMinusSnps ],candsStore[ tblnameLigMinusSnps ]

    nPairs=0
    nPairsSuppressed,nPairsPenalized=0,0

    for ipair,rpair in tblPairsPlus.iterrows():

        rpairMinus = tblPairsMinus.ix[ipair]

        if chromLast != rpair['chrom']:
            chromLast = rpair['chrom']
            seqChromLast = sfGenome.fetch( str( chromLast ) )
        
        armExt = tblExtArms.ix[ rpair['extarmidx'] ]
        armLig = tblLigArms.ix[ rpair['ligarmidx'] ]

        armExtMinus = tblExtArmMinus.ix[ rpairMinus['extarmidx'] ]
        armLigMinus = tblLigArmMinus.ix[ rpairMinus['ligarmidx'] ]

        hitExt = tblExtHits.ix[ rpair['extarmidx'] ]
        snpExt = tblExtSnps.ix[ rpair['extarmidx'] ] 
        hitLig = tblLigHits.ix[ rpair['ligarmidx'] ]
        snpLig = tblLigSnps.ix[ rpair['ligarmidx'] ] 

        hitExtMinus = tblExtMinusHits.ix[ rpairMinus['extarmidx'] ]
        snpExtMinus = tblExtMinusSnps.ix[ rpairMinus['extarmidx'] ] 
        hitLigMinus = tblLigMinusHits.ix[ rpairMinus['ligarmidx'] ]
        snpLigMinus = tblLigMinusSnps.ix[ rpairMinus['ligarmidx'] ] 

        # score pair
        larmRec=[ armExt,armLig,armExtMinus,armLigMinus ]
        lhitRec=[ hitExt,hitLig,hitExtMinus,hitLigMinus  ]
        lsnpRec=[ snpExt,snpLig,snpExtMinus,snpLigMinus ]

        pairScore = rpair['score']

        # pairScore = scoreFxn( larmRec, lhitRec, lsnpRec )
        pairScore = scorer.score_pair( larmRec, lhitRec, lsnpRec )

        suppressPair, scorePen, mipSeq = generateMipSeqFxn( larmRec, lhitRec, lsnpRec )

        if suppressPair:
            print 'warning - supposed to suppress row but didnt ',ipair
            nPairsSuppressed+=1

        print 'old score %.2f new score %.2f'%( rpair['score'], pairScore )

        pairScore += scorePen

        if scorePen>0:
            nPairsPenalized+=1

        tblPairsPlus.ix[ipair,'mip_seq_5to3'] = mipSeq
        tblPairsPlus.ix[ipair,'score'] = pairScore
        tblPairsMinus.ix[ipair,'mip_seq_5to3'] = mipSeq
        tblPairsMinus.ix[ipair,'score'] = pairScore

    print '%d pairs suppressed'%nPairsSuppressed
    print '%d pairs penalized'%nPairsPenalized

    outPairStore[ '%spairsPlus'%o.pairBasename ] = tblPairsPlus
    outPairStore[ '%spairsMinus'%o.pairBasename ] = tblPairsMinus

    copyMetadata( inPairStore, outPairStore, '%spairsPlus'%o.pairBasename,'%spairsPlus'%o.pairBasename)
    copyMetadata( inPairStore, outPairStore, '%spairsMinus'%o.pairBasename ,'%spairsMinus'%o.pairBasename )
