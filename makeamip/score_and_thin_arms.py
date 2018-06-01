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

import primer3.bindings as p3b
# import primer3.wrappers as p3w

import pysam
# import pybedtools as pbt

from makeamip.capcommon import *
#     

from makeamip.capcommon import *

"""

step through arms sorted by their gap-fill proximal coordinate
group by [pickStep] bases in this coordinate
pick the best [pickNbest] arms

"""


if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--genome', dest='genome', default='/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa',
                      help='path to indexed genome file')
    
    opts.add_argument('--genomeDict', dest='genomeDict', default='/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa.dict', 
                      help='path to genome contig file dictionary, eg from Picard CreateSequenceDictionary')

    opts.add_argument('--inArmStore', dest='inArmStore')

    opts.add_argument('--outArmStore', dest='outArmStore')

    opts.add_argument('--armBaseName', default='', dest='armBaseName')
    
    opts.add_argument('--hitTblSufifx', default='', dest='hitTblSufifx')
    opts.add_argument('--snpTblSuffix', default='', dest='snpTblSuffix')

    opts.add_argument('--scoringClass',default='ScoringModelLearned',dest='scoringClass')

    opts.add_argument('--scoreFxnExt',default='score_topext_arm_model1',dest='scoreFxnExt')
    opts.add_argument('--scoreFxnLig',default='score_toplig_arm_model1',dest='scoreFxnLig')

    opts.add_argument('--pickNbest',default=1,type=int,dest='pickNbest')
    opts.add_argument('--pickStep',default=1,type=int,dest='pickStep')

    opts.add_argument('--learningModel', default = 'GradientBoosting', dest = 'clfName', 
        help = 'name a classifier, GradientBoosting or AdaBoost, for machine learning algorithm')
    opts.add_argument('--trainingMode', default = 'preTrained', dest = 'mode', 
        help = 'choose to use a pre-trained model or train model on the fly using new raw data tablee')
    opts.add_argument('--trainSetPath', default = '/nfs/kitzman2/diaorch/trainsets/', dest = 'train_path', 
        help = 'set the path to the raw table and pklFiles folders')
    opts.add_argument('--setPrefix', default = 'trainset', dest = 'set_prefix', 
        help = 'file name prefix for raw table tsv files or pkl files')
    o = opts.parse_args()
    
    scorer = eval(o.scoringClass)( o.genome, o.genomeDict, AFthresh=0.10, **vars(o) )

    candsStore = pd.HDFStore( o.inArmStore, 'r' )
    outCandsStore = pd.HDFStore( o.outArmStore, 'a' )

    outCandsStore[ '%stargets'%o.armBaseName ] = candsStore[ '%stargets'%o.armBaseName ]
    copyMetadata( candsStore, outCandsStore, '%stargets'%o.armBaseName, '%stargets'%o.armBaseName )

    # load plus strand arms & info
    tblnameExtArm,tblnameLigArm='%sextArm'%o.armBaseName,'%sligArm'%o.armBaseName
    tblExtArms = candsStore[ tblnameExtArm ]
    tblLigArms = candsStore[ tblnameLigArm ]

    tblnameExtHits,tblnameLigHits='%s_%s'%(tblnameExtArm,o.hitTblSufifx),'%s_%s'%(tblnameLigArm,o.hitTblSufifx)
    tblnameExtSnps,tblnameLigSnps='%s_%s'%(tblnameExtArm,o.snpTblSuffix),'%s_%s'%(tblnameLigArm,o.snpTblSuffix)

    tblExtHits,tblLigHits=candsStore[ tblnameExtHits ],candsStore[ tblnameLigHits ]
    tblExtSnps,tblLigSnps=candsStore[ tblnameExtSnps ],candsStore[ tblnameLigSnps ]

    tblnameLigArmMinus = candsStore.get_storer( '%sextArm'%o.armBaseName ).attrs['mippipe__dsmip_paired_table']
    tblnameExtArmMinus = candsStore.get_storer( '%sligArm'%o.armBaseName ).attrs['mippipe__dsmip_paired_table']

    # load minus strand arms & info
    tblLigArmMinus = candsStore[ tblnameLigArmMinus ]
    tblExtArmMinus = candsStore[ tblnameExtArmMinus ]

    tblnameExtMinusHits,tblnameLigMinusHits='%s_%s'%(tblnameExtArmMinus,o.hitTblSufifx),'%s_%s'%(tblnameLigArmMinus,o.hitTblSufifx)
    tblnameExtMinusSnps,tblnameLigMinusSnps='%s_%s'%(tblnameExtArmMinus,o.snpTblSuffix),'%s_%s'%(tblnameLigArmMinus,o.snpTblSuffix)

    tblExtMinusHits,tblLigMinusHits=candsStore[ tblnameExtMinusHits ],candsStore[ tblnameLigMinusHits ]
    tblExtMinusSnps,tblLigMinusSnps=candsStore[ tblnameExtMinusSnps ],candsStore[ tblnameLigMinusSnps ]

    liKeepExt = []
    tblExtArmsSortByInterior = tblExtArms.sort_values( by= ['chrom','end'] )
    ii=0
    while ii < tblExtArmsSortByInterior.shape[0]:

        print ii,tblExtArmsSortByInterior.shape[0]
        # scan [here, here+pickStep-1]
        curArmInterior = tblExtArmsSortByInterior.iloc[ ii ]['end']
        liPickFrom = tblExtArmsSortByInterior.index[ 
            (tblExtArmsSortByInterior.end >= curArmInterior) & 
            (tblExtArmsSortByInterior.end <= curArmInterior + o.pickStep - 1) ]

        lScores = []
        for iii in xrange(len(liPickFrom)):
            lArmRecs=[ tblExtArms.ix[ liPickFrom[iii] ], 
                       tblLigArmMinus.ix[ liPickFrom[iii] ] ]

            lHitRecs=[ tblExtHits.ix[ liPickFrom[iii] ], 
                       tblLigMinusHits.ix[ liPickFrom[iii] ] ]

            lHitSnps=[ tblExtSnps.ix[ liPickFrom[iii] ], 
                       tblLigMinusSnps.ix[ liPickFrom[iii] ] ]

            # scCur = score_topext_arm_model1( lArmRecs, lHitRecs, lHitSnps )
            scCur = scorer.score_topext_arms( lArmRecs, lHitRecs, lHitSnps )

            lScores.append( scCur )

        liPicked = [ si[1] for si in sorted( zip(lScores,liPickFrom) )[:o.pickNbest] ]
        liKeepExt += liPicked

        # print ii
        ii += liPickFrom.shape[0]

    #

    print 'pruned ext arms from %d to %d'%( tblExtArms.shape[0], len(liKeepExt) )

    outCandsStore[ tblnameExtArm ] = tblExtArms.ix[ sorted(liKeepExt) ]
    outCandsStore[ tblnameExtHits ] = tblExtHits.ix[ sorted(liKeepExt) ]
    outCandsStore[ tblnameExtSnps ] = tblExtSnps.ix[ sorted(liKeepExt) ]

    copyMetadata( candsStore, outCandsStore, tblnameExtArm, tblnameExtArm )
    copyMetadata( candsStore, outCandsStore, tblnameExtHits, tblnameExtHits )
    copyMetadata( candsStore, outCandsStore, tblnameExtSnps, tblnameExtSnps )

    outCandsStore[ tblnameLigArmMinus ] = tblLigArmMinus.ix[ sorted(liKeepExt) ]
    outCandsStore[ tblnameLigMinusHits ] = tblLigMinusHits.ix[ sorted(liKeepExt) ]
    outCandsStore[ tblnameLigMinusSnps ] = tblLigMinusSnps.ix[ sorted(liKeepExt) ]

    copyMetadata( candsStore, outCandsStore, tblnameLigArmMinus, tblnameLigArmMinus )
    copyMetadata( candsStore, outCandsStore, tblnameLigMinusSnps, tblnameLigMinusSnps )
    copyMetadata( candsStore, outCandsStore, tblnameLigMinusHits, tblnameLigMinusHits )


    ###

    liKeepLig = []
    tblLigArmsSortByInterior = tblLigArms.sort_values( by=['chrom','start'] )
    ii=0
    while ii < tblLigArmsSortByInterior.shape[0]:

        # scan [here, here+pickStep-1]
        curArmInterior = tblLigArmsSortByInterior.iloc[ ii ]['start']
        liPickFrom = tblLigArmsSortByInterior.index[ 
            (tblLigArmsSortByInterior.start >= curArmInterior) & 
            (tblLigArmsSortByInterior.start <= curArmInterior + o.pickStep - 1) ]

        lScores = []
        for iii in xrange(len(liPickFrom)):
            lArmRecs=[ tblLigArms.ix[ liPickFrom[iii] ], 
                       tblExtArmMinus.ix[ liPickFrom[iii] ] ]

            lHitRecs=[ tblLigHits.ix[ liPickFrom[iii] ], 
                       tblExtMinusHits.ix[ liPickFrom[iii] ] ]

            lHitSnps=[ tblLigSnps.ix[ liPickFrom[iii] ], 
                       tblExtMinusSnps.ix[ liPickFrom[iii] ] ]

            # scCur = score_toplig_arm_model1( lArmRecs, lHitRecs, lHitSnps )
            scCur = scorer.score_toplig_arms( lArmRecs, lHitRecs, lHitSnps )

            lScores.append( scCur )

        liPicked = [ si[1] for si in sorted( zip(lScores,liPickFrom) )[:o.pickNbest] ]
        liKeepLig += liPicked

        # print ii
        ii += liPickFrom.shape[0]


    print 'pruned ext arms from %d to %d'%( tblLigArms.shape[0], len(liKeepLig) )

    outCandsStore[ tblnameLigArm ] = tblLigArms.ix[ sorted(liKeepLig) ]
    outCandsStore[ tblnameLigHits ] = tblLigHits.ix[ sorted(liKeepLig) ]
    outCandsStore[ tblnameLigSnps ] = tblLigSnps.ix[ sorted(liKeepLig) ]

    copyMetadata( candsStore, outCandsStore, tblnameLigArm, tblnameLigArm )
    copyMetadata( candsStore, outCandsStore, tblnameLigHits, tblnameLigHits )
    copyMetadata( candsStore, outCandsStore, tblnameLigSnps, tblnameLigSnps )

    outCandsStore[ tblnameExtArmMinus ] = tblExtArmMinus.ix[ sorted(liKeepLig) ]
    outCandsStore[ tblnameExtMinusHits ] = tblExtMinusHits.ix[ sorted(liKeepLig) ]
    outCandsStore[ tblnameExtMinusSnps ] = tblExtMinusSnps.ix[ sorted(liKeepLig) ]

    copyMetadata( candsStore, outCandsStore, tblnameExtArmMinus, tblnameExtArmMinus )
    copyMetadata( candsStore, outCandsStore, tblnameExtMinusSnps, tblnameExtMinusSnps )
    copyMetadata( candsStore, outCandsStore, tblnameExtMinusHits, tblnameExtMinusHits )
