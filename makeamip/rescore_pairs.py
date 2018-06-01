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

from makeamip.capcommon import *
#     

from makeamip.capcommon import *

"""
rescore probes using a different scoring function

"""
    
if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--genome', dest='genome',
                      help='path to indexed genome file')
    
    opts.add_argument('--genomeDict', dest='genomeDict',
                      help='path to genome contig file dictionary, eg from Picard CreateSequenceDictionary')

    opts.add_argument('--inArmStore', dest='inArmStore')
    opts.add_argument('--inPairStore', dest='inPairStore')
    opts.add_argument('--pairsBaseName', default='', dest='pairsBaseName')
    opts.add_argument('--outPairStore', dest='outPairStore')

    opts.add_argument('--hitTblSufifx', default='', dest='hitTblSufifx')
    opts.add_argument('--snpTblSuffix', default='', dest='snpTblSuffix')

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

    scorer = eval(o.scoringClass)( o.genome, o.genomeDict, AFthresh=0.10, **vars(o) )

    pairStore = pd.HDFStore( o.inPairStore, 'r' )
    armStore = pd.HDFStore( o.inArmStore, 'r' )

    tblnamePairsPlus = '%spairsPlus'%o.pairsBaseName
    tblnamePairsMinus = '%spairsMinus'%o.pairsBaseName

    tblPairsPlus = pairStore[ tblnamePairsPlus ]
    tblPairsMinus = pairStore[ tblnamePairsMinus ]

    (tblnameExtArmPlus,tblnameLigArmPlus)=pairStore.get_storer( tblnamePairsPlus ).attrs['mippipe__ext_and_lig_arm_tables']
    (tblnameExtArmMinus,tblnameLigArmMinus)=pairStore.get_storer( tblnamePairsMinus ).attrs['mippipe__ext_and_lig_arm_tables']

    tblExtArmsPlus = armStore[ tblnameExtArmPlus ]  
    tblLigArmsPlus = armStore[ tblnameLigArmPlus ]     
    tblExtArmsMinus = armStore[ tblnameExtArmMinus ]  
    tblLigArmsMinus = armStore[ tblnameLigArmMinus ]  

    tblnameExtHits,tblnameLigHits='%s_%s'%(tblnameExtArmPlus,o.hitTblSufifx),'%s_%s'%(tblnameLigArmPlus,o.hitTblSufifx)
    tblnameExtSnps,tblnameLigSnps='%s_%s'%(tblnameExtArmPlus,o.snpTblSuffix),'%s_%s'%(tblnameLigArmPlus,o.snpTblSuffix)

    tblExtHits,tblLigHits=armStore[ tblnameExtHits ],armStore[ tblnameLigHits ]
    tblExtSnps,tblLigSnps=armStore[ tblnameExtSnps ],armStore[ tblnameLigSnps ]

    tblnameExtMinusHits,tblnameLigMinusHits='%s_%s'%(tblnameExtArmMinus,o.hitTblSufifx),'%s_%s'%(tblnameLigArmMinus,o.hitTblSufifx)
    tblnameExtMinusSnps,tblnameLigMinusSnps='%s_%s'%(tblnameExtArmMinus,o.snpTblSuffix),'%s_%s'%(tblnameLigArmMinus,o.snpTblSuffix)

    tblExtMinusHits,tblLigMinusHits=armStore[ tblnameExtMinusHits ],armStore[ tblnameLigMinusHits ]
    tblExtMinusSnps,tblLigMinusSnps=armStore[ tblnameExtMinusSnps ],armStore[ tblnameLigMinusSnps ]

    ###

    outPairStore = pd.HDFStore( o.outPairStore, 'w', complib='zlib', complevel=9 )

    for ip,rp in tblPairsPlus.iterrows():

        rm = tblPairsMinus.ix[ip]

        armExtTop = tblExtArmsPlus.ix[ rp['extarmidx'] ]
        armLigTop = tblLigArmsPlus.ix[ rp['ligarmidx'] ]

        armExtBtm = tblExtArmsMinus.ix[ rm['extarmidx'] ]
        armLigBtm = tblLigArmsMinus.ix[ rm['ligarmidx'] ]

        hitExtTop = tblExtHits.ix[ rp['extarmidx' ] ]
        snpExtTop = tblExtSnps.ix[ rp['extarmidx' ] ]

        hitLigTop = tblLigHits.ix[ rp['ligarmidx']   ]
        snpLigTop = tblLigSnps.ix[ rp['ligarmidx']   ]

        hitExtBtm = tblExtMinusHits.ix[ rm['extarmidx'] ]
        snpExtBtm = tblExtMinusSnps.ix[ rm['extarmidx'] ]

        hitLigBtm = tblLigMinusHits.ix[ rm['ligarmidx']  ]
        snpLigBtm = tblLigMinusSnps.ix[ rm['ligarmidx']  ]

        larmRec=[ armExtTop,armLigTop,armExtBtm,armLigBtm ]
        lhitRec=[ hitExtTop,hitLigTop,hitExtBtm,hitLigBtm  ]
        lsnpRec=[ snpExtTop,snpLigTop,snpExtBtm,snpLigBtm ]

        pairScore = scorer.score_pair( larmRec, lhitRec, lsnpRec )
              
        oldScorePlus=tblPairsPlus.ix[ip,'score']
        print '%s   %.3f --> %.3f'%( repr(ip), oldScorePlus , pairScore )

        tblPairsPlus.ix[ip,'score'] = pairScore
        tblPairsMinus.ix[ip,'score'] = pairScore

    outPairStore[ '%spairsPlus'%o.pairsBaseName ] = tblPairsPlus
    outPairStore[ '%spairsMinus'%o.pairsBaseName ] = tblPairsMinus

    copyMetadata( pairStore, outPairStore, '%spairsPlus'%o.pairsBaseName, '%spairsPlus'%o.pairsBaseName )
    copyMetadata( pairStore, outPairStore, '%spairsMinus'%o.pairsBaseName, '%spairsMinus'%o.pairsBaseName )
