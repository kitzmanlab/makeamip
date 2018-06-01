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

if __name__ == '__main__':


    """
    update score for arm pairs based upon observed number of reads
    (ie capture efficiency) from a previous experiment

    delete outright probes with no/little coverage


    """
    
if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--hitsTbl', dest='hitsTbl',
                      help='table of MIP probes, eg from probe_hit_count_from_bam.py\n \
                       expected columns:mip_index, probe_hits_OK')

    opts.add_argument('--fxnProbeToDelete',
                      default="lambda r:r['probe_hits_OK']<5",
                      dest='fxnProbeToDelete',
                      help='function( probe row from hits table ) to true/false to delete probe')

    opts.add_argument('--fxnNewScore',
                      default="lambda r,existing_score:existing_score*(0. if r['probe_hits_OK']>10 else 1.)",
                      dest='fxnNewScore',
                      help='function( probe row from hits table, existing_score ) to new score')

    opts.add_argument('--fxnProbeToFlag',
                      default="lambda r:r['probe_hits_OK']>=10",
                      dest='fxnProbeToFlag',
                      help='function( probe row from hits table ) to true/false to flag probe in new column')

    opts.add_argument('--colFlagGoodProbes',
                      default="flagged_good",
                      dest='colFlagGoodProbes',
                      help='column to flag probes')

    opts.add_argument('--inPairStore', dest='inPairStore')
    # opts.add_argument('--inArmStore', dest='inArmStore')
    opts.add_argument('--pairsBaseName', default='', dest='pairsBaseName')

    opts.add_argument('--outPairStore', dest='outPairStore')

    o = opts.parse_args()

    pairStore = pd.HDFStore( o.inPairStore, 'r' )

    # armStore = pd.HDFStore( o.inArmStore, 'r' )

    tblnamePairsPlus = '%spairsPlus'%o.pairsBaseName
    tblnamePairsMinus = '%spairsMinus'%o.pairsBaseName

    tblPairsPlus = pairStore[ tblnamePairsPlus ]
    tblPairsMinus = pairStore[ tblnamePairsMinus ]

    # (tblnameExtArmPlus,tblnameLigArmPlus)=pairStore.get_storer( tblnamePairsPlus ).attrs['mippipe__ext_and_lig_arm_tables']
    # (tblnameExtArmMinus,tblnameLigArmMinus)=pairStore.get_storer( tblnamePairsMinus ).attrs['mippipe__ext_and_lig_arm_tables']

    # tblExtArmsPlus = armStore[ tblnameExtArmPlus ]  
    # tblLigArmsPlus = armStore[ tblnameLigArmPlus ]     
    # tblExtArmsMinus = armStore[ tblnameExtArmMinus ]  
    # tblLigArmsMinus = armStore[ tblnameLigArmMinus ]  

    outPairStore = pd.HDFStore( o.outPairStore, 'w', complib='zlib', complevel=9 )


    fxnProbeToDelete = eval(o.fxnProbeToDelete)
    fxnProbeToFlag = eval(o.fxnProbeToFlag)
    fxnNewScore = eval(o.fxnNewScore)

    tblHits = pd.read_csv( o.hitsTbl, sep='\t' )

    tblPairsPlus['__todel'] = False
    tblPairsPlus[o.colFlagGoodProbes] = False
    tblPairsMinus['__todel'] = False
    tblPairsMinus[o.colFlagGoodProbes] = False

    for _,r in tblHits.iterrows():
        doDelete = fxnProbeToDelete( r )
        doFlag = fxnProbeToFlag( r )

        assert not (doDelete and doFlag),r

        ixmip = r['mip_index']

        rlookedupPlus = tblPairsPlus.ix[ ixmip ]

        assert (r['chrom'],r['start'],r['end']) == ( r['chrom'],r['start'],r['end'] )

        newScore = fxnNewScore( r, rlookedupPlus['score'] )
        tblPairsPlus.ix[ ixmip, 'score' ] = newScore
        tblPairsMinus.ix[ ixmip, 'score' ] = newScore

        if doDelete:
            tblPairsPlus.ix[ ixmip, '__todel' ] = True
            tblPairsMinus.ix[ ixmip, '__todel' ] = True

        if doFlag:
            tblPairsPlus.ix[ ixmip, o.colFlagGoodProbes ] = True
            tblPairsMinus.ix[ ixmip, o.colFlagGoodProbes ] = True

    nInput=tblPairsPlus.shape[0]

    tblPairsPlus = tblPairsPlus.ix[ ~tblPairsPlus['__todel'] ] 
    tblPairsMinus = tblPairsMinus.ix[ ~tblPairsMinus['__todel'] ] 

    print 'input: %d probes'%nInput
    print '%d probes flagged'%tblPairsPlus[o.colFlagGoodProbes].sum()
    nDiff=nInput-tblPairsPlus.shape[0]
    print '%d probes deleted'%nDiff

    outPairStore[ '%spairsPlus'%o.pairsBaseName ] = tblPairsPlus
    outPairStore[ '%spairsMinus'%o.pairsBaseName ] = tblPairsMinus

    copyMetadata( pairStore, outPairStore, '%spairsPlus'%o.pairsBaseName, '%spairsPlus'%o.pairsBaseName )
    copyMetadata( pairStore, outPairStore, '%spairsMinus'%o.pairsBaseName, '%spairsMinus'%o.pairsBaseName )
