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

import Bio.SeqIO

from makeamip.capcommon import *

from makeamip.capcommon import *

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--inPairStore', dest='inPairStore')
    # opts.add_argument('--inArmStore', dest='inArmStore')
    opts.add_argument('--pairsBaseName', default='', dest='pairsBaseName')

    opts.add_argument('--outPairStore', dest='outPairStore')

    opts.add_argument('--mExtraColToFxnSummarize', default="{}", dest='mExtraColToFxnSummarize',
         help='expression of dict: colname --> fxn for extra columns to summarize')

    o = opts.parse_args()

    mExtraColToFxnSummarize = eval(o.mExtraColToFxnSummarize)

    pairStore = pd.HDFStore( o.inPairStore, 'r' )

    tblnamePairsPlus = '%spairsPlus'%o.pairsBaseName
    tblnamePairsMinus = '%spairsMinus'%o.pairsBaseName

    tblPairsPlus = pairStore[ tblnamePairsPlus ]
    tblPairsMinus = pairStore[ tblnamePairsMinus ]

    # lisDupd = tblPairsPlus.duplicated( subset=['chrom','start','end'] )
    lisDupd = tblPairsPlus.duplicated( subset=['extarmidx','ligarmidx'] )

    print 'input: %d pairs'%tblPairsPlus.shape[0]
    print 'found %d dups'%lisDupd.sum()

    liIsDup = tblPairsPlus.ix[lisDupd].index

    tblPairsPlusUnique = tblPairsPlus.ix[ ~lisDupd ].copy()
    tblPairsMinusUnique = tblPairsMinus.ix[ ~lisDupd ].copy()

    ctrdup=0

    for col in mExtraColToFxnSummarize:

        fxnSummarize = mExtraColToFxnSummarize[col]

        for idup in liIsDup:

            ctrdup+=1
            if (ctrdup%1000)==0 : 
                sys.stdout.write('%d..'%ctrdup)
                sys.stdout.flush()

            # lIsDupOfCur = (tblPairsPlus.chrom==tblPairsPlus.chrom[idup])&(tblPairsPlus.start==tblPairsPlus.start[idup])&(tblPairsPlus.end==tblPairsPlus.end[idup])
            # lIsRetainedDup = (tblPairsPlusUnique.chrom==tblPairsPlus.chrom[idup])&(tblPairsPlusUnique.start==tblPairsPlus.start[idup])&(tblPairsPlusUnique.end==tblPairsPlus.end[idup])

            lIsDupOfCur = (tblPairsPlus.extarmidx==tblPairsPlus.extarmidx[idup])&(tblPairsPlus.ligarmidx==tblPairsPlus.ligarmidx[idup])
            lIsRetainedDup = (tblPairsPlusUnique.extarmidx==tblPairsPlus.extarmidx[idup])&(tblPairsPlusUnique.ligarmidx==tblPairsPlus.ligarmidx[idup])

            liIsRetainedDup =  lIsRetainedDup.ix[lIsRetainedDup].index

            assert liIsRetainedDup.shape[0]==1

            ldupvals = tblPairsPlus.ix[ lIsDupOfCur, col ]

            summary = fxnSummarize(ldupvals)

            tblPairsPlusUnique.ix[ liIsRetainedDup[0], col ] = summary
            tblPairsMinusUnique.ix[ liIsRetainedDup[0], col ] = summary


    print '\nwriting %d pairs'%tblPairsPlusUnique.shape[0]

    outPairStore = pd.HDFStore( o.outPairStore, 'w', complib='zlib', complevel=9 )

    outPairStore[ '%spairsPlus'%o.pairsBaseName ] = tblPairsPlusUnique
    outPairStore[ '%spairsMinus'%o.pairsBaseName ] = tblPairsMinusUnique

    copyMetadata( pairStore, outPairStore, '%spairsPlus'%o.pairsBaseName, '%spairsPlus'%o.pairsBaseName )
    copyMetadata( pairStore, outPairStore, '%spairsMinus'%o.pairsBaseName, '%spairsMinus'%o.pairsBaseName )
