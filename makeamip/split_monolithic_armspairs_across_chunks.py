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

import pybedtools as pbt

import Bio.SeqIO

from makeamip.capcommon import *

if __name__ == '__main__':
    
    def intRange(s):
        try:
            lower,upper=map(int, s.split(','))
            return lower,upper
        except:
            raise argparse.ArgumentTypeError('range must be x,y')


    opts = argparse.ArgumentParser()
   
    opts.add_argument('--key', dest='key')

    opts.add_argument('--colnameBed', dest='colnameBed')

    opts.add_argument('--colnameOutArmTbl', dest='colnameOutArmTbl')
    opts.add_argument('--colnameOutPairTbl', dest='colnameOutPairTbl')
    
    opts.add_argument('--inArmStore', dest='inArmStore')
    opts.add_argument('--inPairStore', dest='inPairStore')

    opts.add_argument('--armBasename', default='', dest='armBasename')
    opts.add_argument('--pairBasename', default='', dest='pairBasename')

    # opts.add_argument('--hitTblSuffix', dest='hitTblSuffix')
    # opts.add_argument('--snpTblSuffix', dest='snpTblSuffix')

    o = opts.parse_args()

    tblk = pd.read_csv(o.key,sep='\t')

    print 'loading pairs,arms'

    armStoreIn = pd.HDFStore( o.inArmStore, 'r' )

    tblInArmExtPlus = armStoreIn['%sextArm'%o.armBasename]
    print '%d ext arms (+)'%tblInArmExtPlus.shape[0]
    tblInArmLigPlus = armStoreIn['%sligArm'%o.armBasename]
    print '%d lig arms (+)'%tblInArmLigPlus.shape[0]
    tblInArmExtMinus = armStoreIn['%sextArmMinus'%o.armBasename]
    print '%d ext arms (-)'%tblInArmExtMinus.shape[0]
    tblInArmLigMinus = armStoreIn['%sligArmMinus'%o.armBasename]
    print '%d lig arms (-)'%tblInArmLigMinus.shape[0]

    pairStoreIn = pd.HDFStore( o.inPairStore, 'r' )

    tblInPairsPlus = pairStoreIn['%spairsPlus'%o.pairBasename]
    print '%d pairs (+)'%tblInPairsPlus.shape[0]
    tblInPairsMinus = pairStoreIn['%spairsMinus'%o.pairBasename]
    print '%d pairs (-)'%tblInPairsMinus.shape[0]


    for _,r in tblk.iterrows():
        
        btChunk = pbt.BedTool( r[o.colnameBed] )

        outArmStore = pd.HDFStore( r[o.colnameOutArmTbl], 'w', complib='zlib', complevel=9 )
        outPairStore = pd.HDFStore(  r[o.colnameOutPairTbl], 'w', complib='zlib', complevel=9 )

        tblOutPairsPlus = bed_intersect_dataframe( 
            btChunk,
            tblInPairsPlus )

        tblOutPairsMinus = tblInPairsMinus.ix[ tblOutPairsPlus.index ]

        outPairStore[ '%spairsPlus'%o.pairBasename ] = tblOutPairsPlus
        outPairStore[ '%spairsMinus'%o.pairBasename ] = tblOutPairsMinus

        copyMetadata( pairStoreIn, outPairStore, '%spairsPlus'%o.pairBasename, '%spairsPlus'%o.pairBasename )
        copyMetadata( pairStoreIn, outPairStore, '%spairsMinus'%o.pairBasename, '%spairsMinus'%o.pairBasename )

        print 'chunk:'
        print r
        print '%d,%d pairs'%( tblOutPairsPlus.shape[0], tblOutPairsMinus.shape[0] )

        tblOutTargets = armStoreIn['%stargets'%o.armBasename].ix[ 
            tblOutPairsPlus.targetidx.unique() ]

        outArmStore[ '%stargets'%o.pairBasename  ] = tblOutTargets

        print '%d targets'%tblOutTargets.shape[0]

        tblOutExtArmPlus = tblInArmExtPlus.ix[ 
            tblOutPairsPlus.extarmidx.unique() ]
        tblOutLigArmPlus = tblInArmLigPlus.ix[
            tblOutPairsPlus.ligarmidx.unique() ]

        print '%d, %d ext, lig arms plus'%( tblOutExtArmPlus.shape[0], tblOutLigArmPlus.shape[0] )

        outArmStore[ '%sextArm'%o.armBasename ] = tblOutExtArmPlus
        outArmStore[ '%sligArm'%o.armBasename ] = tblOutLigArmPlus

        copyMetadata( armStoreIn, outArmStore, '%sextArm'%o.armBasename, '%sextArm'%o.armBasename )
        copyMetadata( armStoreIn, outArmStore, '%sligArm'%o.armBasename, '%sligArm'%o.armBasename )

        tblOutExtArmMinus = tblInArmExtMinus.ix[ 
            tblOutPairsMinus.extarmidx.unique() ]
        tblOutLigArmMinus = tblInArmLigMinus.ix[ 
            tblOutPairsMinus.ligarmidx.unique() ]

        outArmStore[ '%sextArmMinus'%o.armBasename ] = tblOutExtArmMinus
        outArmStore[ '%sligArmMinus'%o.armBasename ] = tblOutLigArmMinus

        print '%d, %d ext, lig arms minus'%( tblOutExtArmMinus.shape[0], tblOutLigArmMinus.shape[0] )

        copyMetadata( armStoreIn, outArmStore, '%sextArmMinus'%o.armBasename, '%sextArmMinus'%o.armBasename )
        copyMetadata( armStoreIn, outArmStore, '%sligArmMinus'%o.armBasename, '%sligArmMinus'%o.armBasename )


