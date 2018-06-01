#! /usr/bin/env python

import sys

import os
import os.path
import argparse
from collections import defaultdict, OrderedDict

import numpy as np
import numpy.random as rand
import re

import pybedtools as pbt

import pandas as pd

from makeamip.capcommon import *

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()
   
    opts.add_argument('--pairStore', dest='pairStore')
    opts.add_argument('--armStore', dest='armStore')

    opts.add_argument('--table', dest='table')

    opts.add_argument('--readLenFwd', type=int, default=75, dest='readLenFwd')
    opts.add_argument('--readLenRev', type=int, default=75, dest='readLenRev')
    
    opts.add_argument('--bedOut', dest='bedOut')
        
    o = opts.parse_args()

    pairStore = pd.HDFStore( o.pairStore, 'r' )
        
    assert o.table in pairStore, 'cannot find table %s in %s'%(o.table, o.pairStore)
    tblPairs = pairStore[ o.table ]       
    
    (tblnameExtArm,tblnameLigArm) = pairStore.get_storer( o.table ).attrs['mippipe__ext_and_lig_arm_tables']

    armStore = pd.HDFStore( o.armStore, 'r' )

    tblArmsExt = armStore[tblnameExtArm]
    tblArmsLig = armStore[tblnameLigArm]

    btReadable = []
    
    for i,r in tblPairs.iterrows():

        armExt = tblArmsExt.ix[ r.extarmidx ]
        armLig = tblArmsLig.ix[ r.ligarmidx ]

        nbpExt=armExt.end-armExt.start+1
        nbpLig=armLig.end-armLig.start+1

        corngRead1 = ( r['start'], 
                       min( r['start'] + (o.readLenFwd - nbpExt) - 1,  r['end'] ))

        corngRead2 = ( max( r['start'], 
                            r['end'] - (o.readLenRev - nbpLig) + 1 ),
                       r['end'] )

        btReadable.append( [ r['chrom'], corngRead1[0], corngRead1[1]+1]  )
        btReadable.append( [ r['chrom'], corngRead2[0], corngRead2[1]+1] )

    btReadable = sorted( btReadable )

    btReadable = pbt.BedTool(btReadable)

    btReadableMerge = btReadable.merge()
    btReadableMerge.saveas(o.bedOut)