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

from makeamip.capcommon import *

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()
   
    opts.add_argument('--inStore', dest='inStore')
    opts.add_argument('--table', dest='table')

    opts.add_argument('--fxnDesc',
        default='lambda c:"%s_%s_%d_%d_%s"%(c.name,c["chrom"],c["start"],c["end"]," +-"[c["strand_fill"]])',
        dest='fxnDesc')

    opts.add_argument('--armStore', dest='armStore')

    opts.add_argument('--fxnScore',
        default='lambda c:c.score',
        dest='fxnScore')
    
    opts.add_argument('--bedOut', dest='bedOut')
        
    o = opts.parse_args()

    fxnDesc = eval(o.fxnDesc)
    fxnScore = eval(o.fxnScore)
         
    pairStore = pd.HDFStore( o.inStore, 'r' )
        
    bedOut = open(o.bedOut,'w')

    assert o.table in pairStore, 'cannot find table %s in %s'%(o.table, o.inStore)
    tblPairs = pairStore[ o.table ]       
    
    (tblnameExtArm,tblnameLigArm) = pairStore.get_storer( o.table ).attrs['mippipe__ext_and_lig_arm_tables']

    armStore = pd.HDFStore( o.armStore, 'r' )

    tblArmsExt = armStore[tblnameExtArm]
    tblArmsLig = armStore[tblnameLigArm]
    
    for i,r in tblPairs.iterrows():

        desc = fxnDesc(r)
        score = fxnScore(r)

        armExt = tblArmsExt.ix[ r.extarmidx ]
        armLig = tblArmsLig.ix[ r.ligarmidx ]

        lout=[ r['chrom'], armExt.start, armLig.end+1, desc, score, r['strand_fill'], r['start'], r['end']+1 ]

        nRepeat = 1 if 'picked' not in r else int(r['picked'])

        for _ in xrange(nRepeat):
            bedOut.write( tabnl(lout) )


    