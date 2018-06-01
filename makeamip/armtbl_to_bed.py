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
    opts.add_argument('--hitTable', default=None, dest='hitTable')

    opts.add_argument('--fxnDesc',
        default='lambda c,h:"%s_%s_%d_%d_%s"%(c.name,c["chrom"],c["start"],c["end"]," +-"[c["strand"]])',
        dest='fxnDesc')

    opts.add_argument('--fxnScore',
        default='lambda c,h:1.',
        dest='fxnScore')
    
    opts.add_argument('--bedOut', dest='bedOut')
        
    o = opts.parse_args()

    fxnDesc = eval(o.fxnDesc)
    fxnScore = eval(o.fxnScore)
         
    candsStore = pd.HDFStore( o.inStore, 'r' )
        
    bedOut = open(o.bedOut,'w')

    assert o.table in candsStore, 'cannot find table %s in %s'%(o.table, o.inStore)
    tblCands = candsStore[ o.table ]       
    
    if o.hitTable is not None:
        tblHits = candsStore[ o.hitTable ]
    else:
        tblHits = None

    for i,r in tblCands.iterrows():

        if tblHits is not None:
            rhits = tblHits.ix[i]
            desc = fxnDesc( r, rhits )
            score = fxnDesc( r, rhits )
        else:
            desc = fxnDesc(r,None)
            score = fxnScore(r,None)

        lout=[ r['chrom'], r['start'], r['end']+1, desc, score, r['strand'] ]

        bedOut.write( tabnl(lout) )


    