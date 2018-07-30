# from __future__ import absolute_import

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
   
    opts.add_argument('--inTbl', dest='inTbl')
    opts.add_argument('--colChrom', dest='colChrom')
    opts.add_argument('--colStart', dest='colStart')
    opts.add_argument('--colEnd', dest='colEnd')

    opts.add_argument('--outTbl',dest='outTbl')
    opts.add_argument('--colOutTarget',dest='colOutTarget')

    opts.add_argument('--inBedTargets', dest='inBedTargets')

    o = opts.parse_args()

    tbl = pd.read_csv(o.inTbl,sep='\t')
    tbl[o.colOutTarget]=''

    btTbl = tbl[ ['chrom','start','end'] ].copy()
    btTbl['name'] = btTbl.index
    btTbl = pbt.BedTool.from_dataframe( btTbl )

    btTargets = pbt.BedTool(o.inBedTargets)
    tblTargets = \
        pd.DataFrame( {'chrom': [iv.chrom for iv in btTargets],
                       'start': [iv.start for iv in btTargets],
                       'end': [iv.end for iv in btTargets],
                       'name': range(len(btTargets)) } )
    tblTargets=tblTargets[ ['chrom','start','end','name'] ]
    btTargets = pbt.BedTool.from_dataframe( tblTargets )

    probesToTgts = btTbl.intersect( btTargets, wa=True, wb=True )
    probesToTgts = pd.DataFrame( [iv.fields for iv in probesToTgts ] )
    probesToTgts = probesToTgts.convert_objects(convert_numeric=True)

    # group by 3 = id of probe.  join 7 (id(s) of targets)
    z=probesToTgts[[3,7]].groupby(3).agg( lambda g:str(','.join( [str(v) for v in g[7]] )))
    
    tbl.ix[ z.index, 'cap_tgt' ] = z[7] 

    tbl.ix[ z.index, o.colOutTarget] = z[7]
    tbl.to_csv( o.outTbl, sep='\t', index=False )
