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
        
    o = opts.parse_args()

    pairStore = pd.HDFStore( o.inStore, 'r' )
    tblPairs = pairStore[ o.table ]       

    for i,r in tblPairs.iterrows():
        nRepeat = 1 if 'picked' not in r else int(r['picked'])
        for _ in xrange(nRepeat):
            lout=[i, r['mip_seq_5to3']]
            print '\t'.join( [str(x) for x in lout] )