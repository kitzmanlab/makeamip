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

    opts.add_argument('--pairStore', dest='pairStore')

    opts.add_argument('--pairsBaseName', default='', dest='pairsBaseName')

    opts.add_argument('--colNullToFalse', default=None, dest='colNullToFalse' )

    o = opts.parse_args()

    pairStore = pd.HDFStore( o.pairStore, 'a' )

    tblnamePairsPlus = '%spairsPlus'%o.pairsBaseName
    tblnamePairsMinus = '%spairsMinus'%o.pairsBaseName

    tblPairsPlus = pairStore[ tblnamePairsPlus ]
    tblPairsMinus = pairStore[ tblnamePairsMinus ]

    tblPairsPlus.ix[  pd.isnull(tblPairsPlus[o.colNullToFalse]), o.colNullToFalse ] = False
    tblPairsMinus.ix[  pd.isnull(tblPairsMinus[o.colNullToFalse]), o.colNullToFalse ] = False

# copyMetadataAndReplaceTable( store, tblname, tblold
# copyMetadataAndReplaceTable( candsStore, tblName, tblCands, prefix='mippipe__' )

    copyMetadataAndReplaceTable( pairStore, tblnamePairsPlus, tblPairsPlus )
    copyMetadataAndReplaceTable( pairStore, tblnamePairsMinus, tblPairsMinus )
