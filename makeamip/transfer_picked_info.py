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
transfer picked column from one pairs table to another (can be partially overlapping sets of pairs)
"""
    
if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--pairStoreWithPickedInfo', dest='pairStoreWithPickedInfo')
    opts.add_argument('--pairsBaseName', default='', dest='pairsBaseName')
    opts.add_argument('--pairStoreToWriteTo', dest='pairStoreToWriteTo')

    opts.add_argument('--zeroOut', action='store_true', default=False, dest='zeroOut')
    opts.add_argument('--add', action='store_true', default=False, dest='add')


    # opts.add_argument('--hitTblSufifx', default='', dest='hitTblSufifx')
    # opts.add_argument('--snpTblSuffix', default='', dest='snpTblSuffix')

    o = opts.parse_args()

    pairStoreWithPickedInfo = pd.HDFStore( o.pairStoreWithPickedInfo, 'r' )

    tblnamePairsPlus = '%spairsPlus'%o.pairsBaseName
    tblnamePairsMinus = '%spairsMinus'%o.pairsBaseName

    tblPairsPlus = pairStoreWithPickedInfo[ tblnamePairsPlus ]
    tblPairsMinus = pairStoreWithPickedInfo[ tblnamePairsMinus ]

    pairStoreToWriteTo = pd.HDFStore( o.pairStoreToWriteTo, 'a' )

    tblPairsPlusUpdate = pairStoreToWriteTo[ tblnamePairsPlus ]
    tblPairsMinusUpdate = pairStoreToWriteTo[ tblnamePairsMinus ]

    #

    linclWithInfoPlus = tblPairsPlus.index.isin( tblPairsPlusUpdate.index )
    linclWithInfoMinus = tblPairsMinus.index.isin( tblPairsMinusUpdate.index )

    pickedInfoPlus = tblPairsPlus.ix[ linclWithInfoPlus, 'picked' ]
    pickedInfoMinus = tblPairsPlus.ix[ linclWithInfoPlus, 'picked' ]

    if o.zeroOut:
        tblPairsPlusUpdate[ 'picked' ]=0
        tblPairsMinusUpdate[ 'picked' ]=0

    if o.add:
        tblPairsPlusUpdate.ix[ pickedInfoPlus.index, 'picked' ] += pickedInfoPlus
        tblPairsMinusUpdate.ix[ pickedInfoMinus.index, 'picked' ] += pickedInfoMinus
    else:
        tblPairsPlusUpdate.ix[ pickedInfoPlus.index, 'picked' ] = pickedInfoPlus
        tblPairsMinusUpdate.ix[ pickedInfoMinus.index, 'picked' ] = pickedInfoMinus

    copyMetadataAndReplaceTable( pairStoreToWriteTo, tblnamePairsPlus, tblPairsPlusUpdate )
    copyMetadataAndReplaceTable( pairStoreToWriteTo, tblnamePairsMinus, tblPairsMinusUpdate )
