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
   
    opts.add_argument('--inOligoList', dest='inOligoList')
    opts.add_argument('--outOligoList', dest='outOligoList')

    opts.add_argument('--colIn', dest='colIn')
    opts.add_argument('--colOut', dest='colOut')

    opts.add_argument('--fxnPrefilterSeqIn',default='lambda s:s', dest='fxnPrefilterSeqIn')

    opts.add_argument('--adaptorLeft', dest='adaptorLeft')
    opts.add_argument('--adaptorRight', dest='adaptorRight')
        
    opts.add_argument('--adaptorLeftReplaceWith', dest='adaptorLeftReplaceWith')
    opts.add_argument('--adaptorRightReplaceWith', dest='adaptorRightReplaceWith')

    o = opts.parse_args()

    inOligoList = pd.read_csv(o.inOligoList,sep='\t')

    fxnPrefilterSeqIn = eval(o.fxnPrefilterSeqIn)

    inOligoList[o.colOut] = [ fxnPrefilterSeqIn(s) for s in inOligoList[o.colIn] ]

    for i,r in inOligoList.iterrows():
        ileftStart = inOligoList.ix[i, o.colOut].index( o.adaptorLeft )
        irightStart = inOligoList.ix[i, o.colOut].index( o.adaptorRight )

        seqMid = inOligoList.ix[i, o.colOut][ ileftStart+o.adaptorLeft.__len__() : irightStart ]

        inOligoList.ix[i, o.colOut] = o.adaptorLeftReplaceWith + seqMid + o.adaptorRightReplaceWith

    inOligoList.to_csv(o.outOligoList,sep='\t',index=False)