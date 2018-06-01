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

    opts.add_argument('--inBed', dest='inBed')

    opts.add_argument('--maxbpSingleInterval', default=1e9, type=int, dest='maxbpSingleInterval')

    opts.add_argument('--fxnChunkbedNumToFn', dest='fxnChunkbedNumToFn')
    opts.add_argument('--outKey', dest='outKey')

    opts.add_argument('--targetBpPerChunk', type=int, default=10000, dest='targetBpPerChunk') 

    o = opts.parse_args()
     
    fxnChunkbedNumToFn = eval(o.fxnChunkbedNumToFn)

    keyout=OrderedDict()
    keyout['chunkidx']=[]
    keyout['targetbed']=[]

    tblInBed = pd.read_csv(o.inBed,sep='\t',header=None)
    tblInBedChunked = dict( [(c,[]) for c in tblInBed ] )

    for i,r in tblInBed.iterrows():
        if r[2]-r[1] > o.maxbpSingleInterval:
            cozStart=r[1]
            while cozStart < r[2]:
                tblInBedChunked[0].append( r[0] )
                tblInBedChunked[1].append( cozStart )
                tblInBedChunked[2].append( min(r[2],cozStart+o.maxbpSingleInterval) )
                cozStart+=o.maxbpSingleInterval
        else:
            tblInBedChunked[0].append( r[0] )
            tblInBedChunked[1].append( r[1] )
            tblInBedChunked[2].append( r[2] )

    tblInBed = pd.DataFrame(tblInBedChunked)

    tblInBed[1]=tblInBed[1].astype('int32')
    tblInBed[2]=tblInBed[2].astype('int32')

    tblInBed = tblInBed.sort_values( by= [0,1] )

    arCumlens = np.cumsum(tblInBed[2]-tblInBed[1])

    lirngChunks = []

    istart,iend=0,-1
    chunkCurNbp=0
    while iend < arCumlens.shape[0] - 1:
        iend+=1
        chunkCurNbp+=(tblInBed[2][iend] - tblInBed[1][iend] )
        if chunkCurNbp >= o.targetBpPerChunk:
            lirngChunks.append( [ istart,iend ] )
            istart=iend+1
            chunkCurNbp = 0

    if istart < arCumlens.shape[0]:
        lirngChunks.append( (istart,iend) )

    for ichunk in xrange(len(lirngChunks)):
        # NB - pandas ix is inclusive, iloc is not 
        tblOut = tblInBed.iloc[ lirngChunks[ichunk][0]:lirngChunks[ichunk][1]+1 ]
        fnOut = fxnChunkbedNumToFn( ichunk )
        tblOut.to_csv( fnOut, sep='\t', header=None, index=False )
        keyout['chunkidx'].append(ichunk)
        keyout['targetbed'].append(fnOut)

    pd.DataFrame(keyout).to_csv( o.outKey, sep='\t', index=False )