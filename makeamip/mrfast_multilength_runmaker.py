#! /usr/bin/env python

import sys

import os
import os.path as op
import argparse
from collections import defaultdict, OrderedDict

import numpy as np
import numpy.random as rand
import re

import pandas as pd

import tempfile

import pysam

from makeamip.capcommon import *

if __name__ == '__main__':
        
    opts = argparse.ArgumentParser()
    
    opts.add_argument('--workDir',default='./',dest='workDir' )
    
    opts.add_argument('--outRunScript', dest='outRunScript' )
    opts.add_argument('--outJoinScript',dest='outJoinScript' )
    opts.add_argument('--outCleanupScript',dest='outCleanupScript' )

    opts.add_argument('--inFastq', dest='inFastq')
    opts.add_argument('--table', action='append', dest='table')
    
    opts.add_argument('--mrfastPath', default='mrfast', dest='mrfastPath')
    opts.add_argument('--mrfastOptions', dest='mrfastOptions')
    opts.add_argument('--outBam', dest='outBam')
       
    o = opts.parse_args()


    fnTempBase = tempfile.NamedTemporaryFile( mode='w', dir=o.workDir ).name
        
    filInFastq = open(o.inFastq,'r')
    
    mLenOut = {}
    
    minReadLen,maxReadLen=99999,0
    l1=filInFastq.readline()
    while len(l1)>0:
        l2=filInFastq.readline().rstrip()
        rl=len(l2)
        maxReadLen=max(maxReadLen,rl)
        minReadLen=min(minReadLen,rl)
        if rl not in mLenOut:
            mLenOut[rl] = open( '%s.len%d.fq'%(fnTempBase,rl), 'w' )
        l3=filInFastq.readline()
        l4=filInFastq.readline()

        mLenOut[rl].write( l1 )
        mLenOut[rl].write( '%s\n'%l2 )
        mLenOut[rl].write( l3 )
        mLenOut[rl].write( l4 )

        l1=filInFastq.readline()
                
    loutRunScript=[]
    loutJoinScript=[]
    loutCleanupScript=[]
    
    for rl in xrange(minReadLen, maxReadLen+1):
        loutRunScript += ['%s %s --seq %s.len%d.fq -u /dev/null -o /dev/stdout | python %s/suppress_bad_mrfast_alignments.py | samtools view -b -S - > %s.len%d.bam'%(  
            o.mrfastPath, o.mrfastOptions, fnTempBase, rl, 
            op.dirname( sys.argv[0] ),
            fnTempBase, rl )]
        
    loutJoinScript = ['samtools cat -o %s %s'%( 
                            o.outBam, 
                            ' '.join( [ '%s.len%d.bam'%(fnTempBase,rl) for rl in xrange(minReadLen,maxReadLen+1)] )) ]
    
    loutCleanupScript = ['rm %s'%( ' '.join('%s.len%d.fq %s.len%d.bam'%(fnTempBase,rl,fnTempBase,rl) for rl in xrange(minReadLen,maxReadLen+1)) ) ]
    
    filOutRunScript=open(o.outRunScript,'w')
    filOutRunScript.write('\n'.join(loutRunScript))

    filOutJoinScript=open(o.outJoinScript,'w')
    filOutJoinScript.write('\n'.join(loutJoinScript))

    filOutCleanupScript=open(o.outCleanupScript,'w')
    filOutCleanupScript.write('\n'.join(loutCleanupScript))
    
    
                                                                       
        
