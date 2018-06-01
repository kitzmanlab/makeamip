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

from makeamip.capcommon import *

if __name__ == '__main__':
        
    opts = argparse.ArgumentParser()

    opts.add_argument('--genome', dest='genome',
                      help='path to indexed genome file')
    
    opts.add_argument('--genomeDict', dest='genomeDict',
                      help='path to genome contig file dictionary, eg from Picard CreateSequenceDictionary')
    
    opts.add_argument('--inStore', dest='inStore')
    
    opts.add_argument('--table', action='append', dest='table')
    
    opts.add_argument('--format', default='fastq', choices=['fastq','fasta'], dest='format')
    
    opts.add_argument('-o','--out', dest='fnOut') 
       
    o = opts.parse_args()
     
    sfGenome = pysam.Fastafile( o.genome )
     
    mCtgLen = parseContigLengths( o.genomeDict )

    candsStore = pd.HDFStore( o.inStore, 'r' )
        
    filOut = open(o.fnOut,'w')
        
    if o.format == 'fasta':
        armFormatter = lambda armseq,tblname,ix:'>%s_%d\n%s\n'%( tblname, ix, armseq )
    elif o.format == 'fastq':
        armFormatter = lambda armseq,tblname,ix:'@%s_%d\n%s\n+\n%s\n'%( tblname, ix, armseq, 'I'*len(armseq) )

    chromLast,seqChromLast=None,None
    for tblName in o.table:
        assert tblName in candsStore, 'cannot find table %s in %s'%(tblName, o.inStore)
        tblCands = candsStore[ tblName ]        
                
        for iarm in tblCands.index:
            
            if iarm%1000==0:
                print '%s:%d..'%(tblName,iarm)
                        
            if chromLast != tblCands.chrom[iarm]:
                chromLast = tblCands.chrom[iarm]
                seqChromLast = sfGenome.fetch( str( chromLast ) ).upper()
                
            seqArm = seqChromLast[ tblCands.start[iarm]:tblCands.end[iarm]+1 ]
             
            filOut.write( armFormatter(seqArm,tblName,iarm) )
                                   

    