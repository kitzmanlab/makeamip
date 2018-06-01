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

from subprocess import PIPE,Popen

from makeamip.capcommon import *

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--genome', dest='genome',
                      help='path to indexed genome file')
    
    opts.add_argument('--genomeDict', dest='genomeDict',
                      help='path to genome contig file dictionary, eg from Picard CreateSequenceDictionary')
    
    opts.add_argument('--inStore', dest='inStore')
    
    opts.add_argument('--table', action='append', dest='table')
    
    opts.add_argument('-k','--kmerSize', type=int, dest='kmerSize')
    
    opts.add_argument('-j','--jellyfishFile', dest='jellyfishFile')
    
#     echo -e "CCCACGATTACACAC" | jellyfish query -i hs37d5_k15_canon.jf
    
    o = opts.parse_args()
     
    sfGenome = pysam.Fastafile( o.genome )
     
    mCtgLen = parseContigLengths( o.genomeDict )

    candsStore = pd.HDFStore( o.inStore, 'a' )
    
    procjfish = Popen( ['jellyfish', 'query', '-i', o.jellyfishFile ], stdin=PIPE, stdout=PIPE )
    
    chromLast,seqChromLast=None,None
    for tblName in o.table:
        assert tblName in candsStore, 'cannot find table %s in %s'%(tblName, o.inStore)
        tblCands = candsStore[ tblName ]        
        
        lArmKmerCounts = np.zeros_like( tblCands.arm_tm )
        
        for iarm in tblCands.index:
            
            if iarm%1000==0:
                print '%s:%d..'%(tblName,iarm)
                        
            if chromLast != tblCands.chrom[iarm]:
                chromLast = tblCands.chrom[iarm]
                seqChromLast = sfGenome.fetch( str( chromLast ) ).upper()
                
            seqArm = seqChromLast[ tblCands.start[iarm]:tblCands.end[iarm]+1 ]
            Larm = len(seqArm)
            
                       
            lkmer = [ seqArm[i:i+o.kmerSize] for i in xrange( Larm-o.kmerSize+1 ) ]
            
            procjfish.stdin.write( '%s\n'%'\n'.join(lkmer) )
            procjfish.stdin.flush()
                                   
            lfreqs = np.array( [ int(procjfish.stdout.readline().rstrip()) for _ in xrange(len(lkmer)) ], dtype=np.int32 )
            
            lArmKmerCounts[ iarm ] = lfreqs.mean()           
            

        tblCands['arm_mean_kmer_freq'] = lArmKmerCounts[:]
    
        copyMetadataAndReplaceTable( candsStore, tblName, tblCands, prefix='mippipe__' )
