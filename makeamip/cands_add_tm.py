#! /usr/bin/env python

import sys

import os
import os.path
import argparse
from collections import defaultdict, OrderedDict
import copy

import numpy as np
import numpy.random as rand
import re

import pandas as pd

# import primer3.bindings as p3b
# import primer3.wrappers as p3w

import primer3.thermoanalysis as p3ta

import pysam
# import pybedtools as pbt

from makeamip.capcommon import *

from makeamip.capcommon import *

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--genome', dest='genome',
                      help='path to indexed genome file')
    
    opts.add_argument('--genomeDict', dest='genomeDict',
                      help='path to genome contig file dictionary, eg from Picard CreateSequenceDictionary')
    
    opts.add_argument('--inStore', dest='inStore')
    
    opts.add_argument('--table', action='append', dest='table')
    
    o = opts.parse_args()
     
    sfGenome = pysam.Fastafile( o.genome )
     
    mCtgLen = parseContigLengths( o.genomeDict )

    candsStore = pd.HDFStore( o.inStore, 'a' )
    
    thermo = p3ta.ThermoAnalysis( mv_conc=50., dv_conc=0.2 )
    
    chromLast,seqChromLast=None,None
    for tblName in o.table:
        assert tblName in candsStore, 'cannot find table %s in %s'%(tblName, o.inStore)
        tblCands = candsStore[ tblName ]        
        
        lArmTm = np.zeros_like( tblCands.arm_tm )
        
        for iarm in tblCands.index:
            
            if iarm%10000==0:
                print '%s:%d..'%(tblName,iarm)
            
            if chromLast != tblCands.chrom[iarm]:
                chromLast = tblCands.chrom[iarm]
                seqChromLast = sfGenome.fetch( str( chromLast ) ).upper()
                
            seqArm = seqChromLast[ tblCands.start[iarm]:tblCands.end[iarm]+1 ]
            lArmTm[iarm] = thermo.calcTm( seqArm )
        
        tblCands['arm_tm'] = lArmTm[:]

        copyMetadataAndReplaceTable( candsStore, tblName, tblCands, prefix='mippipe__' )
# 