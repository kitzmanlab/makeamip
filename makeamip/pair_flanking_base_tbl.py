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
   
    opts.add_argument('--inStorePairs', dest='inStorePairs')
    opts.add_argument('--tablePairs', dest='tablePairs')
        
    opts.add_argument('--inStoreArms', dest='inStoreArms')

    opts.add_argument('--mipTableOut', default=None, dest='mipTableOut')

    opts.add_argument('--refFasta', default=None, dest='refFasta')

    o = opts.parse_args()

    pairStore = pd.HDFStore( o.inStorePairs, 'r' )
    tblPairs = pairStore[ o.tablePairs ]       

    sfRef = pysam.Fastafile(o.refFasta)
    mChromSeq={}

    (tblnameExtArm,tblnameLigArm) = pairStore.get_storer( o.tablePairs  ).attrs['mippipe__ext_and_lig_arm_tables']

    armStore = pd.HDFStore( o.inStoreArms, 'r' )
    tblArmExt = armStore[ tblnameExtArm ]
    tblArmLig = armStore[ tblnameLigArm ]

    hdr = ['mip_index']
    hdr += ['extside_probe_terminal_%s'%bp for bp in 'ACGT']
    hdr += ['ligside_probe_terminal_%s'%bp for bp in 'ACGT']

    hdr += ['extside_gf_terminal_%s'%bp  for bp in 'ACGT']
    hdr += ['ligside_gf_terminal_%s'%bp  for bp in 'ACGT']

    hdr += ['lig_junction_%s%s'%(bpi,bpj) for bpi in 'ACGT' for bpj in 'ACGT']

    K=1

    tblOut = OrderedDict( [ (col,[]) for col in hdr ] )

    for i,r in tblPairs.iterrows():
        tblOut['mip_index'].append( int(i) )

        for col in hdr :
            if col!='mip_index':
                tblOut[col].append(0)

        if r.chrom not in mChromSeq:
            print 'fetching %s'%r.chrom
            mChromSeq[r.chrom] = sfRef.fetch( str(r.chrom ) )

        seq_gapfill = mChromSeq[r.chrom][r.start:r.end+1]
        tblOut[ 'extside_gf_terminal_%s'%seq_gapfill[:K] ][-1]=1
        tblOut[ 'ligside_gf_terminal_%s'%seq_gapfill[-K:] ][-1]=1

        ea=tblArmExt.ix[ r.extarmidx ]
        la=tblArmLig.ix[ r.ligarmidx ]

        seq_ext = mChromSeq[r.chrom][ea.start:ea.end+1]
        seq_lig = mChromSeq[r.chrom][la.start:la.end+1]
        tblOut[ 'extside_probe_terminal_%s'%seq_ext[-K:] ][-1]=1
        tblOut[ 'ligside_probe_terminal_%s'%seq_lig[:K] ][-1]=1

        tblOut[ 'lig_junction_%s'%(seq_gapfill[-1] + seq_lig[0]) ][-1]=1

        #


    tblOut = pd.DataFrame( tblOut )
    tblOut.to_csv( o.mipTableOut, sep='\t', index=False )
