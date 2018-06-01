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
    opts.add_argument('--armAlignHitsSuffix', dest='armAlignHitsSuffix')

    opts.add_argument('--mipTableOut', default=None, dest='mipTableOut')

    opts.add_argument('--gcFlankingBp', action='append', type=int, dest='gcFlankingBp')

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

    tblArmExtHits = armStore[ '%s_%s'%(tblnameExtArm,o.armAlignHitsSuffix) ]
    tblArmLigHits = armStore[ '%s_%s'%(tblnameLigArm,o.armAlignHitsSuffix) ]

    hdr = ['mip_index', 'ext_index', 'lig_index', 'chrom','strand_fill','start','end','score',
            'gapfill_len',
            'gapfill_gc' ] +\
           ['flanking_gc_%dbp'%i for i in o.gcFlankingBp]+\
           ['ext_probe_start','ext_probe_end','lig_probe_start','lig_probe_end',
            'lig_arm_len','ext_arm_len',

            'ext_nmerfreq','lig_nmerfreq',
            'ext_tm','lig_tm',
            'ext_gc','lig_gc',
            'ext_mrfAndBwa_exact','lig_mrfAndBwa_exact',
            'ext_mrfAndBwa_close','lig_mrfAndBwa_close',
            'ext_mrfAndBwa_all','lig_mrfAndBwa_all' ]

    tblOut = OrderedDict( [ (col,[]) for col in hdr ] )

    for i,r in tblPairs.iterrows():
        tblOut['mip_index'].append( int(i) )
        tblOut['chrom'].append( r.chrom )

        if r.chrom not in mChromSeq:
            print 'fetching %s'%r.chrom
            mChromSeq[r.chrom] = sfRef.fetch( str( r.chrom ) )

        tblOut['strand_fill'].append( r.strand_fill )
        tblOut['start'].append( r.start )
        tblOut['end'].append( r.end )

        tblOut['gapfill_len'].append( r.end-r.start+1 )
        tblOut['gapfill_gc'].append( fracGc( mChromSeq[r.chrom], (r.start,r.end) ) )

        for flank in o.gcFlankingBp:
            cozCenter = int(0.5*(r.start+r.end))
            corzFlank = ( max(cozCenter - flank,0), min(cozCenter + flank,len( mChromSeq[r.chrom] )-1) )
            tblOut['flanking_gc_%dbp'%flank].append( fracGc( mChromSeq[r.chrom], corzFlank ))

        tblOut['score'].append( r.score )

        ea=tblArmExt.ix[ r.extarmidx ]
        la=tblArmLig.ix[ r.ligarmidx ]

        tblOut['ext_index'].append( int(r.extarmidx) )
        tblOut['lig_index'].append( int(r.ligarmidx) )

        eah=tblArmExtHits.ix[ r.extarmidx ]
        lah=tblArmLigHits.ix[ r.ligarmidx ]

        tblOut['ext_probe_start'].append( ea.start )
        tblOut['ext_probe_end'].append( ea.end )
        tblOut['lig_probe_start'].append( la.start )
        tblOut['lig_probe_end'].append( la.end )

        tblOut['ext_arm_len'].append( abs(ea.end-ea.start)+1 )
        tblOut['lig_arm_len'].append( abs(la.end-la.start)+1 )

        ## 

        tblOut['ext_nmerfreq'].append( ea.arm_mean_kmer_freq )
        tblOut['lig_nmerfreq'].append( la.arm_mean_kmer_freq )

        tblOut['ext_tm'].append( ea.arm_tm )
        tblOut['lig_tm'].append( la.arm_tm )

        tblOut['ext_mrfAndBwa_exact'].append( eah.num_exact_hits )
        tblOut['ext_mrfAndBwa_close'].append( eah.num_close_hits )
        tblOut['ext_mrfAndBwa_all'].append( eah.num_all_hits )

        tblOut['lig_mrfAndBwa_exact'].append( lah.num_exact_hits )
        tblOut['lig_mrfAndBwa_close'].append( lah.num_close_hits )
        tblOut['lig_mrfAndBwa_all'].append( lah.num_all_hits )

        tblOut['ext_gc'].append( ea.arm_gc )
        tblOut['lig_gc'].append( la.arm_gc )



        #


    tblOut = pd.DataFrame( tblOut )
    tblOut.to_csv( o.mipTableOut, sep='\t', index=False )