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

def perPairMetrics( armStore,
                    pairStore,
                    hitTblSuffix,
                    snpTblSuffix,
                    pairtblName ):
    
    pairTbl = pairStore[ pairtblName ]

    (tblnameExtArm,tblnameLigArm)=pairStore.get_storer( pairtblName ).attrs['mippipe__ext_and_lig_arm_tables']
    
    tblnameExtHits,tblnameLigHits='%s_%s'%(tblnameExtArm,hitTblSufifx),'%s_%s'%(tblnameLigArm,hitTblSufifx)
    tblnameExtSnps,tblnameLigSnps='%s_%s'%(tblnameExtArm,snpTblSuffix),'%s_%s'%(tblnameLigArm,snpTblSuffix)

    tblExt = armStore[tblnameExtArm]
    tblLig = armStore[tblnameLigArm]
    tblExtHits = armStore[tblnameExtHits]
    tblLigHits = armStore[tblnameLigHits]
    tblExtSnps = armStore[tblnameExtSnps]
    tblLigSnps = armStore[tblnameLigSnps]

    return perPairMetricsFromTables( pairTbl,
                                    tblExt,
                                    tblLig,
                                    tblExtHits,
                                    tblLigHits,
                                    tblExtSnps,
                                    tblLigSnps )


def perPairMetricsFromTables( pairTbl,
                              tblExt,
                              tblLig,
                              tblExtHits,
                              tblLigHits,
                              tblExtSnps,
                              tblLigSnps):

    pairTbl['len_gap_fill']  = pairTbl.end-pairTbl.start+1
    pairTbl['mean_arm_kmer_freq'] = 0.5*(tblExt.ix[pairTbl.extarmidx]['arm_mean_kmer_freq'].values+
                                         tblLig.ix[pairTbl.ligarmidx]['arm_mean_kmer_freq'].values )
    pairTbl['max_arm_kmer_freq'] = \
      np.c_[ tblExt.ix[pairTbl.extarmidx]['arm_mean_kmer_freq'].values,
             tblLig.ix[pairTbl.ligarmidx]['arm_mean_kmer_freq'].values ].max(1) 

    pairTbl['mean_arm_gc'] = 0.5*(tblExt.ix[pairTbl.extarmidx]['arm_gc'].values+
                                         tblLig.ix[pairTbl.ligarmidx]['arm_gc'].values )
    pairTbl['max_arm_gc'] = \
      np.c_[ tblExt.ix[pairTbl.extarmidx]['arm_gc'].values,
             tblLig.ix[pairTbl.ligarmidx]['arm_gc'].values ].max(1) 

    pairTbl['mean_arm_tm'] = \
      np.c_[ tblExt.ix[pairTbl.extarmidx]['arm_tm'].values,
             tblLig.ix[pairTbl.ligarmidx]['arm_tm'].values ].mean(1)

    pairTbl['arm_tm_diff'] = \
        tblExt.ix[pairTbl.extarmidx]['arm_tm'].values - \
        tblLig.ix[pairTbl.ligarmidx]['arm_tm'].values 

    pairTbl['mean_arm_exact_hits'] = \
      np.c_[ tblExtHits.ix[pairTbl.extarmidx]['num_exact_hits'].values,
             tblLigHits.ix[pairTbl.ligarmidx]['num_exact_hits'].values ].mean(1)

    pairTbl['mean_arm_close_hits'] = \
      np.c_[ tblExtHits.ix[pairTbl.extarmidx]['num_close_hits'].values,
             tblLigHits.ix[pairTbl.ligarmidx]['num_close_hits'].values ].mean(1)

    pairTbl['mean_arm_all_hits'] = \
      np.c_[ tblExtHits.ix[pairTbl.extarmidx]['num_all_hits'].values,
             tblLigHits.ix[pairTbl.ligarmidx]['num_all_hits'].values ].mean(1)

    pairTbl['ext_arm_len'] = \
      tblExt.ix[pairTbl.extarmidx].end.values-tblExt.ix[pairTbl.extarmidx].start.values+1

    pairTbl['lig_arm_len'] = \
      tblLig.ix[pairTbl.ligarmidx].end.values-tblLig.ix[pairTbl.ligarmidx].start.values+1

    pairTbl['mean_arm_len'] = 0.5*(pairTbl['ext_arm_len']+pairTbl['lig_arm_len'])

    return pairTbl
