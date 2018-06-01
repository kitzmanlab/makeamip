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


def hard_filter_2_ds(mtbl):
    lpassAllreqdExt = np.c_[ 
              mtbl[ 'ext' ]['arm_gc']>=10.,
              mtbl[ 'ext' ]['arm_gc']<=90.,
              mtbl[ 'ligminus' ]['arm_gc']>=10.,
              mtbl[ 'ligminus' ]['arm_gc']<=90.,
              mtbl[ 'ext' ]['arm_mean_kmer_freq']<=125,
              mtbl[ 'ligminus' ]['arm_mean_kmer_freq']<=125]

    lpassAllreqdLig = np.c_[ 
              mtbl[ 'lig' ]['arm_gc']>=10.,
              mtbl[ 'lig' ]['arm_gc']<=90.,
              mtbl[ 'extminus' ]['arm_gc']>=10.,
              mtbl[ 'extminus' ]['arm_gc']<=90.,
              mtbl[ 'lig' ]['arm_mean_kmer_freq']<=125,
              mtbl[ 'extminus' ]['arm_mean_kmer_freq']<=125 ]

    lpassExt = np.logical_and.reduce(lpassAllreqdExt,1)
    lpassLig = np.logical_and.reduce(lpassAllreqdLig,1)

    return lpassExt,lpassLig


def hard_filter_uniq_ds_5exactclose( mtbl ):

    lpassAllreqdExt = np.c_[ 
              (mtbl[ 'ext_hits' ]['num_exact_hits']+mtbl[ 'ext_hits' ]['num_close_hits']) < 20,
              (mtbl[ 'ligminus_hits' ]['num_exact_hits']+mtbl[ 'ligminus_hits' ]['num_close_hits']) < 20 
              ]

    lpassAllreqdLig = np.c_[ 
              (mtbl[ 'lig_hits' ]['num_exact_hits']+mtbl[ 'lig_hits' ]['num_close_hits']) < 20,
              (mtbl[ 'extminus_hits' ]['num_exact_hits']+mtbl[ 'extminus_hits' ]['num_close_hits']) < 20 
              ]

    lpassExt = np.logical_and.reduce(lpassAllreqdExt,1)
    lpassLig = np.logical_and.reduce(lpassAllreqdLig,1)

    return lpassExt,lpassLig



if __name__ == '__main__':
    
    def intRange(s):
        try:
            lower,upper=map(int, s.split(','))
            return lower,upper
        except:
            raise argparse.ArgumentTypeError('range must be x,y')
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--inStore', dest='inStore')

    opts.add_argument('--outStore', dest='outStore')

    opts.add_argument('--baseName', default='', dest='baseName')
    opts.add_argument('--hitTblSuffix', dest='hitTblSuffix') 
    opts.add_argument('--snpTblSuffix', dest='snpTblSuffix')

    opts.add_argument('--dsMIPs', action='store_true', default=False, dest='dsMIPs')
    
    # function( {tblname:tbl} ) --> array of pass/fail for extension arms, same for ligation arms
    opts.add_argument('--fxnPass', 
                      default="lambda mt:np.ones( (mt[mt.keys()[0]].shape[0],), dtype='bool' )", 
                      dest='fxnPass')

    o = opts.parse_args()
     
    fxnPass=eval(o.fxnPass)

    # read in store
    stIn = pd.HDFStore( o.inStore, 'r' )
    # output store
    stOut = pd.HDFStore( o.outStore, 'w', complib='zlib', complevel=9 )
    
    mtbls = {}

    tblInExtArm = stIn['%sextArm'%o.baseName]
    tblInLigArm = stIn['%sligArm'%o.baseName]

    tblnameExtArmHits='%s_%s'%( '%sextArm'%o.baseName, o.hitTblSuffix )
    tblnameLigArmHits='%s_%s'%( '%sligArm'%o.baseName, o.hitTblSuffix )

    tblnameExtArmSnps='%s_%s'%( '%sextArm'%o.baseName, o.snpTblSuffix )
    tblnameLigArmSnps='%s_%s'%( '%sligArm'%o.baseName, o.snpTblSuffix )

    hitsFound=True
    if tblnameExtArmHits in stIn:
        tblInExtArmHits = stIn[ tblnameExtArmHits ] 
        tblInLigArmHits = stIn[ tblnameLigArmHits ]
    else:
        print 'hits table records not found, skipping'
        hitsFound = False
        tblInExtArmHits = None 
        tblInLigArmHits = None

    snpsFound=True
    if tblnameExtArmSnps in stIn:
        tblInExtArmSnps = stIn[ tblnameExtArmSnps ] 
        tblInLigArmSnps = stIn[ tblnameLigArmSnps ]
    else:
        print 'snps table records not found, skipping'
        snpsFound=False
        tblInExtArmSnps = None 
        tblInLigArmSnps = None 

    tblnameTargets = stIn.get_storer( '%sextArm'%o.baseName  ).attrs['mippipe__target_table']
    tblTargets = stIn[ tblnameTargets ]

    mtbls['ext']=tblInExtArm
    mtbls['lig']=tblInLigArm
    mtbls['ext_hits']=tblInExtArmHits
    mtbls['lig_hits']=tblInLigArmHits
    mtbls['ext_snps']=tblInExtArmSnps
    mtbls['lig_snps']=tblInLigArmSnps

    if o.dsMIPs:
        tblInExtArmMinus = stIn[ '%sextArmMinus'%o.baseName ]
        tblInLigArmMinus = stIn[ '%sligArmMinus'%o.baseName ]
    
        if hitsFound:
            tblInExtArmMinusHits = stIn[ '%s_%s'%( '%sextArmMinus'%o.baseName, o.hitTblSuffix ) ]
            tblInLigArmMinusHits = stIn[ '%s_%s'%( '%sligArmMinus'%o.baseName, o.hitTblSuffix ) ]
        else:
            tblInExtArmMinusHits=None
            tblInLigArmMinusHits=None

        if snpsFound:
            tblInExtArmMinusSnps = stIn[ '%s_%s'%( '%sextArmMinus'%o.baseName, o.snpTblSuffix ) ]
            tblInLigArmMinusSnps = stIn[ '%s_%s'%( '%sligArmMinus'%o.baseName, o.snpTblSuffix ) ]
        else:
            tblInExtArmMinusSnps=None
            tblInLigArmMinusSnps=None

        mtbls['extminus']=tblInExtArmMinus
        mtbls['ligminus']=tblInLigArmMinus

        mtbls['extminus_hits']=tblInExtArmMinusHits
        mtbls['ligminus_hits']=tblInLigArmMinusHits
        mtbls['extminus_snps']=tblInExtArmMinusSnps
        mtbls['ligminus_snps']=tblInLigArmMinusSnps

    print 'in ext: %d rows'%tblInExtArm.shape[0]
    print 'in lig: %d rows'%tblInLigArm.shape[0]
    
    lidxPassExt,lidxPassLig = fxnPass( mtbls ) 

    print 'out ext : %d rows'%( lidxPassExt.sum() )
    print 'out lig : %d rows'%( lidxPassLig.sum() )

    tblOutExtArm = tblInExtArm[ lidxPassExt ]
    tblOutLigArm = tblInLigArm[ lidxPassLig ]

    if hitsFound:
        tblOutExtArmHits = tblInExtArmHits[ lidxPassExt ]
        tblOutLigArmHits = tblInLigArmHits[ lidxPassLig ]

    if snpsFound:
        tblOutExtArmSnps = tblInExtArmSnps[ lidxPassExt ]
        tblOutLigArmSnps = tblInLigArmSnps[ lidxPassLig ]

    stOut['%sextArm'%o.baseName] = tblOutExtArm
    stOut['%sligArm'%o.baseName] = tblOutLigArm

    if hitsFound:
        stOut[ '%s_%s'%( '%sextArm'%o.baseName, o.hitTblSuffix ) ] = tblOutExtArmHits
        stOut[ '%s_%s'%( '%sligArm'%o.baseName, o.hitTblSuffix ) ] = tblOutLigArmHits

    if snpsFound:
        stOut[ '%s_%s'%( '%sextArm'%o.baseName, o.snpTblSuffix ) ] = tblOutExtArmSnps
        stOut[ '%s_%s'%( '%sligArm'%o.baseName, o.snpTblSuffix ) ] = tblOutLigArmSnps

    if o.dsMIPs:
        tblOutExtArmMinus = tblInExtArmMinus[ lidxPassLig ]
        tblOutLigArmMinus = tblInLigArmMinus[ lidxPassExt ]

        if hitsFound:
            tblOutExtArmMinusHits = tblInExtArmMinusHits[ lidxPassLig ]
            tblOutLigArmMinusHits = tblInLigArmMinusHits[ lidxPassExt ]

        if snpsFound:
            tblOutExtArmMinusSnps = tblInExtArmMinusSnps[ lidxPassLig ]
            tblOutLigArmMinusSnps = tblInLigArmMinusSnps[ lidxPassExt ]

        stOut[ '%sextArmMinus'%o.baseName ] = tblOutExtArmMinus 
        stOut[ '%sligArmMinus'%o.baseName ] = tblOutLigArmMinus 
    
        if hitsFound:
            stOut[ '%s_%s'%( '%sextArmMinus'%o.baseName, o.hitTblSuffix ) ] = tblOutExtArmMinusHits
            stOut[ '%s_%s'%( '%sligArmMinus'%o.baseName, o.hitTblSuffix ) ] = tblOutLigArmMinusHits 

        if snpsFound:
            stOut[ '%s_%s'%( '%sextArmMinus'%o.baseName, o.snpTblSuffix ) ] = tblOutExtArmMinusSnps
            stOut[ '%s_%s'%( '%sligArmMinus'%o.baseName, o.snpTblSuffix ) ] = tblOutLigArmMinusSnps 


    tblnameTargets = stIn.get_storer( '%sextArm'%o.baseName  ).attrs['mippipe__target_table']
    tblTargets = stIn[ tblnameTargets ]

    stOut[tblnameTargets]=tblTargets

    stOut.get_storer( '%sextArm'%o.baseName  ).attrs['mippipe__target_table'] = tblnameTargets
    stOut.get_storer( '%sligArm'%o.baseName  ).attrs['mippipe__target_table'] = tblnameTargets

    if hitsFound:
        stOut.get_storer( '%s_%s'%('%sextArm'%o.baseName,o.hitTblSuffix) ).attrs['mippipe__source_table_name'] = '%sextArm'%o.baseName
        stOut.get_storer( '%s_%s'%('%sligArm'%o.baseName,o.hitTblSuffix) ).attrs['mippipe__source_table_name'] = '%sextArm'%o.baseName

    if snpsFound:
        stOut.get_storer( '%s_%s'%('%sextArm'%o.baseName,o.snpTblSuffix) ).attrs['mippipe__source_table_name'] = '%sextArm'%o.baseName
        stOut.get_storer( '%s_%s'%('%sligArm'%o.baseName,o.snpTblSuffix) ).attrs['mippipe__source_table_name'] = '%sextArm'%o.baseName

    if o.dsMIPs:
        stOut.get_storer( '%sextArmMinus'%o.baseName  ).attrs['mippipe__target_table'] = tblnameTargets
        stOut.get_storer( '%sligArmMinus'%o.baseName  ).attrs['mippipe__target_table'] = tblnameTargets

        if hitsFound:
            stOut.get_storer( '%s_%s'%('%sextArmMinus'%o.baseName,o.hitTblSuffix) ).attrs['mippipe__source_table_name'] = '%sextArmMinus'%o.baseName
            stOut.get_storer( '%s_%s'%('%sligArmMinus'%o.baseName,o.hitTblSuffix) ).attrs['mippipe__source_table_name'] = '%sextArmMinus'%o.baseName

        if snpsFound:
            stOut.get_storer( '%s_%s'%('%sextArmMinus'%o.baseName,o.snpTblSuffix) ).attrs['mippipe__source_table_name'] = '%sextArmMinus'%o.baseName
            stOut.get_storer( '%s_%s'%('%sligArmMinus'%o.baseName,o.snpTblSuffix) ).attrs['mippipe__source_table_name'] = '%sextArmMinus'%o.baseName

        stOut.get_storer( '%sextArm'%o.baseName ).attrs['mippipe__dsmip_paired_table'] = '%sligArmMinus'%o.baseName
        stOut.get_storer( '%sligArm'%o.baseName ).attrs['mippipe__dsmip_paired_table'] = '%sextArmMinus'%o.baseName

        stOut.get_storer( '%sextArmMinus'%o.baseName ).attrs['mippipe__dsmip_paired_table'] = '%sligArm'%o.baseName
        stOut.get_storer( '%sligArmMinus'%o.baseName ).attrs['mippipe__dsmip_paired_table'] = '%sextArm'%o.baseName





    # copy hit numbers tbl
    # copy targets tbl
    #
    # copy attrs


