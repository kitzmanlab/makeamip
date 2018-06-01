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

#import pysam
import vcf
# import pybedtools as pbt

from makeamip.capcommon import *

if __name__ == '__main__':
        
    opts = argparse.ArgumentParser()
   
    opts.add_argument('--inStore', dest='inStore')

    opts.add_argument('--inVcfList', dest='inVcfList')
    
    opts.add_argument('--table', action='append', dest='table')
        
    opts.add_argument('--newtblSuffix', default=None, dest='newtblSuffix')
        
    o = opts.parse_args()
        
    assert o.newtblSuffix is not None
         
    candsStore = pd.HDFStore( o.inStore, 'a' )
        
    def newSnpTbl( candsTbl ):
        tbl = OrderedDict()

        for col in ['most_common_var_af','af_most_common_var_5p_1base','af_most_common_var_3p_1base','af_most_common_var_5p_5base','af_most_common_var_3p_5base']:
            tbl[col]=np.zeros( (candsTbl.shape[0],), dtype=np.float32 )

        for col in ['most_common_var_dist_from_5p','most_common_var_dist_from_3p']:
            tbl[col]=-999*np.ones( (candsTbl.shape[0],), dtype=np.int32 )

        tbl = pd.DataFrame(tbl)
        tbl['aidx']=candsTbl.index[:]
        tbl=tbl.set_index('aidx',drop=True)
    
        return tbl

    lvcf = [ vcf.Reader( filename=fn.rstrip() ) for fn in open(o.inVcfList ,'r') ]

    for tblname in o.table:

        tblArm = candsStore[tblname]
        tblnameTargets = candsStore.get_storer( tblname  ).attrs['mippipe__target_table']
        tblTargets = candsStore[ tblnameTargets ]

        tblSnps = newSnpTbl( tblArm )

        for itgt in tblTargets.index:

            lvl = []
            for vcfFile in lvcf:
                try:
                    lvl += [ l for l in vcfFile.fetch( 
                        str(tblTargets.chrom[itgt]), 
                        tblTargets.start[itgt], tblTargets.end[itgt]+1 ) ]
                except ValueError as e:
                    print 'warning: trying to iterate target %s:%d-%d failed'%(str(tblTargets.chrom[itgt]), 
                        tblTargets.start[itgt], tblTargets.end[itgt]+1 )

            tblVarCurTgt = pd.DataFrame({
                'start':[ vl.affected_start for vl in lvl ],
                'end':[ vl.affected_end-1  for vl in lvl ],
                'af':[ max(vl.INFO.get('AF',0.))  for vl in lvl ] } )

            tblArmsCurTgt = tblArm[ tblArm.targetidx==itgt ]

            print 'table %s target %d/%d : %d arms, %d variants'%( tblname, itgt, tblTargets.shape[0], tblArmsCurTgt.shape[0], tblVarCurTgt.shape[0] )

            if tblArmsCurTgt.shape[0] == 0:
                print 'warning: skipping target %d %s'%( itgt, repr(tblTargets.ix[itgt])) 
                continue


            # generally will be many more arms than variants so let's loop over the variants
            # and vectorize over the arms.
            for ivar,rvar in tblVarCurTgt.iterrows():

                # overlap in arm at all
                arArmsOverlapCurVar = rgOverlap_array_single( 
                    np.c_[ tblArmsCurTgt.start, tblArmsCurTgt.end ],
                    ( rvar.start, rvar.end ))

                arUpdateMostCommon = tblSnps.ix[ tblArmsCurTgt.index[arArmsOverlapCurVar], 'most_common_var_af' ] < rvar.af
                arUpdateMostCommon = arUpdateMostCommon[arUpdateMostCommon].index

                if arUpdateMostCommon.shape[0]>0:
                    tblSnps.ix[ arUpdateMostCommon, 'most_common_var_af' ] = rvar.af
                    tblSnps.ix[ arUpdateMostCommon, 'most_common_var_dist_from_5p' ] = \
                        np.clip( tblArmsCurTgt.start[arArmsOverlapCurVar] - rvar.start, 0, int(1e9) )
                    tblSnps.ix[ arUpdateMostCommon, 'most_common_var_dist_from_3p' ] = \
                        np.clip( tblArmsCurTgt.end[arArmsOverlapCurVar] - rvar.end, 0, int(1e9) )

                # overlap in arm - 5p most base
                arArmsOverlapCurVar = rgOverlap_array_single( 
                    np.c_[ tblArmsCurTgt.start, tblArmsCurTgt.start ],
                    ( rvar.start, rvar.end ))

                arUpdateMostCommon = tblSnps.ix[ tblArmsCurTgt.index[arArmsOverlapCurVar], 'af_most_common_var_5p_1base' ] < rvar.af
                arUpdateMostCommon = arUpdateMostCommon[arUpdateMostCommon].index
                if arUpdateMostCommon.shape[0]>0:
                    tblSnps.ix[ arUpdateMostCommon, 'af_most_common_var_5p_1base' ] = rvar.af

                # overlap in arm - 5p most 5 bases
                arArmsOverlapCurVar = rgOverlap_array_single( 
                    np.c_[ tblArmsCurTgt.start, tblArmsCurTgt.start+4 ],
                    ( rvar.start, rvar.end ))

                arUpdateMostCommon = tblSnps.ix[ tblArmsCurTgt.index[arArmsOverlapCurVar], 'af_most_common_var_5p_5base' ] < rvar.af
                arUpdateMostCommon = arUpdateMostCommon[arUpdateMostCommon].index
                if arUpdateMostCommon.shape[0]>0:
                    tblSnps.ix[ arUpdateMostCommon, 'af_most_common_var_5p_5base' ] = rvar.af

                # overlap in arm - 3p most base
                arArmsOverlapCurVar = rgOverlap_array_single( 
                    np.c_[ tblArmsCurTgt.end, tblArmsCurTgt.end ],
                    ( rvar.start, rvar.end ))

                arUpdateMostCommon = tblSnps.ix[ tblArmsCurTgt.index[arArmsOverlapCurVar], 'af_most_common_var_3p_1base' ] < rvar.af
                arUpdateMostCommon = arUpdateMostCommon[arUpdateMostCommon].index
                if arUpdateMostCommon.shape[0]>0:
                    tblSnps.ix[ arUpdateMostCommon, 'af_most_common_var_3p_1base' ] = rvar.af

                # overlap in arm - 3p most 5 bases
                arArmsOverlapCurVar = rgOverlap_array_single( 
                    np.c_[ tblArmsCurTgt.end-4, tblArmsCurTgt.end ],
                    ( rvar.start, rvar.end ))

                arUpdateMostCommon = tblSnps.ix[ tblArmsCurTgt.index[arArmsOverlapCurVar], 'af_most_common_var_3p_5base' ] < rvar.af
                arUpdateMostCommon = arUpdateMostCommon[arUpdateMostCommon].index

                if arUpdateMostCommon.shape[0]>0:
                    tblSnps.ix[ arUpdateMostCommon, 'af_most_common_var_3p_5base' ] = rvar.af


        candsStore.put( '%s_%s'%(tblname, o.newtblSuffix), tblSnps )
        candsStore.get_storer( '%s_%s'%(tblname,o.newtblSuffix) ).attrs['source_table_name'] = tblname




        # for col in ['af_most_common_var_5p_1base','af_most_common_var_3p_1base',
        #           'af_most_common_var_5p_5base','af_most_common_var_3p_5base']:
        #     tbl[col]=np.zeros( (candsTbl.shape[0],), dtype=np.float32 )




        #     ctrArm=0
        #     for iarm in tblArmsCurTgt.index:
        #         if ctrArm%1000==0:
        #             print '   %d'%ctrArm
        #         ctrArm+=1

        #         rarm = tblArmsCurTgt.ix[iarm]

        #         # any snps in this MIP arm?
        #         arOverlaps = rgOverlap_array_single( 
        #             np.c_[ tblVarCurTgt.start, tblVarCurTgt.end ],
        #             ( rarm.start, rarm.end ))

        #         nOverlapSnp = arOverlaps.sum()

        #         if nOverlapSnp > 0:
        #             tblVarsCurArm = tblVarCurTgt[ arOverlaps ]

        #             varMostCommon = tblVarsCurArm.ix[ tblVarsCurArm.af.argmax() ]

        #             tblSnps.ix[ iarm, 'most_common_var_af' ] = varMostCommon.af
        #             tblSnps.ix[ iarm, 'most_common_var_dist_from_5p' ] = \
        #                 max( 0, varMostCommon.start - rarm.start )

        #             tblSnps.ix[ iarm, 'most_common_var_dist_from_3p' ] = \
        #                 max( 0, rarm.end - varMostCommon.end ) 

        #             arOverlaps_5p_1base = rgOverlap_array_single( 
        #                 np.c_[ tblVarCurTgt.start, tblVarCurTgt.end ],
        #                 ( rarm.start, rarm.start ))
        #             if arOverlaps_5p_1base.sum()>0:                        
        #                 tblSnps.ix[ iarm, 'af_most_common_var_5p_1base' ] = \
        #                     tblVarCurTgt[ arOverlaps_5p_1base ].af.max()

        #             arOverlaps_3p_1base = rgOverlap_array_single( 
        #                 np.c_[ tblVarCurTgt.start, tblVarCurTgt.end ],
        #                 ( rarm.end, rarm.end ))
        #             if arOverlaps_3p_1base.sum()>0:                        
        #                 tblSnps.ix[ iarm, 'af_most_common_var_3p_1base' ] = \
        #                     tblVarCurTgt[ arOverlaps_3p_1base ].af.max()


        #             arOverlaps_5p_5base = rgOverlap_array_single( 
        #                 np.c_[ tblVarCurTgt.start, tblVarCurTgt.end ],
        #                 ( rarm.start, rarm.start+4 ))
        #             if arOverlaps_5p_5base.sum()>0:                        
        #                 tblSnps.ix[ iarm, 'af_most_common_var_5p_5base' ] = \
        #                     tblVarCurTgt[ arOverlaps_5p_5base ].af.max()

        #             arOverlaps_3p_5base = rgOverlap_array_single( 
        #                 np.c_[ tblVarCurTgt.start, tblVarCurTgt.end ],
        #                 ( rarm.end-4, rarm.end ))
        #             if arOverlaps_3p_5base.sum()>0:                        
        #                 tblSnps.ix[ iarm, 'af_most_common_var_3p_5base' ] = \
        #                     tblVarCurTgt[ arOverlaps_3p_5base ].af.max()


        # candsStore.put( '%s_%s'%(tblname, o.newtblSuffix), tblSnps )
        # candsStore.get_storer( '%s_%s'%(tblname,o.newtblSuffix) ).attrs['source_table_name'] = tblname

