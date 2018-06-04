#! /usr/bin/env python

import sys

import os
import os.path
import argparse
from collections import defaultdict, OrderedDict

from itertools import groupby

import numpy as np
import numpy.random as rand
import re

import pandas as pd

import pybedtools as pbt

from makeamip.capcommon import *

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()
   
    opts.add_argument('--inTargetBed', dest='inTargetBed')
    opts.add_argument('--inMipBed', dest='inMipBed')
    opts.add_argument('--out', default='/dev/stdout', dest='out')
    opts.add_argument('--outByTarg', default='/dev/stdout', dest='outByTarg')
    opts.add_argument('--libname', default='', dest='libname')
        
    o = opts.parse_args()

    # we want to measure: 
    #  (1) average coverage of TARGET bases (ie not flanks)
    #  (2) average DISTINCT coverage of TARGET bases (ie not flanks)
    #  (3) percent coverage which is ON-TARGET

    btTargets = pbt.BedTool( o.inTargetBed )
    btMipWholeprobe = pbt.BedTool( o.inMipBed ).sort()

    btMipsGapfill = pbt.BedTool( [pbt.Interval( i.chrom, int(i.fields[6]), int(i.fields[7]), i.name, i.score ) for i in btMipWholeprobe ] )
    btMipsGapfillUnique = pbt.BedTool( [ btMipsGapfill[j] for j in xrange(len(btMipsGapfill)) if j==0 or str(btMipsGapfill[j-1])!=str(btMipsGapfill[j])  ]  )

    
    tblByTarg = pd.DataFrame( {'chrom':[i.chrom for i in btTargets],
                               'start':[i.start for i in btTargets],
                               'end':[i.end for i in btTargets] } ) 
    # tblByTarg['start']=tblByTarg.start.astype('d')
    # tblByTarg['end']=tblByTarg.start.astype('d')
    tblByTarg['target_length']=tblByTarg.end - tblByTarg.start + 1
    tblByTarg['mean_perbp_coverage_by_probes']=0.
    tblByTarg['mean_perbp_coverage_by_distinct_probes']=0.
    tblByTarg['mean_probe_score']=0.
    tblByTarg['num_probes']=0.
    tblByTarg = tblByTarg.set_index( ['chrom','start','end'] )
    
    # btCvgByTarget = btMipsGapfill.coverage( btTargets, hist=True, all=False ).remove_invalid()
    btCvgByTarget = btTargets.coverage( btMipsGapfill, hist=True, all=False ).remove_invalid()

    for r in btCvgByTarget:
        if r.chrom == 'all' : continue
        tblByTarg.loc[ (r.chrom,int(r.start),int(r.end) ), 'mean_perbp_coverage_by_probes' ] += float(r.fields[3]) * float(r.fields[4]) / float(r.fields[5])

    btCvgByTargetDistinct = btTargets.coverage( btMipsGapfillUnique, hist=True, all=False ).remove_invalid()

    for r in btCvgByTargetDistinct:
        if r.chrom == 'all' : continue
        tblByTarg.loc[ (r.chrom,int(r.start),int(r.end) ), 'mean_perbp_coverage_by_distinct_probes' ] += float(r.fields[3]) * float(r.fields[4]) / float(r.fields[5])

    tblByProbe = pd.DataFrame( {'chrom':[i.chrom for i in btMipsGapfill],
                               'start':[i.start for i in btMipsGapfill],
                               'end':[i.end for i in btMipsGapfill] } ) 
    tblByProbe['gapfill_length']=tblByProbe.end - tblByProbe.start + 1
    tblByProbe['frac_off_target']=0
    tblByProbe = tblByProbe.set_index( ['chrom','start','end'] )

    btOnTargByProbe = btTargets.coverage( btMipsGapfill, hist=True ).remove_invalid()
    for r in btOnTargByProbe:
        if r.chrom == 'all': continue
        if int(r.fields[-4])==0:
            tblByProbe.loc[ (r.chrom,int(r.start),int(r.end)), 'frac_off_target' ] += float( r.fields[-1] )


    for targ,ltxp in groupby( btTargets.intersect(btMipsGapfill,wb=True), lambda txp:txp.fields[:3] ):
        lscores = []
        for txp in ltxp:
            lscores.append( float(txp.fields[7]) ) 

        tblByTarg.loc[ (targ[0],int(targ[1]),int(targ[2])), 'mean_probe_score'] = np.mean(lscores)

    tblByTarg.to_csv(o.outByTarg,sep='\t',index=True)

    mRes=pd.DataFrame(
        OrderedDict([
            ('libname',o.libname),
            ('perbase_average_coverage', (tblByTarg.target_length * tblByTarg.mean_perbp_coverage_by_probes).sum() / float(tblByTarg.target_length.sum())),
            ('perbase_average_coverage_distinct_probe', (tblByTarg.target_length * tblByTarg.mean_perbp_coverage_by_distinct_probes).sum() / float(tblByTarg.target_length.sum())),
            ('frac_gapfills_off_target', (tblByProbe.gapfill_length * tblByProbe.frac_off_target).sum() / float(tblByProbe.gapfill_length.sum())),
            ('total_probes', len(btMipsGapfill) ),
            ('total_probes_unique', len(btMipsGapfillUnique) ),
            ('total_probes_repeated', len(btMipsGapfill)-len(btMipsGapfillUnique) ),

        ]), index=[0])

    mRes.to_csv( o.out, sep='\t', header=True, index=False      )
