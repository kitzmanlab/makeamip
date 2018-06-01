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

from makeamip.capcommon import *

if __name__ == '__main__':
       
    opts = argparse.ArgumentParser()
   
    opts.add_argument('--inStore', dest='inStore')
    
    opts.add_argument('--table', action='append', dest='table')
    
    opts.add_argument('--maxedsCloseFrac', type=float, default=0.1, dest='maxedsClose')
    
    opts.add_argument('--newtblSuffix', default=None, dest='newtblSuffix')
    
    opts.add_argument('--inBam', dest='inBam')

    opts.add_argument('--disregardDoubleHitsAtSameLocation', default=False, action='store_true', dest='disregardDoubleHitsAtSameLocation' )
        
    o = opts.parse_args()
        
    assert o.newtblSuffix is not None
         
    candsStore = pd.HDFStore( o.inStore, 'a' )
        
    def newHitsTbl( candsTbl ):
        tbl = OrderedDict()
        for col in ['num_exact_hits','num_close_hits','num_all_hits']: tbl[col]=np.zeros( (candsTbl.shape[0],), dtype=np.int32 )
        return tbl
        
    for tblName in o.table:
        assert tblName in candsStore, 'cannot find table %s in %s'%(tblName, o.inStore)
        tblCands = candsStore[ tblName ]       
        tblHits = newHitsTbl( tblCands )
        tblHits = pd.DataFrame( tblHits )
        tblHits['theindex']=tblCands.index.copy()
        tblHits = tblHits.sort_values(by='theindex')
        tblHits = tblHits.set_index('theindex')
        
        bfIn = pysam.Samfile( o.inBam, 'rb' )
        qlast,nExact,nClose,nAll=None,0,0,0
        qIsInCurTbl = False

        # hits already recorded for name-sorted current query
        sHitAlreadyCurq=set()
        skipDueToRepeat=False

        ctr=0

        for l in bfIn:
            
            if ctr%100000 == 0:
                print '%d (%d %d %d)..'%(ctr, nExact, nClose, nAll )
            ctr+=1
            
            if l.qname!=qlast:
                if (qlast is not None) and qIsInCurTbl:
                    if int(qlast.split('_')[1]) in tblHits.index:
                        # this is apparently VERY slow
                        #tblHits.ix[ int(qlast.split('_')[1]) ] = [nExact,nClose,nAll]  
                        tblHits.ix[ int(qlast.split('_')[1]), 'num_exact_hits'] = nExact
                        tblHits.ix[ int(qlast.split('_')[1]), 'num_close_hits'] = nClose
                        tblHits.ix[ int(qlast.split('_')[1]), 'num_all_hits'] = nAll

                nExact,nClose,nAll=0,0,0
                qlast=l.qname
                qIsInCurTbl = qlast.split('_')[0]==tblName
                sHitAlreadyCurq=set()
            
            if qIsInCurTbl:

                tupCurHit = ( l.reference_id, l.is_reverse, l.reference_start, l.reference_end )
                if o.disregardDoubleHitsAtSameLocation:
                    if tupCurHit not in sHitAlreadyCurq:
                        sHitAlreadyCurq.add( tupCurHit )
                        skipDueToRepeat=False
                    else:
                        skipDueToRepeat=True
                
                if not skipDueToRepeat:
                    nm = dict(l.tags)['NM'] + (l.query_length - l.query_alignment_length)
                    if nm == 0:
                        nExact+=1
                    elif nm <= int( o.maxedsClose * l.query_length ):
                        nClose+=1
                    else:
                        nAll+=1
        
    
        if (qlast is not None) and qIsInCurTbl:
            if int(qlast.split('_')[1]) in tblHits.index:
                tblHits.ix[ int(qlast.split('_')[1]) ] = [nExact,nClose,nAll]

        candsStore.put( '%s_%s'%(tblName,o.newtblSuffix), tblHits )
        candsStore.get_storer( '%s_%s'%(tblName,o.newtblSuffix) ).attrs['mippipe__source_table_name'] = tblName

    