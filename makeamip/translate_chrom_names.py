import sys
import os
from optparse import OptionParser

import copy

from makeamip.capcommon import *

import pysam

if __name__=='__main__':
    opts=OptionParser()
    opts.add_option('','--inSam',dest='inSam')
    opts.add_option('','--outSam',dest='outSam')

    opts.add_option('','--outGenomeDict',default=None,dest='outGenomeDict')
    
    opts.add_option('','--chromMap',dest='chromMap')
    
    (o,args)=opts.parse_args()
    
    assert o.outGenomeDict is not None
    assert o.chromMap is not None

    mChromLenNew = parseContigLengths( o.outGenomeDict )

    mChromOldNew = mapFromCols(o.chromMap,'chrOld','chrNew')

    # for every chrom mapping from old to new, make sure we know the length of the new on.e
    assert all( [ chrom in mChromLenNew for chrom in list(mChromOldNew.values()) ] )

    if o.inSam == '/dev/stdin':
        bi=pysam.Samfile('-','rb')
    else:
        bi=pysam.Samfile(o.inSam,'rb')

    hdrIn = bi.header

    # hdrIn['SQ'] = [ {'LN':ci['LN'],'SN':mChromOldNew[ci['SN']]} for ci in hdrIn['SQ'] if ci['SN'] in mChromOldNew ]

    # build header with new lengths, in same order as provided within genome dictionary
    hdrNewSq = [ {'LN': mChromLenNew[chrom] , 'SN': chrom}  for chrom in mChromLenNew ]

    # hdrOut = copy.deepcopy(hdrIn)
    hdrOut['SQ']=hdrNewSq

    if o.outSam=='/dev/stdout':
        bo=pysam.Samfile('-','wb', header=hdrOut )
    else:
        bo=pysam.Samfile(o.outSam,'wb', header=hdrOut )

    sTidSuppress = set()
    mTidOldNew = {}
    for sqOld in bi.header['SQ']:

        tidOld = bi.gettid( sqOld['SN'] )
        if sqOld['SN'] not in mChromOldNew:
            sTidSuppress.add( tidOld )
        else:
            tidNew = bo.gettid( mChromOldNew[sqOld['SN']] )
            mTidOldNew[ tidOld ] = tidNew
        
    while True:
        try:
            l1=next(bi)
            tid = l1.tid
            if tid==-1:
                bo.write(l1)
            elif tid not in sTidSuppress:
                if (not l1.is_paired) or (l1.mate_is_unmapped):
                    l1.tid = mTidOldNew[ l1.tid ]
                    bo.write(l1)
                elif l1.rnext not in sTidSuppress:
                    l1.rnext = mTidOldNew[ l1.rnext ]
                    bo.write(l1)
        except StopIteration as e :
            break
