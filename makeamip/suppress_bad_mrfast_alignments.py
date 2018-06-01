import sys
import os
from optparse import OptionParser

import tempfile

import pysam 

from makeamip.capcommon import *

def parse_cigar(cigar):
    ops = dict( (i,j) for (j, i) in enumerate('MIDNSHP') )
    result = []
    n = ''
    for c in cigar:
        if c.isdigit():
            n += c
        elif c in ops:
            result.append( (ops[c], int(n)) )
            n = ''
    return result

if __name__=='__main__':
    opts=OptionParser()
    
    opts.add_option('','--inSam',default='/dev/stdin',dest='inSam') 
    opts.add_option('','--outSam',default='/dev/stdout',dest='outSam')

    opts.add_option('','--debugOut',default='/dev/null',dest='debugOut') 
    
    (o,args)=opts.parse_args()
 
    filOut=open(o.outSam,'w')
    filDbg=open(o.debugOut,'w')

    nSuppressedSeqLenNeQualLen, nSuppressedCigLenNeSeqLen, nSuppressedTruncdline, nPassed = 0 , 0, 0, 0

    for l in open(o.inSam,'r'):
        if l[0]=='@': 
            filOut.write( l )
        else:
            _l=l
            l = l.rstrip('\n').split('\t')
            if len(l)<10:
                nSuppressedTruncdline+=1
                filDbg.write('truncated: %s\n'%l[0])
                continue
            cigstr = parse_cigar(l[5])
            ciglen = sum( [ x[1] for x in cigstr if x[0]!=2 ] )
            seqlen = len( l[9] )
            quallen = len( l[10] )

            if seqlen!=quallen:
                nSuppressedSeqLenNeQualLen+=1
                filDbg.write('seqlen != quallen : %s  %d %d \n'%(l[0],seqlen,quallen))
            elif ciglen!=seqlen:
                nSuppressedCigLenNeSeqLen+=1
                filDbg.write('ciglen != seqlen : %s  %d %d \n'%(l[0],ciglen,seqlen))
            else:
                nPassed+=1
                filOut.write(_l)

    sys.stderr.write('%d written, %d suppressed (%d with cigar len!=seq len, %d with truncated line, %d seqlen!=quallen)\n'%( nPassed, nSuppressedTruncdline+nSuppressedCigLenNeSeqLen+nSuppressedSeqLenNeQualLen,nSuppressedCigLenNeSeqLen,nSuppressedTruncdline,nSuppressedSeqLenNeQualLen))


