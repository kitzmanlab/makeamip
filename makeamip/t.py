from __future__ import print_function
from builtins import range
import sys
from optparse import OptionParser
from collections import defaultdict

import os.path
import os.path as op

import random

import copy

from makeamip.capcommon import *

def fl( x, defval=0 ):
    if x=='NA': return float(defval)
    else: return float(x)

if __name__ == '__main__':

    opts=OptionParser()

    opts.add_option('-i','--in',default=None,dest='fnIn')
    opts.add_option('-o','--out',default=None,dest='fnOut')

    opts.add_option('','--pandas',
                    action='store_true',default=False,dest='pandas')

    opts.add_option('','--pandasText',
                    action='store_true',default=False,dest='pandasText')

    ########################################################################

    opts.add_option('','--filter',
                    action='store_true',default=False,dest='filter')

    opts.add_option('','--filterFxn',dest='filterFxn')

    ########################################################################

    opts.add_option('','--fillInCols',
                    action='store_true',default=False,dest='fillInCols')

    opts.add_option('','--colInsBefore',
            default=None,dest='colInsBefore')

    ########################################################################

    opts.add_option('','--expandBy', action='store_true',default=False ,dest='expandBy')

    opts.add_option('','--in2', default=None,dest='fnIn2')

    opts.add_option('','--lcol2', default=None,dest='lcol2')

    opts.add_option('','--fxnSuppressUponExpand', default="lambda l1,l2:False",dest='fxnSuppressUponExpand')  # l1, l2 -> True will suppress that combination

    ########################################################################

    opts.add_option('','--join', action='store_true',default=False ,dest='join')

    opts.add_option('','--joinType', default='left',dest='joinType')
    opts.add_option('','--joinOn1', default=None,dest='joinOn1')
    opts.add_option('','--joinOn2', default=None,dest='joinOn2')

    opts.add_option('','--suff1', default='_1',dest='suff1')
    opts.add_option('','--suff2', default='_2',dest='suff2')

    ########################################################################

    opts.add_option('','--subtract', action='store_true',default=False ,dest='subtract')

    ########################################################################

    opts.add_option('','--chunk', action='store_true',default=False ,dest='chunk')

    opts.add_option('','--chunkSize', default=int(1e6), type=int, dest='chunkSize' )

    opts.add_option('','--chunkInputIsKey', default=False, action='store_true', dest='chunkInputIsKey')
    opts.add_option('','--chunkInputCol', default=None,dest='chunkInputCol')
    opts.add_option('','--chunkOutputCol', default=None,dest='chunkOutputCol')

    # fnIn,(line),idx-->filename
    opts.add_option('','--fxnChunkToOutFn', default=None,dest='fxnChunkToOutFn')

    ########################################################################

    (o, args) = opts.parse_args()

    if o.pandas:
        import pandas as pd

    lTasksOut = None


    if o.chunk:

        assert not o.pandasText
        assert not o.pandas

        if not o.pandas:
            fxnChunkToOutFn = eval(o.fxnChunkToOutFn)
            lkeylineOut=[]
            rk=None
            if o.chunkInputIsKey:
                rk = tabgridwio( open(o.fnIn,'r') )
                lFnIn = [ rl[o.chunkInputCol] for rl in rk ]
                lkeylineIn = [ rl for rl in rk ]
            else:
                lFnIn = [ o.fnIn ]
                lkeylineIn = [ None ]

            for ifn in range(len(lFnIn)):
                fnIn = lFnIn[ifn]
                filIn = pgzopen(o.fnIn,'r')
                hdrIn = filIn.readline().rstrip('\n')
                ifnoutNext = 0
                ioutNext = 0
                for l in filIn:
                    if ioutNext%o.chunkSize == 0:
                        fnOut = fxnChunkToOutFn( fnIn, lkeylineIn[ifn], ifnoutNext )
                        lkeylineOut.append({ '_chunkIdx':ifnoutNext,
                                                 o.chunkOutputCol:fnOut } )
                        filOut = open(fnOut,'w')
                        filOut.write(hdrIn)
                        filOut.write('\n')
                        ifnoutNext+=1
                    filOut.write(l)
                    ioutNext+=1

            rkout = tabgridwio( None )
            rkout.hdr=(rk.hdr if rk is not None else [])+['_chunkIdx',o.chunkOutputCol]
            rkout.ll = lkeylineOut
            rkout.write( open(o.fnOut,'w') )
        else:
            print('not implemented')

    elif o.join:

        assert not o.pandasText
        assert not o.pandas

        import pandas as pd

        if not o.pandas:
            in1 = pd.read_csv( o.fnIn, sep='\t' )
            in2 = pd.read_csv( o.fnIn2, sep='\t' )

            merged = pd.merge( in1, in2, left_on=o.joinOn1.split(','), right_on=o.joinOn2.split(','),
                               sort=False,
                               how=o.joinType,
                               suffixes=(o.suff1,o.suff2), copy=True   )

            merged.to_csv( o.fnOut, header=True, index=False, sep='\t' )

    elif o.subtract:

        assert not o.pandasText
        assert not o.pandas

        import pandas as pd

        if not o.pandas:
            in1 = pd.read_csv( o.fnIn, sep='\t' )
            in2 = pd.read_csv( o.fnIn2, sep='\t' )

            in1 = in1.set_index( o.joinOn1.split(',') )
            in2 = in2.set_index( o.joinOn2.split(',') )

            tblIntersection = pd.merge( in1, in2, left_index=True, right_index=True,
                                           sort=False,
                                           how='inner',
                                           suffixes=(o.suff1,o.suff2), copy=True   )

            tblUnion =  pd.merge( in1, in2,left_index=True, right_index=True,
                                           sort=False,
                                           how='outer',
                                           suffixes=(o.suff1,o.suff2), copy=True   )

            tblJoined = tblUnion[ ~tblUnion.index.isin( tblIntersection.index ) ]

            tblJoined = tblJoined.reset_index()

            tblJoined.to_csv( o.fnOut, header=True, index=False, sep='\t' )
        else:
            assert 1==0,'not implemented'

    elif o.fillInCols:

        assert not o.pandasText

        if not o.pandas:
            rk=tabgridwio(pgzopen(o.fnIn,'r'))

            lNewCols,mNewColFxn=[],{}
            for a in args:
                colname,fxnLineCol=a.split(':')[0],':'.join(a.split(':')[1:])
                fxnLineCol=eval(fxnLineCol)
                lNewCols.append(colname)
                mNewColFxn[colname]=fxnLineCol

            for cn in lNewCols:
                if not cn in rk.hdr:
                    if o.colInsBefore==None:
                        rk.hdr.append(cn)
                    else:
                        rk.hdr[ rk.hdr.index(o.colInsBefore):rk.hdr.index(o.colInsBefore) ] = [cn]

            for l in rk:
                for col in lNewCols:
                    l[col]=mNewColFxn[col](l)

            filOut=pgzopen(o.fnOut,'w') if o.fnOut!=None else sys.stdout

            rk.write(filOut)
        else:
            lNewCols,mNewColFxn=[],{}
            for a in args:
                colname,fxnLineCol=a.split(':')[0],':'.join(a.split(':')[1:])
                fxnLineCol=eval(fxnLineCol)
                lNewCols.append(colname)
                mNewColFxn[colname]=fxnLineCol

            assert len( o.fnIn.split(':')) ==2
            assert len( o.fnOut.split(':'))==2

            fnIn = os.path.abspath( o.fnIn.split(':')[0] )
            fnOut = os.path.abspath( o.fnOut.split(':')[0] )

            gkIn,gkOut=o.fnIn.split(':')[1], o.fnOut.split(':')[1]

            if fnIn==fnOut:
                hsIn = pd.HDFStore( fnIn, 'a' )
                hsOut = hsIn
            else:
                hsIn = pd.HDFStore( fnIn, 'r' )
                hsOut = pd.HDFStore( fnOut, 'w' )

            tblIn = hsIn[ gkIn ]
            lColsOld = [ c for c in tblIn.columns ]

            for cn in lNewCols:
                fxn = mNewColFxn[cn]
                lnewvals = [fxn(l[-1]) for l in tblIn.iterrows()]
                tblIn[cn] = lnewvals

            for cn in lNewCols:
                if not cn in lColsOld:
                    if o.colInsBefore==None:
                        lColsOld.append(cn)
                    else:
                        lColsOld[ lColsOld.index(o.colInsBefore):lColsOld.index(o.colInsBefore) ] = [cn]


            tblIn = tblIn[ lColsOld ]

            hsOut[gkOut] = tblIn

            if fnIn!=fnOut:
                for _gk in list(hsIn.keys()):
                    gk = _gk.split('/')[-1]
                    if gk!=gkIn:
                        hsOut[gk]=hsIn[gk]

            hsOut.close()
            if fnIn!=fnOut:
                hsIn.close()

    elif o.filter:

        assert not o.pandasText

        if not o.pandas:
            filIn=pgzopen(o.fnIn,'r')
            filOut = pgzopen(o.fnOut,'w')

            hdr=parseHdrToMap(filIn.readline().rstrip('\n'))

            hdrOut=sorted([k for k in list(hdr.keys())],key=lambda k:hdr[k])
#            if not o.skipHdrOut:
            filOut.write('%s\n'%'\t'.join(hdrOut))

            filterfxn=eval(o.filterFxn)

            for l in filIn:
                ml = parseLineToMap( l.rstrip('\n'), hdr )

                if filterfxn(ml):
                    filOut.write( l )
            filOut.close()
        else:
            assert len( o.fnIn.split(':')) ==2
            assert len( o.fnOut.split(':'))==2
            fnIn = os.path.abspath( o.fnIn.split(':')[0] )
            fnOut = os.path.abspath( o.fnOut.split(':')[0] )

            gkIn,gkOut=o.fnIn.split(':')[1], o.fnOut.split(':')[1]

            if fnIn==fnOut:
                hsIn = pd.HDFStore( fnIn, 'a' )
                hsOut = hsIn
            else:
                hsIn = pd.HDFStore( fnIn, 'r' )
                hsOut = pd.HDFStore( fnOut, 'w' )

            filterfxn=eval(o.filterFxn)

            tblIn = hsIn[ gkIn ]
            tblIn = tblIn[ [filterfxn(l[-1]) for l in tblIn.iterrows()] ]

            hsOut[gkOut] = tblIn

            if fnIn!=fnOut:
                for _gk in list(hsIn.keys()):
                    gk = _gk.split('/')[-1]
                    if gk!=gkIn:
                        hsOut[gk]=hsIn[gk]

            hsOut.close()
            if fnIn!=fnOut:
                hsIn.close()

    elif o.expandBy:
        assert not o.pandasText
        assert not o.pandas

        rk=tabgridwio(pgzopen(o.fnIn,'r'))
        rk2=tabgridwio(pgzopen(o.fnIn2,'r'))

        fxnSuppressUponExpand = eval(o.fxnSuppressUponExpand)

        for col2 in o.lcol2.split(','):
            assert col2 in rk2.hdr, col2
            assert not col2 in rk.hdr, col2

        rkout = tabgridwio( None )
        rkout.hdr = rk.hdr + o.lcol2.split(',')

        for l1 in rk:
            for l2 in rk2:
                if fxnSuppressUponExpand( l1, l2 ): continue

                lout = dict( [ (k1,l1[k1]) for k1 in l1 if k1!='__lineidx'] )
                for k2 in o.lcol2.split(','):
                    lout[k2]=l2[k2]

                rkout.ll.append( lout )

        filOut=pgzopen(o.fnOut,'w') if o.fnOut!=None else sys.stdout

        rkout.write( filOut )



