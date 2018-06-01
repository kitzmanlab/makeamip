from __future__ import print_function
from builtins import range
import sys
from optparse import OptionParser
from collections import defaultdict

import os.path
import os.path as op

import copy

from makeamip.capcommon import *

#import numpy as np

#import pysam

def isLocalPath(fn):
    return ':' in fn

def hostPath(fn):
    return fn.split(':')[0],fn.split(':')[1]

class RunKey(tabgridwio):
    def __init__(self, keyfn):
        tabgridwio.__init__(self,open(keyfn,'r'))
        
    def attachKnownBarcs(self, lBarcs):
        self.lAllBarcs = lBarcs[:]
        
        for i in range(len(self)):
            if self[i]['barcs_incl']=='ALL':
                self[i]['barcs_incl']=self.lAllBarcs
        
    def write(self, fil):
        fil.write( tabnl( self.hdr ) )
        for i in range(len(self)):
            lout = [ self[i][k] if k in self[i] else '' for k in self.hdr ]
            for j in range(len(lout)):
                if type(lout[j])==list: lout[j]=','.join(lout[j])
            fil.write(tabnl(lout))
        
        fil.flush()
        
    def expandToRunLibKeys( self, 
                           colBcLut,  # column giving bc -> libmapping
                           lNewColNames=[],  
                           lfxnNewColExprs=[],  # list of f(rkline,libname,barcode) -> value
                           lColsPropagate=None,  # if set, assume these are the only ones to propagate
                           lColsExclude=None,   # if set, assume these are the only to exclude
                            ):
        return self.OLD_expandToRunLibKeys( colBcLut,lNewColNames,lfxnNewColExprs,lColsPropagate,lColsExclude )

    # otherTbl is likely bc --> libmapping
    #    many --> 1 is OK
    # 
    #  so then for each line in self:
    #       for each unique value in otherTbl[ otbl_colGroupBy ], 
    #           write out an entry with self cols, otherTbl[ otbl_colGroupBy ], and aggregated list of otherTbl[ colAg ] for colAg in otbl_lColsAggregate
    def expandRunKeyGroupByCol(self,
                               colFnOtherTbl,
                               otbl_colGroupBy,
                               otbl_lColsAggregate=None,
                               lNewColNames=[],  
                               lfxnNewColExprs=[],  # list of f(rkline,otblGroupByVal,listOfOtblLine) -> value
                               lColsPropagate=None,  # if set, assume these are the only ones to propagate
                               lColsExclude=None,   # if set, assume these are the only to exclude
                              ):

        assert lColsPropagate is None or lColsExclude is None

        lRlkOut =  [ tabgridwio(None) for _ in range(len(self)) ]

        for ir in range(len(self)):

            assert os.path.exists( self[ir][colFnOtherTbl] ), self[ir][colFnOtherTbl] 
            otherTbl = tabgridwio( open(self[ir][colFnOtherTbl],'r') )

            assert otbl_colGroupBy in otherTbl.hdr
            assert otbl_lColsAggregate is None or all( [ oc in otherTbl.hdr for oc in otbl_lColsAggregate ])

            mgrpbyvalLiotbl = defaultdict(list)
            for io in range(len(otherTbl)):
                mgrpbyvalLiotbl[ otherTbl[io][otbl_colGroupBy] ].append( io )

            rlkOut=lRlkOut[ir]

            rlkOut.hdr.append( otbl_colGroupBy )
            for ocol in otbl_lColsAggregate:
                rlkOut.hdr.append( ocol )

            mColExpr={}
            if lColsPropagate is not None:
                for cp in lColsPropagate:
                    if cp in self.hdr:
                        rlkOut.hdr.append( cp )
                        mColExpr[ cp ] = lambda lin,otgrpBy,lotlin,_cp=cp:lin[_cp]
            elif lColsExclude is not None:
                for cp in self.hdr:
                    if cp not in lColsExclude:
                        rlkOut.hdr.append( cp )
                        mColExpr[ cp ] = lambda lin,otgrpBy,lotlin,_cp=cp:lin[_cp]

            for icn in range(len(lNewColNames)):
                rlkOut.hdr.append( lNewColNames[icn] )
                mColExpr[ lNewColNames[icn] ] = \
                    lambda lin,otgrpBy,lotlin,_icn=icn:lfxnNewColExprs[_icn](lin,otgrpBy,lotlin)

            for ogv in mgrpbyvalLiotbl:

                mOutCur={}

                mOutCur[otbl_colGroupBy] = ogv
                if otbl_lColsAggregate is not None:
                    for ocol in otbl_lColsAggregate:
                        mOutCur[ocol] = ','.join( [otherTbl[ io ][ ocol ] for io in mgrpbyvalLiotbl[ ogv ]] )

                for col in mColExpr:
                    mOutCur[col] = mColExpr[ col ]( self[ir], ogv, [ otherTbl[io] for io in mgrpbyvalLiotbl[ogv] ] )
                
                rlkOut.ll.append( mOutCur )

        return lRlkOut

    # bc -> libmapping 
    #  this should be encoded by a column
    #  can either be
    #      (default) mapping consisting of [bc:libname],...,[bc:libname],:libname  
    #                                     (last one is default & optional)
    #      (else)    filename of a lut
    # 

    def OLD_expandToRunLibKeys( self, 
                           colBcLut,  # column giving bc -> libmapping
                           lNewColNames=[],  
                           lfxnNewColExprs=[],  # list of f(rkline,libname,barcode) -> value
                           lColsPropagate=None,  # if set, assume these are the only ones to propagate
                           lColsExclude=None,   # if set, assume these are the only to exclude
                            ):
        
        assert lColsPropagate==None or lColsExclude==None
        assert colBcLut in self.hdr
        
        mColExpr={}

        lRlkOut =  [ tabgridwio(None) for _ in range(len(self)) ]
        
        if lColsPropagate:
            for cp in lColsPropagate:
                if cp in self.hdr:
                    for rlk in lRlkOut:
                        rlk.hdr.append( cp )
                    mColExpr[ cp ] = lambda lin,libname,bc,_cp=cp:lin[_cp]
        elif lColsExclude:
            for cp in self.hdr:
                if not cp in lColsExclude:
                    for rlk in lRlkOut:
                        rlk.hdr.append( cp )
                    mColExpr[ cp ] = lambda lin,libname,bc,_cp=cp:lin[_cp]
        
        for rlk in lRlkOut:
            rlk.hdr.append( 'barcname' )
        for rlk in lRlkOut:
            rlk.hdr.append( 'libname' )
        
        for icn in range(len(lNewColNames)):
            for rlk in lRlkOut:
                rlk.hdr.append( lNewColNames[icn] )
            mColExpr[ lNewColNames[icn] ] = \
                lambda lin,libname,bc,_icn=icn:lfxnNewColExprs[_icn](lin,libname,bc)       
        
        # always each bc shuold go to only one libname
        lBcMapping=[]
        for ir in range(len(self)):
            r = self[ir]
            bcLut = r[ colBcLut ]
            if bcLut.find('lambda')!=-1:
                lBcMapping.append({})
                fxnBcToLib=eval(bcLut)
                #for bc in self.lAllBarcs:
                
                assert len(r['barcs_incl'])>0 , 'if using bc->lib fxn then need explicit list of bcs'
                
                for bc in r['barcs_incl'].split(','):
                    lBcMapping[-1][ bc ] = fxnBcToLib(bc)
            elif ':' in bcLut:
                lBcMapping.append({})
                bcLut=[lutitem.strip() for lutitem in bcLut.split(',')]                
                for lutitem in bcLut:
                    assert not lutitem.split(':')[0] in lBcMapping[-1], 'each bc should go to only 1 libname'
                    lBcMapping[-1][ lutitem.split(':')[0] ] = lutitem.split(':')[1]
            else:            
                mBarcName=mapFromCols( bcLut , 'barcname', 'libname' )
                lBcMapping.append( mBarcName )
                        
            #assert len(set(lBcMapping[-1].keys()))==len(set(lBcMapping[-1].values())), 'bc <-> libname should be 1:1' 
            if len(set(lBcMapping[-1].keys()))!=len(set(lBcMapping[-1].values())):

                sys.stderr.write( 'WARNING : bc <-> libname should be 1:1\n' )
                sys.stderr.flush()

                for k in lBcMapping[-1]:
                    for j in lBcMapping[-1]:
                        if j!=k and lBcMapping[-1][j] == lBcMapping[-1][k]:
                            sys.stderr.write( 'bc k, j: %s %s --> %s %s\n'%(j,k,lBcMapping[-1][j],lBcMapping[-1][k]) )
                            sys.stderr.flush()

        
        for ir in range(len(self)):
            mBcLibname=lBcMapping[ir]
            rlkOut=lRlkOut[ir]
            
            # for bc in sorted( mBcLibname,
            #                  key=lambda _bc:dict(zip( self.lAllBarcs, range(len(self.lAllBarcs)) )+[('ALL_OTHER',9999999)])[_bc]) :
                
            for bc in sorted( mBcLibname.keys() ):


                        
                if (self[ir]['barcs_incl'].__len__()>0) and (not bc in self[ir]['barcs_incl'].split(',')): continue
                
                mOutCur={}

                mOutCur['barcname'] = bc
                mOutCur['libname'] = mBcLibname[bc]

                for col in mColExpr:
                    mOutCur[col] = mColExpr[ col ]( self[ir], mOutCur['libname'], mOutCur['barcname'])
                
                rlkOut.ll.append( mOutCur )
        
        return lRlkOut
        
        
def getBarcs(fn):
    filBarcs=open(fn,'r')
    l=filBarcs.readline()
    lb=[]
    if l[0]=='>':
        lb.append(l[1:].rstrip())
        l=filBarcs.readline()
        l=filBarcs.readline()
        while len(l)>0:
            lb.append(l[1:].rstrip())
            l=filBarcs.readline()
            l=filBarcs.readline()
    else:
        lb.append(l.rstrip())
        while len(l)>0:
            lb.append(l.rstrip())
    return lb
       


def tasklistToSimpleSh( lTasks, backgroundify=False, expandify=False ):
    llOut=[]
    for i in range(len(lTasks)):
        if type(lTasks[i])==tuple:
            curTask=lTasks[i][1]
        elif type(lTasks[i])==dict:
            curTask=lTasks[i]['cmd']
        else:
            curTask=lTasks[i]
        llOut.append( '(  %s  )%s'%( os.path.expandvars(os.path.expanduser(curTask)) if expandify else curTask, 
                         '&' if backgroundify else '') )
    return llOut

def tasklistToMakepp( lTasks, backgroundify=False, expandify=False ):
    fxnFilt = lambda t:'%s'%( os.path.expandvars(os.path.expanduser(t)) if expandify else t.replace('$','$$')) 
    
    #llOut=['$(phony all): %s'%(' '.join( ['job%05d'%(i) for i in xrange(len(lTasks))] )) ]
    llOut=[]
    ljNames, mJobnamePhony=[], {}
    
    for i in range(len(lTasks)):
        if isinstance(lTasks[i],basestring):
            llOut.append('job%05d:\n\t%s'%(i,fxnFilt( lTasks[i] ) ))
            ljNames.append( 'job%05d'%i )
            mJobnamePhony[ ljNames[-1] ]=True
        elif type(lTasks[i])==tuple or type(lTasks[i])==list:
            assert ( type(lTasks[i][0])==list or type(lTasks[i][0])==tuple) and isinstance(lTasks[i][1],str)
            assert len(lTasks[i])==2 
            llOut.append('job%05d: %s\n\t%s'%(i, ' '.join(lTasks[i][0]) , fxnFilt( lTasks[i][1] )) )
            ljNames.append( 'job%05d'%i )
            mJobnamePhony[ ljNames[-1] ]=True
        elif type(lTasks[i])==dict:
            if 'output' not in lTasks[i]:
                if 'deps' in lTasks[i]:
                    #llOut.append('$(phony job%05d): %s\n\t%s'%(i, ' '.join(lTasks[i]['deps']) , fxnFilt( lTasks[i]['cmd'] )) )
                    llOut.append('job%05d: %s\n\t%s'%(i, ' '.join(lTasks[i]['deps']) , fxnFilt( lTasks[i]['cmd'] )) )
                    ljNames.append( 'job%05d'%i )
                    mJobnamePhony[ ljNames[-1] ]=True
                else:
                    #llOut.append('$(phony job%05d): %s\n\t%s'%(i, ' '.join(lTasks[i]['deps']) , fxnFilt( lTasks[i]['cmd'] )) )
                    llOut.append('job%05d: %s\n\t%s'%(i, ' '.join(lTasks[i][0]) , fxnFilt( lTasks[i][1] )) )
                    ljNames.append( 'job%05d'%i )
                    mJobnamePhony[ ljNames[-1] ]=True
            else:
                llOut.append('%s: %s # job %05d\n\t%s'%( lTasks[i]['output'],
                                                         ' '.join(lTasks[i]['deps']),
                                                         i,
                                                         fxnFilt(lTasks[i]['cmd']) ) )
                ljNames.append( lTasks[i]['output'] )
                mJobnamePhony[ ljNames[-1] ]=False

    llOut = ['SHELL := /bin/bash\n'
             '.PHONY : all %s'%' '.join( [jn for jn in mJobnamePhony if mJobnamePhony[jn]] ) ,
             'all : %s'%(' '.join( ljNames )) ] + llOut
    
    return llOut
        

"""       
--fillInCols
   
   fills in columns using args
   each arg should be colname:lambda line->value
   
   
--loopCmds 
    
    runs the passed shell snippet once for every item in this key, and 
    sets the env vars to equal the values of the column for that item
    
--taskListToMakefile
    
    convert a tasklist to a makefile
    
""" 

import gc
def disable_stdout_buffering():
    # Appending to gc.garbage is a way to stop an object from being
    # destroyed.  If the old sys.stdout is ever collected, it will
    # close() stdout, which is not good.
    gc.garbage.append(sys.stdout)
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    
disable_stdout_buffering()

       
if __name__ == '__main__':       

    opts=OptionParser()
    
    opts.add_option('','--rkIn',dest='rkIn')
    opts.add_option('','--rkOut',default=None,dest='rkOut')
    
    opts.add_option('','--tasksBackgroundify',action='store_true',default=False,dest='tasksBackgroundify')
    opts.add_option('','--tasksAsMakefile',action='store_true',default=False,dest='tasksAsMakefile')
    
    ########################################################################

    opts.add_option('','--fillInCols',
                    action='store_true',default=False,dest='fillInCols')

    ########################################################################

    opts.add_option('','--loopCmds',
                    action='store_true',default=False,dest='loopCmds')
    
    opts.add_option('','--quiet',default=False,action='store_true',dest='quiet')

    opts.add_option('','--cmdFile',default=None,dest='cmdFile')    

    ########################################################################

    opts.add_option('','--taskListToMakefile',
                action='store_true',default=False,dest='taskListToMakefile')

    opts.add_option('','--makefileFxnLDeps',
                default=None, dest='makefileFxnLDeps')
    opts.add_option('','--makefileFxnOut',
                default=None, dest='makefileFxnOut')

    ########################################################################

    (o, args) = opts.parse_args()

    lTasksOut = None

    if o.fillInCols:
        rk=RunKey(o.rkIn)
        
        lNewCols,mNewColFxn=[],{}
        for a in args:
            colname,fxnLineCol=a.split(':')[0],':'.join(a.split(':')[1:])
            fxnLineCol=eval(fxnLineCol)
            lNewCols.append(colname)
            mNewColFxn[colname]=fxnLineCol

        for cn in lNewCols:            
            if not cn in rk.hdr: 
                rk.hdr.append(cn)
        
        for l in rk:
            for col in lNewCols:
                l[col]=mNewColFxn[col](l)
                
        filOut=open(o.rkOut,'w', 0) if o.rkOut!=None else sys.__stdout__
        
        rk.write(filOut)
               
    elif o.loopCmds:
        if o.makefileFxnLDeps is not None and o.makefileFxnOut is not None:
            makefileFxnLDeps=eval(o.makefileFxnLDeps)
            makefileFxnOut=eval(o.makefileFxnOut)
        else:
            makefileFxnLDeps=makefileFxnOut=None

        lTasksOut=[]
        
        rk=RunKey(o.rkIn)
        
        il=0
        
        if o.cmdFile!=None:
            snippet='; '.join([ l.rstrip() for l in open(o.cmdFile,'r') ])
        elif len(args)>0 and len(args[0])>0:
            snippet = '; '.join([a for a in args if len(a)>0])
        else:
            snippet='; '.join([ l.rstrip() for l in sys.stdin ])
        
        for l in rk:           
            lCmd = 'echo "job %d";'%il if not o.quiet else ''
            for k in l:
                if len( l[k].replace('"','\\"') ) == 0: continue
                cmdPiece='export %s="%s"; '%( k.replace('.','_') ,l[k].replace('"','\\"'))
                lCmd+=cmdPiece
                # print(type(lCmd),type(cmdPiece),type(l[k]),isinstance(l[k],basestring))


            lCmd+=snippet          
        
            if makefileFxnOut is not None:
                lDeps = makefileFxnLDeps( l )
                output = makefileFxnOut( l )
                lTasksOut.append({'deps': lDeps, 'cmd': lCmd, 'output': output } )
            else:
                lTasksOut.append(lCmd)
            
            il+=1
            
    elif o.taskListToMakefile:
        lTasks=[ striplin(t) for t in sys.stdin ]
        sys.stdout.write('\n'.join(tasklistToMakepp( lTasks )))        
        
    if lTasksOut!=None:
        if o.tasksAsMakefile:
            print('\n'.join( tasklistToMakepp(lTasksOut, o.tasksBackgroundify) ))
        else:
            print('\n'.join( tasklistToSimpleSh(lTasksOut, o.tasksBackgroundify) ))
                
