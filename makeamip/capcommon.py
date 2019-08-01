import random 
import numpy as np

from past.utils import old_div

from collections import OrderedDict, defaultdict
import gzip
import sys
import subprocess

######################################################################
# from capcommon


_aaAbbr=[ ('A','Ala'), ('R','Arg'), ('N','Asn'), ('D','Asp'), ('C','Cys'),
          ('Q','Gln'), ('E','Glu'), ('G','Gly'), ('H','His'), ('I','Ile'),
          ('L','Leu'),
           ('K','Lys'), ('M','Met'), ('F','Phe'), ('P','Pro'), ('S','Ser'),
           ('T','Thr'), ('W','Trp'), ('Y','Tyr'), ('V','Val'), ('*','End') ]
mAaAbbr = dict( _aaAbbr )
mAbbrAa = dict( [(x[1],x[0]) for x in _aaAbbr ] )

rcMapRNA = {'A':'U', 'C': 'G', 'G': 'C', 'U': 'A', 'T': 'A', 'N': 'N', '-': '-', ' ': ' ',
            'a':'u', 'c': 'g', 'g': 'c', 'u': 'a', 't': 'a', 'n': 'n'}
rcMapDNA = {'A':'T', 'C': 'G', 'G': 'C', 'U': 'A', 'T': 'A', 'N': 'N', '-': '-', ' ': ' ',
            'a':'t', 'c': 'g', 'g': 'c', 'u': 'a', 't': 'a', 'n': 'n'}

rcMapDNAAmbigs = { 'A':'T', 'C':'G', 'G':'C', 'T':'A',
                   'M':'K', 'R':'Y', 'W':'W', 'S':'S',
                   'Y':'R', 'K':'M', 'V':'B', 'H':'D',
                   'D':'H', 'B':'V', 'N':'N' }

mDisambig= {'A':'A', 'C':'C', 'G':'G', 'U':'U', 'T':'T', 'R':'N', 'Y':'N', 'M':'N', \
            'K':'N', 'W':'N', 'S':'N', 'B':'N', 'D':'N', 'H':'N', 'V':'N', 'N':'N' }

#    Code Ambiguity  Complement Mnemonic
#A    A          T          adenine
#B    CGT        V          not_adenine
#C    C          G          cytosine
#D    AGT        H          not_cytosine
#G    G          C          guanine
#H    ACT        D          not_guanine
#K    GT         M          keto_base
#M    AC         K          amino_base
#N    ACGT       N          any_base
#R    AG         Y          purine_base
#S    CG         S          strong_bond
#T    T          A          thymine
#U    T          A          uracil
#V    ACG        B          not_thymine/uracil
#W    AT         W          weak_bond
#X    ACGT       X          unknown
#Y    CT         R          pyrimidine

mAmbigBases = { 'A':'A', 'B':'CGT', 'C':'C', 'D':'AGT', 'G':'G',
               'H':'ACT', 'K':'GT', 'M':'AC', 'N':'ACGT', 'R':'AG',
               'S':'CG', 'T':'T', 'V':'ACG', 'W': 'AT', 'X': 'ACGT',
               'Y': 'CT' }
mBasesAmbig = dict( [ (mAmbigBases[k],k) for k in mAmbigBases ] )

mDisambigBothCases = {}

for k in mDisambig: 
    mDisambigBothCases[k]=mDisambig[k]
    mDisambigBothCases[k.lower()]=mDisambig[k].lower()    

def disambigToUCase(seq):
    return ''.join( [mDisambig[b] for b in seq.upper()] )            

def disambig(seq):
    return ''.join( [mDisambigBothCases[b] for b in seq] )            

def revComp(seq, useRNA=False):
    return ''.join( [useRNA and rcMapRNA[b] or rcMapDNA[b] for b in reversed(seq)])

def revCompDnaWAmbigs(seq):
    return ''.join( [ rcMapDNAAmbigs[b] for b in reversed(seq)] )

def sort_dict_keys_by_value( dict ):
    return sorted( list(dict.keys()), key = lambda x: dict[x] )

# pygr seqdb's are [00).  
# intvCorrection here is the correction applied to the input to get pygr coordinates.
# python in general is [00) so this is a general tool.
# def pygrCoordCorrection( intvType='[00]' ):
#     intvCorrection=[0,0]
#     intvCorrection[0] = ( intvType[0]=='(' and 1 or
#                           intvType[0]=='[' and 0 or
#                           0 ) +\
#                         ( intvType[1]=='0' and 0 or
#                           intvType[1]=='1' and -1 or
#                           0 ) 
#     intvCorrection[1] = ( intvType[3]==')' and 0 or
#                           intvType[3]==']' and 1 or
#                           0 ) +\
#                         ( intvType[1]=='0' and 0 or
#                           intvType[1]=='1' and -1 or
#                           0 ) 
#     return tuple(intvCorrection)

def pygrCoordCorrection( intvType='[00]' ):
    return coordConversion( typeFrom=intvType, typeTo='[00)' )

def coordConversion( typeFrom='[01]', typeTo='[00]' ):
    intvCorrectionFrom=[0,0]
    intvCorrectionFrom[0] = ( typeFrom[0]=='(' and 1 or
                          typeFrom[0]=='[' and 0 or
                          0 ) +\
                        ( typeFrom[1]=='0' and 0 or
                          typeFrom[1]=='1' and -1 or
                          0 ) 
    intvCorrectionFrom[1] = ( typeFrom[3]==')' and 0 or
                          typeFrom[3]==']' and 1 or
                          0 ) +\
                        ( typeFrom[1]=='0' and 0 or
                          typeFrom[1]=='1' and -1 or
                          0 ) 

    intvCorrectionTo=[0,0]
    intvCorrectionTo[0] = ( typeTo[0]=='(' and 1 or
                          typeTo[0]=='[' and 0 or
                          0 ) +\
                        ( typeTo[1]=='0' and 0 or
                          typeTo[1]=='1' and -1 or
                          0 ) 
    intvCorrectionTo[1] = ( typeTo[3]==')' and 0 or
                          typeTo[3]==']' and 1 or
                          0 ) +\
                        ( typeTo[1]=='0' and 0 or
                          typeTo[1]=='1' and -1 or
                          0 ) 

    return ( intvCorrectionFrom[0]-intvCorrectionTo[0], 
             intvCorrectionFrom[1]-intvCorrectionTo[1] )


coordCorrectionToPyStyle = pygrCoordCorrection

tabnl=lambda l:'\t'.join([str(x) for x in l])+'\n'

def parseHdrToMap(lineHdr, delim='\t'):
    hdr = lineHdr.split(delim)
    return dict(list(zip(hdr,list(range(len(hdr))))))    

def parseLineToMap(line, mHdr, delim='\t'):
    m={}
    l=line.split(delim)
    for k in mHdr:
        i=mHdr[k]
        if i<len(l): m[k]=l[i]
    return m

def tabDelimFileToColList( fil ):
    hdr=fil.readline().rstrip()
    mhdr=parseHdrToMap(hdr)
    mColLL=defaultdict(list)
    for l in fil:
        l=parseLineToMap(l.rstrip(),mhdr)
        for k in l: mColLL[k].append(l[k])
    return mColLL

def mapFromCols(fn,colFrom,colTo):
    fil=pgzopen(fn,'r',nThreads=1)
    mf=tabDelimFileToColList(fil)
    lk,lv=mf[colFrom],mf[colTo]
    assert len(lk)==len(set(lk)), 'key column must be unique'
    return dict(list(zip(lk,lv)))

def tabDelimFileHdrRowsColsStatic( fil ):
    hdr=fil.readline()
    mhdr=parseHdrToMap( hdr.rstrip() )
    ll=[ parseLineToMap(l.rstrip(),mhdr) for l in fil ]
    lc=dict( [(col,[ll[i][col] for i in range(len(ll))]) for col in mhdr ] )
    return hdr.rstrip().split('\t'),ll,lc

class tabgrid(object):
    
    class bycolv(object):
        def __init__(self, tginst, col):
            self.tginst=tginst
            self.col=col
        def __getitem__(self,k):
            if type(k)==slice:
                ixs = k.indices( len(self.tginst.ll ) )
                return [ self.tginst.ll[i][self.col] for i in range(ixs[0],ixs[1],ixs[2]) ]
            else:
                assert type(k)==int
                return self.tginst.ll[k][self.col]
        def __setitem__(self,k,v):
            if type(k)==slice:
                ixs = k.indices( len(self.tginst.ll ) )
                j=0
                for i in range(ixs[0],ixs[1],ixs[2]):
                    self.tginst.ll[i][self.col]=v[j]
                    j+=1
            else:
                assert type(k)==int
                self.tginst.ll[k][self.col]=v
            
    
    def __init__(self,ll):
        self.ll=ll
        
    def __getitem__(self,k):
        if type(k) in (slice,int):
            return self.ll[k]
        else:
            assert type(k)==str
            return tabgrid.bycolv(self,k)    
        
    def __len__(self):
        return len(self.ll)
    
    def __iter__(self):
        return self.ll.__iter__()

def snl(x):
    if len(x)==0: return x
    if x[-1]=='\n': return x[:-1]
    else: return x

class tabgridwio(tabgrid):
    def __init__(self, fil=None):
        if fil is not None:
            self.hdr=fil.readline().rstrip('\n').lstrip('#')
            mhdr=parseHdrToMap( self.hdr )
            self.hdr=self.hdr.split('\t')
            #ll=[ parseLineToMap(l.rstrip(), mhdr) for l in fil ]
            ll=[ parseLineToMap(snl(l), mhdr) for l in fil ]

            for i in range(len(ll)):
                ll[i]['__lineidx'] = str(i)
        else:
            ll=[]
            self.hdr=[]

        tabgrid.__init__( self, ll )
        
    def write(self, fil):
        fil.write( tabnl( self.hdr ) )
        for i in range(len(self)):
            lout = [ self[i][k] for k in self.hdr ]
            fil.write(tabnl(lout))
            
        
def tabDelimFileHdrRowsCols( fil ):
    hdr=fil.readline()
    mhdr=parseHdrToMap( hdr.rstrip('\n') )
    ll=[ parseLineToMap(l.rstrip('\n'),mhdr) for l in fil ]
    
    return hdr.rstrip('\n').split('\t'), tabgrid(ll)
    
#    lc=dict( [(col,[ll[i][col] for i in xrange(len(ll))]) for col in mhdr ] )
#    return hdr.rstrip().split('\t'),ll,lc


def striplin(l) : return l.replace('\n','').strip()

def spanIncl(rng): return abs(rng[1]-rng[0]+1)

def spanInclAllowNone(rng): 
    if rng==None: return 0
    else: return abs(rng[1]-rng[0]+1) 

def zeroPadChrom(chrom):
    matches = re.compile(r'(?P<chr>[^0-9]+)(?P<num>[0-9]+)').match(chrom)
    return matches and \
        '%s%02d'%(matches.group('chr'),int(matches.group('num'))) or\
        chrom

def deZeroPadChrom(chrom):
    matches = re.compile(r'(?P<chr>[^0-9]+)(?P<num>[0-9]+)').match(chrom)
    return matches and \
        '%s%d'%(matches.group('chr'),int(matches.group('num'))) or\
        chrom
    
def chrStartStop(s):
    return s.split(':')[0],int(s.split(':')[1].split('-')[0]),int(s.split(':')[1].split('-')[1])
    
def chrTupStartStop(s):
    return s.split(':')[0],(int(s.split(':')[1].split('-')[0]),int(s.split(':')[1].split('-')[1]))
    
def splitRejoinSlice(sep, x, slc):
    return sep.join(x.split(sep)[slc])

def splitjoin(sep, x, slc):
    return sep.join(x.split(sep)[:slc])

def pick5pMost(cooRange,strand):
    return max(cooRange) if strand==-1 else min(cooRange)
def pick3pMost(cooRange,strand):
    return min(cooRange) if strand==-1 else max(cooRange)
def is5pOf(coord,cooRange,strand):
    if strand==1: return coord < min(cooRange)
    if strand==-1: return coord > max(cooRange)
    assert 1==-1, '%s %s'%(repr(strand),type(strand))
def is3pOf(coord,cooRange,strand):
    if strand==1: return coord > max(cooRange)
    if strand==-1: return coord < min(cooRange)
    assert 1==-1, '%s %s'%(repr(strand),type(strand))
def isWithin(cooCoord, cooRange, fIsInclusive=True):
    if fIsInclusive: return (cooCoord <= max(cooRange) and cooCoord >= min(cooRange))
    else:            return (cooCoord <  max(cooRange) and cooCoord >  min(cooRange))
def rgOverlap(bounds1, bounds2):
    if (min(bounds1) >= min(bounds2) and min(bounds1) <= max(bounds2) ) or \
       (max(bounds1) >= min(bounds2) and max(bounds1) <= max(bounds2) ) or \
       (min(bounds2) >= min(bounds1) and min(bounds2) <= max(bounds1) ) or \
       (max(bounds2) >= min(bounds1) and max(bounds2) <= max(bounds1) ):
        return 1
    else:
        return 0
    
def rgOverlap_array_single(arbounds1, bounds2):

    if arbounds1.shape[0]==0: return np.zeros_like( arbounds1 )

    sbounds2=(min(bounds2),max(bounds2))
    arsbounds1=np.sort(arbounds1,1)

    return ( (( arsbounds1[:,0] >= sbounds2[0] ) & ( arsbounds1[:,0] <= sbounds2[1] )) | \
             (( arsbounds1[:,1] >= sbounds2[0] ) & ( arsbounds1[:,1] <= sbounds2[1] )) | \
             (( sbounds2[0] >= arsbounds1[:,0] ) & ( sbounds2[0] <= arsbounds1[:,1] )) | \
             (( sbounds2[1] >= arsbounds1[:,0] ) & ( sbounds2[1] <= arsbounds1[:,1] )) )


def ovlFraction( corngA, corngB ):
    corngIntAB=rngIntersect(corngA,corngB)
    if corngIntAB==None: return 0.
    lcorngUnionAB=unionBounds([corngA],[corngB])
    assert len(lcorngUnionAB)==1
    return old_div(float(spanInclAllowNone( corngIntAB )),float(spanInclAllowNone( lcorngUnionAB[0] )))        
    
def pad_centered_corng( targLen, curStart, curEnd ):
    curLen = curEnd-curStart+1
    if curLen >= targLen: return (curStart,curEnd)
    padEach = old_div((targLen-curLen),2)
    padL = min(curStart,padEach)
    padR = (targLen - curLen) - padL
    return (max(curStart-padL,0), curEnd+padR)
    
def rngIntersect(bounds1, bounds2): 
    if bounds1==None or bounds2==None: return None
    elif rgOverlap(bounds1,bounds2) : return tuple(sorted(list(bounds1)+list(bounds2))[1:3])
    else: return None
    
def rgSubtract(rgA,rgB):
    rgi=rngIntersect(rgA,rgB)
    if rgi==None:
        return [rgA]
    else:
        lRet = [(rgA[0], rgi[0]-1),(rgi[1]+1, rgA[1])]
        lRet = [rng for rng in lRet if rng[0]<=rng[1] ]                
        return lRet 
    
def rgSubtract53ret53(lco1,lco2,strand):
    lrg = rgSubtract((min(lco1),max(lco1)),
                     (min(lco2),max(lco2)) )
    if strand==1:
        return sorted( [ (min(rng),max(rng)) for rng in lrg ], key=lambda x:x[0] )
    else:
        return sorted( [ (max(rng),min(rng)) for rng in lrg ], key=lambda x:-x[1] )

def unionBounds( srclbA, srclbB, doRetMap=False ):
    mUnionBlkToSrcBlk = [{}, {}]

    lblks = sorted( srclbA+srclbB, key=min )

    i = 0
    while i < len( lblks ) - 1:
        if rgOverlap( (lblks[i][0], lblks[i][1]+1), lblks[i+1] ):
            unionBlks = (min(min(lblks[i]),min(lblks[i+1])),
                         max(max(lblks[i]),max(lblks[i+1])) )

            lblks[i:i+2] = [ unionBlks ]

            continue

        i += 1

    if doRetMap:
        for blk in lblks:
            if not blk in mUnionBlkToSrcBlk[0]: mUnionBlkToSrcBlk[0][blk] = []
            if not blk in mUnionBlkToSrcBlk[1]: mUnionBlkToSrcBlk[1][blk] = []
    
            for iA in range( len(srclbA) ):
                if rgOverlap( srclbA[iA], blk ):
                    mUnionBlkToSrcBlk[0][ blk ].append( srclbA[iA] )
            for iB in range( len(srclbB) ):
                if rgOverlap( srclbB[iB], blk ):
                    mUnionBlkToSrcBlk[1][ blk ].append( srclbB[iB] )\
    
        return (lblks, mUnionBlkToSrcBlk[0], mUnionBlkToSrcBlk[1] )
    else: return lblks

# for an array, return inclusive ranges which are continuously
# equal to a single value (excluding ranges equal to zeroDefVal)    
def arrayToNonZeroContigRanges( A, zeroDefVal=0 ):
    dA=np.diff(A)
    N=A.shape[0]
    nzdA,=np.nonzero(dA)
    liRangesAll=np.c_[ np.r_[0,nzdA+1],
                       np.r_[nzdA,N-1] ]
    iiR = A.take(liRangesAll[:,0]) != zeroDefVal
    return liRangesAll[iiR,:]

#
def parseContigLengths(fnContigLengths):
    mCn=OrderedDict()
    filContigLengths = open(fnContigLengths,'r')
    l=filContigLengths.readline()
    if l[0]=='@':
        while len(l)>0:
            l=l.replace('\n','').split('\t')
            if l[0]=='@SQ':
                m={}
                for kv in l[1:]: m[kv.split(':')[0]]=kv.split(':')[1]
                if 'SN' in m and 'LN' in m:
                    mCn[m['SN']]=int(m['LN'])
            l=filContigLengths.readline()            
    else:
        while len(l)>0:
            l=l.replace('\n','').split('\t')
            mCn[l[0]]=int(l[1])
            l=filContigLengths.readline()
    return mCn


# read through a file and rather than returning one line
# return list of lines grouped by a column
class FileGroupize(object):
    def __init__( self, fnIn, colName=None, colIdx=None, skipN=0 ):
        self.ctr = -1
        self.fil=pgzopen(fnIn,'r')
        for i in range(skipN): 
            self.fil.readline()
            self.ctr+=1
        if colIdx!=None:
            self.col=colIdx
        else:
            self.hdrInOrder=self.fil.readline().rstrip('\n')
            self.hdr=parseHdrToMap( self.hdrInOrder )
            self.hdrInOrder=self.hdrInOrder.split('\t')
            self.ctr+=1
            self.col=self.hdr[colName]
        self.lastK,self.lastL=None,None
        
    def __iter__(self): return self
        
    def __next__(self):
        l=self.fil.readline().rstrip('\n')
        self.ctr+=1
        if len(l)==0: 
            if self.lastK!=None:
                self.lastK,ret,self.lastL=None,[self.lastL],None
                return ret
            else:
                raise StopIteration
        
        if self.lastK!=None:
            ll=[self.lastL]
            curK=self.lastK
        else:
            ll=[]
            curK=l.rstrip('\n').split('\t')[self.col]
        
        while len(l)>0 and l.rstrip('\n').split('\t')[self.col]==curK:
            ll.append(l)
            l=self.fil.readline().rstrip('\n')
            self.ctr+=1
        
        if len(l)>0:
            self.lastK,self.lastL=l.rstrip('\n').split('\t')[self.col],l
        else:
            self.lastK,self.lastL=None,None
        
        return ll
    
def vcfGetHdr(fn):
    if fn.endswith('gz'):
        fil=gzip.GzipFile(fn,'r')
    else:
        fil=open(fn,'r')
    for l in fil:
        l=l.rstrip()
        if l[0]=='#' and not l.startswith('##'):
            return l[1:]
        elif l[0]!='#':return None

def openSamAutoMode(fn,mode):
    return '%sb'%mode if fn.endswith('.bam') else mode  

def zopen(fn,mode):
    return gzip.GzipFile(fn,mode) if fn.endswith('gz') else open(fn,mode)

def pgzopenAndProc(fn,mode,nThreads=8):
    if mode=='w':
        return pgz_openwrite_and_proc(fn,nThreads=nThreads) if fn.endswith('gz') else (open(fn,mode),None)
    elif mode=='r':
        return pgz_openread_and_proc(fn,nThreads=nThreads) if fn.endswith('gz') else (open(fn,mode),None)


def pgzopen(fn,mode,nThreads=1):
    if nThreads>1:
        if mode=='w':
            return pgz_openwrite(fn,nThreads=nThreads) if fn.endswith('gz') else open(fn,mode)
        elif mode=='r':
            return pgz_openread(fn,nThreads=nThreads) if fn.endswith('gz') else open(fn,mode)
    else:
        return gzip.GzipFile( fn, mode ) if fn.endswith('gz') else open(fn,mode)

def pgz_openwrite(fn, nThreads=8):
    proc_bufsize=2**18
    #pgbufsize= (proc_bufsize / nThreads)/2**10
    filOut=open(fn,'w')
    prcPigz=subprocess.Popen(['pigz','-p','%d'%nThreads],#,'-b','%d'%pgbufsize],
                             stdin=subprocess.PIPE, stdout=filOut, stderr=sys.stderr, close_fds=True, bufsize=proc_bufsize ) # 16 Mbp
    return prcPigz.stdin

def pgz_openread(fn, nThreads=8):
    proc_bufsize=2**18
#    filIn=open(fn,'r')
    prcPigz=subprocess.Popen(['pigz','-c','-d','-p','%d'%nThreads,fn],#,'-b','%d'%pgbufsize],
                             stdout=subprocess.PIPE, stderr=sys.stderr, close_fds=True, bufsize=proc_bufsize ) # 16 Mbp
    return prcPigz.stdout

def pgz_openread_and_proc(fn, nThreads=8):
    proc_bufsize=2**18
#    filIn=open(fn,'r')
    prcPigz=subprocess.Popen(['pigz','-c','-d','-p','%d'%nThreads,fn],#,'-b','%d'%pgbufsize],
                             stdout=subprocess.PIPE, stderr=sys.stderr, close_fds=True, bufsize=proc_bufsize ) # 16 Mbp
    return prcPigz.stdout,prcPigz


def pgz_openwrite_and_proc(fn, nThreads=8):
    proc_bufsize=2**18
    #pgbufsize= (proc_bufsize / nThreads)/2**10
    filOut=open(fn,'w')
    prcPigz=subprocess.Popen(['pigz','-p','%d'%nThreads],#,'-b','%d'%pgbufsize],
                             stdin=subprocess.PIPE, stdout=filOut, stderr=sys.stderr, close_fds=True, bufsize=proc_bufsize ) # 16 Mbp
    return prcPigz.stdin,prcPigz

    
# eg gte3 --> lambda x:x>=3
def comparisonToLambda(comparstr,dtype=int):
    if comparstr.startswith('gte'):
        val=dtype(comparstr[3:])
        return lambda x,_val=val:x>=_val
    elif comparstr.startswith('lte'):
        val=dtype(comparstr[3:])
        return lambda x,_val=val:x<=_val
    elif comparstr.startswith('eq'):   
        val=dtype(comparstr[2:])
        return lambda x,_val=val:x==_val
    elif comparstr.startswith('gt'):   
        val=dtype(comparstr[2:])
        return lambda x,_val=val:x>_val
    elif comparstr.startswith('lt'):   
        val=dtype(comparstr[2:])
        return lambda x,_val=val:x<_val
    else:
        raise Error('dont understand that expression')


def stripAndKeepExtensions(_fn, lext):
    fn=_fn
    lextRet=[]
    while True:
        if not any( [ fn.endswith(ext) for ext in lext ] ):
            break
        nextExt=[ ext for ext in lext if  fn.endswith(ext)][0]
        fn=fn[:-len(nextExt)]            
        lextRet[0:0]=[nextExt]
    return fn, ''.join(lextRet)

import itertools

def ddply(rows, cols, op=lambda group_rows: group_rows):
    def group_key(row):
        return (row[col] for col in cols)
    rows = sorted(rows, key=group_key)
    return itertools.chain.from_iterable( op(group_rows) for k,group_rows in itertools.groupby(rows, key=group_key))

def flatten( listIn=[], listOut = [] ):
    ''' returns flattened listIn '''
    for element in listIn:
        if getattr(element, '__iter__', False): listOut += flatten( element, [] )
        else: listOut.append( element)
    return listOut

def dictToOrderedDict( din, lcol ):
    return OrderedDict( [ (k,din[k]) for k in lcol ] )




def fasta_bed_fetch_perregion_seq( bed, fasta ):

    import pybedtools as pbt
    import pysam 

    if isinstance( bed, str ):
        btIn = pbt.BedTool(open(bed,'r')).saveas()
    elif isinstance( bed, pbt.bedtool.BedTool ):
        btIn = bed
    else:
        raise ValueError("bed should be a path to a .BED file, or a pybedtools BedTool instance")

    if isinstance( fasta, str ):
        sfIn = pysam.FastaFile( fasta )
    elif isinstance( fasta, pysam.cfaidx.FastaFile ):
        sfIn = fasta
    else:
        raise ValueError("fasta should be a path to a FASTA file, or a pysam FastaFile instance")

    lczOfs = [ iv.start for iv in btIn ]
    lchromName = [ iv.chrom for iv in btIn ]
    lseq = [ sfIn.fetch( iv.chrom, iv.start, iv.end ) for iv in btIn ]

    return lseq 


def fasta_fetch_oneseq( fasta, record_idx=0 ):

    import Bio.SeqIO 

    sfIn = Bio.SeqIO.parse( open(o.inFasta,'r'), 'fasta' )

    for j in range(record_idx):
        sfIn.next()

    seq = sfIn.next()

    return ( seq.name, str(seq.seq) )


######################################################################


def fracGc( seqChrom, corz ):
    seq = seqChrom[corz[0]:corz[1]+1].upper()
    return float( seq.count('C')+seq.count('G') )/float(len(seq))


def copyMetadataAndReplaceTable( store, tblname, tblold, prefix='mippipe__' ):
    attrsOld = dict( [ 
        (k, store.get_storer( tblname ).attrs[k])
        for k in store.get_storer( tblname ).attrs._v_attrnames 
        if k.startswith(prefix) ] )

    store.put( tblname,tblold )

    for k in attrsOld:
        if k not in store.get_storer( tblname ).attrs._v_attrnames :
            store.get_storer( tblname ).attrs[k] = attrsOld[k]

    del attrsOld


def copyMetadata( storeOrig, storeNew, tblnameOrig, tblnameNew, prefix='mippipe__' ):
    attrsOld = dict( [ 
        (k, storeOrig.get_storer( tblnameOrig ).attrs[k])
        for k in storeOrig.get_storer( tblnameOrig ).attrs._v_attrnames 
        if k.startswith(prefix) ] )

    for k in attrsOld:
        if k not in storeNew.get_storer( tblnameNew ).attrs._v_attrnames :
            storeNew.get_storer( tblnameNew ).attrs[k] = attrsOld[k]

    del attrsOld

def intRange(s):
    try:
        lower,upper=map(int, s.strip().split(','))
        return lower,upper
    except:
        raise argparse.ArgumentTypeError('range must be x,y')

    
def floatRange(s):
    try:
        lower,upper=map(float, s.strip().split(','))
        return lower,upper
    except:
        raise argparse.ArgumentTypeError('range must be x,y')


from scipy.special import expit

# def ilogit(p):
#     return np.exp(p) / (1 + np.exp(p)) 


yoonOlEvenF   = 'ggtagcaaagtgcagatgtgctcttc'
yoonOlEvenRrc = 'gaagagttcagagcagtgtgagttca'

yoonOlEvenFshort    = 'caaagtgcagatgtgctcttc'
yoonOlEvenRrcshort = 'gaagagttcagagcagtgtga'

yoonOlOddF    = 'ctatgagcatgttcttcaggctcttc'
yoonOlOddRrc  = 'gaagagcagcaagcgtaattaactgc'

yoonOlOddFshort    = 'agcatgttcttcaggctcttc'
yoonOlOddRrcshort  = 'gaagagcagcaagcgtaatta'

yoonBackbone  = 'GTTGGAGGCTCATCGTTCCCAGATGTTATCGAGGTCCGAC'

mipAmp3F = 'cacctatcgttccgtccaggctcttc'  
mipAmp3Fshort='atcgttccgtccaggctcttc'
mipAmp3Rrc = 'gaagagcacatgttcacggctgctta'
mipAmp3Rrcshort='gaagagcacatgttcacggct'

mipAmp4F = 'aagctactcgacagcagtcgctcttc'
mipAmp4Fshort='actcgacagcagtcgctcttc' 
mipAmp4Rrc = 'gaagagaacgtcgggctgttcagaat'    
mipAmp4Rrcshort='gaagagaacgtcgggctgttc'

mipAmp5F = 'tgagtgatatggttcccgcgctcttc'  
mipAmp5Fshort='gatatggttcccgcgctcttc'
mipAmp5Rrc = 'gaagagatgacgtcgtgggattctct'   
mipAmp5Rrcshort='gaagagatgacgtcgtgggat'


### 2019 MIPS

new2019backbone='AGATCGGAAGAGCACACGTCAGCGTCAGATGTGTATAAGAGACAG'

new2019mipsPrep1f='GAGATCGGCGCGTTAGAAGAC'
new2019mipsPrep1r='TGCAGGATCTAGGGCGAAGAC'
new2019mipsPrep1r_rc=revComp(new2019mipsPrep1r) #GTCTTCGCCCTAGATCCTGCA

new2019mipsPrep2f='GATTGTGGTCGCCTGGAAGAC'
new2019mipsPrep2r='ACTACCGCGTACCTCGAAGAC'
new2019mipsPrep2r_rc=revComp(new2019mipsPrep2r)

new2019mipsPrep3f='GGGCGTAGTTGTTGGGAAGAC'
new2019mipsPrep3r='CTCGCTGGTGTGTCAGAAGAC'
new2019mipsPrep3r_rc=revComp(new2019mipsPrep3r)

new2019mipsPrep4f='CTCGTCAACGACACGGAAGAC'
new2019mipsPrep4r='TTATCTCGTCGCGCAGAAGAC'
new2019mipsPrep4r_rc=revComp(new2019mipsPrep4r)




class ScoringModelBase:

    def __init__(self, fastaFn, fastaDictFn):
        self.chromLast = None
        self.chromSeqLast = None
        self.sfGenome = pysam.Fastafile( fastaFn )    
        self.mCtgLen = parseContigLengths( fastaDictFn )

    def update_chrom( self, chrom ):
        if chrom != self.chromLast:
            self.chromLast = chrom
            self.chromSeqLast = self.sfGenome.fetch( chrom.encode('ascii') )



    # larmRecords: list of arms for [top/ext, top/lig, btm/ext, btm/lig]
    # normalize length - make all designed sequences this long
    # tag len - randomer tag length
    # maxLigArmLen 
    # maxExtArmLen 
    def makeseq_dsmip_bbsi(self,
                           larmRecords,
                           lhitRecords,
                           lsnpRecords,
                           normalizeLength=160,
                           tagLen=10,
                           maxLigArmLen=30,
                           disallowIn3pNbases=8,
                           scorePenaltyToMutate=0.25,
                           adaptorLeft=yoonOlEvenF,
                           adaptorRight=yoonOlEvenRrc ):

        self.update_chrom( larmRecords[0].chrom )

        lenExtTop = larmRecords[0].end - larmRecords[0].start + 1 + 6
        lenLigTop = larmRecords[1].end - larmRecords[1].start + 1 + 2 
        
        # assert lenLigTop <= maxLigArmLen
        if lenLigTop > maxLigArmLen:
            return (True, -1, '')

        # for BbsI, pad with the in-between two bases equaling the genome at this site
        # add LIG arm
        seqLigArm = self.chromSeqLast[ larmRecords[1].start-2:larmRecords[1].start ]
        seqLigArm += self.chromSeqLast[ larmRecords[1].start:larmRecords[1].end+1 ]
            
        # add EXT arm
        seqExtArm = self.chromSeqLast[ larmRecords[0].start:larmRecords[0].end+1 ]
        # for EarI, pad with the in-between *six* bases equaling genome @ this site
        seqExtArm += self.chromSeqLast[ larmRecords[0].end+1:larmRecords[0].end+7 ]

        seqLigArm,seqExtArm=seqLigArm.upper(), seqExtArm.upper()


        penalty=0
        suppress=False

        # for arms w/ EarI sites, mutate them if far enough away
        # from 3' ends and add a penalty, or indicate that they should be suppressed entirely.
        if seqLigArm.count('GAAGAC') > 0:
            if seqLigArm.index('GAAGAC') +1 >= disallowIn3pNbases:
                seqLigArm = seqLigArm.replace('GAAGAC','GAAGAg')
                penalty += scorePenaltyToMutate
            else:
                suppress = True

        if seqLigArm.endswith('GAAGA'):
            # prevent GAAGAn (ie the beginning degen base of tag)
            if len(seqLigArm)-5 >= disallowIn3pNbases:
                seqLigArm = seqLigArm[:-5]+'GAAGt'
                penalty += scorePenaltyToMutate
            else:
                suppress = True

        # on neg strand:
        if seqLigArm.count('GTCTTC') > 0:
            if seqLigArm.index('GTCTTC') +1 >= disallowIn3pNbases:
                seqLigArm = seqLigArm.replace('GTCTTC','GTCTTg')
                penalty += scorePenaltyToMutate
            else:
                suppress = True

        if seqLigArm.endswith('GAAGA'):
            # prevent GAAGAN (ie the beginning degen base of tag)
            if len(seqLigArm)-5 >= disallowIn3pNbases:
                seqLigArm = seqLigArm[:-5]+'GAAGt'
                penalty += scorePenaltyToMutate
            else:
                suppress = True


        if seqExtArm.count('GAAGAC') > 0:
            if seqExtArm.index('GAAGAC') + 5 >= len(seqExtArm) - disallowIn3pNbases:
                seqExtArm = seqExtArm.replace('GAAGAC','GAAGAg')
                penalty += scorePenaltyToMutate
            else:
                suppress = True

        # prev base is G
        if seqExtArm.startswith('AAGAC'):
            # prevent gAAGAC (ie from last base of linker)
            seqExtArm='c'+seqExtArm[1:]
            penalty += scorePenaltyToMutate

        # on neg strand:
        if seqExtArm.count('GTCTTC') > 0:
            if seqExtArm.index('GTCTTC') + 5 >= len(seqExtArm) - disallowIn3pNbases:
                seqExtArm = seqExtArm.replace('GTCTTC','GTCTTg')
                penalty += scorePenaltyToMutate
            else:
                suppress = True


        seqMip = adaptorLeft
        seqMip += seqLigArm
        seqMip += 'N'*tagLen
        seqMip += new2019backbone
        seqMip += seqExtArm
        seqMip += adaptorRight

        # need to pad the "L" side and clip the "R" side to keep the Ns in the same cycles for 
        # synthesis.  Can adjust if need be.
        designedLen = len(seqMip)
        # [lenOuterF] - [<= maxArmLen]

        nPrepad = maxLigArmLen - lenLigTop
        nCrop = max( 0, (designedLen + nPrepad) - normalizeLength )

        seqMipPrepad = ''.join(['ACGT'[random.randint(0,3)] for _ in xrange(nPrepad) ]) + seqMip
        seqMipCrop = seqMipPrepad[ :normalizeLength ]

        # print seqMipCrop
        return (suppress, penalty, seqMipCrop)

    # larmRecords: list of arms for [top/ext, top/lig, btm/ext, btm/lig]
    # normalize length - make all designed sequences this long
    # tag len - randomer tag length
    # maxLigArmLen 
    # maxExtArmLen 
    def makeseq_dsmip_eari(self,
                           larmRecords,
                           lhitRecords,
                           lsnpRecords,
                           normalizeLength=160,
                           tagLen=10,
                           maxLigArmLen=30,
                           disallowIn3pNbases=8,
                           scorePenaltyToMutate=0.25,
                           adaptorLeft=yoonOlEvenF,
                           adaptorRight=yoonOlEvenRrc ):

        self.update_chrom( larmRecords[0].chrom )

        lenExtTop = larmRecords[0].end - larmRecords[0].start + 1 + 3
        lenLigTop = larmRecords[1].end - larmRecords[1].start + 1 + 1 
        
        # assert lenLigTop <= maxLigArmLen
        if lenLigTop > maxLigArmLen:
            return (True, -1, '')

        # for EarI, pad with the in-between base equaling the genome at this site
        # add LIG arm
        seqLigArm = self.chromSeqLast[ larmRecords[1].start-1 ].lower()
        seqLigArm += self.chromSeqLast[ larmRecords[1].start:larmRecords[1].end+1 ]
            
        # add EXT arm
        seqExtArm = self.chromSeqLast[ larmRecords[0].start:larmRecords[0].end+1 ]
        # for EarI, pad with the in-between *three* bases equaling genome @ this site
        seqExtArm += self.chromSeqLast[ larmRecords[0].end+1:larmRecords[0].end+3 ].lower()

        penalty=0
        suppress=False

        # for arms w/ EarI sites, mutate them if far enough away
        # from 3' ends and add a penalty, or indicate that they should be suppressed entirely.
        if seqLigArm.count('CTCTTC') > 0:
            if seqLigArm.index('CTCTTC') +1 >= disallowIn3pNbases:
                seqLigArm = seqLigArm.replace('CTCTTC','CTCTTa')
                penalty += scorePenaltyToMutate
            else:
                suppress = True

        if seqLigArm.endswith('CTCTT'):
            # prevent CTCTTN (ie the beginning degen base of tag)
            if len(seqLigArm)-5 >= disallowIn3pNbases:
                seqLigArm = seqLigArm[:-5]+'CTCTa'
                penalty += scorePenaltyToMutate
            else:
                suppress = True

        # on neg strand:
        if seqLigArm.count('GAAGAG') > 0:
            if seqLigArm.index('GAAGAG') +1 >= disallowIn3pNbases:
                seqLigArm = seqLigArm.replace('GAAGAG','GAAGAc')
                penalty += scorePenaltyToMutate
            else:
                suppress = True

        if seqLigArm.endswith('CTCTT'):
            # prevent GAAGAN (ie the beginning degen base of tag)
            if len(seqLigArm)-5 >= disallowIn3pNbases:
                seqLigArm = seqLigArm[:-5]+'CTCTa'
                penalty += scorePenaltyToMutate
            else:
                suppress = True


        if seqExtArm.count('CTCTTC') > 0:
            if seqExtArm.index('CTCTTC') + 5 >= len(seqExtArm) - disallowIn3pNbases:
                seqExtArm = seqExtArm.replace('CTCTTC','CTCTTa')
                penalty += scorePenaltyToMutate
            else:
                suppress = True

        # prev base is C
        if seqExtArm.startswith('TCTTC'):
            # prevent cTCTTC (ie from last base of linker)
            seqExtArm='g'+seqExtArm[1:]
            penalty += scorePenaltyToMutate

        # on neg strand:
        if seqExtArm.count('GAAGAG') > 0:
            if seqExtArm.index('GAAGAG') + 5 >= len(seqExtArm) - disallowIn3pNbases:
                seqExtArm = seqExtArm.replace('GAAGAG','GAAGAc')
                penalty += scorePenaltyToMutate
            else:
                suppress = True


        seqMip = adaptorLeft
        seqMip += seqLigArm
        seqMip += 'N'*tagLen
        seqMip += yoonBackbone
        seqMip += seqExtArm
        seqMip += adaptorRight

        # need to pad the "L" side and clip the "R" side to keep the Ns in the same cycles for 
        # synthesis.  Can adjust if need be.
        designedLen = len(seqMip)
        # [lenOuterF] - [<= maxArmLen]

        nPrepad = maxLigArmLen - lenLigTop
        nCrop = max( 0, (designedLen + nPrepad) - normalizeLength )

        seqMipPrepad = ''.join(['ACGT'[random.randint(0,3)] for _ in xrange(nPrepad) ]) + seqMip
        seqMipCrop = seqMipPrepad[ :normalizeLength ]

        # print seqMipCrop
        return (suppress, penalty, seqMipCrop)


class ScoringModelV1(ScoringModelBase):

    def __init__( self, fastaFn, fastaDictFn, AFthresh, **kwargs ):

        ScoringModelBase.__init__(self, fastaFn, fastaDictFn)

        self.AFthresh = AFthresh
        self.endSnpPenalty = 0.25

        self.model1_coefs=np.array(
        [[-3.70825266,-0.17230565,-0.14256875,0.05550825,0.03905947,-0.20969247,
           0.00571007, 0.01107943],
         [-4.95934269, 0.28154416, 0.28032762,-0.00660438, 0.44924173, 0.38232691,
          -0.01736603,-0.0126063 ]])

        """
        [1] "(Intercept)"                           
        [2] "arm_tm_ext_diffMean"                   
        [3] "num_exactclosehits_cumul_ext_isgt5TRUE"
        [4] "arm_mean_kmer_freq_extarm_clip"        
        """
        self.model1_coefs_extarm=np.array(
            [[-3.85545654,-0.12251115,-0.22278188, 0.00626786],
             [-3.8672664,  0.32301799, 0.45956199,-0.02436648]] )

        """
        [1] "(Intercept)"                           
        [2] "arm_tm_lig_diffMean"                   
        [3] "num_exactclosehits_cumul_lig_isgt5TRUE"
        [4] "arm_mean_kmer_freq_ligarm_clip"    
        """
        self.model1_coefs_ligarm=np.array(
            [[-3.89274488,-0.15036587, 0.02659756,  0.01116147],
             [-3.97947798, 0.32784149, 0.53942639,-0.02066444]] )

    # larmRecords: list of arms for [top/ext, btm/lig]
    def score_topext_arms(self, 
                          larmRecords,
                          lhitRecords,
                          lsnpRecords ):

        arm_tm_ext_diffMean_top = abs( larmRecords[0].arm_tm - 56.34 )
        arm_tm_lig_diffMean_btm = abs( larmRecords[1].arm_tm - 56.34 )

        num_exactclosehits_gt5_ext_top = 1 if (lhitRecords[0].num_exact_hits + lhitRecords[0].num_close_hits)>5 else 0
        num_exactclosehits_gt5_lig_btm = 1 if (lhitRecords[1].num_exact_hits + lhitRecords[1].num_close_hits)>5 else 0

        arm_mean_kmer_freq_extarm_top = np.clip( larmRecords[0].arm_mean_kmer_freq ,0.,100. )
        arm_mean_kmer_freq_ligarm_btm = np.clip( larmRecords[1].arm_mean_kmer_freq ,0.,100. )

        x_top = np.array( [1., arm_tm_ext_diffMean_top,
                           num_exactclosehits_gt5_ext_top,
                           arm_mean_kmer_freq_extarm_top ])

        y_top_toohigh = np.dot(x_top,self.model1_coefs_extarm[0,:])
        y_top_toolow = np.dot(x_top,self.model1_coefs_extarm[1,:])

        x_btm = np.array( [1., arm_tm_lig_diffMean_btm, 
                           num_exactclosehits_gt5_lig_btm,
                           arm_mean_kmer_freq_ligarm_btm ] )

        y_btm_toohigh = np.dot(x_btm,self.model1_coefs_ligarm[0,:])
        y_btm_toolow = np.dot(x_btm,self.model1_coefs_ligarm[1,:])

        p_topbot_highlow = np.array( 
            [ [ expit(y_top_toohigh),expit(y_top_toolow) ],
              [ expit(y_btm_toohigh),expit(y_btm_toolow) ] ] )

        # penalize strongly end-SNPs
        #  top ext,   ------->*
        if lsnpRecords[0].af_most_common_var_3p_5base >= self.AFthresh:
            p_topbot_highlow[0,1] += self.endSnpPenalty

        #  bottom lig,   <-------*
        if lsnpRecords[1].af_most_common_var_3p_5base >= self.AFthresh:
            p_topbot_highlow[1,1] += self.endSnpPenalty

        score = max( p_topbot_highlow[0,1], p_topbot_highlow[1,1] )

        return score

    # larmRecords: list of arms for [top/lig, btm/ext]
    def score_toplig_arms( self,
                           larmRecords,
                           lhitRecords,
                           lsnpRecords ):

        arm_tm_lig_diffMean_top = abs( larmRecords[0].arm_tm - 56.34 )
        arm_tm_ext_diffMean_btm = abs( larmRecords[1].arm_tm - 56.34 )

        num_exactclosehits_gt5_lig_top = 1 if (lhitRecords[0].num_exact_hits + lhitRecords[0].num_close_hits)>5 else 0
        num_exactclosehits_gt5_ext_btm = 1 if (lhitRecords[1].num_exact_hits + lhitRecords[1].num_close_hits)>5 else 0

        arm_mean_kmer_freq_ligarm_top = np.clip( larmRecords[0].arm_mean_kmer_freq ,0.,100. )
        arm_mean_kmer_freq_extarm_btm = np.clip( larmRecords[1].arm_mean_kmer_freq ,0.,100. )

        x_top = np.array( [1., arm_tm_lig_diffMean_top,
                           num_exactclosehits_gt5_lig_top,
                           arm_mean_kmer_freq_ligarm_top ])

        y_top_toohigh = np.dot(x_top,self.model1_coefs_ligarm[0,:])
        y_top_toolow = np.dot(x_top,self.model1_coefs_ligarm[1,:])

        x_btm = np.array( [1., arm_tm_ext_diffMean_btm, 
                           num_exactclosehits_gt5_ext_btm,
                           arm_mean_kmer_freq_extarm_btm ] )

        y_btm_toohigh = np.dot(x_btm,self.model1_coefs_extarm[0,:])
        y_btm_toolow = np.dot(x_btm,self.model1_coefs_extarm[1,:])

        p_topbot_highlow = np.array( 
            [ [ expit(y_top_toohigh),expit(y_top_toolow) ],
              [ expit(y_btm_toohigh),expit(y_btm_toolow) ] ] )

        # penalize strongly end-SNPs

        #  top lig,   *------->
        if lsnpRecords[0].af_most_common_var_5p_5base >= self.AFthresh:
            p_topbot_highlow[0,1] += self.endSnpPenalty

        #  bottom ext,   *<-------
        if lsnpRecords[1].af_most_common_var_5p_5base >= self.AFthresh:
            p_topbot_highlow[1,1] += self.endSnpPenalty

        score = max( p_topbot_highlow[0,1], p_topbot_highlow[1,1] )

        return score


    # larmRecords: list of arms for [top/ext, top/lig, btm/ext, btm/lig]
    def score_pair( self,
                   larmRecords,
                   lhitRecords,
                   lsnpRecords ):

        arm_tm_ext_diffMean_top = abs( larmRecords[0].arm_tm - 56.34 )
        arm_tm_lig_diffMean_top = abs( larmRecords[1].arm_tm - 56.34 )
        arm_tm_ext_diffMean_btm = abs( larmRecords[2].arm_tm - 56.34 )
        arm_tm_lig_diffMean_btm = abs( larmRecords[3].arm_tm - 56.34 )

        arm_abs_tmdiff_top = abs( larmRecords[0].arm_tm - larmRecords[1].arm_tm )
        arm_abs_tmdiff_btm = abs( larmRecords[2].arm_tm - larmRecords[3].arm_tm )

        num_exactclosehits_gt5_lig_top = 1 if (lhitRecords[1].num_exact_hits + lhitRecords[1].num_close_hits)>5 else 0
        num_exactclosehits_gt5_ext_top = 1 if (lhitRecords[0].num_exact_hits + lhitRecords[0].num_close_hits)>5 else 0

        num_exactclosehits_gt5_lig_btm = 1 if (lhitRecords[3].num_exact_hits + lhitRecords[3].num_close_hits)>5 else 0
        num_exactclosehits_gt5_ext_btm = 1 if (lhitRecords[2].num_exact_hits + lhitRecords[2].num_close_hits)>5 else 0

        arm_mean_kmer_freq_ligarm_top = np.clip( larmRecords[1].arm_mean_kmer_freq ,0.,100. )
        arm_mean_kmer_freq_extarm_top = np.clip( larmRecords[0].arm_mean_kmer_freq ,0.,100. )

        arm_mean_kmer_freq_ligarm_btm = np.clip( larmRecords[3].arm_mean_kmer_freq ,0.,100. )
        arm_mean_kmer_freq_extarm_btm = np.clip( larmRecords[2].arm_mean_kmer_freq ,0.,100. )

        x_top = np.array( [1., arm_tm_lig_diffMean_top, arm_tm_ext_diffMean_top,
                           arm_abs_tmdiff_top, num_exactclosehits_gt5_lig_top,
                           num_exactclosehits_gt5_ext_top, arm_mean_kmer_freq_extarm_top,
                           arm_mean_kmer_freq_ligarm_top ] )

        y_top_toohigh = np.dot(x_top,self.model1_coefs[0,:])
        y_top_toolow = np.dot(x_top,self.model1_coefs[1,:])

        x_btm = np.array( [1., arm_tm_lig_diffMean_btm, arm_tm_ext_diffMean_btm,
                           arm_abs_tmdiff_btm, num_exactclosehits_gt5_lig_btm,
                           num_exactclosehits_gt5_ext_btm, arm_mean_kmer_freq_extarm_btm,
                           arm_mean_kmer_freq_ligarm_btm ] )

        y_btm_toohigh = np.dot(x_btm,self.model1_coefs[0,:])
        y_btm_toolow = np.dot(x_btm,self.model1_coefs[1,:])

        p_topbot_highlow = np.array( 
            [ [ expit(y_top_toohigh),expit(y_top_toolow) ],
              [ expit(y_btm_toohigh),expit(y_btm_toolow) ] ] )

        # penalize strongly end-SNPs
        #  top ext,   ------->*
        if lsnpRecords[0].af_most_common_var_3p_5base >= self.AFthresh:
            p_topbot_highlow[0,1] += 5

        #  top lig,   *------->
        if lsnpRecords[1].af_most_common_var_5p_5base >= self.AFthresh:
            p_topbot_highlow[0,1] += 5

        #  bottom ext,   *<-------
        if lsnpRecords[2].af_most_common_var_5p_5base >= self.AFthresh:
            p_topbot_highlow[1,1] += 5

        #  bottom lig,   <-------*
        if lsnpRecords[3].af_most_common_var_3p_5base >= self.AFthresh:
            p_topbot_highlow[1,1] += 5

        score = max( p_topbot_highlow[0,1], p_topbot_highlow[1,1] )

        return score



# score arm pairs using machine learning method
# Rucheng Diao, diaorch@umich.edu
# Winter, 2016

import sklearn

# import numpy as np
import pandas as pd
import pysam

from sklearn.externals import joblib

class ScoringModelLearned( ScoringModelV1 ):

    # mode: use pre-trained model or train on-the-fly
    # train_path: folder path holding both raw table and pickel files
    def __init__(self,
                 fastaFn, 
                 fastaDictFn,
                 AFthresh,   
                 mode = None, 
                 clfName = None, 
                 train_path = None, 
                 set_prefix = None,
                 **kwargs ):

        ScoringModelV1.__init__(self, fastaFn, fastaDictFn, AFthresh, **kwargs)

        assert mode in ('preTrained','onthefly')

        # HARD CODED TRANSFORMATION PARAMETERS
        ## gapfill gc content, z-transfromation parameters, on data v2 and v3
        ##    mu: 0.530444078166
        ## sigma: 0.115187356795
        self.gapfill_gc_mu, self.gapfill_gc_sigma = 0.530444078166, 0.115187356795
        ## arm gc content difference, z-trans params, on data v2 and v3
        ##    mu: 1.13702717852
        ## sigma: 12.0168607408
        self.arm_gc_diff_mu, self.arm_gc_diff_sigma = 1.13702717852, 12.0168607408
        ## variance stablizing transformation parameters, Poisson distribution
        ## modeled by statmodel, Python
        ## arm hit sum, poisson param, sqrt lambda, on data v2 and v3
        self.arm_hit_sqLmbd = 3.0676275468939727
        ## arm nmerfreq sum, poisson param, sqrt lambda, on data v2 and v3
        self.arm_nmerfreq_sqLmbd = 6.5946973963705577

        if mode == 'preTrained':
            # load pre-trained models
            print '=====loading pre-trained models====='
            self.clf = self.get_pretrained_classifier(clfName, train_path, set_prefix)
        elif mode == 'onthefly':
            # load pre-trained models
            print '=====training models====='
            self.clf = self.train_scoring_model_onthefly(clfName, train_path, set_prefix = None)

    def get_pretrained_classifier(
        self,
        clfName = None, 
        train_path = None, 
        set_prefix = None):
        '''
        load pre-trained models from pkl files
        '''
        if clfName == None:
            clfName = 'GradientBoosting'
        # todo - remove this harcoded path
        if train_path == None:
            train_path = '/nfs/kitzman2/diaorch/trainsets/'
        if set_prefix == None:
            set_prefix = 'trainset'
        
        pklFileName = train_path + 'pklFiles/clf' + clfName + set_prefix + '.pkl'
        print 'loading pre-trained set: ' + pklFileName
        clf = joblib.load(pklFileName)
        return clf

    ## the function to train the machine learning model on the fly
    ## using the data identified by train_path and set_prefix
    ## and classifier identified by clfName
    def train_scoring_model_onthefly(
        self,
        clfName = None, 
        train_path = None, 
        set_prefix = None):
        '''
        check, transform and train model on the new training 
        set defined as new background or user specific bg, 
        return a model and write a pickle file for model 
        '''
        # read train data
        trainTblName = train_path + 'table/clf' + clfName + set_prefix + '.tsv'
        trainTbl = pd.read_csv(trainTblName, sep = '\t')
        # define training features and target
        features = ['arm_gc_diff', 'arm_hit','arm_nmerfreq', 'flanking_gc_1kbp', 'gapfill_gc', 'gapfill_len']
        # features = ['arm_gc_diff', 'arm_hit', 'arm_len', 'arm_nmerfreq', 'flanking_gc_1kbp', 'gapfill_gc', 'gapfill_len']
        target = ['flag']
        trainX = np.asarray(trainTbl[features])
        trainy = np.asarray(trainTbl[target].values).ravel()
        # train
        if clfName == None:
            clfName = 'GrandientBoosting'
        # mode of Gradient Boosting
        if clfName == 'GrandientBoosting':
            from sklearn.ensemble import GradientBoostingClassifier
            clf = GradientBoostingClassifier() #default n_estimators = 100
            clf.fit(trainX, trainy)
        #mode of AdaBoost
        elif clfName == 'AdaBoost':
            from sklearn.ensemble import AdaBoostClassifier
            clf = AdaBoostClassifier() #default n_estimators = 50
            clf.fit(trainX, trainy)
        return clf

    def score_topext_arms(self, 
                          larmRecords,
                          lhitRecords,
                          lsnpRecords ):

        self.update_chrom( larmRecords[0].chrom )

        # categorize and transform features
        # print '=====transforming features====='
        tbl = self.categorize_transform_features_single(larmRecords, lhitRecords, lsnpRecords )
        # select model and train set
        # predict
        # print '=====predicting score====='
        scoresTopBtm = 1.0 - expit( 50.0 * self.clf.decision_function(tbl) )
        # scoresTopBtm = self.clf.predict_proba( tbl )[1]
        # penalize snp
        # print '=====penalize SNPs====='
        snpPens = self.penalize_snp([lsnpRecords[0],None,lsnpRecords[1],None] )
        scoresTopBtm[0] += snpPens[0]*self.endSnpPenalty
        scoresTopBtm[1] += snpPens[1]*self.endSnpPenalty
        # select score, choose the worse
        # print '=====select top/btm score====='
        score = max(scoresTopBtm)+snpPens[0]
        # return the score
        return score

    def score_toplig_arms(self, 
                          larmRecords,
                          lhitRecords,
                          lsnpRecords ):

        self.update_chrom( larmRecords[0].chrom )

        # categorize and transform features
        # print '=====transforming features====='
        tbl = self.categorize_transform_features_single(larmRecords, lhitRecords, lsnpRecords )
        # select model and train set
        # predict
        # print '=====predicting score====='
        scoresTopBtm = 1.0 - expit( 50.0 * self.clf.decision_function(tbl) )
        # scoresTopBtm = self.clf.predict_proba( tbl )[1]
        # penalize snp
        # print '=====penalize SNPs====='
        snpPens = self.penalize_snp([None,lsnpRecords[0],None,lsnpRecords[1]] )
        scoresTopBtm[0] += snpPens[0]*self.endSnpPenalty
        scoresTopBtm[1] += snpPens[1]*self.endSnpPenalty
        # select score, choose the worse
        # print '=====select top/btm score====='
        score = max(scoresTopBtm)+snpPens[1]
        # return the score
        return score


    # larmRecords: list of arms for [top/ext, top/lig, btm/ext, btm/lig]
    def score_pair( self, 
                   larmRecords,
                   lhitRecords,
                   lsnpRecords ):

        self.update_chrom( larmRecords[0].chrom )

        # categorize and transform features
        # print '=====transforming features====='
        tbl = self.categorize_transform_features_pair(larmRecords, lhitRecords, lsnpRecords )
        # select model and train set
        # predict
        # print '=====predicting score====='
        scoresTopBtm = 1.0 - expit( 50.0 * self.clf.decision_function(tbl) )
        # scoresTopBtm = self.clf.predict_proba( tbl )[1]
        # penalize snp
        # print '=====penalize SNPs====='
        snpPens = self.penalize_snp(lsnpRecords)
        scoresTopBtm[0] += snpPens[0]*self.endSnpPenalty
        scoresTopBtm[1] += snpPens[1]*self.endSnpPenalty
        # select score, choose the worse
        # print '=====select top/btm score====='
        score = max(scoresTopBtm)
        # return the score
        return score



    def categorize_transform_features_pair( self, larmRecords, lhitRecords, lsnpRecords ):
        '''
        categorize and tranform features, takes inputs as 
        seperated list and returns a merged pd.dataframe
        '''
        # larmRecords: list of arms for [top/ext, top/lig, btm/ext, btm/lig]

        gapfill_len_top = abs(larmRecords[1].start - larmRecords[0].end) + 1
        gapfill_len_btm = abs(larmRecords[2].start - larmRecords[3].end) + 1

        # hard coded genome path
        # genomePath = '/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa'
        # sfRef = pysam.Fastafile(genomePath)
        
        # gap_fill_gc, and flanking_gc_1kbp
        # arm_gc
        flankLen = 1000
        
        cozCenterTop = int(0.5 * (larmRecords[1].start + larmRecords[0].end))
        corzFlankTop = (max(cozCenterTop - flankLen, 0), min(cozCenterTop + flankLen, len(self.chromSeqLast) - 1))
        
        flanking_gc_1kbp_top = fracGc(self.chromSeqLast, corzFlankTop)
        
        gapfill_gc_top = fracGc(self.chromSeqLast, (larmRecords[0].end, larmRecords[1].start))
        # arm_gc_top_ext = fracGc(self.chromSeqLast, (larmRecords[0].start, larmRecords[0].end))
        # arm_gc_top_lig = fracGc(self.chromSeqLast, (larmRecords[1].start, larmRecords[1].end))

        cozCenterBtm = int(0.5 * (larmRecords[3].start + larmRecords[2].end))
        corzFlankBtm = (max(cozCenterBtm - flankLen, 0), min(cozCenterBtm + flankLen, len(self.chromSeqLast) - 1))
        flanking_gc_1kbp_btm = fracGc(self.chromSeqLast, corzFlankBtm)
        
        gapfill_gc_btm = fracGc(self.chromSeqLast, (larmRecords[3].end, larmRecords[2].start))
        arm_gc_btm_ext = fracGc(self.chromSeqLast, (larmRecords[2].start, larmRecords[2].end))
        arm_gc_btm_lig = fracGc(self.chromSeqLast, (larmRecords[3].start, larmRecords[3].end))

        # arm_gc_diff_top = abs(arm_gc_top_ext - arm_gc_top_lig)
        # arm_gc_diff_btm = abs(arm_gc_btm_ext - arm_gc_btm_lig)

        # this does not match how data were scaled to train classifier
        #  there, difference wasn't abs'd
        #  there, values used were PCTs not fracs
        # 
        #  so - don't abs()
        #  and use 'arm_gc' which is percent
        arm_gc_diff_top = larmRecords[0]['arm_gc']-larmRecords[1]['arm_gc']
        arm_gc_diff_btm = larmRecords[2]['arm_gc']-larmRecords[3]['arm_gc']

        arm_len_top = (larmRecords[0].end - larmRecords[0].start + 1) + (larmRecords[1].end - larmRecords[1].start + 1)
        arm_len_btm = (larmRecords[2].end - larmRecords[2].start + 1) + (larmRecords[3].end - larmRecords[3].start + 1)

        # average kmer frequence of ext/lig arm
        # this does not match how data were scaled to train classifier
        # arm_nmerfreq_top = (np.clip(larmRecords[0].arm_mean_kmer_freq, 0., 100.) + np.clip( larmRecords[1].arm_mean_kmer_freq, 0., 100. )) / 2
        # arm_nmerfreq_btm = (np.clip(larmRecords[2].arm_mean_kmer_freq, 0., 100.) + np.clip( larmRecords[3].arm_mean_kmer_freq, 0., 100. )) / 2
        arm_nmerfreq_top = larmRecords[0].arm_mean_kmer_freq+larmRecords[1].arm_mean_kmer_freq
        arm_nmerfreq_btm = larmRecords[2].arm_mean_kmer_freq+larmRecords[3].arm_mean_kmer_freq


        # arm hit, mrfAndBwa
        # this does not match how data were scaled to train classifier
        # arm_hit_top = lhitRecords[0].num_exact_hits + lhitRecords[0].num_close_hits + lhitRecords[1].num_exact_hits + lhitRecords[1].num_close_hits
        # arm_hit_btm = lhitRecords[2].num_exact_hits + lhitRecords[2].num_close_hits + lhitRecords[3].num_exact_hits + lhitRecords[3].num_close_hits



        arm_hit_top = lhitRecords[0].num_all_hits + lhitRecords[1].num_all_hits
        arm_hit_btm = lhitRecords[2].num_all_hits + lhitRecords[3].num_all_hits
        



        def z_trans(x, mu, sigma):
            return (x - mu)/sigma
        def vst_trans(x, paraSqLmbd):
            return 2*(np.sqrt(x) - paraSqLmbd)
        # z-trans
        gapfill_gc_z_top = z_trans(gapfill_gc_top, self.gapfill_gc_mu, self.gapfill_gc_sigma)
        gapfill_gc_z_btm = z_trans(gapfill_gc_btm, self.gapfill_gc_mu, self.gapfill_gc_sigma)
        arm_gc_diff_z_top = z_trans(arm_gc_diff_top, self.arm_gc_diff_mu, self.arm_gc_diff_sigma)
        arm_gc_diff_z_btm = z_trans(arm_gc_diff_btm, self.arm_gc_diff_mu, self.arm_gc_diff_sigma)

        # variance stablizing transformation
        arm_hit_vst_top = vst_trans(arm_hit_top, self.arm_hit_sqLmbd)
        arm_hit_vst_btm = vst_trans(arm_hit_btm, self.arm_hit_sqLmbd)
        arm_nmerfreq_vst_top = vst_trans(arm_nmerfreq_top, self.arm_nmerfreq_sqLmbd)
        arm_nmerfreq_vst_btm = vst_trans(arm_nmerfreq_btm, self.arm_nmerfreq_sqLmbd)

        # features should be in the same order as for training ...
        transTbl = pd.DataFrame(OrderedDict([
                ('arm_gc_diff',[arm_gc_diff_z_top, arm_gc_diff_z_btm]),
                ('arm_hit', [arm_hit_vst_top, arm_hit_vst_btm]),
                ('arm_len', [arm_len_top, arm_len_btm]),
                ('arm_nmerfreq', [arm_nmerfreq_vst_top, arm_nmerfreq_vst_btm]),
                ('flanking_gc_1kbp', [flanking_gc_1kbp_top, flanking_gc_1kbp_btm]),
                ('gapfill_gc', [gapfill_gc_z_top, gapfill_gc_z_btm]),
                ('gapfill_len',[gapfill_len_top, gapfill_len_btm])
            ]),
            index = ['top', 'btm'])
        # pass a data frame on 
        # print transTbl
        return transTbl


    def categorize_transform_features_single( self, larmRecords, lhitRecords, lsnpRecords):
        
        # hard code
        gapfill_len = 100
        flankLen = 1000
        cozCenterTop = int(0.5 * (larmRecords[0].start + larmRecords[0].end))
        corzFlankTop = (max(cozCenterTop - flankLen, 0), min(cozCenterTop + flankLen, len(self.chromSeqLast) - 1))
        
        flanking_gc_1kbp = fracGc(self.chromSeqLast, corzFlankTop)
        
        gapfill_gc = fracGc(self.chromSeqLast, (larmRecords[0].end,larmRecords[0].end+gapfill_len))

        arm_gc_diff = 0

        arm_len_top = (larmRecords[0].end - larmRecords[0].start + 1)
        arm_len_btm = (larmRecords[1].end - larmRecords[1].start + 1)

        # average kmer frequence of ext/lig arm
        arm_nmerfreq_top = larmRecords[0].arm_mean_kmer_freq
        arm_nmerfreq_btm = larmRecords[1].arm_mean_kmer_freq

        arm_hit_top = lhitRecords[0].num_all_hits
        arm_hit_btm = lhitRecords[1].num_all_hits
        
        def z_trans(x, mu, sigma):
            return (x - mu)/sigma
        def vst_trans(x, paraSqLmbd):
            return 2*(np.sqrt(x) - paraSqLmbd)
        # z-trans
        gapfill_gc_z_top = z_trans(gapfill_len, self.gapfill_gc_mu, self.gapfill_gc_sigma)
        gapfill_gc_z_btm = z_trans(gapfill_len, self.gapfill_gc_mu, self.gapfill_gc_sigma)
        arm_gc_diff_z_top = z_trans(arm_gc_diff, self.arm_gc_diff_mu, self.arm_gc_diff_sigma)
        arm_gc_diff_z_btm = z_trans(arm_gc_diff, self.arm_gc_diff_mu, self.arm_gc_diff_sigma)

        # variance stablizing transformation
        arm_hit_vst_top = vst_trans(arm_hit_top, self.arm_hit_sqLmbd)
        arm_hit_vst_btm = vst_trans(arm_hit_btm, self.arm_hit_sqLmbd)
        arm_nmerfreq_vst_top = vst_trans(arm_nmerfreq_top, self.arm_nmerfreq_sqLmbd)
        arm_nmerfreq_vst_btm = vst_trans(arm_nmerfreq_btm, self.arm_nmerfreq_sqLmbd)

        # features should be in the same order as for training ...
        transTbl = pd.DataFrame(OrderedDict([
                ('arm_gc_diff',[arm_gc_diff_z_top, arm_gc_diff_z_btm]),
                ('arm_hit', [arm_hit_vst_top, arm_hit_vst_btm]),
                ('arm_len', [arm_len_top, arm_len_btm]),
                ('arm_nmerfreq', [arm_nmerfreq_vst_top, arm_nmerfreq_vst_btm]),
                ('flanking_gc_1kbp', [flanking_gc_1kbp, flanking_gc_1kbp]),
                ('gapfill_gc', [gapfill_gc_z_top, gapfill_gc_z_btm]),
                ('gapfill_len',[gapfill_len, gapfill_len])
            ]),
            index = ['top', 'btm'])
        # pass a data frame on 
        # print transTbl
        return transTbl

    def penalize_snp(self, lsnpRecords):
        
        pen=[0.,0.]

        # penalize strongly end-SNPs
        #  top ext,   ------->*
        if lsnpRecords[0] is not None and lsnpRecords[0].af_most_common_var_3p_5base >= self.AFthresh:
            pen[0] += 1

        #  top lig,   *------->
        if lsnpRecords[1] is not None and lsnpRecords[1].af_most_common_var_5p_5base >= self.AFthresh:
            pen[0] += 1

        #  bottom ext,   *<-------
        if lsnpRecords[2] is not None and lsnpRecords[2].af_most_common_var_5p_5base >= self.AFthresh:
            pen[1] += 1

        #  bottom lig,   <-------*
        if lsnpRecords[3] is not None and lsnpRecords[3].af_most_common_var_3p_5base >= self.AFthresh:
            pen[1] += 1

        return pen



###
# from common.utils.pbtutils

def data_frame_to_pbt_nameidx( df, 
                               col_chrom='chrom', 
                               col_start='start', 
                               col_end='end',
                               col_idx=None,
                               use_index=True,
                               coords='[00]' ):

    import pybedtools as pbt

    if col_idx is None and use_index:
        df_tobt = df[ [ col_chrom, col_start, col_end ] ].copy()
        df_tobt['name'] = df.index
    elif col_idx is None and not use_index:
        df_tobt = df[ [ col_chrom, col_start, col_end ] ].copy()
        df_tobt['name'] = np.arange( df.shape[0] )
    else:
        df_tobt = df[ [ col_chrom, col_start, col_end, col_idx ] ].copy()

    df_tobt.columns = ['chrom','start','end','name']
        
    coord_offset = coordConversion( coords, '[00)' )
    if coord_offset!=(0,0):
        df_tobt['start']+=coord_offset[0]
        df_tobt['end']+=coord_offset[1]

    outbt = pbt.BedTool.from_dataframe( df_tobt )

    return outbt


# TODO add optional param to this and data_frame_to_pbt_nameidx to keep strand info
def bed_intersect_dataframe( bed,
                             df_table, 
                             col_chrom='chrom', 
                             col_start='start', 
                             col_end='end',
                             col_idx=None,
                             use_index=True,
                             coords='[00]',
                             *intersect_args ):
  
    bt_tbl = data_frame_to_pbt_nameidx( df_table, col_chrom, col_start, col_end, col_idx, use_index, coords )

    bt_tbl_ovl = bt_tbl.intersect( bed, *intersect_args )

    # get indices from bed intersect results
    # cast back to proper type to index this table.
    lidcs = np.unique( np.array( [ iv.fields[3] for iv in bt_tbl_ovl ] ).astype( df_table.index.dtype ) )

    df_tbl_out = df_table.ix[ lidcs ]

    return df_tbl_out
