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

import primer3.bindings as p3b
# import primer3.wrappers as p3w

import pysam
import pybedtools as pbt

from makeamip.capcommon import *

if __name__ == '__main__':
    
    def intRange(s):
        try:
            lower,upper=map(int, s.split(','))
            return lower,upper
        except:
            raise argparse.ArgumentTypeError('range must be x,y')
    
    def floatRange(s):
        try:
            lower,upper=map(float, s.split(','))
            return lower,upper
        except:
            raise argparse.ArgumentTypeError('range must be x,y')
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--genome', dest='genome',
                      help='path to indexed genome file')
    
    opts.add_argument('--genomeDict', dest='genomeDict',
                      help='path to genome contig file dictionary, eg from Picard CreateSequenceDictionary')
    
    opts.add_argument('--targetBed', dest='targetBed',
                      help='bed file of targets, chr start stop name')

    opts.add_argument('--padBy', type=int, dest='padBy',
                      help='pad targets by this amount')

    opts.add_argument('--outTable', dest='outTable')
    
    opts.add_argument('--baseName', default='', dest='baseName')    
    
    opts.add_argument('--ligArmLengthRange', default='20,25', type=intRange, dest='ligArmLengthRange',
                      help='ligation arm length range eg 20,25')

    opts.add_argument('--extArmLengthRange', default='20,25', type=intRange, dest='extArmLengthRange',
                      help='extension arm length range eg 20,25')    
    
    opts.add_argument('--gapFillRange', default='100,120', type=intRange, dest='gapFillRange',
                      help='gap fill range, eg 100,120')    

    #     
    opts.add_argument('--gcRange', default='10,90', type=intRange, dest='gcRange',
                      help='range of arm GC percentages')    

    opts.add_argument('--maxHomopolLen', default=4, type=int, dest='maxHomopolLen',
                      help='maximum length of homopolymer')    



    #     
    opts.add_argument('--generateDsMIPs', action='store_true', default=False, dest='generateDsMIPs')
    opts.add_argument('--enforceTopStrandOnly', action='store_true', default=False, dest='enforceTopStrandOnly')
    opts.add_argument('--dsmipLeftIndent', type=int, default=0, dest='dsmipLeftIndent' )
    opts.add_argument('--dsmipRightIndent', type=int, default=0, dest='dsmipRightIndent' )
        
    o = opts.parse_args()
    
    sfGenome = pysam.Fastafile( o.genome )
    
    mCtgLen = parseContigLengths( o.genomeDict )
    
    bedTargs = pbt.BedTool( o.targetBed ) 

    candsStore = pd.HDFStore( o.outTable, 'w', complib='zlib', complevel=9 )
    
    tblTargets = OrderedDict(
        [('chrom',[iv.chrom for iv in bedTargs]),
         ('start',[iv.start for iv in bedTargs]),
         ('end',[iv.end for iv in bedTargs]),
         ('name',[iv.name for iv in bedTargs])] )
    
    tblTargets=pd.DataFrame(tblTargets)
    tblTargets=tblTargets.sort_values( by= ['chrom','start'] )
    tblTargets=tblTargets.reset_index()
       
    candsStore[ '%stargets'%o.baseName ] = tblTargets

    chromLast,seqChromLast=None,None

    def newCandsTbl():
        tbl = OrderedDict()
        for col in ['targetidx','chrom','start','end','strand','arm_gc']: tbl[col]=[]
        return tbl
    
    def addAddlCols(tbl):
        tbl['arm_tm']=-1.
        tbl['arm_mean_kmer_freq']=0.
        
    lHpsBad=[ b*(o.maxHomopolLen+1) for b in 'ACTG' ]
    def candCheck(seq):
        candGc = 100.*(seqCand.count('G')+seqCand.count('C')) / float(len(seqCand))
        numN = seqCand.count('N')
        hasBadHp = any( [ seq.count(hpb)>0 for hpb in lHpsBad ] )
        return ((not hasBadHp) and candGc>=o.gcRange[0] and candGc<=o.gcRange[1] and numN==0,
                candGc)        
        
    tblCandsExt = newCandsTbl()
    tblCandsLig = newCandsTbl()
    
    if o.generateDsMIPs:
        tblCandsExt_minus = newCandsTbl()
        tblCandsLig_minus = newCandsTbl()
        
    for itgt in tblTargets.index:

        if itgt%100==0 : print 'target:%d..'%itgt

        if chromLast != tblTargets.chrom[itgt]:
            chromLast = tblTargets.chrom[itgt]
            seqChromLast = sfGenome.fetch( str( chromLast) )
        
        # generate MIPs on top-strand
        
        # pad targets by indicated amount or to minimum gap fill size


        
        # padBy = o.padBy if tblTargets.end[itgt]-tblTargets.start[itgt]+1 >= o.gapFillRange[0] \
        #         else ( tblTargets.end[itgt]-tblTargets.start[itgt]+1 - o.gapFillRange[0] + 1 ) / 2
                
        corzTargAndPad = ( max(0, tblTargets.start[itgt] - o.padBy ),
                           min(tblTargets.end[itgt] + o.padBy, mCtgLen[tblTargets.chrom[itgt]]-1 ) )
            
        
        corzExtArmStarts = ( corzTargAndPad[0] - o.extArmLengthRange[1],
                           corzTargAndPad[1] - o.gapFillRange[0] + 1 )
                
        corzLigArmStarts = ( corzTargAndPad[0] - o.ligArmLengthRange[0] + o.gapFillRange[0] + 1,
                             corzTargAndPad[1] + 1 )
               
        if o.generateDsMIPs:
            for coz in xrange( corzExtArmStarts[0], corzExtArmStarts[1]+1 ):
                for armLen in xrange( o.extArmLengthRange[0], o.extArmLengthRange[1]+1 ):
                    seqCand = seqChromLast[ coz:coz+armLen ].upper()
                    
                    # now have to check the resulting ligation arm
                    # cozLig = ( coz + o.dsmipLeftIndent, coz + armLen - 1 )
                    
                    # cozLig = ( coz + o.dsmipLeftIndent, coz + armLen - 1 )
                    cozLig = ( coz, coz+armLen - 1 + o.dsmipLeftIndent )
                    if cozLig[1]-cozLig[0]+1 < o.ligArmLengthRange[0] or \
                       cozLig[1]-cozLig[0]+1 > o.ligArmLengthRange[1]: continue
                    
                    seqCandLig = revComp( seqChromLast[ cozLig[0]:cozLig[1]+1 ].upper() )
                    
                    isOKAndGC=candCheck(seqCand)
                    isOKAndGCLig=(True,candCheck(seqCandLig)[1]) if o.enforceTopStrandOnly else candCheck(seqCandLig)                                                            
                                        
                    if not (isOKAndGC[0] and isOKAndGCLig[0]): continue
                            
                    # todo - should clip to chrom length.
                    # base off ext(+), don't recalc since will be virtually the same for lig(-)
                    # corzExtFlanking1kbp = ( max(0, coz-1000), coz+armLen+1000 )
                    # fracGcExtFlanking1kbp = fracGc( seqChromLast, corzExtFlanking1kbp )

                    tblCandsExt['targetidx'].append(itgt)
                    tblCandsExt['chrom'].append( chromLast )
                    tblCandsExt['start'].append( coz )
                    tblCandsExt['end'].append( coz + armLen - 1 )
                    tblCandsExt['strand'].append( 1 )
                    tblCandsExt['arm_gc'].append( isOKAndGC[1] )
                    # tblCandsExt['arm_1kbflank_gc'].append( fracGcExtFlanking1kbp )
            
                    tblCandsLig_minus['targetidx'].append(itgt)
                    tblCandsLig_minus['chrom'].append( chromLast )
                    tblCandsLig_minus['start'].append( cozLig[0] )
                    tblCandsLig_minus['end'].append( cozLig[1] )
                    tblCandsLig_minus['strand'].append( -1 )
                    tblCandsLig_minus['arm_gc'].append( isOKAndGCLig[1] )
                    # tblCandsLig_minus['arm_1kbflank_gc'].append( fracGcExtFlanking1kbp )
            
            for coz in xrange( corzLigArmStarts[0], corzLigArmStarts[1]+1 ):
                for armLen in xrange( o.ligArmLengthRange[0], o.ligArmLengthRange[1]+1 ):
                    
                    seqCand = seqChromLast[ coz:coz+armLen ].upper()

                    # now have to check the resulting extension arm
                    # cozExt = ( coz, coz + armLen - 1 + o.dsmipRightIndent )

                    cozExt = ( coz + o.dsmipRightIndent, coz+armLen-1 )
                    if cozExt[1]-cozExt[0]+1 < o.extArmLengthRange[0] or \
                       cozExt[1]-cozExt[0]+1 > o.extArmLengthRange[1]: continue
                    
                    seqCandExt = revComp( seqChromLast[ cozExt[0]:cozExt[1]+1 ].upper() )
                    
                    isOKAndGC=candCheck(seqCand)
                    isOKAndGCExt=(True,candCheck(seqCandExt)[1]) if o.enforceTopStrandOnly else candCheck(seqCandExt)                                                            
                                        
                    if not (isOKAndGC[0] and isOKAndGCExt[0]): continue

                    # todo - should clip to chrom length.
                    # base off ext(+), don't recalc since will be virtually the same for lig(-)
                    # corzLigFlanking1kbp = ( max(0, coz-1000), coz+armLen+1000 )
                    # fracGcLigFlanking1kbp = fracGc( seqChromLast, corzLigFlanking1kbp )
              
                    tblCandsLig['targetidx'].append(itgt)
                    tblCandsLig['chrom'].append( chromLast )
                    tblCandsLig['start'].append( coz )
                    tblCandsLig['end'].append( coz + armLen - 1 )
                    tblCandsLig['strand'].append( 1 )
                    tblCandsLig['arm_gc'].append( isOKAndGC[1] )
                    # tblCandsLig['arm_1kbflank_gc'].append( fracGcLigFlanking1kbp )
            
                    tblCandsExt_minus['targetidx'].append(itgt)
                    tblCandsExt_minus['chrom'].append( chromLast )
                    tblCandsExt_minus['start'].append( cozExt[0] )
                    tblCandsExt_minus['end'].append( cozExt[1] )
                    tblCandsExt_minus['strand'].append( -1 )
                    tblCandsExt_minus['arm_gc'].append( isOKAndGCExt[1] )
                    # tblCandsExt_minus['arm_1kbflank_gc'].append( fracGcLigFlanking1kbp )
        else:
            for coz in xrange( corzExtArmStarts[0], corzExtArmStarts[1]+1 ):
                for armLen in xrange( o.extArmLengthRange[0], o.extArmLengthRange[1]+1 ):
                    seqCand = seqChromLast[ coz:coz+armLen ].upper()
                    isOKAndGC=candCheck(seqCand)
                    if not isOKAndGC[0]:
                        continue      

                    # todo - should clip to chrom length.
                    # base off ext(+), don't recalc since will be virtually the same for lig(-)
                    # corzExtFlanking1kbp = ( max(0, coz-1000), coz+armLen+1000 )
                    # fracGcExtFlanking1kbp = fracGc( seqChromLast, corzExtFlanking1kbp )

                    tblCandsExt['targetidx'].append(itgt)
                    tblCandsExt['chrom'].append( chromLast )
                    tblCandsExt['start'].append( coz )
                    tblCandsExt['end'].append( coz + armLen - 1 )
                    tblCandsExt['strand'].append( 1 )
                    tblCandsExt['arm_gc'].append( isOKAndGC[1] )
                    # tblCandsExt['arm_1kbflank_gc'].append( fracGcExtFlanking1kbp )
                            
            for coz in xrange( corzLigArmStarts[0], corzLigArmStarts[1]+1 ):
                for armLen in xrange( o.ligArmLengthRange[0], o.ligArmLengthRange[1]+1 ):
                    seqCand = seqChromLast[ coz:coz+armLen ].upper()
                    isOKAndGC=candCheck(seqCand)
                    if not isOKAndGC[0]:
                        continue   

                    # todo - should clip to chrom length.
                    # base off ext(+), don't recalc since will be virtually the same for lig(-)
                    # corzLigFlanking1kbp = ( max(0, coz-1000), coz+armLen+1000 )
                    # fracGcLigFlanking1kbp = fracGc( seqChromLast, corzLigFlanking1kbp )

                    tblCandsLig['targetidx'].append(itgt)
                    tblCandsLig['chrom'].append( chromLast )
                    tblCandsLig['start'].append( coz )
                    tblCandsLig['end'].append( coz + armLen - 1 )
                    tblCandsLig['strand'].append( 1 )
                    tblCandsLig['arm_gc'].append( isOKAndGC[1] )
                    # tblCandsLig['arm_1kbflank_gc'].append( fracGcLigFlanking1kbp )

                           
    tblCandsExt = pd.DataFrame( tblCandsExt )
    tblCandsLig = pd.DataFrame( tblCandsLig )    

    addAddlCols(tblCandsExt)
    addAddlCols(tblCandsLig)
    
    candsStore[ '%sextArm'%o.baseName ] = tblCandsExt
    candsStore[ '%sligArm'%o.baseName ] = tblCandsLig
       
    print 'made table %sextArm with %d entries'%(o.baseName,tblCandsExt.shape[0])
    print 'made table %sligArm with %d entries'%(o.baseName,tblCandsLig.shape[0])

    candsStore.get_storer( '%sextArm'%o.baseName ).attrs['mippipe__target_table'] = '%stargets'%o.baseName
    candsStore.get_storer( '%sligArm'%o.baseName ).attrs['mippipe__target_table'] = '%stargets'%o.baseName

    if o.generateDsMIPs:
        tblCandsExt_minus = pd.DataFrame( tblCandsExt_minus )
        tblCandsLig_minus = pd.DataFrame( tblCandsLig_minus )    
    
        addAddlCols(tblCandsExt_minus)
        addAddlCols(tblCandsLig_minus)
        
        candsStore[ '%sligArmMinus'%o.baseName ] = tblCandsLig_minus
        candsStore[ '%sextArmMinus'%o.baseName ] = tblCandsExt_minus
           
        print 'made table %sextArmMinus with %d entries'%(o.baseName,tblCandsExt_minus.shape[0])
        print 'made table %sligArmMinus with %d entries'%(o.baseName,tblCandsLig_minus.shape[0])

        candsStore.get_storer( '%sextArm'%o.baseName ).attrs['mippipe__dsmip_paired_table'] = '%sligArmMinus'%o.baseName
        candsStore.get_storer( '%sligArm'%o.baseName ).attrs['mippipe__dsmip_paired_table'] = '%sextArmMinus'%o.baseName

        candsStore.get_storer( '%sextArmMinus'%o.baseName ).attrs['mippipe__dsmip_paired_table'] = '%sligArm'%o.baseName
        candsStore.get_storer( '%sligArmMinus'%o.baseName ).attrs['mippipe__dsmip_paired_table'] = '%sextArm'%o.baseName

        candsStore.get_storer( '%sextArmMinus'%o.baseName ).attrs['mippipe__target_table'] = '%stargets'%o.baseName
        candsStore.get_storer( '%sligArmMinus'%o.baseName ).attrs['mippipe__target_table'] = '%stargets'%o.baseName

       
    sys.exit(0)
    
        