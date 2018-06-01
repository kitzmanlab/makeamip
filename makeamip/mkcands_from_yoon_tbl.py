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

import pysam
import pybedtools as pbt

from makeamip.capcommon import *

from makeamip.capcommon import *

if __name__ == '__main__':
      
    opts = argparse.ArgumentParser()

    opts.add_argument('--genome', dest='genome',
                      help='path to indexed genome file')
    
    opts.add_argument('--genomeDict', dest='genomeDict',
                      help='path to genome contig file dictionary, eg from Picard CreateSequenceDictionary')
    
    opts.add_argument('--inTable', dest='inTable')
    
    opts.add_argument('--outTable', dest='outTable')
    
    opts.add_argument('--baseName', default='', dest='baseName')    
    
    opts.add_argument('--fxnArmCoordsToMipseq', 
                      default='lambda corzExt,corzLig,corzExtMinus,corzLigMinus,seqChrom:"%sNNN%s"%( seqChrom[corzExt[0]:corzExt[1]+1], seqChrom[corzLig[0]:corzLig[1]+1] )', 
                      help='fxn( coordrange0i extension, ligation, [extminus], [ligminus], chromSeq ) -> mip seq',
                      dest='fxnArmCoordsToMipseq')

    #     
    opts.add_argument('--dsmipLeftIndent', type=int, default=0, dest='dsmipLeftIndent' )
    opts.add_argument('--dsmipRightIndent', type=int, default=0, dest='dsmipRightIndent' )
        
    o = opts.parse_args()

    sfGenome = pysam.Fastafile( o.genome )

    fxnArmCoordsToMipseq =eval(o.fxnArmCoordsToMipseq)

    #'the yoon data are ds only'
    
    mCtgLen = parseContigLengths( o.genomeDict )
        
    candsStore = pd.HDFStore( o.outTable, 'w' )

    tblIn = pd.read_csv(o.inTable,sep='\t')
    tblIn = tblIn.sort_values(by=['chrom','target_start'])
    tblIn = tblIn.set_index( np.arange( tblIn.shape[0] ) )

    # targets are defined as gap-fill intervals, instead of
    # trying to go backwards to the design targets

    tblTargets=pd.DataFrame()
    tblTargets['chrom']=tblIn['chrom']
    tblTargets['start']=tblIn['target_start']
    tblTargets['end']=tblIn['target_end']-1
    tblTargets['name']=tblIn['probeid']
    
    candsStore[ '%stargets'%o.baseName ] = tblTargets
    
    chromLast,seqChromLast=None,None

    def newCandsTbl():
        tbl = OrderedDict()
        for col in ['targetidx','chrom','start','end','strand','arm_gc']: tbl[col]=[]
        return tbl
    
    def addAddlCols(tbl):
        tbl['arm_tm']=-1.
        tbl['arm_mean_kmer_freq']=0.
                
    def countGC(seq):
        candGC = 100.*(seq.count('G')+seq.count('C')) / float(len(seq))
        return candGC                       

    def newPairsTbl():
        tbl = OrderedDict()
        for col in ['targetidx','extarmidx','ligarmidx','chrom','start','end','strand_fill','mip_seq_5to3']:
            tbl[col]=[]
        return tbl


    tblCandsExt = newCandsTbl()
    tblCandsLig = newCandsTbl()

    tblCandsExt_minus = newCandsTbl()
    tblCandsLig_minus = newCandsTbl()
        
    tblPairsPlus, tblPairsMinus = newPairsTbl(), newPairsTbl()

    for iincomingMip in tblIn.index:
        
        itgt = iincomingMip # mip<->target 1:1 now
        
        if chromLast != tblTargets.chrom[itgt]:
            chromLast = tblTargets.chrom[itgt]
            seqChromLast = sfGenome.fetch( str( chromLast ) ).upper()

# browser details extarm            21     1    21    21 100.0%     9   +  133589678 133589698     21
# browser details ligarm            22     1    22    22 100.0%     9   +  133589810 133589831     22
# target_start: 133589698       target_end:133589809
# ---> 
#         [89677...89697][00]                                      [89809...89830][00]
#                                 [133589698 .. 133589809) [00)


        # target_start --> target_end seems to be defined as [00)

        corzExt = ( tblIn.target_start[iincomingMip] - tblIn.extarm_len[iincomingMip],
                    tblIn.target_start[iincomingMip] - 1 )
        
        corzLig = ( tblIn.target_end[iincomingMip], 
                    tblIn.target_end[iincomingMip] + tblIn.ligarm_len[iincomingMip] - 1 )
        
        seqArmExt = seqChromLast[ corzExt[0]:corzExt[1]+1 ]
        seqArmLig = seqChromLast[ corzLig[0]:corzLig[1]+1 ]
        
        assert seqArmExt == tblIn.extarm_seq[ iincomingMip ], \
            ('extarm does not match stated %s!=%s'%(seqArmExt,tblIn.extarm_seq[ iincomingMip ]), \
                repr(tblIn.ix[iincomingMip]))
        assert seqArmLig == tblIn.ligarm_seq[ iincomingMip ], \
            ('ligarm does not match stated %s!=%s'%(seqArmLig,tblIn.ligarm_seq[ iincomingMip ]), \
                repr(tblIn.ix[iincomingMip]))
        
        corzLigMinus = ( corzExt[0] + o.dsmipLeftIndent, corzExt[1] )
        
        corzExtMinus = ( corzLig[0], corzLig[1] + o.dsmipLeftIndent )

        seqArmLigMinus = seqChromLast[ corzLigMinus[0]:corzLigMinus[1]+1 ]
        seqArmExtMinus = seqChromLast[ corzExtMinus[0]:corzExtMinus[1]+1 ]

        fracGcExt,fracGcLig,fracGcLigMinus,fracGcExtMinus = \
            countGC( seqArmExt ), countGC( seqArmLig ),\
            countGC( seqArmLigMinus ), countGC( seqArmExtMinus )

        tblCandsExt['targetidx'].append(itgt)        
        tblCandsExt['chrom'].append( chromLast )
        tblCandsExt['start'].append( corzExt[0] )
        tblCandsExt['end'].append( corzExt[1] ) 
        tblCandsExt['strand'].append( 1 )
        tblCandsExt['arm_gc'].append( fracGcExt )

        tblCandsLig_minus['targetidx'].append(itgt)
        tblCandsLig_minus['chrom'].append( chromLast )
        tblCandsLig_minus['start'].append( corzLigMinus[0] )
        tblCandsLig_minus['end'].append( corzLigMinus[1] )
        tblCandsLig_minus['strand'].append( -1 )
        tblCandsLig_minus['arm_gc'].append( fracGcLigMinus )

        tblCandsExt_minus['targetidx'].append(itgt)        
        tblCandsExt_minus['chrom'].append( chromLast )
        tblCandsExt_minus['start'].append( corzExtMinus[0] )
        tblCandsExt_minus['end'].append( corzExtMinus[1] ) 
        tblCandsExt_minus['strand'].append( -1 )
        tblCandsExt_minus['arm_gc'].append( fracGcExtMinus )

        tblCandsLig['targetidx'].append(itgt)
        tblCandsLig['chrom'].append( chromLast )
        tblCandsLig['start'].append( corzLig[0] )
        tblCandsLig['end'].append( corzLig[1] )
        tblCandsLig['strand'].append( 1 )
        tblCandsLig['arm_gc'].append( fracGcLig )

        # now make the pairings. these are defined for us so only one per mip x strand.
        tblPairsPlus[ 'targetidx' ].append(itgt)
        tblPairsPlus[ 'extarmidx' ].append(iincomingMip)
        tblPairsPlus[ 'ligarmidx' ].append(iincomingMip)
        tblPairsPlus[ 'chrom' ].append(chromLast)
        tblPairsPlus[ 'start' ].append( corzExt[1] + 1 )
        tblPairsPlus[ 'end' ].append( corzLig[0] - 1 )
        tblPairsPlus[ 'strand_fill' ].append( 1 )

        seqCurMipTop = fxnArmCoordsToMipseq( corzExt,corzLig,corzExtMinus,corzLigMinus,seqChromLast )

        tblPairsPlus[ 'mip_seq_5to3' ].append( seqCurMipTop )

        tblPairsMinus[ 'targetidx' ].append(itgt)
        tblPairsMinus[ 'extarmidx' ].append(iincomingMip)
        tblPairsMinus[ 'ligarmidx' ].append(iincomingMip)
        tblPairsMinus[ 'chrom' ].append(chromLast)
        tblPairsMinus[ 'start' ].append( corzLigMinus[1] + 1 )
        tblPairsMinus[ 'end' ].append( corzExtMinus[0] - 1 )
        tblPairsMinus[ 'strand_fill' ].append( -1 )

        # just use the top strand as the sequence.

        tblPairsMinus[ 'mip_seq_5to3' ].append( seqCurMipTop )



    ##### 
                          
    tblCandsExt = pd.DataFrame( tblCandsExt )
    tblCandsLig = pd.DataFrame( tblCandsLig )    

    addAddlCols(tblCandsExt)
    addAddlCols(tblCandsLig)
    
    candsStore[ '%sextArm'%o.baseName ] = tblCandsExt
    candsStore[ '%sligArm'%o.baseName ] = tblCandsLig
       
    print 'made table %sextArm with %d entries'%(o.baseName,tblCandsExt.shape[0])
    print 'made table %sligArm with %d entries'%(o.baseName,tblCandsLig.shape[0])


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

    tblPairsPlus = pd.DataFrame( tblPairsPlus )
    tblPairsMinus = pd.DataFrame( tblPairsMinus )

    candsStore[ '%spairsPlus'%o.baseName ] = tblPairsPlus
    candsStore[ '%spairsMinus'%o.baseName ] = tblPairsMinus

    print 'made table %spairsPlus with %d entries'%(o.baseName,tblPairsPlus.shape[0])
    print 'made table %spairsMinus with %d entries'%(o.baseName,tblPairsMinus.shape[0])

    candsStore.get_storer( '%sextArm'%o.baseName  ).attrs['mippipe__target_table'] = '%stargets'%o.baseName
    candsStore.get_storer( '%sligArm'%o.baseName ).attrs['mippipe__target_table'] = '%stargets'%o.baseName
    candsStore.get_storer( '%sextArmMinus'%o.baseName  ).attrs['mippipe__target_table'] = '%stargets'%o.baseName
    candsStore.get_storer( '%sligArmMinus'%o.baseName ).attrs['mippipe__target_table'] = '%stargets'%o.baseName

    candsStore.get_storer( '%spairsPlus'%o.baseName  ).attrs['mippipe__target_table'] = '%stargets'%o.baseName
    candsStore.get_storer( '%spairsMinus'%o.baseName ).attrs['mippipe__target_table'] = '%stargets'%o.baseName

    candsStore.get_storer( '%spairsPlus'%o.baseName  ).attrs['mippipe__ext_and_lig_arm_tables'] = ('%sextArm'%o.baseName,'%sligArm'%o.baseName)
    candsStore.get_storer( '%spairsMinus'%o.baseName  ).attrs['mippipe__ext_and_lig_arm_tables'] = ('%sextArmMinus'%o.baseName,'%sligArmMinus'%o.baseName)
    
    candsStore.get_storer( '%spairsPlus'%o.baseName  ).attrs['mippipe__dsmip_paired_pair_table'] = '%spairsMinus'%o.baseName
    candsStore.get_storer( '%spairsMinus'%o.baseName  ).attrs['mippipe__dsmip_paired_pair_table'] = '%spairsPlus'%o.baseName

