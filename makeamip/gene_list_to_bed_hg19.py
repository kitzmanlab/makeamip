import sys
import argparse
from collections import defaultdict,OrderedDict

import pandas as pd
import numpy as np

import pybedtools as pbt
from pybedtools.featurefuncs import extend_fields

import itertools

def rgOverlap_array_single(arbounds1, bounds2):

    if arbounds1.shape[0]==0: return np.zeros_like( arbounds1 )

    sbounds2=(min(bounds2),max(bounds2))
    arsbounds1=np.sort(arbounds1,1)

    return ( (( arsbounds1[:,0] >= sbounds2[0] ) & ( arsbounds1[:,0] <= sbounds2[1] )) | \
             (( arsbounds1[:,1] >= sbounds2[0] ) & ( arsbounds1[:,1] <= sbounds2[1] )) | \
             (( sbounds2[0] >= arsbounds1[:,0] ) & ( sbounds2[0] <= arsbounds1[:,1] )) | \
             (( sbounds2[1] >= arsbounds1[:,0] ) & ( sbounds2[1] <= arsbounds1[:,1] )) )

if __name__=='__main__':

    opts = argparse.ArgumentParser()
    
    opts.add_argument('--inGeneList',dest='inGeneList')

    opts.add_argument('--knownGeneTbl',dest='knownGeneTbl')
    opts.add_argument('--refGeneTbl',dest='refGeneTbl')
    opts.add_argument('--gencodeTbl',dest='gencodeTbl')

    opts.add_argument('--outBed',dest='outBed')

    opts.add_argument('--inclExclusivelyNoncoding',action='store_true',dest='inclExlcusivelyNoncoding')

    o = opts.parse_args()

    print 'only taking exons from coding gene models'
    print 'and clipping to cds'

    # was going to use cruzdb module but could not find
    # a way to get gene symbol w/ genes

    # instead download these tables through table browser

    # coordinates from these tables are [0,0)  (0-based, right-open)
    # pybedtools deals in [0,0)  (pythonically)
    # our output bed should be [1,1) (but this should happen automatically)

    tblGeneList = pd.read_csv( o.inGeneList, header=None )
    tblGeneList[0] = [nm.strip() for nm in tblGeneList[0]]

    tblKnownGene = pd.read_csv( o.knownGeneTbl, sep='\t', compression='gzip' )
    tblRefGene = pd.read_csv( o.refGeneTbl, sep='\t', compression='gzip' )
    tblGencode = pd.read_csv( o.gencodeTbl, sep='\t', compression='gzip' )

    tblRefGene = tblRefGene.set_index( 'name2' )
    tblKnownGene = tblKnownGene.set_index( 'hg19.kgXref.geneSymbol' )
    tblGencode = tblGencode.set_index( 'hg19.wgEncodeGencodeAttrsV24lift37.geneName' )

    tblExons = OrderedDict( [(k,[]) for k in ['chrom','start','end','gene'] ])
    tblNcExons = OrderedDict( [(k,[]) for k in  ['chrom','start','end','gene'] ])    
    tblCds = OrderedDict( [(k,[]) for k in  ['chrom','start','end','gene'] ])

    for gn in tblGeneList[0]:

        ctrAddedCdsCurGene  =  0

        tblExonsCur = OrderedDict( [(k,[]) for k in 
            ['chrom','start','end','gene'] ])
        tblNcExonsCur = OrderedDict( [(k,[]) for k in 
            ['chrom','start','end','gene'] ])
        tblCdsCur = OrderedDict( [(k,[]) for k in 
            ['chrom','start','end','gene'] ])

        # table, prefix
        for tblname,tbl,pfx in \
            zip( ['refgene','knowngene','gencode'],
                 [tblRefGene,tblKnownGene,tblGencode],
                 ['','hg19.knownGene.','hg19.wgEncodeGencodeCompV24lift37.'] ):

            # pass [gn] to preven the result from being a 
            # series (if only one hit)

            if gn not in tbl.index:
                print 'no genes named %s in %s'%( gn, tblname )
                continue

            tblMatches = tbl.loc[ [gn] ]

            ctrAddedCdsCurGeneCurTbl  =  0

            if tblMatches.shape[0]>0:
                
                for _,r in tblMatches.iterrows():

                    chrom = r['%schrom'%pfx]

                    isCds = False

                    # only take coding exons.
                    if r['%scdsEnd'%pfx]>r['%scdsStart'%pfx]:

                        ctrAddedCdsCurGene += 1
                        ctrAddedCdsCurGeneCurTbl += 1

                        tblCdsCur['chrom'].append( chrom )
                        tblCdsCur['start'].append( r['%scdsStart'%pfx] )
                        tblCdsCur['end'].append( r['%scdsEnd'%pfx] )
                        tblCdsCur['gene'].append(gn)

                        isCds = True

                    if len(r['%sexonStarts'%pfx]) > 0:
                        # print r['%sexonStarts'%pfx]
                        lexStarts = np.array( [int(x) for x in r['%sexonStarts'%pfx].split(',') if len(x)>0 ], int )
                        lexEnds = np.array( [int(x) for x in r['%sexonEnds'%pfx].split(',') if len(x)>0 ], int )

                        if isCds:
                            lexIsCoding = rgOverlap_array_single( np.c_[ lexStarts,lexEnds ], 
                                                                  (tblCdsCur['start'][-1], tblCdsCur['end'][-1]) )

                            lexStarts,lexEnds = lexStarts.clip(tblCdsCur['start'][-1], tblCdsCur['end'][-1] ),\
                                                lexEnds.clip(tblCdsCur['start'][-1], tblCdsCur['end'][-1] )

                            for iex in range(len(lexStarts)):

                                if lexIsCoding[iex]:
                                    tblExonsCur['chrom'].append( chrom )
                                    tblExonsCur['start'].append( lexStarts[iex] )
                                    tblExonsCur['end'].append( lexEnds[iex] )
                                    tblExonsCur['gene'].append(gn)
                                    assert lexEnds[iex]>=lexStarts[iex]
                        else:
                            for iex in range(len(lexStarts)):
                                tblNcExonsCur['chrom'].append( chrom )
                                tblNcExonsCur['start'].append( lexStarts[iex] )
                                tblNcExonsCur['end'].append( lexEnds[iex] )
                                tblNcExonsCur['gene'].append(gn)
                                assert lexEnds[iex]>=lexStarts[iex]
                            

            if ctrAddedCdsCurGeneCurTbl  == 0 :
                print ' in %s - no coding exons found for %s'%(tblname,gn)


        if ctrAddedCdsCurGene==0:
            if len(tblNcExonsCur['chrom'])>0:
                print 'adding noncodingexons found for %s'%gn
                
                for k in tblNcExons: tblNcExons[k]+=tblNcExonsCur[k]
            else:
                print 'WARNING - nothing added for %s'%gn

        for k in tblCds: tblCds[k]+=tblCdsCur[k]
        for k in tblExons: tblExons[k]+=tblExonsCur[k]


    tblExons = pd.DataFrame( tblExons )
    tblExons['chrom'] = tblExons['chrom'].astype('str')
    tblExons = tblExons.sort_values( by=['chrom','start'] ).reset_index(drop=True)

    # tblCds = pd.DataFrame( tblCds )
    # tblCds['chrom'] = tblCds['chrom'].astype('str')
    # tblCds = tblCds.sort_values( by= ['chrom','start'] ).reset_index(drop=True)

    tblNcExons = pd.DataFrame( tblNcExons )
    tblNcExons['chrom'] = tblNcExons['chrom'].astype('str')
    tblNcExons = tblNcExons.sort_values( by= ['chrom','start'] ).reset_index(drop=True)

    btCodingExons = pbt.BedTool.from_dataframe( tblExons )
    btCodingExonsMrg = btCodingExons.merge().saveas()

    mrgCexToGene = btCodingExonsMrg.intersect( b=btCodingExons, wo=True )
    m_cexintv_lgene = dict([ (tuple(g), ','.join(list(set([i.fields[6] for i in li])))) for  g,li in  itertools.groupby( mrgCexToGene, lambda x:x.fields[:3] ) ])


    def genenameadder(bt):
        for iv in bt:
            feat=extend_fields(iv,4)
            if tuple(iv.fields[:3]) in m_cexintv_lgene:
                feat.name=m_cexintv_lgene[tuple(iv.fields[:3])]
            yield feat

    btCodingExonsMrg_wgn = pbt.BedTool( genenameadder(btCodingExonsMrg) )

    # btCds = pbt.BedTool.from_dataframe( tblCds ).merge().saveas()
    # btCodingExons = btExons.intersect( btCds )

    if o.inclExlcusivelyNoncoding and tblNcExons.shape[0]>0:
        btNcExons = pbt.BedTool.from_dataframe( tblNcExons ).merge().saveas()
        btTarget =  btCodingExons.cat( btNcExons ).saveas().sort().merge()
        btTarget.saveas(o.outBed)
    else:
        btCodingExonsMrg_wgn.saveas(o.outBed)