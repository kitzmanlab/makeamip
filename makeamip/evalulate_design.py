#! /usr/bin/env python

import sys

import os
import os.path
import argparse
from collections import defaultdict, OrderedDict

import numpy as np
import re

import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
import seaborn as sns

import pandas as pd

from common.utils import pbtutils

from makeamip.capcommon import *

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()
   
    opts.add_argument('--inStorePairs', dest='inStorePairs')
    opts.add_argument('--inStorePairsPicked', dest='inStorePairsPicked')

    opts.add_argument('--pairsBaseName', default='', dest='pairsBaseName')

    opts.add_argument('--outPlotBase', dest='outPlotBase')

    opts.add_argument('--pairsOutWithOtherScores', default=None, dest='pairsOutWithOtherScores')

    opts.add_argument('--scoreClipRange', default='0,1', dest='scoreClipRange')

    opts.add_argument('--considerUnpickedFracOvl', default=0.50, type=float, dest='considerUnpickedFracOvl')

    o = opts.parse_args()

    scoreClipRange=( float(o.scoreClipRange.split(',')[0]), 
                     float(o.scoreClipRange.split(',')[1]) )

    pairStoreAll = pd.HDFStore( o.inStorePairs, 'r' )
    tblPairsAll = pairStoreAll[ '%spairsPlus'%o.pairsBaseName ]       

    pairStorePicked = pd.HDFStore( o.inStorePairsPicked, 'r' )
    tblPairsPicked = pairStorePicked[ '%spairsPlus'%o.pairsBaseName ]       

    scoresPairsAll = tblPairsAll['score']
    scoresPairsPicked = tblPairsPicked['score']

    weightsPairsPicked = tblPairsPicked['picked']

    # generate overall histograms of pair scores.

    f,ax=plt.subplots(1,1,figsize=(8,8))

    bins = np.linspace( 
            scoresPairsAll.clip(scoreClipRange[0],scoreClipRange[1]).min(),
            scoresPairsAll.clip(scoreClipRange[0],scoreClipRange[1]).max(),
            100 )

    plt.ylim(0.,1.)
        
    # ax.hist( scoresPairsAll.clip(scoreClipRange[0],scoreClipRange[1]),
    #          bins=bins,
    #          label='all pairs',
    #          normed=True,
    #          histtype='step',
    #          lw=2 )

    # ax.hist( scoresPairsPicked.clip(scoreClipRange[0],scoreClipRange[1]),
    #          bins=bins,
    #          label='picked pairs',
    #          normed=True,
    #          histtype='step',
    #          lw=2 )

    countsAll,binedges=np.histogram( scoresPairsAll.clip(scoreClipRange[0],scoreClipRange[1]), bins=bins, density=False )
    cdfAll = np.cumsum(countsAll/sum(countsAll).astype('f'))
    plt.step(binedges[1:],cdfAll,lw=2,label='all pairs' )
    # print (lPicked[i],(mIxPickedScores[i]>0.5).mean(),(mIxPickedScores[i]).mean(),(mIxPickedScores[i]).max())

    countsPicked,binedges=np.histogram( scoresPairsPicked.clip(scoreClipRange[0],scoreClipRange[1]), weights=weightsPairsPicked, bins=bins, density=False )
    cdfPicked = np.cumsum(countsPicked/sum(countsPicked).astype('f'))
    plt.step(binedges[1:],cdfPicked,lw=2,label='picked pairs')

    countsPicked,binedges=np.histogram( scoresPairsPicked.clip(scoreClipRange[0],scoreClipRange[1]), bins=bins, density=False )
    cdfPicked = np.cumsum(countsPicked/sum(countsPicked).astype('f'))
    plt.step(binedges[1:],cdfPicked,lw=0.5,label='picked pairs (unique probes)')

    plt.xlabel('pair score')

    ax.legend()

    f.savefig('%s_scorehist.png'%o.outPlotBase)

    # intersect each pair w/ all other pairs

    pbtPairsAll = pbtutils.data_frame_to_pbt_nameidx( tblPairsAll )
    pbtPairsPicked = pbtutils.data_frame_to_pbt_nameidx( tblPairsPicked )

    mPbidxLotherScore=defaultdict(list)

    for pairs_isect in pbtPairsPicked.intersect( pbtPairsAll, f=o.considerUnpickedFracOvl, wao=True ):
        
        nbpOverlap = int(pairs_isect.fields[8])
        lenPbPicked = int(pairs_isect.fields[2])-int(pairs_isect.fields[1])
        fracOverlap = float(nbpOverlap)/float(lenPbPicked)

        ixpbPicked,ixpbOther=float(pairs_isect.fields[3]), float(pairs_isect.fields[7])
                
        # check that this other probe was not actually picked.
        if ixpbOther not in tblPairsPicked.index:
    
            rpbOther = tblPairsAll.loc[ ixpbOther ]

            mPbidxLotherScore[ ixpbPicked ].append( rpbOther['score'] )


    nHaveNoOverlappedUnpicked=0

    tblPairsPicked['min_other_score']=None

    for ixpbPicked,rpbPicked in tblPairsPicked.iterrows():

        if len(mPbidxLotherScore[ixpbPicked])==0:
            nHaveNoOverlappedUnpicked+=1
        else:
            cur_score = rpbPicked['score']

            min_other_score = min( mPbidxLotherScore[ixpbPicked] )

            tblPairsPicked.loc[ixpbPicked,'min_other_score'] = min_other_score

    if o.pairsOutWithOtherScores is not None:
        tblPairsPicked.to_csv(o.pairsOutWithOtherScores,sep='\t',index=False)

    # generate histogram of (picked score - other score)

    f,ax=plt.subplots(1,1,figsize=(8,8))

    tblPairsPickedAndOthers = tblPairsPicked[ ~pd.isnull(tblPairsPicked['min_other_score']) ]

    ax.hist(
        tblPairsPickedAndOthers['score'] - tblPairsPickedAndOthers['min_other_score'],
        bins=100,
        normed=True,
        histtype='step',
        lw=2)

    plt.title( 'picked probe score minus min overlapping other probe score (>%.2f overlap)\n%d probes of %d picked total'%(
        o.considerUnpickedFracOvl, 
        tblPairsPickedAndOthers.shape[0],
        tblPairsPicked.shape[0] ))

    f.savefig('%s_pickedLessOvlScoreDiff.png'%o.outPlotBase)
