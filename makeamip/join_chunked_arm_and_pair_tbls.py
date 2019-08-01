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

import Bio.SeqIO

from makeamip.capcommon import *

def remap( lstores, tablename ):
    lnRemap = {'ichunk':[], 'oldidx':[], 'newidx':[]}
    lnRemap = pd.DataFrame( lnRemap )
    lnRemap = lnRemap.astype(int)
    
    itargNext=0
    for ifn in xrange( len(lstores) ):
        tblcurchunkTarg = lstores[ ifn ][ tablename ]
        lnRemapNew = pd.DataFrame( {'newidx':np.arange( itargNext, itargNext+tblcurchunkTarg.shape[0] ) } )
        lnRemapNew['ichunk'] = ifn
        lnRemapNew['oldidx'] = tblcurchunkTarg.index[:]
        lnRemap = pd.concat( [lnRemap, lnRemapNew], sort=False )
        itargNext += tblcurchunkTarg.shape[0]

    lnRemap = lnRemap.sort_values( by= ['ichunk','oldidx'] ).set_index( ['ichunk','oldidx'] )

    return lnRemap

if __name__ == '__main__':
    
    def intRange(s):
        try:
            lower,upper=map(int, s.split(','))
            return lower,upper
        except:
            raise argparse.ArgumentTypeError('range must be x,y')


    opts = argparse.ArgumentParser()
   
    opts.add_argument('--key', dest='key')

    opts.add_argument('--colnameArmTbl', dest='colnameArmTbl')
    opts.add_argument('--colnamePairTbl', dest='colnamePairTbl')
    
    opts.add_argument('--armBasename', default='', dest='armBasename')
    opts.add_argument('--pairBasename', default='', dest='pairBasename')

    opts.add_argument('--hitTblSuffix', dest='hitTblSuffix')
    opts.add_argument('--snpTblSuffix', dest='snpTblSuffix')

    opts.add_argument('--outArmStore', dest='outArmStore')
    opts.add_argument('--outPairStore', dest='outPairStore')

    opts.add_argument('--skipMergeArmTbls', default=False, action='store_true', dest='skipMergeArmTbls')

    o = opts.parse_args()

    tblk = pd.read_csv(o.key,sep='\t')

    if not o.skipMergeArmTbls:
        lArmStore=[ pd.HDFStore( r[o.colnameArmTbl], 'r' ) for _,r in tblk.iterrows() ]

    lPairStore=[ pd.HDFStore( r[o.colnamePairTbl], 'r' ) for _,r in tblk.iterrows() ]

    if not o.skipMergeArmTbls:
        outArmStore = pd.HDFStore( o.outArmStore, 'w', complib='zlib', complevel=9 )

    outPairStore = pd.HDFStore( o.outPairStore, 'a', complib='zlib', complevel=9 )

    # REMAP target index
    
    if not o.skipMergeArmTbls:
        tblRemapTargidx = remap( lArmStore, '%stargets'%o.armBasename )

        # REMAP ARM idcs
        tblRemapExtarmidx = remap( lArmStore, '%sextArm'%o.armBasename )
        tblRemapLigarmidx = remap( lArmStore, '%sligArm'%o.armBasename )
        tblRemapExtarmMinusidx = remap( lArmStore, '%sextArmMinus'%o.armBasename )
        tblRemapLigarmMinusidx = remap( lArmStore, '%sligArmMinus'%o.armBasename )
        
    # REMAP PAIR idcs
    tblRemapPairPlusidx = remap( lPairStore, '%spairsPlus'%o.armBasename )
    tblRemapPairMinusidx = remap( lPairStore, '%spairsMinus'%o.armBasename )

    if not o.skipMergeArmTbls:
        lTargets = [ astore['%stargets'%o.armBasename] for astore in lArmStore ]
        for ichunk in xrange(len(lTargets)):
        
            # deal with chunks with zero length but warn
            if ichunk not in tblRemapTargidx.index:
                print 'warning: chunk %d has zero targets'%(ichunk)
                continue

            newidx = tblRemapTargidx.ix[ ichunk ].ix[ lTargets[ichunk].index ]['newidx']

            lTargets[ichunk]['index'] = newidx
            lTargets[ichunk] = lTargets[ichunk].set_index( 'index', drop=True )
            # del lTargets[ichunk]['index'] # get rid of old index.

        tblJoinedTargets = pd.concat( lTargets )

        outArmStore[ '%stargets'%o.armBasename ] = tblJoinedTargets

        print 'joined to table of %d target intervals '%tblJoinedTargets.shape[0]

    ####### 

    if not o.skipMergeArmTbls:

        lAstoreTblnameRemap = ['%sextArm'%o.armBasename,
                               '%sligArm'%o.armBasename,
                               '%sextArmMinus'%o.armBasename,
                               '%sligArmMinus'%o.armBasename]  

        if o.hitTblSuffix is not None: 
            lAstoreTblnameRemap+=[
                               '%sextArm_%s'%(o.armBasename, o.hitTblSuffix),
                               '%sligArm_%s'%(o.armBasename, o.hitTblSuffix),
                               '%sextArmMinus_%s'%(o.armBasename, o.hitTblSuffix),
                               '%sligArmMinus_%s'%(o.armBasename, o.hitTblSuffix)]

        if o.snpTblSuffix is not None:
            lAstoreTblnameRemap+=[
                               '%sextArm_%s'%(o.armBasename,o.snpTblSuffix),
                               '%sligArm_%s'%(o.armBasename,o.snpTblSuffix),
                               '%sextArmMinus_%s'%(o.armBasename,o.snpTblSuffix),
                               '%sligArmMinus_%s'%(o.armBasename,o.snpTblSuffix)]

        mAstoreTblnameRemapUsing = { '%sextArm'%o.armBasename: tblRemapExtarmidx,
                           '%sligArm'%o.armBasename: tblRemapLigarmidx,
                           '%sextArmMinus'%o.armBasename: tblRemapExtarmMinusidx,
                           '%sligArmMinus'%o.armBasename: tblRemapLigarmMinusidx,
                           '%sextArm_%s'%(o.armBasename, o.hitTblSuffix):tblRemapExtarmidx,
                           '%sligArm_%s'%(o.armBasename, o.hitTblSuffix):tblRemapLigarmidx,
                           '%sextArmMinus_%s'%(o.armBasename, o.hitTblSuffix):tblRemapExtarmMinusidx,
                           '%sligArmMinus_%s'%(o.armBasename, o.hitTblSuffix):tblRemapLigarmMinusidx,
                           '%sextArm_%s'%(o.armBasename,o.snpTblSuffix):tblRemapExtarmidx,
                           '%sligArm_%s'%(o.armBasename,o.snpTblSuffix):tblRemapLigarmidx,
                           '%sextArmMinus_%s'%(o.armBasename,o.snpTblSuffix):tblRemapExtarmMinusidx,
                           '%sligArmMinus_%s'%(o.armBasename,o.snpTblSuffix):tblRemapLigarmMinusidx }

        mAstoreTblnameJoined = {}

        for atblname in lAstoreTblnameRemap:
            lChunks = [ astore[atblname] for astore in lArmStore ]
            for ichunk in xrange(len(lChunks)):

                # deal with chunks with zero length but warn
                if ichunk not in mAstoreTblnameRemapUsing[atblname].index:
                    print 'warning: chunk %d has zero %s'%(ichunk, atblname)
                    continue

                newidx = mAstoreTblnameRemapUsing[atblname].ix[ ichunk ].ix[ lChunks[ichunk].index ]['newidx']
                lChunks[ichunk]['newidx'] = newidx
                lChunks[ichunk] = lChunks[ichunk].set_index( 'newidx', drop=True )
            
            mAstoreTblnameJoined[atblname] = pd.concat( lChunks )

            print 'joined %s --> %d items '%(atblname, mAstoreTblnameJoined[atblname].shape[0])

            outArmStore[ atblname ] = mAstoreTblnameJoined[atblname]

    ####### 

    lPstoreTblnameRemap = ['%spairsPlus'%o.pairBasename,
                           '%spairsMinus'%o.pairBasename ]

    mPstoreTblnameRemapUsing = { '%spairsPlus'%o.pairBasename: tblRemapPairPlusidx,
                       '%spairsMinus'%o.pairBasename: tblRemapPairMinusidx}
    mPstoreTblnameJoined = {}

    ptblname='%spairsPlus'%o.pairBasename
    lChunks = [ pstore[ptblname] for pstore in lPairStore ]
    for ichunk in xrange(len(lChunks)):
        
        # deal with chunks with zero length but warn
        if ichunk not in mPstoreTblnameRemapUsing[ptblname].index:
            print 'warning: chunk %d has zero %s'%(ichunk, ptblname)
            continue

        newidx = mPstoreTblnameRemapUsing[ptblname].ix[ ichunk ].ix[ lChunks[ichunk].index ]['newidx']

        lChunks[ichunk]['newidx'] = newidx

        if not o.skipMergeArmTbls:
            lChunks[ichunk]['targetidx'] = tblRemapTargidx.ix[ ichunk ].ix[ lChunks[ichunk].targetidx ]['newidx'].values
            lChunks[ichunk]['extarmidx'] = tblRemapExtarmidx.ix[ ichunk ].ix[ lChunks[ichunk].extarmidx ]['newidx'].values
            lChunks[ichunk]['ligarmidx'] = tblRemapLigarmidx.ix[ ichunk ].ix[ lChunks[ichunk].ligarmidx ]['newidx'].values

        lChunks[ichunk] = lChunks[ichunk].set_index( 'newidx', drop=True )

    mPstoreTblnameJoined[ptblname] = pd.concat( lChunks )

    print 'joined %s --> %d items '%(ptblname, mPstoreTblnameJoined[ptblname].shape[0])

    outPairStore[ ptblname ] = mPstoreTblnameJoined[ ptblname ]
    outPairStore.get_storer( ptblname  ).attrs['mippipe__ext_and_lig_arm_tables'] = ('%sextArm'%o.armBasename,'%sligArm'%o.armBasename)

    ptblname='%spairsMinus'%o.pairBasename
    lChunks = [ pstore[ptblname] for pstore in lPairStore ]
    for ichunk in xrange(len(lChunks)):

        # deal with chunks with zero length but warn
        if ichunk not in mPstoreTblnameRemapUsing[ptblname].index:
            print 'warning: chunk %d has zero %s'%(ichunk, ptblname)
            continue

        newidx = mPstoreTblnameRemapUsing[ptblname].ix[ ichunk ].ix[ lChunks[ichunk].index ]['newidx']

        lChunks[ichunk]['newidx'] = newidx
        
        if not o.skipMergeArmTbls:
            lChunks[ichunk]['targetidx'] = tblRemapTargidx.ix[ ichunk ].ix[ lChunks[ichunk].targetidx ]['newidx'].values
            lChunks[ichunk]['extarmidx'] = tblRemapExtarmMinusidx.ix[ ichunk ].ix[ lChunks[ichunk].extarmidx ]['newidx'].values
            lChunks[ichunk]['ligarmidx'] = tblRemapLigarmMinusidx.ix[ ichunk ].ix[ lChunks[ichunk].ligarmidx ]['newidx'].values

        lChunks[ichunk] = lChunks[ichunk].set_index( 'newidx', drop=True )


    mPstoreTblnameJoined[ptblname] = pd.concat( lChunks )

    print 'joined %s --> %d items '%(ptblname, mPstoreTblnameJoined[ptblname].shape[0])

    outPairStore[ ptblname ] = mPstoreTblnameJoined[ ptblname ]
    outPairStore.get_storer( ptblname  ).attrs['mippipe__ext_and_lig_arm_tables'] = ('%sextArmMinus'%o.armBasename,'%sligArmMinus'%o.armBasename)
