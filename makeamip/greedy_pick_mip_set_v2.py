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

from common.utils.pbtutils import *

from makeamip.capcommon import *
#     

from makeamip.capcommon import *


# 
#  ref 5 -------> 3
#          3 <--------- 5  query
#            ????  

def spanOverlap_3pends_rev_array_vs_single_fwd(arbounds1, bounds2):

    if arbounds1.shape[0]==0: return np.zeros_like( arbounds1 )

    sbounds2=(min(bounds2),max(bounds2))
    arsbounds1=np.sort(arbounds1,1)

    doesOvl = (arsbounds1[:,0] >= bounds2[0]) & (arsbounds1[:,0] <= bounds2[1])

    spanOvl = (bounds2[1] - arsbounds1[:,0]).clip(0,1e9)

    return (spanOvl * doesOvl).squeeze()

# 
#  query 5 -------> 3
#          3 <--------- 5  ref
#            ????  

def spanOverlap_3pends_fwd_array_vs_single_rev(arbounds1, bounds2):

    if arbounds1.shape[0]==0: return np.zeros_like( arbounds1 )

    sbounds2=(min(bounds2),max(bounds2))
    arsbounds1=np.sort(arbounds1,1)

    doesOvl = (arsbounds1[:,1] >= bounds2[0]) & (arsbounds1[:,1] <= bounds2[1])

    spanOvl = (bounds2[0] - arsbounds1[:,1]).clip(0,1e9)

    return (spanOvl * doesOvl).squeeze()





def spanOverlap_array_single(arbounds1, bounds2):

    if arbounds1.shape[0]==0: return np.zeros_like( arbounds1 )

    sbounds2=(min(bounds2),max(bounds2))
    arsbounds1=np.sort(arbounds1,1)

    arsboundsj=np.sort( np.c_[ arsbounds1, 
                               np.repeat(sbounds2[0], arsbounds1.shape[0]),
                               np.repeat(sbounds2[1], arsbounds1.shape[0]) ], 
                        1 )

    spanovlAny = ( arsboundsj[:,2]-arsboundsj[:,1]+1 ) 

    doesOvl = (  (( arsbounds1[:,0] >= sbounds2[0] ) & ( arsbounds1[:,0] <= sbounds2[1] )) |
            (( arsbounds1[:,1] >= sbounds2[0] ) & ( arsbounds1[:,1] <= sbounds2[1] )) |
             (( sbounds2[0] >= arsbounds1[:,0] ) & ( sbounds2[0] <= arsbounds1[:,1] )) | \
             (( sbounds2[1] >= arsbounds1[:,0] ) & ( sbounds2[1] <= arsbounds1[:,1] )) ).astype('d')

    return (spanovlAny * doesOvl).squeeze()

# returns (array of does-overlap, array of overlap coordinates, and, -1, 0 )
def rngIntersect_array_single(arbounds1, bounds2):

    if arbounds1.shape[0]==0: return np.zeros_like( arbounds1 )

    sbounds2=(min(bounds2),max(bounds2))
    arsbounds1=np.sort(arbounds1,1)

    arsboundsj=np.sort( np.c_[ arsbounds1, 
                               np.repeat(sbounds2[0], arsbounds1.shape[0]),
                               np.repeat(sbounds2[1], arsbounds1.shape[0]) ], 
                        1 )

    # spanovlAny = ( arsboundsj[:,2]-arsboundsj[:,1]+1 ) 
    corngOvl = arsboundsj[:,1:3]

    doesOvl = (  (( arsbounds1[:,0] >= sbounds2[0] ) & ( arsbounds1[:,0] <= sbounds2[1] )) |
            (( arsbounds1[:,1] >= sbounds2[0] ) & ( arsbounds1[:,1] <= sbounds2[1] )) |
             (( sbounds2[0] >= arsbounds1[:,0] ) & ( sbounds2[0] <= arsbounds1[:,1] )) | 
             (( sbounds2[1] >= arsbounds1[:,0] ) & ( sbounds2[1] <= arsbounds1[:,1] )) )

    corngOvl[ ~doesOvl, 0 ] = -1
    corngOvl[ ~doesOvl, 1 ] = 0

    return (doesOvl.squeeze(), corngOvl)


def oneTargetTallyCvg(
    chrom, 
    corzTarget,
    tblPairs,
    fixedValue=False ):

    arCoverage = np.zeros( (corzTarget[1]-corzTarget[0]+1), dtype=np.int32 )

    # only using coverage from the top strand.
    ielig = tblPairs.index[  
               (tblPairs.chrom==chrom) &
               rgOverlap_array_single( np.c_[ tblPairs.start, tblPairs.end ],
                                       (corzTarget[0],corzTarget[1]) ) ]

    tblPairsOvl = tblPairs.ix[ ielig ]

    for _,rp in tblPairsOvl.iterrows():
        corzCvgToAdd = rngIntersect( (rp.start, rp.end), corzTarget )
        if not fixedValue:
            arCoverage[ corzCvgToAdd[0] - corzTarget[0] : corzCvgToAdd[1] - corzTarget[0] + 1 ] += rp.picked
        else:
            arCoverage[ corzCvgToAdd[0] - corzTarget[0] : corzCvgToAdd[1] - corzTarget[0] + 1 ] += 1

    return arCoverage


# have a hard max of arm-arm overlaps (of EA[i] overlapping with LAminus[i+1] )
# and a contribution to the score.

def greedySetMulCov( 
    chrom,     
    corzTarget,
    tblPairsPlus,
    tblPairsMinus,
    tblExtArmsPlus,
    tblLigArmsPlus,
    tblExtArmsMinus,
    tblLigArmsMinus,
    targetCov,
    scorePenaltyforOverlapExisting, 
    scorePenaltyForNonOverlapTarget,
    maxBpDaisychain,
    colPickPriority=None,
    useAlreadyPickedInfo=False ):

    if useAlreadyPickedInfo:
        print 'comparing against existing already-picked probes for coverage purposes'
    
    # count probes already picked at this target against coverage
    if useAlreadyPickedInfo:
        arCoverage = oneTargetTallyCvg( chrom, corzTarget, tblPairsPlus, fixedValue=False )
    else:
        arCoverage = np.zeros( (corzTarget[1]-corzTarget[0]+1), dtype=np.int32 ) 

    # # only using coverage from the top strand.
    # ielig = tblPairsPlus.index[  
    #            (tblPairsPlus.chrom==chrom) &
    #            rgOverlap_array_single( np.c_[ tblPairsPlus.start, tblPairsPlus.end ],
    #                                    (corzTarget[0],corzTarget[1]) ) ]

    # include probes which don't cover this target at all 
    # but which may, if already picked, disqualify existing probes

    # step out +/- 200 bp 
    ielig = tblPairsPlus.index[  
               (tblPairsPlus.chrom==chrom) &
               rgOverlap_array_single( np.c_[ tblPairsPlus.start, tblPairsPlus.end ],
                                       (corzTarget[0] - 300,corzTarget[1] + 300) ) ]


    tblPairsPlusElig = tblPairsPlus.ix[ ielig ].copy()
    tblPairsMinusElig = tblPairsMinus.ix[ ielig ].copy()

    # count regions that CAN be covered by probes that we can pick.    
    arPotentialCoverage = oneTargetTallyCvg( chrom, corzTarget, tblPairsPlusElig, fixedValue=True )

    corzPairsPlusElig = np.c_[ 
        np.array(tblPairsPlusElig.start,dtype=np.int32),
        np.array(tblPairsPlusElig.end,dtype=np.int32) ]

    tblPairsPlusElig_extArms = tblExtArmsPlus.ix[ tblPairsPlusElig.extarmidx ]
    tblPairsMinusElig_extArms = tblExtArmsMinus.ix[ tblPairsMinusElig.extarmidx ]

    corzPairsPlusEligExtArms = np.c_[
        np.array(tblPairsPlusElig_extArms.start,dtype=np.int32),
        np.array(tblPairsPlusElig_extArms.end,dtype=np.int32) ]

    corzPairsMinusEligExtArms = np.c_[
        np.array(tblPairsMinusElig_extArms.start,dtype=np.int32),
        np.array(tblPairsMinusElig_extArms.end,dtype=np.int32) ]


    pairIsDisallowed = pd.Series( index=tblPairsPlusElig.index )
    pairIsDisallowed[:] = False

    penaltyFactor = np.zeros( (tblPairsPlusElig.shape[0],), dtype='float32' )

    inearButNotOverlap = tblPairsPlusElig.index[  
               ((tblPairsPlusElig.chrom==chrom) &
               ~(rgOverlap_array_single( np.c_[ tblPairsPlusElig.start, tblPairsPlusElig.end ],
                                       (corzTarget[0],corzTarget[1])))) ]

    pairIsDisallowed[ inearButNotOverlap ] = True 

    prevPicked = tblPairsPlusElig.query('picked>0')

    if prevPicked.shape[0]>0:
        'screening for overlap vs previously picked'

    for ip, rpchosen in prevPicked.iterrows():

        # compute overlap penalty based on already selected probes, if specified    
        if useAlreadyPickedInfo:
            ioverlapCurrentChosen = tblPairsPlusElig.index[  
                   rgOverlap_array_single( np.c_[ tblPairsPlusElig.start, tblPairsPlusElig.end ],
                                           ( rpchosen.start, rpchosen.end ) ) ]

            for iother in ioverlapCurrentChosen:
                iiother = tblPairsPlusElig.index.get_loc( iother )
                penaltyFactor[ iiother ] = max( penaltyFactor[ iiother ],
                                            ovlFraction( ( rpchosen.start, rpchosen.end ),
                                                          ( tblPairsPlusElig.start[iother], tblPairsPlusElig.end[iother] ) ) )

        # dq other probes based on already selected probes.
        rpminus = tblPairsMinusElig.ix[ip]

        eaplus = tblExtArmsPlus.ix[ rpchosen.extarmidx ]
        eaminus = tblExtArmsMinus.ix[ rpminus.extarmidx ]

        nbpovl_EAprime_of_curEA = \
            spanOverlap_array_single( corzPairsMinusEligExtArms, (eaplus.start,eaplus.end) )
            # spanOverlap_3pends_rev_array_vs_single_fwd( corzPairsMinusEligExtArms, (eaplus.start,eaplus.end) )

        nbpovl_EA_of_curEAprime = \
            spanOverlap_array_single( corzPairsPlusEligExtArms, (eaminus.start,eaminus.end) )
            # spanOverlap_3pends_fwd_array_vs_single_rev( corzPairsPlusEligExtArms, (eaminus.start,eaminus.end) )
            
            

        nbpovl_EAprime_of_curEA = pd.Series( nbpovl_EAprime_of_curEA, index=tblPairsMinusElig.index )
        nbpovl_EA_of_curEAprime = pd.Series( nbpovl_EA_of_curEAprime, index=tblPairsPlusElig.index )

        liDisallow1=nbpovl_EAprime_of_curEA[ nbpovl_EAprime_of_curEA>maxBpDaisychain ].index 
        liDisallow2=nbpovl_EAprime_of_curEA[ nbpovl_EA_of_curEAprime>maxBpDaisychain ].index

        for ipdq in set(list(liDisallow1)+list(liDisallow2)):
            if not pairIsDisallowed.ix[ipdq]:
                corzCvgToRmv = rngIntersect( (tblPairsPlusElig.start[ipdq],
                                              tblPairsPlusElig.end[ipdq]),
                                            corzTarget )
                arPotentialCoverage[ corzCvgToRmv[0] - corzTarget[0] : corzCvgToRmv[1] - corzTarget[0] + 1 ] -= 1
                assert not np.any(arPotentialCoverage<0)

        pairIsDisallowed.ix[ liDisallow1 ] = True 
        pairIsDisallowed.ix[ liDisallow2 ] = True 

        print '%d/%d Hard Filtered'%(pairIsDisallowed.sum(),pairIsDisallowed.shape[0])


    # reset out current picked count to 0.  
    tblPairsPlusElig['picked'] = 0

    # while np.any( tblCoverage.coverage < targetCov ):
    while np.any( (arCoverage < targetCov) & (arPotentialCoverage>0) ):


        #jkv2
        #arCostSel = tblPairsPlusElig.score+tblPairsMinusElig.score

        # things start to break if the p_fail=1 (score=1) are given exp value of 0.
        # so we'll clip the pr to 0.8 and then add another 0.1 if it's above that.
        arCostSel = np.c_[tblPairsPlusElig.score, tblPairsMinusElig.score].max(axis=1).clip(0,0.8)+ \
                    (0.1*(np.c_[tblPairsPlusElig.score, tblPairsMinusElig.score].max(axis=1)>0.8))
        

        # only using coverage based on top...
        # NB - pandas indexing IS inclusive.  
        # arNewlyCovered = np.array( [ (tblCoverage.ix[ pair.start:pair.end ]['coverage'] < targetCov).sum() for _,pair in tblPairsPlusElig.iterrows() ] )


        #<- BUG - NEED TO REPLACE THIS WITH INTERSECT ON ARRAY
        # right now these arent being clipped
        # arNewlyCovered = np.array(
        #     [ (arCoverage[ corzPairsPlusElig[ip,0] - corzTarget[0] :
        #                    corzPairsPlusElig[ip,1] + 1 - corzTarget[0] ]
        #         <targetCov).sum()
        #        for ip in xrange(corzPairsPlusElig.shape[0]) ],
        #     dtype=np.int32 )

    
        arDoesCoverAny,arCorzNewlyCovered = rngIntersect_array_single( corzPairsPlusElig, corzTarget )
        arNewlyCovered = np.array(
                            [ (arCoverage[arCorzNewlyCovered[ii,0] - corzTarget[0] :
                                          arCorzNewlyCovered[ii,1] - corzTarget[0] +1]<targetCov).sum()
                               if arDoesCoverAny[ii] else 0
                               for ii in xrange(corzPairsPlusElig.shape[0] ) ], dtype=np.int32 )

        if (tblPairsPlusElig.ix[(~pairIsDisallowed).values].shape[0] == 0):
            print 'warning - all remaining arms are excluded by daisychain potential'
            break

        if np.all(arNewlyCovered[(~pairIsDisallowed).values]==0):
            print 'WARNING - ALL REMAINING ARMS ARE EXCLUDED BY DAISYCHAIN POTENTIAL OR ARE REDUNDANT'
            break

        
 


        # COULD instead make this not binary but instead partial relative to the 'amount' of additional coverage added vs required

        # dont fully recalculate inside of loop (expensive!)
        # penaltyFactor = np.zeros_like(arNewlyCovered)
        # for ipairAlreadyPicked in tblPairsPlusElig[ tblPairsPlusElig.picked>0 ].index:
        #     corzAlreadyPicked = (tblPairsPlusElig.start[ipairAlreadyPicked],tblPairsPlusElig.end[ipairAlreadyPicked])
        #     for jipair in xrange(tblPairsPlusElig.shape[0]):


        #         penaltyFactor[ jipair ] = max( penaltyFactor[ jipair ], 
        #                                        ovlFraction( corzAlreadyPicked,
        #                                                     (tblPairsPlusElig.iloc[jipair].start,
        #                                                      tblPairsPlusElig.iloc[jipair].end) ) )
        
        # nonOverlapPenaltyFactor = np.array([
        #     1.0 - 
        #         spanInclAllowNone(rngIntersect( (pair.start,pair.end), (corzTarget[0],corzTarget[1]) )) / 
        #         float(pair.end-pair.start+1)  for 
        #     _,pair in tblPairsPlusElig.iterrows() ] )

        nonOverlapPenaltyFactor = \
            1.0 - (spanOverlap_array_single( corzPairsPlusElig, corzTarget ) /
                  ( corzPairsPlusElig[:,1]-corzPairsPlusElig[:,0]+1 ).astype('f'))
   
        nonOverlapPenaltyFactor = nonOverlapPenaltyFactor.squeeze()

        # arCostAdj = arCostSel * (scorePenaltyforOverlapExisting * penaltyFactor) / arNewlyCovered
        #jkv2
        # arCostAdj = (arCostSel + scorePenaltyForNonOverlapTarget * nonOverlapPenaltyFactor + scorePenaltyforOverlapExisting * penaltyFactor) / arNewlyCovered
        # expected yield
        arExpYield = arNewlyCovered * (1.-arCostSel)

        wpenNonovl = (scorePenaltyForNonOverlapTarget * nonOverlapPenaltyFactor)
        wpenOvlexist = (scorePenaltyforOverlapExisting * penaltyFactor)

        arAdjExpYield = (1. - (wpenNonovl + wpenOvlexist).clip(0.,0.95)) * arExpYield
        arAdjExpYield = arAdjExpYield+0.1

        arCostAdj = 1.0 / arAdjExpYield

        arCostAdj = pd.Series( arCostAdj, index=tblPairsPlusElig.index )

        # pd.merge( pd.DataFrame({'costadj':arCostAdj}), pd.DataFrame({'disallow':pairIsDisallowed}), how='inner', left_index=True, right_index=True ).to_csv('./debug.txt',sep='\t')

        assert not np.any( np.isneginf(arCostAdj.ix[~pairIsDisallowed]) )

        if np.all(np.isneginf(arCostAdj.ix[~pairIsDisallowed])|np.isinf(arCostAdj.ix[~pairIsDisallowed])):
            print 'warning - all remaining arms are excluded by daisychain potential'
            break

        if colPickPriority is None:
            imincost = arCostAdj.ix[ ~pairIsDisallowed ].argmin()
        else:
            arCostPriority = arCostAdj.ix[ (~pairIsDisallowed)&(tblPairsPlusElig[colPickPriority] == True) ]
            #if arCostPriority.shape[0]>0:
            if np.any( ~(np.isinf( arCostPriority )|np.isneginf( arCostPriority )) ):
                imincost = arCostPriority.argmin()
            else:
                imincost = arCostAdj.ix[ ~pairIsDisallowed ].argmin()   

        pairplusMinCost = tblPairsPlusElig.ix[imincost]
        pairminusMinCost = tblPairsMinusElig.ix[imincost]

        eaplusMinCost = tblExtArmsPlus.ix[ pairplusMinCost.extarmidx ]
        eaminusMinCost = tblExtArmsMinus.ix[ pairminusMinCost.extarmidx ]

        # now flag any others that are disallowed based upon probe-probe daisychain potential
        
        # have a hard max of arm-arm overlaps (of EA[i] overlapping with EAminus[i+1] )
        # and a contribution to the score.

        # TODO - above, if using already picked info, need to also do this. 
        # TODO - here we need to remove those probes that are disallowed.
        #  AND - we need to see if AFTER removal of those probes, any other probes overlap the region to be covred.

        nbpovl_EAprime_of_curEA = \
            spanOverlap_array_single( corzPairsMinusEligExtArms, (eaplusMinCost.start,eaplusMinCost.end) )
            # spanOverlap_3pends_rev_array_vs_single_fwd( corzPairsMinusEligExtArms, (eaplusMinCost.start,eaplusMinCost.end) )
            

        nbpovl_EA_of_curEAprime = \
            spanOverlap_array_single( corzPairsPlusEligExtArms, (eaminusMinCost.start,eaminusMinCost.end) )
            # spanOverlap_3pends_fwd_array_vs_single_rev( corzPairsPlusEligExtArms, (eaminusMinCost.start,eaminusMinCost.end) )
            

        nbpovl_EAprime_of_curEA = pd.Series( nbpovl_EAprime_of_curEA, index=tblPairsMinusElig.index )
        nbpovl_EA_of_curEAprime = pd.Series( nbpovl_EA_of_curEAprime, index=tblPairsPlusElig.index ) 

        liDisallow1=nbpovl_EAprime_of_curEA[ nbpovl_EAprime_of_curEA>maxBpDaisychain ].index 
        liDisallow2=nbpovl_EAprime_of_curEA[ nbpovl_EA_of_curEAprime>maxBpDaisychain ].index

        for ipdq in set(list(liDisallow1)+list(liDisallow2)):
            if not pairIsDisallowed.ix[ipdq]:
                corzCvgToRmv = rngIntersect( (tblPairsPlusElig.start[ipdq],
                                              tblPairsPlusElig.end[ipdq]),
                                            corzTarget )
                arPotentialCoverage[ corzCvgToRmv[0] - corzTarget[0] : corzCvgToRmv[1] - corzTarget[0] + 1 ] -= 1
                assert not np.any(arPotentialCoverage<0)

        pairIsDisallowed.ix[ liDisallow1 ] = True 
        pairIsDisallowed.ix[ liDisallow2 ] = True 

        print '%d/%d hard filtered'%(pairIsDisallowed.sum(),pairIsDisallowed.shape[0])

        tblPairsPlusElig.ix[imincost,'picked']+=1
        # tblCoverage.ix[ tblPairsPlusElig.start[imincost]:tblPairsPlusElig.end[imincost], 'coverage'] += 1   

        corzCvgToAdd = rngIntersect( (tblPairsPlusElig.start[imincost],
                                      tblPairsPlusElig.end[imincost]),
                                    corzTarget )
        arCoverage[ corzCvgToAdd[0] - corzTarget[0] : corzCvgToAdd[1] - corzTarget[0] + 1 ] += 1

        ioverlapCurrentChosen = tblPairsPlusElig.index[  
               rgOverlap_array_single( np.c_[ tblPairsPlusElig.start, tblPairsPlusElig.end ],
                                       ( tblPairsPlusElig.start[imincost], tblPairsPlusElig.end[imincost] ) ) ]

        for iother in ioverlapCurrentChosen:
            iiother = tblPairsPlusElig.index.get_loc( iother )
            penaltyFactor[ iiother ] = max( penaltyFactor[ iiother ],
                                        ovlFraction( ( tblPairsPlusElig.start[imincost], tblPairsPlusElig.end[imincost] ),
                                                      ( tblPairsPlusElig.start[iother], tblPairsPlusElig.end[iother] ) ) )

        print 'picked %d / %d eligible, coverage mean=%.2f, ncov@tgt=%d/%d    %s %s:%d-%d (%d x)'%( 
            tblPairsPlusElig.picked.sum(),
            tblPairsPlusElig.shape[0],
            arCoverage.mean(),
            (arCoverage>=targetCov).sum(),
            # tblCoverage.coverage.mean(),
            # (tblCoverage.coverage>=targetCov).sum(),
            arCoverage.shape[0],
            str(imincost),
            tblPairsPlusElig.chrom[imincost],
            tblPairsPlusElig.start[imincost],
            tblPairsPlusElig.end[imincost],
            tblPairsPlusElig.ix[imincost,'picked'] )



    return tblPairsPlusElig[ tblPairsPlusElig.picked>0 ].picked

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--targetBed', dest='targetBed',
                      help='bed file of targets, chr start stop name')

    opts.add_argument('--coverage', type=int, dest='coverage',
                      help='desired coverage')

    opts.add_argument('--scorePenaltyForOverlapExisting', 
                      type=float, 
                      default=0.1,
                      dest='scorePenaltyForOverlapExisting',
                      help='x*(max frac overlapping existing single MIP)')

    opts.add_argument('--scorePenaltyForOverlapExistingPass2', 
                      type=float, 
                      default=0.1,
                      dest='scorePenaltyForOverlapExistingPass2',
                      help='x*(max frac overlapping existing single MIP)')

    opts.add_argument('--scorePenaltyForNonOverlapTarget',
                      type=float,
                      default=0.0,
                      dest='scorePenaltyForNonOverlapTarget',
                      help='x*(fraction of gap fill NOT overlapping target)')

    opts.add_argument('--maxBpDaisychain',
                      type=int,
                      default=14,
                      dest='maxBpDaisychain',
                      help='maximum #bp overlap between ext arm of a given probe and minus ext arm of another, to prevent daisychainning')

    opts.add_argument('--inPairStore', dest='inPairStore')
    opts.add_argument('--inArmStore', dest='inArmStore')
    opts.add_argument('--pairsBaseName', default='', dest='pairsBaseName')
   
    opts.add_argument('--colPickPriority', default=None,
                      dest='colPickPriority',
                      help='')

    opts.add_argument('--targetNumProbes', default=12000,
                      type=int,
                      dest='targetNumProbes',
                      help='target number of probes to synthesize')

    opts.add_argument('--useAlreadyPickedInfoForCoverage', 
        default=False,
        action='store_true',
        help='look for a picked column indicating probe has already been selected',
        dest='useAlreadyPickedInfoForCoverage')

    opts.add_argument('--useAlreadyPickedInfoForDaisyChain', 
        default=False,
        action='store_true',
        help='look for a picked column indicating probe has already been selected',
        dest='useAlreadyPickedInfoForDaisyChain')

    opts.add_argument('--outPairStore', dest='outPairStore')

    o = opts.parse_args()

    bedTargs = pbt.BedTool( o.targetBed ).merge()
    nbpTargs = sum([ (iv.end-iv.start+1) for iv in bedTargs ])

    print '%d bases in target'%nbpTargs

    pairStore = pd.HDFStore( o.inPairStore, 'r' )

    armStore = pd.HDFStore( o.inArmStore, 'r' )

    tblnamePairsPlus = '%spairsPlus'%o.pairsBaseName
    tblnamePairsMinus = '%spairsMinus'%o.pairsBaseName

    tblPairsPlus = pairStore[ tblnamePairsPlus ]
    tblPairsMinus = pairStore[ tblnamePairsMinus ]

    (tblnameExtArmPlus,tblnameLigArmPlus)=pairStore.get_storer( tblnamePairsPlus ).attrs['mippipe__ext_and_lig_arm_tables']
    (tblnameExtArmMinus,tblnameLigArmMinus)=pairStore.get_storer( tblnamePairsMinus ).attrs['mippipe__ext_and_lig_arm_tables']

    tblExtArmsPlus = armStore[ tblnameExtArmPlus ]  
    tblLigArmsPlus = armStore[ tblnameLigArmPlus ]     
    tblExtArmsMinus = armStore[ tblnameExtArmMinus ]  
    tblLigArmsMinus = armStore[ tblnameLigArmMinus ]  

    outPairStore = pd.HDFStore( o.outPairStore, 'w', complib='zlib', complevel=9 )

    bedUnionGapfills=[]
    for _,r in tblPairsPlus.sort_values(by=['chrom','start']).iterrows():
        bedUnionGapfills.append( [ r.chrom, r.start, r.end+1 ] )
    # for _,r in tblPairsMinus.sort(['chrom','start']).iterrows():
    #     bedUnionGapfills.append( [ r.chrom, r.start, r.end+1 ] )

    bedUnionGapfills=sorted(bedUnionGapfills)

    bedUnionGapfills = pbt.BedTool(bedUnionGapfills).merge()

    bedTargetable = bedTargs.intersect( bedUnionGapfills )
    nbpTargetable = sum([ (iv.end-iv.start+1) for iv in bedTargetable ])

    assert nbpTargetable > 0

    print '%d targetable bases in probe'%nbpTargetable

    print 'picking params %s %s %s %s'%(str(o.coverage),
                                    str(o.scorePenaltyForOverlapExisting),
                                    str(o.scorePenaltyForNonOverlapTarget),
                                    str(o.maxBpDaisychain))

    if 'picked' not in tblPairsPlus:
        tblPairsPlus['picked']=0
        tblPairsMinus['picked']=0
    else:
        if not (o.useAlreadyPickedInfoForCoverage or o.useAlreadyPickedInfoForDaisyChain):
            tblPairsPlus['picked']=0
            tblPairsMinus['picked']=0
        else:
            print 'NOT overwriting picked counts'

    arPickedThisRun = pd.Series(  index=tblPairsPlus.index, dtype=np.int32 )
    arPickedThisRun[:]=0

    iivtgt=-1
    for ivtgt in bedTargetable:

        iivtgt+=1

        print 'consol target %d/%d: %s'%(iivtgt, len(bedTargetable), repr(ivtgt))

        # if iivtgt<20:
        #     continue

        lPicked = greedySetMulCov( 
            ivtgt.chrom,
            (ivtgt.start,ivtgt.end-1),
            tblPairsPlus,
            tblPairsMinus,
            tblExtArmsPlus,
            tblLigArmsPlus, 
            tblExtArmsMinus,
            tblLigArmsMinus,
            o.coverage,
            o.scorePenaltyForOverlapExisting,
            o.scorePenaltyForNonOverlapTarget,
            o.maxBpDaisychain,
            o.colPickPriority,
            o.useAlreadyPickedInfoForCoverage )

        tblPairsPlus.ix[ lPicked.index, 'picked' ] += lPicked
        tblPairsMinus.ix[ lPicked.index, 'picked' ] += lPicked

        arPickedThisRun[ lPicked.index ] += lPicked

        # if iivtgt==0:
        # # if iivtgt==20:
        #     tblPickedPlus = tblPairsPlus.ix[ lPicked.index ].copy()
        #     tblPickedPlus['picked'] = np.array(lPicked,np.int32)
        #     tblPickedMinus = tblPairsMinus.ix[ lPicked.index ].copy()
        #     tblPickedMinus['picked'] = np.array(lPicked,np.int32)
        # else:
        #     lii = tblPickedPlus.index.intersection( lPicked.index )
        #     tblPickedPlus.ix[ lii,'picked' ] += lPicked.ix[lii]
        #     tblPickedMinus.ix[ lii,'picked' ] += lPicked.ix[lii]
            
        #     lio = lPicked.index.difference( tblPickedPlus.index )

        #     tblPickedPlusNew = tblPairsPlus.ix[lio].copy()
        #     tblPickedPlusNew['picked'] = lPicked.ix[lio]
        #     tblPickedMinusNew = tblPairsMinus.ix[lio].copy()
        #     tblPickedMinusNew['picked'] = lPicked.ix[lio]

        #     tblPickedPlus = pd.concat( [tblPickedPlus, tblPickedPlusNew ] )
        #     tblPickedMinus = pd.concat( [tblPickedMinus, tblPickedMinusNew ] )
        
        # if iivtgt==100:
        #     break

    print 'PASS 1: picked %d probes, %d unique'%( tblPairsPlus.picked.sum(), tblPairsPlus[tblPairsPlus.picked>0].shape[0] )

    # WE CAN TAKE AN ADDITIONAL N probes
    # (1) count gapfill depth of coverage at every target base
    # (2) aggregate per target the number of bases with coverage @ min
    # (3) call greedy pick
    #     if it picked any more probes, update coverage at this target
    #     if not then mark this target as done.

    allTargetsDone=False

    # target index --> array of cvg
    mIivtArcvg = {}

    # target index --> #nbp @ min cvg.
    mIivtNbpAtMinCvg = {}

    # per target, have we exhausted all elig probes here?
    arTargIsDone = np.zeros( (len(bedTargetable),), dtype=bool )

    iivtgt=-1
    for ivtgt in bedTargetable:
        iivtgt+=1
        mIivtArcvg[ iivtgt ] = \
            oneTargetTallyCvg( ivtgt.chrom, (ivtgt.start, ivtgt.end-1 ), tblPairsPlus.query('picked>0') )


    while (not allTargetsDone) and tblPairsPlus.picked.sum() < o.targetNumProbes: 

        # current min cvg
        minBpCvg = 9e9

        iivtgt=-1
        for ivtgt in bedTargetable:
            iivtgt+=1
            if not arTargIsDone[iivtgt]:
                minBpCvg = min( minBpCvg, mIivtArcvg[iivtgt].min() )

        iivtgt=-1
        for ivtgt in bedTargetable:
            iivtgt+=1
            if not arTargIsDone[iivtgt]:
                mIivtNbpAtMinCvg[iivtgt] = (mIivtArcvg[ iivtgt ]<=minBpCvg).sum()

        # which target to pick?

        iivtgtCur = sorted( [ iivtgt  for iivtgt in mIivtNbpAtMinCvg.keys() if not arTargIsDone[iivtgt] ], key=lambda iivtgt:-mIivtNbpAtMinCvg[iivtgt] )[0]
        # dbg:[ (iivtgt,-mIivtNbpAtMinCvg[iivtgt]  )for iivtgt in mIivtNbpAtMinCvg.keys() if not arTargIsDone[iivtgt] ]


        print 'pass 2: target %d, %d bp @ cvg=%d'%( iivtgtCur, mIivtNbpAtMinCvg[iivtgtCur], minBpCvg )

        ivtgt = bedTargetable[iivtgtCur]

        lPicked = greedySetMulCov( 
            ivtgt.chrom,
            (ivtgt.start,ivtgt.end-1),
            tblPairsPlus,
            tblPairsMinus,
            tblExtArmsPlus,
            tblLigArmsPlus, 
            tblExtArmsMinus,
            tblLigArmsMinus,
            minBpCvg+1,  # attempting to improve coverage by 1 at this target
            o.scorePenaltyForOverlapExistingPass2,
            o.scorePenaltyForNonOverlapTarget,
            o.maxBpDaisychain,
            o.colPickPriority,
            True )

        if lPicked.shape[0]==0:
            # nothing picked
            arTargIsDone[ iivtgtCur ]=True

            allTargetsDone = np.all(arTargIsDone)
        else:

            arOldCvg = mIivtArcvg[ iivtgtCur ]

            tblPairsPlusNewlyPicked = tblPairsPlus.ix[ lPicked.index ].copy()
            tblPairsPlusNewlyPicked.ix[ lPicked.index, 'picked' ] = lPicked
 
            # per-bp coverage by newly picked probes
            arOldCvg = mIivtArcvg[ iivtgtCur ]
            arNewPbsCvg = oneTargetTallyCvg( ivtgt.chrom, (ivtgt.start, ivtgt.end-1 ), tblPairsPlusNewlyPicked )
            arNewCvg = arOldCvg + arNewPbsCvg           

            if (arNewCvg<=minBpCvg).sum() >= (arOldCvg<=minBpCvg).sum():
                # no improvement
                # blacklist this target
                arTargIsDone[ iivtgtCur ]=True

                allTargetsDone = np.all(arTargIsDone)
            else:
                print ' --> @ %d cvg, %d bp --> %d bp'%( minBpCvg, (arOldCvg<=minBpCvg).sum(), (arNewCvg<=minBpCvg).sum())

                # lii = tblPickedPlus.index.intersection( lPicked.index )
                # tblPickedPlus.ix[ lii,'picked' ] += lPicked.ix[lii]
                # tblPickedMinus.ix[ lii,'picked' ] += lPicked.ix[lii]
                
                # lio = lPicked.index.difference( tblPickedPlus.index )

                # tblPickedPlusNew = tblPairsPlus.ix[lio].copy()
                # tblPickedPlusNew['picked'] = lPicked.ix[lio]
                # tblPickedMinusNew = tblPairsMinus.ix[lio].copy()
                # tblPickedMinusNew['picked'] = lPicked.ix[lio]

                # tblPickedPlus = pd.concat( [tblPickedPlus, tblPickedPlusNew ] )
                # tblPickedMinus = pd.concat( [tblPickedMinus, tblPickedMinusNew ] )

                # update coverage
                mIivtArcvg[ iivtgtCur ] = arNewCvg

                # need to update master pairs list w/ #picked.
                # tblPairsPlus.ix[ tblPickedPlus.index, 'picked' ] = tblPickedPlus.picked
                # tblPairsMinus.ix[ tblPickedMinus.index, 'picked' ] = tblPickedMinus.picked
                tblPairsPlus.ix[ lPicked.index, 'picked' ] += lPicked
                tblPairsMinus.ix[ lPicked.index, 'picked' ] += lPicked

                arPickedThisRun[ lPicked.index ] += lPicked

    print 'PASS 2: picked %d probes, %d unique'%( tblPairsPlus.picked.sum(), tblPairsPlus[tblPairsPlus.picked>0].shape[0] )

    # tblPickedPlus = tblPairsPlus.query('picked>0')
    tblPickedPlus = tblPairsPlus.ix[ arPickedThisRun>0 ].copy()
    tblPickedMinus = tblPairsMinus.ix[ tblPickedPlus.index ].copy()

    print '(picked this run: %d unique probes, %d total)'%( (arPickedThisRun>0).sum(), arPickedThisRun.sum() )

    tblPickedPlus['picked' ] = arPickedThisRun
    tblPickedMinus['picked' ] = arPickedThisRun    

    outPairStore[ '%spairsPlus'%o.pairsBaseName ] = tblPickedPlus
    outPairStore[ '%spairsMinus'%o.pairsBaseName ] = tblPickedMinus

    copyMetadata( pairStore, outPairStore, '%spairsPlus'%o.pairsBaseName, '%spairsPlus'%o.pairsBaseName )
    copyMetadata( pairStore, outPairStore, '%spairsMinus'%o.pairsBaseName, '%spairsMinus'%o.pairsBaseName )
