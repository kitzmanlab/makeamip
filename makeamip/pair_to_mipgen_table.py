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

from makeamip.capcommon import *

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()
   
    opts.add_argument('--inStorePairs', dest='inStorePairs')
    opts.add_argument('--tablePairs', dest='tablePairs')
        
    opts.add_argument('--inStoreArms', dest='inStoreArms')        

    opts.add_argument('--mipGenOut', default=None, dest='mipGenOut')

    o = opts.parse_args()

    pairStore = pd.HDFStore( o.inStorePairs, 'r' )
    tblPairs = pairStore[ o.tablePairs ]       

    (tblnameExtArm,tblnameLigArm) = pairStore.get_storer( o.tablePairs  ).attrs['mippipe__ext_and_lig_arm_tables']

    armStore = pd.HDFStore( o.inStoreArms, 'r' )
    tblArmExt = armStore[ tblnameExtArm ]
    tblArmLig = armStore[ tblnameLigArm ]

    tblOut = OrderedDict( [ (col,[]) 
                            for col in ['mip_key','chr','ext_probe_start','ext_probe_stop','lig_probe_start','lig_probe_stop'] ] )

    for i,r in tblPairs.iterrows():
        tblOut['mip_key'].append( int(i) )
        tblOut['chr'].append( r.chrom )
        tblOut['ext_probe_start'].append( tblArmExt.ix[ r.extarmidx, 'start' ] )
        tblOut['ext_probe_stop'].append( tblArmExt.ix[ r.extarmidx, 'end' ] )
        tblOut['lig_probe_start'].append( tblArmLig.ix[ r.ligarmidx, 'start' ] )
        tblOut['lig_probe_stop'].append( tblArmLig.ix[ r.ligarmidx, 'end' ] )

    tblOut = pd.DataFrame( tblOut )
    tblOut.to_csv( o.mipGenOut, sep='\t', index=False )