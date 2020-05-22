# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:27:25 2020

@author: willi
"""

import numpy as np
import rSNAPsim 
rss = rSNAPsim.rSNAPsim()

#lets load in the sequence for the human insulin receptor INSR protein
try:
    rss.open_seq_file('gene_files/HUMINSR.gb')
except:
    rss.get_gb_file('M10051')
    

def get_binned_k(k,bins):
    '''
    evenly bins elongation rates as best it can.
    
    '''
    binsize = int(np.floor(len(k)/bins))
    binned_ks = []
    
    k_binned = np.zeros(bins)
    k_lens = np.ones(bins)*binsize
    
    to_redistribute = len(k)%bins

    k_lens[-to_redistribute:] = binsize+1
    
    inds = np.hstack(([0.], np.cumsum(k_lens))).astype(int)     
    print(inds)

    
    for i in range(0,bins):
        binned_ks = binned_ks +   [k[inds[i]:inds[i+1]].tolist(),]  
     
    for i in range(0,bins):            
        k_binned[i] = 1/np.sum(1/np.array(binned_ks[i]))

    return k_binned,k_lens


rss.run_default()

ke = rss.get_k(rss.pois_seq[0],.033,10 )[1:-1]
pv,pl = rss.get_probvec() 
binned_k = get_binned_k(np.array(ke),100)

ssa_obj = rss.ssa_solver()


pv_binned, pl_binned = binned_probe = rss.get_binned_probe_vec(pl,100)


t, mI, soln,m  = rss.build_ODE(binned_k[0],np.linspace(0,100,101),.033, pl_binned)


