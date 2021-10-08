# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 15:42:06 2021

@author: willi
"""
import os
cwd = os.getcwd()
os.chdir('../../..')

import rsnapsim as rss
from rsnapsim import seqmanip

import numpy as np
import time
import matplotlib.pyplot as plt

os.chdir(cwd)




'''
Cases

 Solve SSA under:
     protein given
     protein ungiven with blank probes
     just kelong given (no probes)
     
     leaky probes
     FRAP
     Harringtonine
     
     low memory
     record stats
     
     binned with different footprint
     
     

 Solve ODE:
    w/correlation
    w/o correlation
    binned

 Solve Ballistic: 
     given just L and kelong
     given protein oobject


'''



aa_seq, poi_objs, tagged_pois, raw_seq = rss.seqmanip.open_seq_file(
    './Bactin_WithTags.txt')

bactin = poi_objs['1'][0]

bactin_nt = bactin.nt_seq[337*3:]
forward_rates = rss.propf.get_k(bactin_nt, .033,10,10)[1:-1]


t_array = np.linspace(0,2000,2001)
n_traj = 10

#blank probe
probe_loc = np.zeros((1,len(forward_rates))).astype(int)
probe_vec = np.zeros((1,len(forward_rates))).astype(int)

probe_loc[0,0] =0
probe_vec[0,0] =0
 
rss.solver.solve_ssa(forward_rates, t_array,
                     ki=.033, kt=10, n_traj=n_traj, probe_loc=probe_loc,
                     probe_vec=probe_vec, record_stats=False)

