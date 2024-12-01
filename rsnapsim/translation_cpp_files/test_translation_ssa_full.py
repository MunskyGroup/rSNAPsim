# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 21:54:50 2021

@author: willi
"""

import numpy as np
import ssa_translation_lowmem
import matplotlib.pyplot as plt
import time

# load the elongation 
kelong = np.loadtxt('elongationrates.txt')
kbind = kelong[0]
kcompl = kelong[-1]
kelong = kelong[1:-1]

ncolor = 1
t_array = np.array([0,10,20,30,50,100,250,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,1000,1000,dtype=np.float64)
N_rib = 20
result = np.zeros((len(t_array)*ncolor),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 100
start = time.time()








#preallocated arrays here

all_ribosome_locations = np.zeros((n_trajectories,len(t_array),N_rib),dtype=np.int32)

all_results = np.zeros((n_trajectories,len(t_array),ncolor),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)

all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories,dtype=np.int32)
x0 = np.zeros((N_rib),dtype=np.int32)

pl = np.zeros((1, len(kelong)), dtype=np.int32)

pl[:,0] = 1

pv = pl.astype(np.int32)
#pl[ [10,140],1  ] = 1

pl = np.cumsum(pl,axis=1)
pl = pl.astype(np.int32)
print(x0.shape)
print(all_results.shape)
print(all_frapresults.shape)
print(all_ribtimes.shape)
print(all_coltimes.shape)

all_col_points = []



kon = 3
koff = 3

k_probe = np.array([.8,.1],dtype=np.float64)

#record locations
flags = np.array([0,0,0],dtype=np.int32)
start = time.time()


for i in range(n_trajectories):
    
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    ribloc = np.zeros((N_rib,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)

    print(result.dtype)
    print(ribtimes.dtype)
    print(coltimes.dtype)
    print(colpointsx.dtype)
    print(colpointst.dtype)
    print(kelong.dtype)
    print(frapresult.dtype)
    print(t_array.dtype)
    print(type(.03))
    print(type(kcompl))
    print(type(300))
    print(type(seeds[i]))
    print(nribs.dtype)
    print(x0.dtype)
    print(pl.dtype)
    print(type(kon))
    print(type(koff))
    print(k_probe.dtype)
    print(pv.dtype)
    print(flags.dtype)
    print(type(N_rib))

    ssa_translation_lowmem.run_SSA_full( ribloc, result,ribtimes,coltimes,colpointsx,colpointst, kelong,t_array,.03,kcompl, 1,0,200,250, seeds[i],nribs,x0,9, pl, ncolor, kon, koff,k_probe,pv,flags)
    
    all_results[i,:,:] = result.T
    all_ribosome_locations[i,:,:] = ribloc.T
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)

full_rec_t = time.time()-start
print(full_rec_t/n_trajectories)


