# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 16:01:24 2020

@author: willi
"""

import numpy as np
import ssa_translation_lowmem_leaky
import matplotlib.pyplot as plt
import time

# load the elongation 
kelong = np.loadtxt('elongationrates.txt')
kbind = kelong[0]
kcompl = kelong[-1]
kelong = kelong[1:-1]

t_array = np.array([0,10,20,30,50,100,250,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,1000,1000,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 10
start = time.time()

#preallocated arrays here
all_results = np.zeros((n_trajectories,len(t_array)),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)

all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)
x0 = np.zeros((N_rib),dtype=np.int32)

pl = np.zeros((len(kelong)+1), dtype=np.int32)

pl[ [10,20,30,100,120,140]  ] = 1

pl = np.cumsum(pl)
print(x0.shape)
print(all_results.shape)
print(all_frapresults.shape)
print(all_ribtimes.shape)
print(all_coltimes.shape)

all_col_points = []


k_probe = .01
probe_loc = np.array([1,2,3,4,5,6,7,8,9,10] ,dtype=np.int32)

for i in range(n_trajectories):
    result = np.zeros((len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)
    
    ssa_translation_lowmem_leaky.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,x0,9, pl,k_probe,probe_loc )
    
    plt.plot(result)
  
    all_results[i,:] = result
    
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
print('time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))
#plt.hist(result[result>0])
#plt.show()
#traj = result.reshape((N_rib,len(t_array))).T
##print('The result is \n {0}'.format(result.reshape((N_rib,len(t_array))).T))
#plt.plot(traj[-1,:])
#plt.show()

# map to fluorescence.
ntimes = len(t_array)
intensity_vec = np.zeros(ntimes)
pv = np.loadtxt('probe_design.txt')
tstart = 0
I = all_results


# Plotting
all_traj = np.loadtxt('ivec_1000t')
#all_traj = ssa_obj_01.intensity_vec
f,ax = plt.subplots(2,1)
ax[0].plot(I[0:50,-500:].T)
ax[1].plot(all_traj[0:50,-500:].T)
f2,ax2 = plt.subplots(1,2)
ax2[0].hist(I[:,-500:].ravel())
ax2[1].hist(all_traj[:100,-500:].ravel())
plt.show()

