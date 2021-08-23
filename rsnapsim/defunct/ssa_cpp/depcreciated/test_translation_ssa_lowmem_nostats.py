# -*- coding: utf-8 -*-
"""
Created on Thu May 21 15:31:44 2020

@author: willi
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 16:01:24 2020

@author: willi
"""

import numpy as np
import ssa_translation_lowmem_nostats
import matplotlib.pyplot as plt
import time

# load the elongation 
kelong = np.loadtxt('elongationrates.txt')
kbind = kelong[0]
kcompl = kelong[-1]
kelong = kelong[1:-1]
ncolor = 2
t_array = np.array([0,10,20,30,50,100,250,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,1000,1000,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*ncolor),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)


n_trajectories = 1000
start = time.time()

#preallocated arrays here
all_results = np.zeros((n_trajectories,len(t_array),ncolor),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)


seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)
x0 = np.zeros((N_rib),dtype=np.int32)

pl = np.zeros((len(kelong)-2,ncolor), dtype=np.int32)

pl[ [10,20,30,100,120,140],0  ] = 1
pl[ [10,140],1  ] = 1

pl = np.cumsum(pl,axis=0)
pl = pl.T.copy(order='C')


print(x0.shape)
print(all_results.shape)
print(all_frapresults.shape)


for i in range(n_trajectories):
    print(seeds[i])
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    
    ssa_translation_lowmem_nostats.run_SSA(result, kelong,frapresult,t_array,.53,kcompl, 1,0,300, seeds[i],x0,9, pl,2)
    

  
    all_results[i,:,:] = result.T
    
    all_frapresults[i,:] = frapresult


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

