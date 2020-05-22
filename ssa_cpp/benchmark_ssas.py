# -*- coding: utf-8 -*-
"""
Created on Thu May 21 16:02:42 2020

@author: willi
"""

import numpy as np
import ssa_translation_lowmem
import ssa_translation_lowmem_leaky
import ssa_translation_lowmem_nostats
import ssa_translation
import matplotlib.pyplot as plt
import time
import os

os.chdir('..')
from rss import ProbeVectorFactory as pvf
from rss import PropensityFactory as pff
os.chdir('ssa_cpp')


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
n_trajectories = 100


#preallocated arrays here
all_results = np.zeros((n_trajectories,len(t_array),ncolor),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)

all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)
x0 = np.zeros((N_rib),dtype=np.int32)

pl = np.zeros((len(kelong),ncolor), dtype=np.int32)

pl[ [10,20,30,100,120,140],0  ] = 1
#pl[ [10,140],1  ] = 1

pl = np.cumsum(pl,axis=0)
pl = pl.T.copy(order='C')

print('-----------------------')
print('GENERATING REPORT')
print('-----------------------')
print('1 cpu core, 2 color')
print('{0} base pairs'.format(len(kelong)-2))
print('-----------------------')



all_col_points = []
start = time.time()
for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)
    
    ssa_translation_lowmem.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,x0,9, pl,2)
    

  
    all_results[i,:,:] = result.T
    
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
print('low memory w/recording_stats: time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))
#plt.hist(result[result>0])
#plt.show()
#traj = result.reshape((N_rib,len(t_array))).T
##print('The result is \n {0}'.format(result.reshape((N_rib,len(t_array))).T))
#plt.plot(traj[-1,:])
#plt.show()

plt.plot(result.T,'--')

start = time.time()
for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    
    ssa_translation_lowmem_nostats.run_SSA(result, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],x0,9, pl,2)
    

  
    all_results[i,:,:] = result.T
    
    all_frapresults[i,:] = frapresult


print('Low memory w/o recording stats: time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))




plt.plot(result.T)




pl_2color = np.atleast_2d(pl)
probe_loc = (pl_2color[:,1:]-pl_2color[:,:-1] > 0).astype(int)
inds = pff.intellegent_bin(np.atleast_2d(probe_loc),100)
bpv,bpl = pvf.bin_probe_vecs(probe_loc,inds)

kelong_2color = kelong
k_bin = pff.bin_k(kelong, inds)




start = time.time()
for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    
    ssa_translation_lowmem_nostats.run_SSA(result, k_bin,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],x0,9, bpl,2)
    

  
    #all_results[i,:,:] = result.T
    
    #all_frapresults[i,:] = frapresult


print('Low memory 100 bins: time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))



plt.plot(result.T,'.')



plt.legend(['w/ stats color 1','w/ stats color 2','w/o stats color 1','w/o stats color 2'])
plt.xlabel('time')
plt.ylabel('intensity')


print('-----------------------')
print('1 cpu core, 1 color')
print('-----------------------')

all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)

all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)
x0 = np.zeros((N_rib),dtype=np.int32)

pv = np.loadtxt('probe_design.txt')
all_col_points = []
start = time.time()



for i in range(n_trajectories):
    result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)
    
    ssa_translation.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,x0,9,N_rib)
    
    

    all_results[i,:] = result

    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
    
ntimes = len(t_array)
intensity_vec = np.zeros(ntimes)

tstart = 0
I = np.zeros((n_trajectories,ntimes-tstart))
for i in range(n_trajectories):
    traj = all_results[i,:].reshape((N_rib,len(t_array))).T
    for j in range(tstart,ntimes):
        temp_output = traj[j,:]
        I[i,j] = np.sum(pv[temp_output[temp_output>0]-1])
        
print('Full SSA with Recording: time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))
plt.figure()
plt.plot(I[-1,:],'x')



pl = np.atleast_2d(pv.astype(int))
ncolor=1

all_results = np.zeros((n_trajectories,len(t_array),ncolor),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)

all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
#seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)
x0 = np.zeros((N_rib),dtype=np.int32)


all_col_points = []
start = time.time()
for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)
    
    ssa_translation_lowmem.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,x0,9, pl,1)
    

  
    all_results[i,:] = result.T
    
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
print('Low memory w/recording_stats: time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))
#plt.hist(result[result>0])
#plt.show()
#traj = result.reshape((N_rib,len(t_array))).T
##print('The result is \n {0}'.format(result.reshape((N_rib,len(t_array))).T))
#plt.plot(traj[-1,:])
#plt.show()
plt.plot(result.T,'o')


start = time.time()
for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    
    ssa_translation_lowmem_nostats.run_SSA(result, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],x0,9, pl,1)
    

  
    all_results[i,:] = result.T
    
    all_frapresults[i,:] = frapresult


print('Low memory w/o recording stats: time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))


plt.plot(result.T)

pl = pl.flatten()
probe_loc = (np.where(pv[1:]-pv[:-1] > 0)[0]+1).astype(np.int32)
k_probe = .2

all_results = np.zeros((n_trajectories,len(t_array)),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)

all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)
x0 = np.zeros((N_rib),dtype=np.int32)

start = time.time()
for i in range(n_trajectories):
    result = np.zeros((len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)
    
    ssa_translation_lowmem_leaky.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,x0,9, pl,k_probe,probe_loc )
    
  
    all_results[i,:] = result
    
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
print('Low memory w/ leaky probes: time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))

plt.plot(result.T)
plt.legend(['full ssa','lowmem w/ stats','lowmem w/o stats','lowmem leaky'])
plt.xlabel('time')
plt.ylabel('intensity')