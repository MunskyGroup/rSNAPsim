# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 16:01:24 2020

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
ncolor = 2
t_array = np.array([0,10,20,30,50,100,250,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,1000,1000,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*ncolor),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 1
start = time.time()

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

pl = np.zeros((2, len(kelong)), dtype=np.int32)

pl[:,[5,8,10,20,30,100,120,140]] = 1

pv = pl
#pl[ [10,140],1  ] = 1

pl = np.cumsum(pl,axis=1)

print(x0.shape)
print(all_results.shape)
print(all_frapresults.shape)
print(all_ribtimes.shape)
print(all_coltimes.shape)

all_col_points = []



kon = 3
koff = 3

k_probe = np.array([.8,.1],dtype=np.float64)


flags = np.array([0,0,0])
start = time.time()
for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)

    ssa_translation_lowmem.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,x0,9, pl, 2, kon, koff,k_probe,pv,flags)
    

  
    all_results[i,:,:] = result.T
    
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
    
print('time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))

plt.figure(dpi=300)
plt.plot(all_results[0])
plt.xlabel('time')
plt.ylabel('intensity')
plt.title('original lowmemory')
plt.legend(['c1','c2'])



# bursting, leaky, stats
flags = np.array([1,0,0])
start = time.time()
for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)

    ssa_translation_lowmem.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,x0,9, pl, 2, kon, koff,k_probe,pv,flags)
    

  
    all_results[i,:,:] = result.T
    
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
    
print('time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))

plt.figure(dpi=300)
plt.plot(all_results[0])
plt.xlabel('time')
plt.ylabel('intensity')
plt.title('bursting kon/koff = 3/3,  lowmemory')
plt.legend(['c1','c2'])




flags = np.array([0,1,0])
start = time.time()
for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)

    ssa_translation_lowmem.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,x0,9, pl, 2, kon, koff,k_probe,pv,flags)
    

  
    all_results[i,:,:] = result.T
    
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
    
print('time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))

plt.figure(dpi=300)
plt.plot(all_results[0])
plt.xlabel('time')
plt.ylabel('intensity')
plt.title('leaky,  lowmemory')
plt.legend(['c1','c2'])








flags = np.array([1,1,0])
start = time.time()
for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)

    ssa_translation_lowmem.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,x0,9, pl, 2, kon, koff,k_probe,pv,flags)
    

  
    all_results[i,:,:] = result.T
    
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
    
print('time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))

plt.figure(dpi=300)
plt.plot(all_results[0])
plt.xlabel('time')
plt.ylabel('intensity')
plt.title('leaky and bursting,  lowmemory')
plt.legend(['c1','c2'])











flags = np.array([0,0,1])
start = time.time()
for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)

    ssa_translation_lowmem.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,x0,9, pl, 2, kon, koff,k_probe,pv,flags)
    

  
    all_results[i,:,:] = result.T
    
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
    
print('time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))

plt.figure(dpi=300)
plt.hist(ribtimes[ribtimes > 0])
plt.xlabel('time')
plt.ylabel('intensity')
plt.title('ribosometimes,  lowmemory')
plt.legend(['c1','c2'])

