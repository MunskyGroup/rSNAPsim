# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 18:31:02 2020

@author: willi
"""

import numpy as np
import ssa_trna
import matplotlib.pyplot as plt
import time

# load the elongation 
kelong = np.loadtxt('elongationrates.txt')
kbind = kelong[0]
kcompl = kelong[-1]
kelong = kelong[1:-1]

t_array = np.array([0,100,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,1000,1000,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 25
start = time.time()





strGeneCopy = {'TTT': 17.6, 'TCT': 15.2, 'TAT': 12.2, 'TGT': 10.6, 'TTC': 20.3,
                            'TCC': 17.7, 'TAC': 15.3, 'TGC': 12.6, 'TTA': 7.7, 'TCA': 12.2,
                            'TAA': 1.0, 'TGA': 1.6, 'TTG': 12.9, 'TCG':  4.4, 'TAG': 0.8,
                            'TGG': 13.2, 'CTT': 13.2, 'CCT': 17.5, 'CAT': 10.9, 'CGT': 4.5,
                            'CTC': 19.6, 'CCC': 19.8, 'CAC': 15.1, 'CGC': 10.4, 'CTA':  7.2,
                            'CCA': 16.9, 'CAA': 12.3, 'CGA':  6.2, 'CTG': 39.6, 'CCG':  6.9,
                            'CAG': 34.2, 'CGG': 11.4, 'ATT': 16.0, 'ACT': 13.1, 'AAT': 17.0,
                            'AGT': 12.1, 'ATC': 20.8, 'ACC': 18.9, 'AAC': 19.1, 'AGC': 19.5,
                            'ATA':  7.5, 'ACA': 15.1, 'AAA': 24.4, 'AGA': 12.2, 'ATG': 22.0,
                            'ACG': 6.1, 'AAG': 31.9, 'AGG': 12.0, 'GTT': 11.0, 'GCT': 18.4,
                            'GAT': 21.8, 'GGT': 10.8, 'GTC': 14.5, 'GCC': 27.7, 'GAC': 25.1,
                            'GGC': 22.2, 'GTA':  7.1, 'GCA': 15.8, 'GAA': 29.0, 'GGA': 16.5,
                            'GTG': 28.1, 'GCG': 7.4, 'GAG': 39.6, 'GGG': 16.5}

strGeneCopy.pop('TAG')
strGeneCopy.pop('TAA')
strGeneCopy.pop('TGA')
k_trna = np.array(list(strGeneCopy.values()))*.5


#preallocated arrays here
all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)

all_trna_results = np.zeros((n_trajectories,61*len(t_array)),dtype=np.int32)

lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)

all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
np.random.seed(0)
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)
x0 = np.zeros((N_rib),dtype=np.int32)


print(x0.shape)
print(all_results.shape)
print(all_frapresults.shape)
print(all_ribtimes.shape)
print(all_coltimes.shape)

all_col_points = []
k_diffusion = 1

np.random.seed(0)
k_index = np.random.randint(0,61,len(kelong))
for i in range(n_trajectories):
    result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
    trna_result = np.zeros((len(t_array)*61),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)
    
    ssa_trna.run_SSA(result,trna_result,ribtimes,coltimes,colpointsx,colpointst, k_index,k_trna,k_diffusion,frapresult,t_array,.1,kcompl, 0,0,0, seeds[i],nribs,x0,10)
    
    

    all_results[i,:] = result
    print(all_results.shape)
    print(result.shape)
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]
    all_trna_results[i,:] = trna_result

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
I = np.zeros((n_trajectories,ntimes-tstart))
plt.figure()
for i in range(n_trajectories):
    traj = all_results[i,:].reshape((N_rib,len(t_array))).T
    trna_pool_traj = all_trna_results[i,:].reshape((61,len(t_array))).T
    plt.plot(trna_pool_traj)
    for j in range(tstart,ntimes):
        temp_output = traj[j,:]
        I[i,j] = np.sum(pv[temp_output[temp_output>0]-1])

plt.figure(dpi=300)
plt.plot(traj)
plt.ylabel('rib position')
plt.xlabel('time')

plt.figure(dpi=300)
plt.plot(trna_pool_traj)

plt.ylabel('trna_levels')
plt.xlabel('time')
#plt.legend(np.linspace(0,60,61).tolist())
# # Plotting
# all_traj = np.loadtxt('ivec_1000t')
# #all_traj = ssa_obj_01.intensity_vec
# f,ax = plt.subplots(2,1)
# ax[0].plot(I[0:50,-500:].T)
# ax[1].plot(all_traj[0:50,-500:].T)
# f2,ax2 = plt.subplots(1,2)
# ax2[0].hist(I[:,-500:].ravel())
# ax2[1].hist(all_traj[:100,-500:].ravel())
# plt.show()

