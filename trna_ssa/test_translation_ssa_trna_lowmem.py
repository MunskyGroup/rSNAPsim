# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 16:20:47 2020

@author: willi
"""


import numpy as np
import ssa_trna_lowmem
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
n_trajectories = 1
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
k_trna = np.array(list(strGeneCopy.values()))*.05 #D


#preallocated arrays here
ncolor = 2
all_results = np.zeros((n_trajectories,len(t_array),ncolor),dtype=np.int32)

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

ke = 10

ncolor = 2
pl = np.zeros((len(kelong),ncolor), dtype=np.int32)

pl[ [10,20,30,100,120,140],0  ] = 1
pl[ [10,140],1  ] = 1

pl = np.cumsum(pl,axis=0)
pl = pl.T.copy(order='C')

np.random.seed(0)



k_index = np.random.randint(0,61,len(kelong))

k_index = np.sort(k_index)



for i in range(n_trajectories):
    result = np.zeros((ncolor,len(t_array)),dtype=np.int32)    
    trna_result = np.zeros((len(t_array)*61),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    

    

    
    ssa_trna_lowmem.run_SSA(result,trna_result, k_index,k_trna,k_diffusion,ke,frapresult,t_array,kbind,kcompl, 0,0,0, seeds[i],x0,9 , pl, ncolor)
    
    

    all_results[i,:,:] = result.T
    print(all_results.shape)
    print(result.shape)
    all_frapresults[i,:] = frapresult

    all_ribs[i,:] = nribs[0]
    all_trna_results[i,:] = trna_result

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
    
    trna_pool_traj = all_trna_results[i,:].reshape((61,len(t_array))).T
    plt.plot(trna_pool_traj)

plt.figure(dpi=300)
plt.plot(all_results[0])
plt.ylabel('intensities')
plt.xlabel('time')


def moving_average1d(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def moving_average(a, n=3) :
    ret = np.cumsum(a,axis=1, dtype=float)
    ret[:,n:] = ret[:,n:] - ret[:,:-n]
    return ret[:,n - 1:] / n

import matplotlib.cm as cm
colors = cm.viridis(np.linspace(0, 1, 61))

plt.figure(dpi=500)
for i in range(0,61):
    plt.plot(moving_average1d(trna_pool_traj[:,i],n=50),alpha=.4,color=colors[i])

for i in range(0,61):
    plt.plot([0,1000],[k_trna[i], k_trna[i]],color=colors[i],alpha=.5)

plt.title('simulated tRNA levels per time (moving average of 50)')
plt.ylabel('trna_levels')
plt.xlabel('time')






plt.figure(dpi=300)
plt.bar(np.linspace(0,60,61).tolist(),np.mean(trna_pool_traj,axis=0))
plt.scatter(np.linspace(0,60,61).tolist(), k_trna)
plt.xlabel('tRNA ID')
plt.ylabel('mean levels')
plt.legend(['copy number','simulated level'])





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

