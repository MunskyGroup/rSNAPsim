# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 18:17:09 2020

@author: willi
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:07:30 2020

@author: willi
"""

#Test file for generalized SSA

#import rSNAPsim as rss
import numpy as np
import time
import matplotlib.pyplot as plt

import ssa_translation_generic

def generate_additional_ks(enters,pauses,jumps,stops,L):
    
    def frame_check_1(L,arr):        
        return (L- arr[:,1]+1)*(arr[:,1]>0) + L*(arr[:,1]>1)
    
    def frame_check_3(L,arr):        
        return (L- arr[:,3]+1)*(arr[:,3]>0) + L*(arr[:,3]>1)            
                
    def gen_ks_1_loc(L,arr):
        arr[:,0] = arr[:,0]+frame_check_1(L,arr)
        arr[:,1] = arr[:,2]    
        arr = arr[:,0:2]
        max_arr = np.max( arr[:,0])     
        return arr,max_arr
    
    def gen_ks_3_loc(L,arr):
        arr[:,0] = arr[:,0]+ frame_check_1(L,arr)     
        arr[:,1] = arr[:,2]+ frame_check_3(L,arr)
        arr[:,2] = arr[:,4]
        arr = arr[:,0:3]
        max_arr = max([np.max( arr[:,0]),np.max( arr[:,1])])
        return arr,max_arr

    max_enter = 0
    max_pause = 0
    max_stop = 0
    max_jump = 0
    k_jumps = np.copy(jumps)
    k_pauses = np.copy(pauses)
    k_stops = np.copy(stops)
    k_enters = np.copy(enters)
        
    
    if len(k_enters) != 0:
        k_enters,max_enter = gen_ks_1_loc(L,k_enters)

    if len(k_pauses) != 0:
        k_pauses,max_pause = gen_ks_1_loc(L,k_pauses)

    if len(k_stops) != 0:
        k_stops,max_stop = gen_ks_1_loc(L,k_stops)
    
    if len(k_jumps) != 0:
        k_jumps,max_jump = gen_ks_3_loc(L,k_jumps)
        
    max_loc = max(max_jump,max_stop,max_pause,max_enter)
    
    if max_loc <=L: 
        frames_used = 0
    if max_loc > L:
        frames_used = 1
    if max_loc > 2*L-1 :
        frames_used = 2
    
    return k_enters, k_pauses, k_stops, k_jumps, frames_used



#rsnap = rss.rSNAPsim()
#rsnap.open_seq_file('gene_files/H2B_withTags.txt')
#rsnap.run_default()

k = np.ones((1,300)).flatten()

kelong = k[1:-1]
kelong[49] = 0
kelong[149]= 0

kelong[248] = 0

#k_fss = np.array([[200,0,200,1,.3]])
k_pause = np.array([[30,0,.01]])
k_enters = np.array([[10,0,.02],[10,2,.04]],dtype=np.float64)
k_stops = np.array([[50,0,10],[50,1,10],[50,2,10]],dtype=np.float64)
k_fss = np.array([[20,0,20,1,1]],dtype=np.float64)
#k_pause =  np.array([[30,2,100],[40,2,100]],dtype=np.float64)

k_enters,k_pauses,k_stops,k_jumps,frames_used = generate_additional_ks(k_enters,[],k_fss,k_stops,100)




t_array = np.array([0,100,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,400,400,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 1
start = time.time()
all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)




all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)


k_add = np.hstack((k_enters.flatten(),k_pauses.flatten(),k_stops.flatten(),k_jumps.flatten() ))
t_array_copy = np.copy(t_array)
while t_array_copy.shape[0] != 200:
    t_array_copy = np.vstack((t_array_copy,t_array))
 


for i in range(n_trajectories):
    result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    ssa_translation_generic.run_SSA_generic(result,ribtimes,coltimes, kelong,frapresult,t_array, np.array([0,0,0],dtype=np.float64), seeds[i],nribs, k_add.flatten()  ,2,0,3,1 )
    all_results[i,:] = result
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    

traj = all_results[0,:].reshape((N_rib,len(t_array))).T

f,ax = plt.subplots(2,1)    

ax[0].set_ylim([0,300])
ax[0].fill_between([0,400],[100,100],color='red',alpha=.2)
ax[0].fill_between([0,400],[200,200],color='green',alpha=.2)
ax[0].fill_between([0,400],[300,300],color='blue',alpha=.2)

ax[0].plot(traj,'.')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Ribosome Location')
ax[0].set_title(' 100 codons,  enters: 0,10 and +2,10  FSS: 0,20 to +1,20 Stops: 50 0,1,2'  )

spatial_x = (traj + (traj > 100) + (traj > 199))%100

ax[1].set_ylim([0,100])


#ax[1].plot(t_array,spatial_x,'.')

ax[1].plot(t_array_copy.T[traj<=100],spatial_x[traj <= 100],'r.')
ax[1].plot(t_array_copy.T[traj>100],spatial_x[traj > 100],'g.')

ax[1].plot(t_array_copy.T[traj>199],spatial_x[traj > 199],'b.')

ax[1].set_xlabel('Time')
ax[1].set_ylabel('Ribosome Location')
ax[1].set_title(' spatial location '  )
ax[1].legend(['0','+1','+2'])












###################################################################

k = np.ones((1,300)).flatten()

kelong = k[1:-1]
kelong[49] = 3
kelong[79] = 0

k_enters = np.array([[10,0,.04]],dtype=np.float64)
k_stops = np.array([[50,0,10],[80,0,10]],dtype=np.float64)
k_fss = []
k_pause = []
#k_pause =  np.array([[30,2,100],[40,2,100]],dtype=np.float64)

k_enters,k_pauses,k_stops,k_jumps,frames_used = generate_additional_ks(k_enters,k_pause,k_fss,k_stops,100)




t_array = np.array([0,100,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,400,400,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 1
start = time.time()
all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)




all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)


k_add = np.hstack((k_enters.flatten(),k_pauses.flatten(),k_stops.flatten(),k_jumps.flatten() ))
t_array_copy = np.copy(t_array)
while t_array_copy.shape[0] != 200:
    t_array_copy = np.vstack((t_array_copy,t_array))
 


for i in range(n_trajectories):
    result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    ssa_translation_generic.run_SSA_generic(result,ribtimes,coltimes, kelong,frapresult,t_array, np.array([0,0,0],dtype=np.float64), seeds[i],nribs, k_add.flatten()  ,len(k_enters),len(k_pauses),len(k_stops),len(k_jumps) )
    all_results[i,:] = result
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    

traj = all_results[0,:].reshape((N_rib,len(t_array))).T

f,ax = plt.subplots(2,1)    

ax[0].set_ylim([0,300])
ax[0].fill_between([0,400],[100,100],color='red',alpha=.2)
ax[0].fill_between([0,400],[200,200],color='green',alpha=.2)
ax[0].fill_between([0,400],[300,300],color='blue',alpha=.2)

ax[0].plot(traj,'.')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Ribosome Location')
ax[0].set_title(' 100 codons,  enters: 0,10 stops: 0,50 and 0,80'  )

spatial_x = (traj + (traj > 100) + (traj > 199))%100

ax[1].set_ylim([0,100])


#ax[1].plot(t_array,spatial_x,'.')

ax[1].plot(t_array_copy.T[traj<=100],spatial_x[traj <= 100],'r.')
ax[1].plot(t_array_copy.T[traj>100],spatial_x[traj > 100],'g.')

ax[1].plot(t_array_copy.T[traj>199],spatial_x[traj > 199],'b.')

ax[1].set_xlabel('Time')
ax[1].set_ylabel('Ribosome Location')
ax[1].set_title(' spatial location '  )
ax[1].legend(['0','+1','+2'])




#####################################################
k = np.ones((1,300)).flatten()

kelong = k[1:-1]

kelong[49] = 0
kelong[179] = 0

k_enters = np.array([[10,0,.04]],dtype=np.float64)
k_stops = np.array([[50,0,10],[80,1,10]],dtype=np.float64)
k_fss = np.array([[30,0,30,1,1]],dtype=np.float64)
k_pause = []
#k_pause =  np.array([[30,2,100],[40,2,100]],dtype=np.float64)

k_enters,k_pauses,k_stops,k_jumps,frames_used = generate_additional_ks(k_enters,k_pause,k_fss,k_stops,100)




t_array = np.array([0,100,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,400,400,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 1
start = time.time()
all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)




all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)


k_add = np.hstack((k_enters.flatten(),k_pauses.flatten(),k_stops.flatten(),k_jumps.flatten() ))
t_array_copy = np.copy(t_array)
while t_array_copy.shape[0] != 200:
    t_array_copy = np.vstack((t_array_copy,t_array))
 


for i in range(n_trajectories):
    result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    ssa_translation_generic.run_SSA_generic(result,ribtimes,coltimes, kelong,frapresult,t_array, np.array([0,0,0],dtype=np.float64), seeds[i],nribs, k_add.flatten()  ,1,0,2,1 )
    all_results[i,:] = result
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    

traj = all_results[0,:].reshape((N_rib,len(t_array))).T

f,ax = plt.subplots(2,1)    

ax[0].set_ylim([0,300])
ax[0].fill_between([0,400],[100,100],color='red',alpha=.2)
ax[0].fill_between([0,400],[200,200],color='green',alpha=.2)
ax[0].fill_between([0,400],[300,300],color='blue',alpha=.2)

ax[0].plot(traj,'.')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Ribosome Location')
ax[0].set_title(' 100 codons,  enters: 0,10 stops: 0,50 and 1,80  FSS: 0,30 to 1,30'  )

spatial_x = (traj + (traj > 100) + (traj > 199))%100

ax[1].set_ylim([0,100])


#ax[1].plot(t_array,spatial_x,'.')

ax[1].plot(t_array_copy.T[traj<=100],spatial_x[traj <= 100],'r.')
ax[1].plot(t_array_copy.T[traj>100],spatial_x[traj > 100],'g.')

ax[1].plot(t_array_copy.T[traj>199],spatial_x[traj > 199],'b.')

ax[1].set_xlabel('Time')
ax[1].set_ylabel('Ribosome Location')
ax[1].set_title(' spatial location '  )
ax[1].legend(['0','+1','+2'])




######################

k = np.ones((1,300)).flatten()

kelong = k[1:-1]

kelong[49] = 0
kelong[278] = 0

k_enters = np.array([[10,0,.04],[10,2,.02]],dtype=np.float64)
k_stops = np.array([[50,0,10],[80,2,10]],dtype=np.float64)
k_fss = []
k_pause = []
#k_pause =  np.array([[30,2,100],[40,2,100]],dtype=np.float64)

k_enters,k_pauses,k_stops,k_jumps,frames_used = generate_additional_ks(k_enters,k_pause,k_fss,k_stops,100)




t_array = np.array([0,100,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,400,400,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 1
start = time.time()
all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)




all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)


k_add = np.hstack((k_enters.flatten(),k_pauses.flatten(),k_stops.flatten(),k_jumps.flatten() ))
t_array_copy = np.copy(t_array)
while t_array_copy.shape[0] != 200:
    t_array_copy = np.vstack((t_array_copy,t_array))
 


for i in range(n_trajectories):
    result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    ssa_translation_generic.run_SSA_generic(result,ribtimes,coltimes, kelong,frapresult,t_array, np.array([0,0,0],dtype=np.float64), seeds[i],nribs, k_add.flatten()  ,2,0,2,0 )
    all_results[i,:] = result
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    

traj = all_results[0,:].reshape((N_rib,len(t_array))).T

f,ax = plt.subplots(2,1)    

ax[0].set_ylim([0,300])
ax[0].fill_between([0,400],[100,100],color='red',alpha=.2)
ax[0].fill_between([0,400],[200,200],color='green',alpha=.2)
ax[0].fill_between([0,400],[300,300],color='blue',alpha=.2)

ax[0].plot(traj,'.')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Ribosome Location')
ax[0].set_title(' 100 codons,  enters: 0,10 2,20 stops: 0,50 and 2,80'  )

spatial_x = (traj + (traj > 100) + (traj > 199))%100

ax[1].set_ylim([0,100])


#ax[1].plot(t_array,spatial_x,'.')

ax[1].plot(t_array_copy.T[traj<=100],spatial_x[traj <= 100],'r.')
ax[1].plot(t_array_copy.T[traj>100],spatial_x[traj > 100],'g.')

ax[1].plot(t_array_copy.T[traj>199],spatial_x[traj > 199],'b.')

ax[1].set_xlabel('Time')
ax[1].set_ylabel('Ribosome Location')
ax[1].set_title(' spatial location '  )
ax[1].legend(['0','+1','+2'])








###########

k = np.ones((1,300)).flatten()

kelong = k[1:-1]

kelong[49] = 0

kelong[39] = 0.1
kelong[278] = 0

k_enters = np.array([[10,0,.04],[10,2,.02]],dtype=np.float64)
k_stops = np.array([[50,0,10],[80,2,10]],dtype=np.float64)
k_fss = []
k_pause =  np.array([[40,0,100]],dtype=np.float64)
#k_pause =  np.array([[30,2,100],[40,2,100]],dtype=np.float64)

k_enters,k_pauses,k_stops,k_jumps,frames_used = generate_additional_ks(k_enters,k_pause,k_fss,k_stops,100)




t_array = np.array([0,100,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,400,400,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 1
start = time.time()
all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)




all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)


k_add = np.hstack((k_enters.flatten(),k_pauses.flatten(),k_stops.flatten(),k_jumps.flatten() ))
t_array_copy = np.copy(t_array)
while t_array_copy.shape[0] != 200:
    t_array_copy = np.vstack((t_array_copy,t_array))
 


for i in range(n_trajectories):
    result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    ssa_translation_generic.run_SSA_generic(result,ribtimes,coltimes, kelong,frapresult,t_array, np.array([0,0,0],dtype=np.float64), seeds[i],nribs, k_add.flatten()  ,len(k_enters),len(k_pauses),len(k_stops),len(k_jumps))
    all_results[i,:] = result
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    

traj = all_results[0,:].reshape((N_rib,len(t_array))).T

f,ax = plt.subplots(2,1)    

ax[0].set_ylim([0,300])
ax[0].fill_between([0,400],[100,100],color='red',alpha=.2)
ax[0].fill_between([0,400],[200,200],color='green',alpha=.2)
ax[0].fill_between([0,400],[300,300],color='blue',alpha=.2)

ax[0].plot(traj,'.')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Ribosome Location')
ax[0].set_title(' 100 codons,  enters: 0,10 2,20 stops: 0,50 and 2,80'  )

spatial_x = (traj + (traj > 100) + (traj > 199))%100

ax[1].set_ylim([0,100])


#ax[1].plot(t_array,spatial_x,'.')

ax[1].plot(t_array_copy.T[traj<=100],spatial_x[traj <= 100],'r.')
ax[1].plot(t_array_copy.T[traj>100],spatial_x[traj > 100],'g.')

ax[1].plot(t_array_copy.T[traj>199],spatial_x[traj > 199],'b.')

ax[1].set_xlabel('Time')
ax[1].set_ylabel('Ribosome Location')
ax[1].set_title(' spatial location '  )
ax[1].legend(['0','+1','+2'])


