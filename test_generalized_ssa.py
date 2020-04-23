# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:07:30 2020

@author: willi
"""

#Test file for generalized SSA

import rSNAPsim as rss
import numpy as np
import time


def generate_additional_ks(enters,pauses,jumps,stops,L):
    
    def frame_check_1(L,arr):        
        return (L- (arr[:,1]-1))*(arr[:,1]>0)
    
    def frame_check_3(L,arr):        
        return (L- (arr[:,3]-1))*(arr[:,3]>0)                
                
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
        
    
    if k_enters != []:
        k_enters,max_enter = gen_ks_1_loc(L,k_enters)

    if k_pauses != []:
        k_pauses,max_pause = gen_ks_1_loc(L,k_pauses)

    if k_stops != []:
        k_stops,max_stop = gen_ks_1_loc(L,k_stops)
    
    if k_jumps !=[]:
        k_jumps,max_jump = gen_ks_3_loc(L,k_jumps)
        
    max_loc = max(max_jump,max_stop,max_pause,max_enter)
    
    if max_loc <=L: 
        frames_used = 0
    if max_loc > L:
        frames_used = 1
    if max_loc > 2*L-1 :
        frames_used = 2
    
    return k_enters, k_pauses, k_stops, k_jumps, frames_used



rsnap = rss.rSNAPsim()
rsnap.open_seq_file('gene_files/H2B_withTags.txt')
rsnap.run_default()

k = rsnap.get_k(rsnap.sequence_str,.02,10)

kelong = k[1:-1]


k_fss = np.array([[200,0,200,1,.3]])
k_pause = np.array([[300,1,.5]])

k_enters,k_pauses,k_stops,k_jumps,frames_used = generate_additional_ks([],k_pause,k_fss,[],len(kelong))




t_array = np.array([0,100,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,1000,1000,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 100
start = time.time()
all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)




all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)


for i in range(n_trajectories):
    result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    ssa_translation_generic.run_SSA(result,ribtimes,coltimes, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs,200,.03)
    all_results[i,:] = result
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    
    
    
    

