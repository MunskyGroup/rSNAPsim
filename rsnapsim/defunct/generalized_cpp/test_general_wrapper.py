# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 15:15:25 2020

@author: willi
"""

import matplotlib.pyplot as plt
import time
import numpy as np
kelong = np.loadtxt('elongationrates.txt')
kbind = kelong[0]
kcompl = kelong[-1]
kelong = (np.ones((464*2,1))*3.9).flatten()

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




kenters = np.array([[0,0,kbind],[600,1,kbind]])
kjumps = np.array([[100,0, 200,1, 10]])



def generate_additional_ks(k_enters,k_pauses,k_jumps,k_stops,L):

    max_enter = 0
    max_pause = 0
    max_stop = 0
    max_jump = 0
    
    if k_enters != []:
        k_enters[:,0] = k_enters[:,0]+L*k_enters[:,1]
        k_enters[:,1] = k_enters[:,2]    
        k_enters = k_enters[:,0:2]
        max_enter = np.max( k_enters[:,0])

    if k_pauses != []:
        k_pauses[:,0] = k_pauses[:,0]+ L*k_pauses[:,1]
        k_pauses[:,1] = k_pauses[:,2]
        k_pauses = k_pauses[:,0:2]
        max_pause = np.max( k_pauses[:,0])

    if k_stops != []:
        k_stops[:,0] = k_stops[:,0]+L*k_stops[:,1]
        k_stops[:,1] = k_stops[:,2]    
        k_stops = k_stops[:,0:2]
        max_stop = np.max( k_stops[:,0])
    
    if k_jumps != []:
        k_jumps[:,0] = k_jumps[:,0]+ L*k_jumps[:,1]
        
        k_jumps[:,1] = k_jumps[:,2]+ L*k_jumps[:,3]
        k_jumps[:,2] = k_jumps[:,4]
        k_jumps = k_jumps[:,0:3]
        
        max_jump = max([np.max( k_jumps[:,0]),np.max( k_jumps[:,1])])
        
    max_loc = max(max_jump,max_stop,max_pause,max_enter)
    
    if max_loc <=L: 
        frames_used = 0
    if max_loc > L:
        frames_used = 1
    if max_loc > 2*L :
        frames_used = 1
    
    return k_enters, k_pauses, k_stops, k_jumps, frames_used
    

