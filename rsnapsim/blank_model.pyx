# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 14:26:26 2021

@author: wsraymon
"""

# cython: infer_types=True
# distutils: language = c++

import numpy as np
cimport numpy as np
cimport cython


#min length goes here 


#parsed rules go here


#original rules go here


#cdef goes here
    void generic_ssa_cpp(int* result, int* intensity, int* states, int* Stoich_states, int* Stoich_lattice, double* forward_rates,
                         double* parameters, int* xi_lattice, int* xi_state, double* time_vector,
                         double tf, int seed, int Nt,
                         int n_rxns, int n_total_rxns, int n_states, int length,
                         int max_particles, int used_frames,
                         int* probe_location_matrix,   int Ncolors);
    
def run_ssa_cpp( np.ndarray[int, ndim=2, mode="c"] result not None,
                 np.ndarray[int, ndim=2, mode="c"] intensity not None,
                 np.ndarray[int, ndim=2, mode="c"] states not None,
                 np.ndarray[int, ndim=2, mode="c"] stoich_states not None,
                 np.ndarray[int, ndim=2, mode="c"] stoich_lattice not None,
                 np.ndarray[double, ndim=1, mode="c"] parameters not None,
                 np.ndarray[double, ndim=1, mode="c"] kelong not None, 
                 np.ndarray[double, ndim=1, mode="c"] t_array not None,
                 np.ndarray[int,ndim=2, mode="c"] xi_lattice not None,
                 np.ndarray[int,ndim=2, mode="c"] xi_states not None,
                 np.ndarray[int,ndim=2, mode="c"] probe_locations not None,
                 int length,
                 unsigned int seed,
                 int n_total_reactions):
    
    
    cdef int used_frames = 1
    
    if len(xi_states) > 1.1*length:
        used_frames = 2
        if len(xi_states) > 2*length:
            used_frames = 3
            
    #print(used_frames)
    
    cdef int max_particles = result.shape[1]
    #print(max_particles)
    
    cdef int ncolors = probe_locations.shape[0]

    cdef int n_states = xi_states.shape[1]
    cdef int n_reaction_states = stoich_states.shape[0]
    cdef int nt = t_array.shape[0]
    #print(nt)
    cdef double tf = t_array[nt-1]

    generic_ssa_cpp(&result[0,0], &intensity[0,0], &states[0,0], &stoich_states[0,0], &stoich_lattice[0,0], &kelong[0], &parameters[0],
                    &xi_lattice[0,0], &xi_states[0,0], &t_array[0], tf, seed, nt, n_reaction_states,
                    n_total_reactions, n_states, length, max_particles, used_frames,&probe_locations[0,0], ncolors)
    

def __original_rules():
    return original_rules_str
    
def __parsed_rules():
    return rules_str

def __min_length():
    return min_length
    
    
    
    
    
    