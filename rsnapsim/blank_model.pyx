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

#model id goes here

#parsed rules go here


#original rules go here


#cdef goes here
    void generic_ssa_cpp(int* particle_array, int* state_array, int* resource_array,  # arrays to fill
                    int* particle_x0, int* state_x0,  int* resource_x0, int* rxns, int* rxnids, double* elongs, int* probes,
                    double* parameters, int npars,       # initial simulation state
                    int max_rib, int n_constant_reactions, int n_ribosome_reactions,  # constants of the simulation
                    int n_colors, int n_states, int n_resources, int Nt, double* time_vector, double tf, double burnin, int seed, int L, int error_code);
    
    
    

    
def run_ssa_cpp( np.ndarray[int, ndim=2, mode="c"] particle_array not None,
                 np.ndarray[int, ndim=2, mode="c"] state_array not None,
                 np.ndarray[int, ndim=2, mode="c"] resource_array not None,
                 np.ndarray[int, ndim=2, mode="c"] particle_x0 not None,
                 np.ndarray[int, ndim=1, mode="c"] state_x0 not None,
                 np.ndarray[int, ndim=1, mode="c"] resource_x0 not None,
                 np.ndarray[int, ndim=2, mode="c"] rxn_mat not None,
                 np.ndarray[int, ndim=1, mode="c"] rxnids not None,
                 np.ndarray[double, ndim=2, mode="c"] elong_mat not None,
                 np.ndarray[int, ndim=2, mode="c"] probes not None,
                 np.ndarray[double, ndim=1, mode="c"] parameters not None,
                 np.ndarray[double, ndim=1, mode="c"] time_vector not None,
                 int max_rib, int n_states, int n_resources, int n_constant_reactions, int n_ribosome_reactions, 
                 double burnin, int seed):
    
    
    cdef int n_colors = probes.max()
    cdef int n_parameters = len(parameters)
    cdef int Nt = time_vector.shape[0]
    cdef double tf = time_vector[Nt-1]
    cdef int n_rxns = n_constant_reactions + n_ribosome_reactions
    cdef int L = int(elong_mat.shape[1])
    cdef int error_code = 0
    

    generic_ssa_cpp(&particle_array[0,0], &state_array[0,0], &resource_array[0,0], &particle_x0[0,0], &state_x0[0], &resource_x0[0],
                    &rxn_mat[0,0], &rxnids[0], &elong_mat[0,0], &probes[0,0], &parameters[0],  n_parameters,
                    max_rib, n_constant_reactions, n_ribosome_reactions, n_colors, n_states, n_resources, Nt,  &time_vector[0],
                    tf, burnin, seed, L, error_code)
    return error_code


def __original_rules():
    return original_rules_str
    
def __parsed_rules():
    return rules_str

def __model_id():
    return model_id
    
    
    
    
    