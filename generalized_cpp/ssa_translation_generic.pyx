# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 13:14:02 2020

@author: willi
"""

# distutils: language=c++

import cython
cimport numpy as np
import numpy as np

# make the c type
#cdef extern void translationSSA(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int *SSA_result)

cdef extern from "ssa_translation_generic_c_w.h":
    void translationSSA_generic(double *kelong, double *t_array, int *SSA_result, int N,int Nt, double* inhibitors, int seed, double *SSA_ribtimes, int *nribs,int ribtimesize, int fNt, int* frap_result, int cNt, int* col_result, double* k_add, int n_enters,int n_pauses,int n_stops, int n_jumps  )

@cython.boundscheck(False)
@cython.wraparound(False)

def run_SSA_generic(np.ndarray[int, ndim=1, mode="c"] result not None, np.ndarray[double, ndim=1, mode="c"] ribtimes not None, np.ndarray[int, ndim=1, mode="c"] coltimes not None,np.ndarray[double, ndim=1, mode="c"] kelong not None, np.ndarray[int, ndim=1, mode="c"] frap_result not None, np.ndarray[double, ndim=1, mode="c"] t_array not None, np.ndarray[double, ndim=1, mode="c"] inhibitors not None, int seed, np.ndarray[int, ndim=1, mode="c"] nribs not None,np.ndarray[double, ndim=1, mode="c"] k_add not None, int n_enters,int n_pauses, int n_stops, int n_jumps ):
    """
    I need to write this. 
    """
    cdef int Nt
    cdef int N
    

    # subtract 2. 
    N = len(kelong)/3
   
    Nt = t_array.shape[0]
    ribtimesize = ribtimes.shape[0]
    coltimesize = coltimes.shape[0]
    fNt = t_array.shape[0]
    cNt = t_array.shape[0]
    #inhibitors = [0,0,0]
    
    translationSSA_generic (&kelong[0], &t_array[0], &result[0], N, Nt, &inhibitors[0], seed, &ribtimes[0],&nribs[0], ribtimesize,fNt,&frap_result[0],coltimesize,&coltimes[0],  &k_add[0], n_enters, n_pauses,  n_stops, n_jumps)

    return None

## Numpy must be initialized. When using numpy from C or Cython you must
## _always_ do that, or you will have segfaults
#np.import_array()
#
#
#cpdef translation_ssa():
#    cdef np.ndarray[int, ndim=1, mode='c'] a
#
#    a = np.zeros((4,), dtype=np.int32)
#    set_integer_arr_ptr(&a[0])
#    return a
