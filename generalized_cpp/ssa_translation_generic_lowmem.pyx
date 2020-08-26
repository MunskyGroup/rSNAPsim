# distutils: language=c++

import cython
cimport numpy as np
import numpy as np

# make the c type
#cdef extern void translationSSA(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int *SSA_result)

cdef extern from "ssa_translation_generic_lowmem_c_w.h":
    void translationSSA_generic(double* kelong, double* t_array, int* SSA_result, int N, int Nt, double* inhibitors, int seed, int fNt, int* frap_result, double* k_add, int n_enters,int n_pauses,int n_stops, int n_jumps,int* SSA_probe, int Ncolor, int Nlocs, int watched_ribs );

@cython.boundscheck(False)
@cython.wraparound(False)

def run_SSA_generic(np.ndarray[int, ndim=2, mode="c"] result not None, np.ndarray[double, ndim=1, mode="c"] kelong not None, np.ndarray[int, ndim=1, mode="c"] frap_result not None, np.ndarray[double, ndim=1, mode="c"] t_array not None, np.ndarray[double, ndim=1, mode="c"] inhibitors not None, int seed, np.ndarray[double, ndim=1, mode="c"] k_add not None, int n_enters,int n_pauses, int n_stops, int n_jumps, np.ndarray[int, ndim=2, mode="c"] probevec not None, int Ncolor, int watched_ribs ):
    """
    I need to write this. 
    """
    cdef int Nt
    cdef int N
    
    
    # subtract 2. 
    N = len(kelong+2)/3
    Nlocs = len(kelong)
    Nt = t_array.shape[0]

    fNt = t_array.shape[0]

    #inhibitors = [0,0,0]
 
    translationSSA_generic (&kelong[0], &t_array[0], &result[0,0], N, Nt, &inhibitors[0], seed,fNt,&frap_result[0],  &k_add[0], n_enters, n_pauses,  n_stops, n_jumps, &probevec[0,0], Ncolor,Nlocs, watched_ribs)

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
