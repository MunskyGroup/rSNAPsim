# distutils: language=c++

import cython
cimport numpy as np
import numpy as np

# make the c type
#cdef extern void translationSSA(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int *SSA_result)

cdef extern from "ssa_translation_c_w_lowmem_nostats.h":
    void translationSSA(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int* SSA_intensity, int N, int FRAP, int Inhibitor, double inhibit_time, int seed,  int fNt, int* frap_result, int* x0, int r_footprint, int* SSA_probe,int Ncolor)

@cython.boundscheck(False)
@cython.wraparound(False)

def run_SSA(np.ndarray[int, ndim=2, mode="c"] intensity not None, np.ndarray[double, ndim=1, mode="c"] kelong not None, np.ndarray[int, ndim=1, mode="c"] frap_result not None, np.ndarray[double, ndim=1, mode="c"] t_array not None, double kbind, double kcompl, int FRAP, int Inhibitor, double inhibit_time, int seed,  np.ndarray[int, ndim=1, mode="c"] x0 not None, int r_footprint, np.ndarray[int, ndim=2, mode="c"] probevec not None, int Ncolor):
    """
    I need to write this. 
    """
    cdef int Nt
    cdef int N
   
    
    # subtract 2. 
    N = len(kelong)
    Nt = t_array.shape[0]

    fNt = t_array.shape[0]

    translationSSA (&kelong[0], &t_array[0], Nt, kbind, kcompl, &intensity[0,0], N, FRAP, Inhibitor, inhibit_time,seed,fNt,&frap_result[0], &x0[0], r_footprint, &probevec[0,0],Ncolor)

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
