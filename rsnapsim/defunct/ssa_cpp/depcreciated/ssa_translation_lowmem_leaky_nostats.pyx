# distutils: language=c++

import cython
cimport numpy as np
import numpy as np

# make the c type
#cdef extern void translationSSA(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int *SSA_result)

cdef extern from "ssa_translation_c_w_lowmem_leaky_nostats.h":
    void translationSSA(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int* SSA_intensity, int N, int FRAP, int Inhibitor, double inhibit_time, int seed, int fNt, int* frap_result, int* x0, int r_footprint, int* SSA_probe,double *k_probe,int* probe_loc,int* n_probes, int Ncolor)

@cython.boundscheck(False)
@cython.wraparound(False)

def run_SSA(np.ndarray[int, ndim=2, mode="c"] intensity not None, np.ndarray[double, ndim=1, mode="c"] kelong not None, np.ndarray[int, ndim=1, mode="c"] frap_result not None, np.ndarray[double, ndim=1, mode="c"] t_array not None, double kbind, double kcompl, int FRAP, int Inhibitor, double inhibit_time, int seed, np.ndarray[int, ndim=1, mode="c"] x0 not None,int r_footprint, np.ndarray[int, ndim=2, mode="c"] probevec not None, np.ndarray[double, ndim=1, mode="c"] k_probe not None,np.ndarray[int, ndim=2, mode="c"] probeloc not None, int Ncolor ):
    """
    I need to write this. 
    """
    cdef int Nt
    cdef int N
    #cdef list n_probes
 
    
    cdef np.ndarray[int, ndim=1, mode="c"] n_probes = np.zeros([Ncolor],dtype=np.int32)
    
    # subtract 2. 
    N = len(kelong)
    Nt = t_array.shape[0]
    fNt = t_array.shape[0]

    #n_probes = []
    for i in range(0,Ncolor):
      
        n_probes[i] = len(np.where(probeloc[i] == 1)[0]) 
        
    pind = np.array(np.where(probeloc== 1)).T.astype(np.int32)
    
    cdef np.ndarray[int, ndim=2, mode="c"] probeloc_indexes = np.zeros(list(pind.shape),dtype=np.int32)
 
    for i in range(0,pind.shape[0]):
        for j in range(0,pind.shape[1]):
            probeloc_indexes[i,j] = pind[i,j]
    

    translationSSA (&kelong[0], &t_array[0], Nt, kbind, kcompl, &intensity[0,0], N, FRAP, Inhibitor, inhibit_time,seed,fNt,&frap_result[0], &x0[0], r_footprint, &probevec[0,0],&k_probe[0],&probeloc_indexes[0,0],&n_probes[0],Ncolor)

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
