# distutils: language=c++

import cython
cimport numpy as np
import numpy as np

# make the c type
#cdef extern void translationSSA(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int *SSA_result)

cdef extern from "ssa_translation_c_w_lowmem_leaky.h":
    void translationSSA(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int* SSA_intensity, int N, int FRAP, int Inhibitor, double inhibit_time, int seed, double *SSA_ribtimes, int *nribs,int ribtimesize, int fNt, int* frap_result, int cNt, int* col_result, double* col_t, int* col_x, int colNp, int* x0, int r_footprint, int* SSA_probe,double k_probe,int* probe_loc,int n_probes)

@cython.boundscheck(False)
@cython.wraparound(False)

def run_SSA(np.ndarray[int, ndim=1, mode="c"] intensity not None, np.ndarray[double, ndim=1, mode="c"] ribtimes not None, np.ndarray[int, ndim=1, mode="c"] coltimes not None, np.ndarray[int, ndim=1, mode="c"] col_x not None, np.ndarray[double, ndim=1, mode="c"] col_t not None ,np.ndarray[double, ndim=1, mode="c"] kelong not None, np.ndarray[int, ndim=1, mode="c"] frap_result not None, np.ndarray[double, ndim=1, mode="c"] t_array not None, double kbind, double kcompl, int FRAP, int Inhibitor, double inhibit_time, int seed, np.ndarray[int, ndim=1, mode="c"] nribs not None,  np.ndarray[int, ndim=1, mode="c"] x0 not None,int r_footprint, np.ndarray[int, ndim=1, mode="c"] probevec not None,double k_probe,np.ndarray[int, ndim=1, mode="c"] probeloc not None ):
    """
    I need to write this. 
    """
    cdef int Nt
    cdef int N
    cdef int n_probes
    
    # subtract 2. 
    N = len(kelong)
    Nt = t_array.shape[0]
    ribtimesize = ribtimes.shape[0]
    coltimesize = coltimes.shape[0]
    fNt = t_array.shape[0]
    cNt = t_array.shape[0]
    colNp = col_x.shape[0]
    n_probes = len(probeloc)
    
    translationSSA (&kelong[0], &t_array[0], Nt, kbind, kcompl, &intensity[0], N, FRAP, Inhibitor, inhibit_time,seed,&ribtimes[0],&nribs[0], ribtimesize,fNt,&frap_result[0],coltimesize,&coltimes[0], &col_t[0],&col_x[0],colNp, &x0[0], r_footprint, &probevec[0],k_probe,&probeloc[0],n_probes)

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
