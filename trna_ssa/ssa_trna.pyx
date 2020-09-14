# distutils: language=c++

import cython
cimport numpy as np
import numpy as np

# make the c type
#cdef extern void translationSSA(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int *SSA_result)

cdef extern from "ssa_trna_c_w.h":
    void translationSSA_trna(int *k_index, double *k_trna, double k_diffusion, double *t_array, int Nt, double kbind, double kcompl, int *SSA_result,int *trna_result, int N, int FRAP, int Inhibitor, double inhibit_time, int seed, double *SSA_ribtimes, int *nribs,int ribtimesize, int fNt, int* frap_result, int cNt, int* col_result, double* col_t, int* col_x, int colNp, int* x0, double kelong)

@cython.boundscheck(False)
@cython.wraparound(False)

def run_SSA(np.ndarray[int, ndim=1, mode="c"] result not None, np.ndarray[int, ndim=1, mode="c"] trna_result not None, np.ndarray[double, ndim=1, mode="c"] ribtimes not None, np.ndarray[int, ndim=1, mode="c"] coltimes not None, np.ndarray[int, ndim=1, mode="c"] col_x not None, np.ndarray[double, ndim=1, mode="c"] col_t not None ,np.ndarray[int, ndim=1, mode="c"] k_index not None, np.ndarray[double, ndim=1, mode="c"] k_trna not None, double k_diffusion, np.ndarray[int, ndim=1, mode="c"] frap_result not None, np.ndarray[double, ndim=1, mode="c"] t_array not None, double kbind, double kcompl, int FRAP, int Inhibitor, double inhibit_time, int seed, np.ndarray[int, ndim=1, mode="c"] nribs not None,  np.ndarray[int, ndim=1, mode="c"] x0 not None, double kelong):
    """
    I need to write this. 
    """
    cdef int Nt
    cdef int N
    
    # subtract 2. 
    N = len(k_index)
    Nt = t_array.shape[0]
    ribtimesize = int(ribtimes.shape[0])
    coltimesize = int(coltimes.shape[0])
    fNt = int(t_array.shape[0])
    cNt = int(t_array.shape[0])
    colNp = int(col_x.shape[0])
    
    
    translationSSA_trna(&k_index[0],&k_trna[0],k_diffusion, &t_array[0], Nt, kbind, kcompl, &result[0], &trna_result[0], N, FRAP, Inhibitor, inhibit_time,seed,&ribtimes[0],&nribs[0], ribtimesize,fNt,&frap_result[0],coltimesize,&coltimes[0], &col_t[0],&col_x[0],colNp, &x0[0], kelong)

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
