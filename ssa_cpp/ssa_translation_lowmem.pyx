# distutils: language=c++

import cython
cimport numpy as np
import numpy as np


# make the c type

#cdef extern from "ssa_translation_c_w.h":
    #void translationSSA(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int *SSA_result, int N, int FRAP, int Inhibitor, double inhibit_time, int seed, double *SSA_ribtimes, int *nribs,int ribtimesize, int fNt, int* frap_result, int cNt, int* col_result, double* col_t, int* col_x, int colNp, int* x0, int r_footprint, int rib_max )

cdef extern from "ssa_translation_c_w_lowmem.h":
    void translationSSA_lowmem(double *kelong, double *t_array, int Nt, double kbind, double kcompl, int* SSA_intensity, int N, int FRAP, int Inhibitor, double inhibit_time, int seed, double *SSA_ribtimes, int *nribs,int ribtimesize, int fNt, int* frap_result, int cNt, int* col_result, double* col_t, int* col_x, int colNp, int* x0, int r_footprint, int* SSA_probe,int Ncolor, int* flags, double kon, double koff,double* k_probe,int* probe_loc,int* n_probes)


@cython.boundscheck(False)
@cython.wraparound(False)


# cdef class PySSA:
    
#     cdef int c_leaky
#     cdef int c_bursting
#     cdef int c_lowmem
#     cdef int c_stats
    
#     cdef int footprint
#     cdef int rib_max
    
    
#     def __cinit__(self):
#         self.c_leaky = 0
#         self.c_bursting = 0
#         self.c_stats = 0
#         self.c_lowmem = 1
#         self.footprint = 9
#         self.rib_max = 200
        

#     def set_flags(self,leaky, bursting, stats, lowmem):
#         self.c_leaky = leaky
#         self.c_bursting = bursting
#         self.c_stats = stats
#         self.c_lowmem = lowmem
        
#     def set_footprints(self,r_footprint,rib_max):
#         self.footprint = r_footprint
#         self.rib_max = rib_max


def run_SSA(np.ndarray[int, ndim=2, mode="c"] intensity not None, np.ndarray[double, ndim=1, mode="c"] ribtimes not None, np.ndarray[int, ndim=1, mode="c"] coltimes not None, np.ndarray[int, ndim=1, mode="c"] col_x not None, np.ndarray[double, ndim=1, mode="c"] col_t not None ,np.ndarray[double, ndim=1, mode="c"] kelong not None, np.ndarray[int, ndim=1, mode="c"] frap_result not None, np.ndarray[double, ndim=1, mode="c"] t_array not None, double kbind, double kcompl, int FRAP, int Inhibitor, double inhibit_time, int seed, np.ndarray[int, ndim=1, mode="c"] nribs not None,  np.ndarray[int, ndim=1, mode="c"] x0 not None,int r_footprint, np.ndarray[int, ndim=2, mode="c"] probevec not None, int Ncolor, double kon, double koff, np.ndarray[double, ndim=1, mode="c"] k_probe not None,np.ndarray[int, ndim=2, mode="c"] probeloc not None, np.ndarray[int, ndim=1, mode="c"] flags not None ):
    """
    I need to write this. 
    """       
    
    cdef int N
    cdef int Nt
    
    cdef np.ndarray[int, ndim=1, mode="c"] n_probes = np.zeros([Ncolor],dtype=np.int32)
    # subtract 2. 
    N = len(kelong)
    Nt = t_array.shape[0]
    ribtimesize = ribtimes.shape[0]
    coltimesize = coltimes.shape[0]
    fNt = t_array.shape[0]
    cNt = t_array.shape[0]
    colNp = col_x.shape[0]
    
    rib_max = 200
    for i in range(0,Ncolor):
      
        n_probes[i] = len(np.where(probeloc[i] == 1)[0]) 
    
    
    pind = np.array(np.where(probeloc== 1)).T.astype(np.int32)
    
    cdef np.ndarray[int, ndim=2, mode="c"] probeloc_indexes = np.zeros(list(pind.shape),dtype=np.int32)
 
    for i in range(0,pind.shape[0]):
        for j in range(0,pind.shape[1]):
            probeloc_indexes[i,j] = pind[i,j]
    
    translationSSA_lowmem (&kelong[0], &t_array[0], Nt, kbind, kcompl, &intensity[0,0], N, FRAP, Inhibitor, inhibit_time,seed,&ribtimes[0],&nribs[0], ribtimesize,fNt,&frap_result[0],coltimesize,&coltimes[0], &col_t[0],&col_x[0],colNp, &x0[0], r_footprint, &probevec[0,0],Ncolor, &flags[0], kon, koff, &k_probe[0], &probeloc_indexes[0,0],&n_probes[0])

    return None


# def run_SSA(np.ndarray[int, ndim=1, mode="c"] result not None, np.ndarray[double, ndim=1, mode="c"] ribtimes not None, np.ndarray[int, ndim=1, mode="c"] coltimes not None, np.ndarray[int, ndim=1, mode="c"] col_x not None, np.ndarray[double, ndim=1, mode="c"] col_t not None ,np.ndarray[double, ndim=1, mode="c"] kelong not None, np.ndarray[int, ndim=1, mode="c"] frap_result not None, np.ndarray[double, ndim=1, mode="c"] t_array not None, double kbind, double kcompl, int FRAP, int Inhibitor, double inhibit_time, int seed, np.ndarray[int, ndim=1, mode="c"] nribs not None,  np.ndarray[int, ndim=1, mode="c"] x0 not None,int r_footprint, int rib_max ):
#     """
#     I need to write this. 
#     """
#     cdef int Nt
#     cdef int N
    
#     # subtract 2. 
#     N = len(kelong)
#     Nt = t_array.shape[0]
#     ribtimesize = ribtimes.shape[0]
#     coltimesize = coltimes.shape[0]
#     fNt = t_array.shape[0]
#     cNt = t_array.shape[0]
#     colNp = col_x.shape[0]
    
    
#     translationSSA (&kelong[0], &t_array[0], Nt, kbind, kcompl, &result[0], N, FRAP, Inhibitor, inhibit_time,seed,&ribtimes[0],&nribs[0], ribtimesize,fNt,&frap_result[0],coltimesize,&coltimes[0], &col_t[0],&col_x[0],colNp, &x0[0], r_footprint, rib_max)

#     return None

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
