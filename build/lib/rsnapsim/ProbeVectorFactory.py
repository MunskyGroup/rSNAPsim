# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:47:28 2020

@author: willi
"""

import numpy as np
import warnings
    
class ProbeVectorFactory():
    def __init__(self):
        pass        
    
    def get_probe_vec(self, tag_epitope_dict, length):
        
        pv = np.zeros( (len(list(tag_epitope_dict)), length))
        for i in range(len(list(tag_epitope_dict))):
            pv[i,[tag_epitope_dict[list(tag_epitope_dict.keys())[i]]]] = 1
        pv = np.cumsum(pv,axis=1)        
        return pv

    def get_probe_loc(self, tag_epitope_dict, length):
        
        pv = np.zeros( (len(list(tag_epitope_dict)), length))
        for i in range(len(list(tag_epitope_dict))):
            pv[i,[tag_epitope_dict[list(tag_epitope_dict.keys())[i]]]] = 1
        return pv

    @staticmethod
    def bin_probe_vecs(probe_loc,inds):

        probeloc_binned = np.zeros((probe_loc.shape[0],   len(inds)-1 ) )
        for i in range(0,len(inds)-1):
            probeloc_binned[:,i] = np.sum(probe_loc[:,inds[i]:inds[i+1]],axis=1)
            
            
        #probeloc_binned[:,-1] = np.sum(probe_loc[:,inds[-1]:],axis=1)
        probevec_binned = np.cumsum(probeloc_binned,axis=1)
        
        return probeloc_binned.astype(int),probevec_binned.astype(int)
        
    

    @staticmethod
    def intellegent_bin(pl,nbins ,min_bin = 1):
        '''
        Function to do intellegent binning, focuses resolution on the areas
        defined in the probe location vector
        
        Note if you pass it a minium bin that when min_bin*nbins > length of your sequence
        this function will fail
        '''
        
        if min_bin*nbins > pl.shape[1]:
            warnings.warn('Desired minimum binsize and desired number of bins is not possible with the length of the probe vector, returning best guess')

        pl_inds = np.where(pl == 1)[1]
        if 0 not in pl_inds:
            pl_inds = np.hstack((np.array([0]), pl_inds))
        if pl.shape[1] not in pl_inds:
            pl_inds = np.hstack(( pl_inds,np.array(pl.shape[1]) ))
     
        used_bins = len(pl_inds)
        k = len(pl_inds)-1
        j = 0
        to_add = []
        while used_bins < nbins+1:
            if j == k:
                
                j = 0

                prev_pl_inds = pl_inds
                pl_inds = np.hstack(( pl_inds,np.array(to_add)) )
                pl_inds = np.sort(pl_inds)
             
                if np.array_equal(prev_pl_inds,pl_inds):
                    break
                k =  pl.shape[1]-1
                to_add = []
                
            newbin = int( pl_inds[j]   +  (pl_inds[j+1]-pl_inds[j]  )/2)

            if newbin not in pl_inds:
                if not (np.abs(pl_inds - newbin) <= min_bin).any():
                    to_add.append(newbin)
                    used_bins+=1
            j+=1    
            
        pl_inds = np.hstack(( pl_inds,np.array(to_add)) )
        pl_inds = np.sort(pl_inds) 
        
        return pl_inds
        
    @staticmethod
    def even_bin(length,nbins):
        '''
        Parameters
        ----------
        length : int
            Length of the vector to bin.
        nbins : int
            Number of desired bins.

        Returns
        -------
        inds : TYPE
            n bin locations over the vector length given.

        '''
        
        binsize = int(np.floor(length/nbins))
        
        k_lens = np.ones(nbins)*binsize
        
        to_redistribute = length%nbins

        k_lens[-to_redistribute:] = binsize+1
        
        inds = np.hstack(([0.], np.cumsum(k_lens))).astype(int)     
        
        return inds            
