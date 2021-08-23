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
        '''
        Generate a probe vector given a tag epitope dictionary and the Length of the transcript:
            
            Ex:   get_probe_bec({Tag_1:[10,15]},20)
            returns:
                pv = np.array([[0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]])
        

        Parameters
        ----------
        tag_epitope_dict : dictionary
            dictionary of tag epitopes, can be accessed from a POI.tag_epitopes.
        length : int
            Length of the transcript.

        Returns
        -------
        probe_vec : Ncolor x L Numpy array
            returns the cumulative sum probe vector per location of a ribosome, meaning: A ribosome at x location
            will have pv[x] intensity.

        '''
        pv = np.zeros( (len(list(tag_epitope_dict)), length))
        for i in range(len(list(tag_epitope_dict))):
            pv[i,[tag_epitope_dict[list(tag_epitope_dict.keys())[i]]]] = 1
        pv = np.cumsum(pv,axis=1)        
        return pv

    def get_probe_loc(self, tag_epitope_dict, length):
        '''
        Get 1, 0 probe location from a tag epitope dictionary

            Ex:   get_probe_bec({Tag_1:[10,15]},20)
            returns:
                pv = np.array([[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0]])
        

        Parameters
        ----------
        tag_epitope_dict : dictionary
            dictionary of tag epitopes, can be accessed from a POI.tag_epitopes.
        length : int
            Length of the transcript.

        Returns
        -------
        probe_loc : Ncolor x L Numpy array
            returns the binary probe per location of an epitope,
        '''
        
        pv = np.zeros( (len(list(tag_epitope_dict)), length))
        for i in range(len(list(tag_epitope_dict))):
            pv[i,[tag_epitope_dict[list(tag_epitope_dict.keys())[i]]]] = 1
        return pv

    @staticmethod
    def bin_probe_vecs(probe_loc,inds):
        '''
        given a probe location vector and an desired binning indices, return pv and pl binned

        Parameters
        ----------
        probe_loc : Ncolor by L binary numpy array
            locations of tag epitopes.
        inds : 1xL ind locations sizes
            binning strategy, the sum of this vector should = L:.

        Returns
        -------
        Probe_loc_binned: Ncolor x L numpy array
            the binned probe locations.
        Probe_vec_binned: Ncolor x L numpy array
            the binned probe intensity vector.
        '''
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

        Parameters
        ----------
        pl : numpy array
            numpy array of Ncolor x Length 0's 1's for probe locations.
        nbins : int
            number of desired bins.
        min_bin : int, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        pl_inds : numpy array
            n bin locations over the vector length given.

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
        inds : numpy array
            n bin locations over the vector length given.

        '''
        
        binsize = int(np.floor(length/nbins))
        
        k_lens = np.ones(nbins)*binsize
        
        to_redistribute = length%nbins

        k_lens[-to_redistribute:] = binsize+1
        
        inds = np.hstack(([0.], np.cumsum(k_lens))).astype(int)     
        
        return inds            
