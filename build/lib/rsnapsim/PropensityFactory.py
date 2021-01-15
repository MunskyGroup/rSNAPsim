# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:45:31 2020

@author: willi
"""

from . import CodonDictionaries
import numpy as np
import warnings

class PropensityFactory():
    '''
    factory class for the k's
    '''
    def __init__(self):
        self.codon_dicts = CodonDictionaries.CodonDictionaries()
        pass    

    @staticmethod
    def bin_k(k,inds):
        try: 
            k = k.tolist()
        except:
            pass
        k_binned = np.zeros(len(inds)-1)
        binned_ks = []
        for i in range(0,len(inds)-1):
            binned_ks = binned_ks +   [k[inds[i]:inds[i+1]],]  
        
         
        for i in range(0,len(inds)-1):            
            k_binned[i] = 1/ np.sum( 1/np.array(binned_ks[i]))
        return k_binned


    def get_trna_ids(self,nt_seq):
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] 
        return [self.codon_dicts.trna_dict[x] for x in seperated_codons]
        
        

    def get_k(self, nt_seq, k_init, k_elong_mean,k_end):
        '''
        returns all propensities for a given nucleotide sequence

        *args*

            **nt_seq**,  nucleotide sequence as a string

            **k_initiation**, initiation rate of ribosome binding

            **k_elong_mean**, average rate of elgonation experimentally found


        '''
        codons = nt_seq.upper()
        
        genelength = int(len(codons)/3)
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
        k_elongation = np.zeros((1, genelength))
        tRNA_copynumber = np.zeros((1, genelength))

     
        for i in range(len(seperated_codons)):
            tRNA_copynumber[0, i] = self.codon_dicts.strGeneCopy[seperated_codons[i]]

        mean_tRNA_copynumber = np.mean(list(self.codon_dicts.strGeneCopy_single.values()))

        k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean
        all_k = [k_init] + k_elongation.flatten().tolist() + [k_end]
        
        return all_k
    
    
    def get_k_3_frame(self,nt_seq,k_elong_mean):
        '''
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.
        k_elong_mean : float
            mean elongation rate aa/s.

        Returns
        -------
        kelongs : list
            elongation rates for all 3 frames of a nucleotide sequence stacked in:
                [0+frame L, 1+frame L-1, 2+frame L-1] .

        '''
        kelongs = []
        
        for n in range(3):
            if n !=0:
                codons = nt_seq[n:-(3-n)]
            else:
                codons = nt_seq
            genelength = int(len(codons)/3)
            seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
            k_elongation = np.zeros((1, genelength))
            tRNA_copynumber = np.zeros((1, genelength))

     
            for i in range(len(seperated_codons)):
                tRNA_copynumber[0, i] = self.codon_dicts.strGeneCopy[seperated_codons[i]]
    
            mean_tRNA_copynumber = np.mean(list(self.codon_dicts.strGeneCopy_single.values()))
    
            k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean
        
           # k_elongation.flatten().tolist()[:-1]
        
            kelongs = kelongs + k_elongation.flatten().tolist()
        
        return kelongs
        
    
    @staticmethod
    def get_binned_k(k,bins):
        '''
        evenly bins elongation rates as best it can.
        
        '''
        binsize = int(np.floor(len(k)/bins))
        binned_ks = []
        
        k_binned = np.zeros(bins)
        k_lens = np.ones(bins)*binsize
        
        to_redistribute = len(k)%bins

        k_lens[-to_redistribute:] = binsize+1
        
        inds = np.hstack(([0.], np.cumsum(k_lens))).astype(int)     

        
        for i in range(0,bins):
            binned_ks = binned_ks +   [k[inds[i]:inds[i+1]].tolist(),]  
         
        for i in range(0,bins):            
            k_binned[i] = np.mean(binned_ks[i])/len(binned_ks[i])

        return k_binned,k_lens
    
    
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
        if pl.shape[1]+1 not in pl_inds:
            pl_inds = np.hstack(( pl_inds,np.array(pl.shape[1]+1) ))
     
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
                k =  len(pl_inds)-1
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
        evenly bins a length over a given amount of bins as best it can
        
        '''
        binsize = int(np.floor(length/nbins))
        
        k_lens = np.ones(nbins)*binsize
        
        to_redistribute = length%nbins

        k_lens[-to_redistribute:] = binsize+1
        
        inds = np.hstack(([0.], np.cumsum(k_lens))).astype(int)     
        
        return inds