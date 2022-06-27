# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:45:31 2020

@author: willi
"""

import warnings
import numpy as np
from . import CodonDictionaries
from . import SequenceManipMethods
from . import custom_errors as custom_err

class PropensityFactory():
    '''
    factory class for the k's
    '''
    def __init__(self):
        self.codon_dicts = CodonDictionaries.CodonDictionaries()

    @staticmethod
    def bin_k(k, inds):
        '''
        Generate a k_binned from a given indexing vector

        Parameters
        ----------
        k : array like
            propensity list or ndarray.
        inds : ndarray
            numpy array of index locations, can be generated from even_bin()
            or intellegent_bin().

        Returns
        -------
        k_binned : list
            propensity rates for the proposed binning strategy.

        '''
        try:
            k = k.tolist()
        except:
            pass
        k_binned = np.zeros(len(inds)-1)
        binned_ks = []
        for i in range(0, len(inds)-1):
            binned_ks = binned_ks + [k[inds[i]:inds[i+1]],]


        for i in range(0, len(inds)-1):
            k_binned[i] = 1/ np.sum(1/np.array(binned_ks[i]))
        return k_binned


    def get_trna_ids(self, nt_seq):
        '''
        Return tRNA ids 0-60 for a given nt_seq

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.

        Returns
        -------
        list
            list of trna IDs 0-60 (one for each codon).

        '''

        if '*' in SequenceManipMethods.SequenceManipMethods().nt2aa(nt_seq)[-1]:
            nt_seq = nt_seq[:-3]

        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)]
        return [self.codon_dicts.trna_dict[x] for x in seperated_codons]


    def get_k_PA(self, nt_seq, k_init, k_elong_mean, k_term, PA_dict=None):
        '''
        Convert a list of Peptide, Acceptor duplet combo dwell times to 
        an elongation vector. This will lower the "Length" of the mRNA lattice
        by one codons since it does not account for the last 3 nucleotides.
        
        Example
        -------
        
        given the sequence AAACCCGGGUUU and the dictionary:
        
        {AAACCC : 1.0
         CCCGGG : 1.5
         GGGUUU : 2.0}
        
        returns: [1.0, 1.5, 2.0]
        
        
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence to convert to values.
        k_init : float
            initiation rate.
        k_elong_mean : float
            average mean elongation rate to scale the construct too.
        k_term : float
            termination rate.
        PA_dict : dict, optional
            Provided dictionary that maps nt triplets to values. 
            The default is the rates from Gobet 2020.
            [https://www.pnas.org/doi/full/10.1073/pnas.1918145117#supplementary-materials]

        Raises
        ------
        UnrecognizedCodonError 
            Error if there is a missing codon-to-value entry, unregonized substring in the nucleotide sequence.

        Returns
        -------
        all_k : list
            A list of elongation rates of size (L + 1) (initation, L-1 elongation rates, termination)

        '''        



        if PA_dict is None:
            codon_dict = self.codon_dicts.Gobet2020_PA_rates
        else:
            codon_dict = PA_dict
        codons = nt_seq.upper()

        genelength = int(len(codons)/3)
        seperated_codons = [codons[i:i+6] for i in range(0, len(codons), 3)][:-1] #split codons by 3
        k_elongation = np.zeros((1, genelength))
        codon_freq = np.zeros((1, genelength))
        

        for i in range(len(seperated_codons)):
            try:
                codon_freq[0, i] = codon_dict[seperated_codons[i]]
            except KeyError:
                msg = 'Unrecognized codon to optimize: "%s", at index %i provide an entry for this'\
                    ' codon and its amino acid in '\
                    'the codon to value dictionary.'%(seperated_codons[i],i) 
                raise custom_err.UnrecognizedCodonError(msg)
            

        mean_value = np.mean(list(codon_dict.values()))

        k_elongation = (codon_freq / mean_value) * k_elong_mean
        all_k = [k_init] + k_elongation.flatten().tolist() + [k_term]

        return all_k


    def get_k_EPA(self, nt_seq, k_init, k_elong_mean, k_term, EPA_dict=None):
        '''
        Convert a list of Exit, Peptide, Acceptor triplet combo dwell times to 
        an elongation vector. This will lower the "Length" of the mRNA lattice
        by two codons since it does not account for the last 6 nucleotides.
        
        Example
        -------
        
        given the sequence AAACCCGGGUUU and the dictionary:
        
        {AAACCCGGG : 1.0
         CCCGGGUUU : 1.5}
        
        returns: [1.0, 1.5]
        
        
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence to convert to values.
        k_init : float
            initiation rate.
        k_elong_mean : float
            average mean elongation rate to scale the construct too.
        k_term : float
            termination rate.
        EPA_dict : dict, optional
            Provided dictionary that maps nt triplets to values. The default is Blank.

        Raises
        ------
        UnrecognizedCodonError 
            Error if there is a missing codon-to-value entry, unregonized substring in the nucleotide sequence.

        Returns
        -------
        all_k : list
            list of elongation rates of size (L) (initation, L-2 elongation rates, termination)

        '''

        if EPA_dict is None:
            codon_dict = None
        else:
            codon_dict = EPA_dict
        codons = nt_seq.upper()

        genelength = int(len(codons)/3)
        seperated_codons = [codons[i:i+9] for i in range(0, len(codons), 3)][:-2] #split codons by 3
        k_elongation = np.zeros((1, genelength))
        codon_freq = np.zeros((1, genelength))
        

        for i in range(len(seperated_codons)):
            try:
                codon_freq[0, i] = codon_dict[seperated_codons[i]]
            except KeyError:
                msg = 'Unrecognized codon to optimize: "%s", at index %i provide an entry for this'\
                    ' codon and its amino acid in '\
                    'the codon to value dictionary.'%(seperated_codons[i],i) 
                raise custom_err.UnrecognizedCodonError(msg)
            

        mean_value = np.mean(list(codon_dict.values()))

        k_elongation = (codon_freq / mean_value) * k_elong_mean
        all_k = [k_init] + k_elongation.flatten().tolist() + [k_term]

    def get_k(self, nt_seq, k_init, k_elong_mean,
              k_term, codon_freq_bias_dict=None):
        '''
        Returns all propensity functions for a given kinit, kelong_mu, and k termination

        Parameters
        ----------
        nt_seq : string
            nucleotide sequence.
        k_init : float
            initation rate.
        k_elong_mean : float
            average kelongation rate for the transcript.
        k_term : float
            termination rate.
            
        Raises
        ------
        UnrecognizedCodonError 
            Error if there is a missing codon-to-value entry, unregonized substring in the nucleotide sequence.


        Returns
        -------
        all_k : list
            list of propensity values per location + initiation and termination at respective ends.

        '''

        if codon_freq_bias_dict is None:
            codon_dict = self.codon_dicts.human_codon_frequency_bias_nakamura
        else:
            codon_dict = codon_freq_bias_dict
        codons = nt_seq.upper()

        genelength = int(len(codons)/3)
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3

        
        k_elongation = np.zeros((1, genelength))
        codon_freq = np.zeros((1, genelength))
        

        for i in range(len(seperated_codons)):
            try:
                codon_freq[0, i] = codon_dict[seperated_codons[i]]
            except KeyError:
                msg = 'Unrecognized codon to optimize: "%s", at index %i provide an entry for this'\
                    ' codon and its amino acid in '\
                    'the codon to value dictionary.'%(seperated_codons[i],i) 
                raise custom_err.UnrecognizedCodonError(msg)
            

        mean_value = np.mean(list(codon_dict.values()))

        k_elongation = (codon_freq / mean_value) * k_elong_mean
        all_k = [k_init] + k_elongation.flatten().tolist() + [k_term]

        return all_k


    def get_k_3_frame(self, nt_seq, k_elong_mean,
                      codon_freq_bias_dict=None):
        '''
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.
        k_elong_mean : float
            mean elongation rate aa/time.

        Returns
        -------
        kelongs : list
            elongation rates for all 3 frames of a nucleotide sequence stacked in:
                [0+frame L, 1+frame L-1, 2+frame L-1] .

        '''
        kelongs = []

        if codon_freq_bias_dict is None:
            codon_dict = self.codon_dicts.human_codon_frequency_bias_nakamura

        for n in range(3):
            if n != 0:
                codons = nt_seq[n:-(3-n)]
            else:
                codons = nt_seq
            genelength = int(len(codons)/3)
            seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
            k_elongation = np.zeros((1, genelength))
            tRNA_copynumber = np.zeros((1, genelength))


            for i in range(len(seperated_codons)):
                tRNA_copynumber[0, i] = codon_dict[seperated_codons[i]]

            mean_tRNA_copynumber = np.mean(list(codon_dict.values()))

            k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean

           # k_elongation.flatten().tolist()[:-1]

            kelongs = kelongs + k_elongation.flatten().tolist()

        return kelongs


    @staticmethod
    def get_binned_k(k, bins):
        '''
        given a propensity vector and an desired binning indices, return pv and pl binned

        Parameters
        ----------
        k : 1 x L  numpy array
            Propensity vector of rates.
        bins : 1xL ind locations sizes
            binning strategy, the sum of this vector should = L.

        Returns
        -------
        propensity_binned: 1 x Nbins numpy array
            the binned propensities.
        bin_sizes: 1 x Nbins numpy array
            list of bin sizes.

        '''
        binsize = int(np.floor(len(k)/bins))
        binned_ks = []

        k_binned = np.zeros(bins)
        k_lens = np.ones(bins)*binsize

        to_redistribute = len(k)%bins

        k_lens[-to_redistribute:] = binsize+1

        inds = np.hstack(([0.], np.cumsum(k_lens))).astype(int)


        for i in range(0, bins):
            binned_ks = binned_ks + [k[inds[i]:inds[i+1]].tolist(),]

        for i in range(0, bins):
            k_binned[i] = np.mean(binned_ks[i])/len(binned_ks[i])

        return k_binned, k_lens


    @staticmethod
    def intellegent_bin(pl, nbins, min_bin=1):
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
            warnings.warn('Desired minimum binsize and desired number'\
                          ' of bins is not possible with the length of'\
                              ' the probe vector, returning best guess')

        pl_inds = np.where(pl == 1)[1]

        if 0 not in pl_inds:
            pl_inds = np.hstack((np.array([0]), pl_inds))
        if pl.shape[1]+1 not in pl_inds:
            pl_inds = np.hstack((pl_inds, np.array(pl.shape[1]+1)))

        used_bins = len(pl_inds)
        k = len(pl_inds)-1
        j = 0
        to_add = []
        while used_bins < nbins+1:
            if j == k:

                j = 0

                prev_pl_inds = pl_inds
                pl_inds = np.hstack((pl_inds, np.array(to_add)))
                pl_inds = np.sort(pl_inds)

                if np.array_equal(prev_pl_inds, pl_inds):
                    break
                k = len(pl_inds)-1
                to_add = []

            newbin = int(pl_inds[j] + (pl_inds[j+1]-pl_inds[j])/2)

            if newbin not in pl_inds:
                if not (np.abs(pl_inds - newbin) <= min_bin).any():
                    to_add.append(newbin)
                    used_bins += 1
            j += 1

        pl_inds = np.hstack((pl_inds, np.array(to_add)))
        pl_inds = np.sort(pl_inds)

        return pl_inds

    @staticmethod
    def even_bin(length, nbins):
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
    