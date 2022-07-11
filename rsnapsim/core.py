# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:35:43 2022

@author: willi
"""
import re
import time
import os
import numpy as np
import itertools as it
from . import CodonDictionaries
from . import AuxCodonDicts
from . import FileParser
from . import custom_errors as custom_err
import warnings

try:
    from Bio import SeqIO
    from Bio import Entrez
except:

    print('BioPython is not installed, polling genbank will not be possible')
    pass



# Any function that does analysis on a given sequence string
class SequenceCore:
    '''
    class that handles manipulation methods dealing with sequences
    '''
    def __init__(self, sequence=''):
        self.sequence = sequence
        self.codon_dicts = CodonDictionaries.CodonDictionaries()
        self.aux_dicts = AuxCodonDicts.AuxCodonDicts()
        self.previous_files = {}
        
        self.letter_dict = ['a', 'c', 'u', 'g']
        # IPUAC SUBSTITIONS - They need to  adjusted /
        #extended for mrna modifications?
        self.sub_dict = {'m':'a', 'r':'a', 'y':'t', 'k':'g', 's':'g',
                         'w':'a', 'h':'a', 'n':'a', 'Ψ':'u'}
        #get the codon dictionaries


    def optimize_ntseq(self, nt_seq, opt_dict=None, nt2aa_dict=None):
        '''
        Optimizes a nucleotide sequence

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence string
        opt_dict : dictionary, optional
            a user defined dictionary to optimize over. The default is None.
        nt2aa_dict : dictionary, optional
            a user defined dictionary to convert nucleotide codons to amino acids. The default is None.
        Returns
        -------
        opt_seq : str
            Optimized NT sequenced based on a given dictionary of rates

        '''

        if opt_dict is None:
            opt_dict = self.codon_dicts.human_codon_frequency_bias_nakamura
        if nt2aa_dict is None:
            nt2aa_dict = self.codon_dicts.aa_table
            aa2nt_dict = self.codon_dicts.aa_table_r

        else:
            aa2nt_dict = {v: k for k, v in nt2aa_dict.items()}

        
        self.__check_sequence_multiple_of_3(nt_seq)
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
        
        aa = []
        for i in range(len(seperated_codons)):
            try:
                aa += [nt2aa_dict[seperated_codons[i]],]
            except KeyError:
                msg = 'Unrecognized codon to amino acid: "%s", at index %i provide an entry for this'\
                    ' codon and its amino acid in '\
                    'the nt2aa dictionary.'%(seperated_codons[i],i)
                raise custom_err.UnrecognizedAAError(msg)
                
            try:
                opt_dict[seperated_codons[i]]
            except KeyError:
                msg = 'Unrecognized codon to optimize: "%s", at index %i provide an entry for this'\
                    ' codon and its amino acid in '\
                    'the optimization dictionary.'%(seperated_codons[i],i) 
                raise custom_err.UnrecognizedCodonError(msg)
        
        opt_seq = ''
        for i in range(0, len(aa)):
            
            ind = np.argmax([opt_dict[x] for x in aa2nt_dict[aa[i]]])
            opt_codon = aa2nt_dict[aa[i]][ind]
            opt_seq = opt_seq + opt_codon
        return opt_seq



    def __get_rate_uniform_opt_dicts(self, nt_seq, codon_dict=None, stop_codons=['UAG','UAA','UGA']):
        if codon_dict is None:
            codon_dict = self.codon_dicts.human_codon_frequency_bias_nakamura
            
        
        only_u_dict = {}
        n_keys = 0
        for key in codon_dict.keys():
            if 'T' not in key.upper():
                only_u_dict[key.upper()] = codon_dict[key]
                n_keys += 1
                
        
        n = 1/n_keys
        
        #worst_norm = eqiv to [60/61 + 1/61 * 60], 
        #a vector of N length with 1 codon making up 100% of the sequence and the rest 0,
        #the least possible even percentage
        worst_norm = (n_keys-1)/n_keys + (1/n_keys)*(n_keys-1)

        for x in stop_codons: #remove stop codons
            only_u_dict.pop(x)

        valid_key = list(only_u_dict.keys())
        
        test_seq = self.clean_seq(nt_seq)
        codons = [test_seq[i:i+3] for i in range(0, len(test_seq), 3)] 
        codonkey, codoncount = np.unique(codons,return_counts=True)
    
        keys_to_add = np.array([ x for  x in valid_key if x not in codonkey])
        vals_to_add = np.zeros(len(keys_to_add))
        
        split_str = '.'.join(codons)
        a_table = self.codon_dicts.aa_table_r
        
        opt_dict = {}
        seq_dict = {}
        ccc_dict = {}
        worst_dict = {}
        opt_dict_aa = {}
        listsum = 0
        for key in a_table.keys():
            
            if key != '*':
                codon_list = [x for x in a_table[key] if 'T' not in x]
                
                values = [only_u_dict[x] for x in codon_list]
                codon_cnts = [split_str.upper().count(x) for x in codon_list ]
                listsum += len(codon_cnts)
                
                porportion = np.array(values)/sum(values)
                worst_count = np.zeros(len(codon_cnts)).astype(int)
                worst_count[np.argmin(values)] = int(np.sum(codon_cnts))
                
                guessed_dist = porportion*np.sum(codon_cnts)
                new_codons = np.ceil(np.array(values)/sum(values)*np.sum(codon_cnts))
    
                if sum(new_codons) != sum(codon_cnts):
                    diff = sum(new_codons) - sum(codon_cnts)
                    order = np.argsort(new_codons)[::-1] 
                    
                    while diff > 0:
                        new_codons[order[0]] -= 1
                        diff -= 1
                        
                new_codons = new_codons.astype(int).tolist()
                
                for i in range(len(codon_list)):
                    opt_dict[codon_list[i]] = new_codons[i] 
                    
                for i in range(len(codon_list)):
                    seq_dict[codon_list[i]] = codon_cnts[i] 
                    
                for i in range(len(codon_list)):
                    worst_dict[codon_list[i]] = worst_count[i] 
                            
                ccc_dict[key] = codon_list
                opt_dict_aa[key] = new_codons
        # if len(keys_to_add) !=0:
        #     codonkey = np.append(codonkey,keys_to_add)
        #     codoncount = np.append(codoncount,vals_to_add)
    
        # codon_per = codoncount/(len(test_seq)/3)
    
        # test_metric = 1  - np.sum(np.abs((codon_per-n))) / worst_norm
        
        
        optkeys = list(opt_dict.keys())
        
        seq_arr = np.array([seq_dict[x] for x in optkeys])
        opt_arr = np.array([opt_dict[x] for x in optkeys])
        worst_arr =  np.array([worst_dict[x] for x in optkeys])
        
        worst_err = np.sum(np.abs(opt_arr - worst_arr))
        seq_err = np.sum(np.abs(opt_arr - seq_arr))
        
             
        
        return opt_dict_aa, ccc_dict, seq_dict, opt_dict, worst_dict, opt_dict

    def get_rate_uniform_metric(self, nt_seq, codon_dict=None, stop_codons=['UAG','UAA','UGA']):
        '''
        Rate Uniform Metric:
            
        This sets the metric between 0 and 1, where 1 is the theoretical evenest usage
        of codons porportional to rates of a given sequence.     

        Example:
            
            for each codon:

                Rate uniform optimization, given a set of codons of usage:
                    [10,5,13,1] sum = 23
                with rate values of:
                    [1, 5, 15, 2]  
                set the relative porpotions of usage to:
                    [1/23, 5/23, 13/23, 2/23]
                Theoretical evenest usage of codons is:
                    [1,5,13,2]       
                Theoretical worst usage of codons is:
                    [23,0,0,0]          
                    
                rank proprotional distance of given sequence vs the worst possible usage
                1 - abs(sequence_counts - evenest_counts) / abs(worst_counts - evenest_counts)
                
                1 - ([10,5,13,1] - [1,5,13,2])/([23,0,0,0] - [1,5,13,2]) 
                                                     
        
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence string.
        codon_dict : dict, optional
            Dictionary of codon with its scaling rate. The default is Human nakumura codon frequency (CAI).
        stop_codons : list, optional
            List of stop codons to ignore. The default is ['UAG', 'UAA','UGA'].

        Returns
        -------
        codon_uniform_metric : float
            0-1 uniform metric.

        '''
        opt_dict_aa, ccc_dict, seq_dict, opt_dict, worst_dict, opt_dict = self.__get_rate_uniform_opt_dicts(nt_seq,codon_dict=codon_dict,
                                                          stop_codons=stop_codons)
        
        optkeys = list(opt_dict.keys())
        
        seq_arr = np.array([seq_dict[x] for x in optkeys])
        opt_arr = np.array([opt_dict[x] for x in optkeys])
        worst_arr =  np.array([worst_dict[x] for x in optkeys])
        
        worst_err = np.sum(np.abs(opt_arr - worst_arr))
        seq_err = np.sum(np.abs(opt_arr - seq_arr))
        
        metric = 1 - seq_err/worst_err    
    
        return metric



    def optimize_ntseq_rate_uniform(self, nt_seq, codon_dict=None, random_placement=True, stop_codons =['UAG','UAA','UGA']):
        '''
        Rate uniform optimization, given a set of codons of usage:
            [10,5,13,1]
        with rate values of:
            [1, 5, 15, 2]  
        set the relative porpotions of usage to:
            [1/23, 5/23, 13/23, 2/23]
        for a new codon usage of :
            [1,5,13,2]
        

        Parameters
        ----------
        nt_seq : TYPE
            DESCRIPTION.
        codon_dict : TYPE, optional
            DESCRIPTION. The default is None.
        random_placement : TYPE, optional
            DESCRIPTION. The default is True.
        stop_codons : TYPE, optional
            DESCRIPTION. The default is ['UAG','UAA','UGA'].

        Returns
        -------
        opt_seq : TYPE
            DESCRIPTION.

        '''
        
        
        test_seq = self.clean_seq(nt_seq).upper().replace('t','u')
        
        
        self.get_codon_uniform_metric(nt_seq,codon_dict=codon_dict, stop_codons=stop_codons )
        
        
        opt_dict_aa, ccc_dict, _,_,_,_, = self.__get_rate_uniform_opt_dicts(nt_seq,codon_dict=codon_dict,
                                                                  stop_codons=stop_codons)
        
        cd = [test_seq[i:i+3] for i in range(0, len(test_seq), 3)]
        new_cd = []
        placement_dict = {}
        keys = list(opt_dict_aa.keys())
        for i in range(len(keys)):
            placement_dict[keys[i]] = -1
        
        for i in range(len(cd)):
            
            aa = self.codon_dicts.aa_table[cd[i]]
            
            if aa != '*':
                if random_placement: # randomly place the new codons
                    
                    val_left = 0
                    while val_left == 0: 
                        ind = np.random.randint(len(opt_dict_aa[aa]))
                        val_left = opt_dict_aa[aa][ind]
                        
                    codons_left = opt_dict_aa[aa]
                    codons_left[ind] = codons_left[ind] - 1
                    opt_dict_aa[aa] = codons_left
                    new_cd = new_cd + [ccc_dict[aa][ind], ]

                else:  # cycle through new codons and place them 
                    codons_left = opt_dict_aa[aa]
                    codon_L = len(opt_dict_aa[aa])
                    val_left = 0
                    ind = placement_dict[aa]
                    while val_left == 0:
                        ind+=1
                        ind %= codon_L
                        val_left = opt_dict_aa[aa][ind]
                        
                    
                    codons_left[ind] = codons_left[ind] - 1
                    opt_dict_aa[aa] = codons_left
                    new_cd = new_cd + [ccc_dict[aa][ind], ]
                    placement_dict[aa] = ind 
                    
            else:
                new_cd = new_cd + [cd[i]]
        opt_seq = ''.join(new_cd)
        
        return opt_seq



    def __get_codon_uniform_opt_dicts(self, nt_seq, codon_dict=None, stop_codons=['UAG','UAA','UGA']):
        if codon_dict is None:
            codon_dict = self.codon_dicts.human_codon_frequency_bias_nakamura
        only_u_dict = {}
        n_keys = 0
        for key in codon_dict.keys():
            if 'T' not in key.upper():
                only_u_dict[key.upper()] = codon_dict[key]
                n_keys += 1
                
       
        n = 1/n_keys
        
        #worst_norm = eqiv to [60/61 + 1/61 * 60], 
        #a vector of N length with 1 codon making up 100% of the sequence and the rest 0,
        #the least possible even percentage
        worst_norm = (n_keys-1)/n_keys + (1/n_keys)*(n_keys-1)

        for x in stop_codons: #remove stop codons

            only_u_dict.pop(x)

        valid_key = list(only_u_dict.keys())
        
        test_seq = self.clean_seq(nt_seq).upper()
        codons = [test_seq[i:i+3] for i in range(0, len(test_seq), 3)] 
        codonkey, codoncount = np.unique(codons,return_counts=True)
    
        keys_to_add = np.array([ x for  x in valid_key if x not in codonkey])
        vals_to_add = np.zeros(len(keys_to_add))
    
        if len(keys_to_add) !=0:
            codonkey = np.append(codonkey,keys_to_add)
            codoncount = np.append(codoncount,vals_to_add)
    
        codon_per = codoncount/(len(test_seq)/3)
        

        new_codon_dict = dict(zip([x.upper() for x  in codonkey], [int(x) for x in codoncount]))
        
        split_str = '.'.join(codons)
        
        a_table = self.codon_dicts.aa_table_r
        
        seq_dict = {}
        ccc_dict = {}
        
        opt_dict_aa = {}
        listsum = 0
        for key in a_table.keys():
            
            if key != '*':
                codon_list = [x for x in a_table[key] if 'T' not in x]
                
                values = [only_u_dict[x] for x in codon_list]
                codon_cnts = [split_str.upper().count(x) for x in codon_list ]
                    
                for i in range(len(codon_list)):
                    seq_dict[codon_list[i]] = codon_cnts[i] 
                
                #create the flattest distribution of codons here
                flat = int(np.floor(sum(codon_cnts) / len(codon_cnts)))
                remainder  = sum(codon_cnts) % len(codon_cnts)
                
                new_codons = [flat,]*len(codon_cnts)
                inds = np.argsort(np.abs(values - np.mean(values))).tolist()
                for i in range(remainder):
                    new_codons[inds[i]] +=1
                    
                ccc_dict[key] = codon_list
                opt_dict_aa[key] = new_codons        
        

        return codon_per, new_codon_dict, worst_norm, n, ccc_dict, opt_dict_aa


    def get_codon_uniform_metric(self,nt_seq, codon_dict=None, stop_codons=['UAG', 'UAA','UGA']):
        '''
        Codon Uniform Metric:
            
        .. math::
            UM =  1 - \frac{\sum_{i=0}^{N codons} | \frac{# codon_i}{L} - \frac{1}{N codons} |}{ \frac{N codons-1}{N codons} - \frac{1}{N codons}*(N codons - 1) }

        This sets the metric between 0 and 1, where 1 is the theoretical evenest usage
        of codons of a given sequence.                                                  
        
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence string.
        codon_dict : dict, optional
            Dictionary of codon with its scaling rate. The default is Human nakumura codon frequency (CAI).
        stop_codons : list, optional
            List of stop codons to ignore. The default is ['UAG', 'UAA','UGA'].

        Returns
        -------
        codon_uniform_metric : float
            0-1 uniform metric.

        '''
        codon_per, _,worst_norm, n, _,_, = self.__get_codon_uniform_opt_dicts(nt_seq, codon_dict=codon_dict, stop_codons=stop_codons)
        test_metric = 1  - np.sum(np.abs((codon_per-n))) / worst_norm
        return test_metric

    def optimize_ntseq_codon_uniform(self, nt_seq, codon_dict=None, random_placement=True, stop_codons=['UAG', 'UAA','UGA']):
        '''
        Optimize a sequence so its codons are used as uniformly as possible ie:
        
        an original sequence has [4,10,6,1] usage of synonomous codons for a given AA
        this will set that usage to [6,5,5,5], remainders are distributed preferntially to the 
        closest values to that amino acids average rate in the codon dictionary.

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.
        codon_dict : dict, optional
            codon to rate dictionary to optimize over. The default is None.
        random_placement : bool, optional
            randomly place the new codons or cyclically place them, 
            ie we have a new set of . The default is True.
        stop_codons : list, optional
            list of stop codons to ignore during optimization. The default is ['UAG', 'UAA','UGA'].

        Returns
        -------
        opt_seq : str
            optomized nucleotide sequence.

        '''
        test_seq = self.clean_seq(nt_seq).upper().replace('t','u')
        _, codon_dict, _, _, ccc_dict, opt_dict_aa = self.__get_codon_uniform_opt_dicts(nt_seq, codon_dict=codon_dict, stop_codons=stop_codons)
        

        a_table = self.codon_dicts.aa_table

        cd = [test_seq[i:i+3] for i in range(0, len(test_seq), 3)]
        new_cd = []
        
        placement_dict = {}
        keys = list(opt_dict_aa.keys())
        for i in range(len(keys)):
            placement_dict[keys[i]] = -1
        
        for i in range(len(cd)):
            
            aa = self.codon_dicts.aa_table[cd[i]]
            
            if aa != '*':
                if random_placement: # randomly place the new codons
                    
                    val_left = 0
                    while val_left == 0: 
                        ind = np.random.randint(len(opt_dict_aa[aa]))
                        val_left = opt_dict_aa[aa][ind]
                        
                    codons_left = opt_dict_aa[aa]
                    codons_left[ind] = codons_left[ind] - 1
                    opt_dict_aa[aa] = codons_left
                    new_cd = new_cd + [ccc_dict[aa][ind], ]

                else:  # cycle through new codons and place them 
                    codons_left = opt_dict_aa[aa]
                    codon_L = len(opt_dict_aa[aa])
                    val_left = 0
                    ind = placement_dict[aa]
                    while val_left == 0:
                        ind+=1
                        ind %= codon_L
                        val_left = opt_dict_aa[aa][ind]
                        
                    
                    codons_left[ind] = codons_left[ind] - 1
                    opt_dict_aa[aa] = codons_left
                    new_cd = new_cd + [ccc_dict[aa][ind], ]
                    placement_dict[aa] = ind 
                    
            else:
                new_cd = new_cd + [cd[i]]

        opt_seq = ''.join(new_cd)
        
        return opt_seq

    def deoptimize_ntseq(self, nt_seq, deopt_dict=None, nt2aa_dict=None):

        '''
        Optimizes a nucleotide sequence

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence string
        deopt_dict : dictionary, optional
            a user defined dictionary to deoptimize over. The default is None.

        Returns
        -------
        deopt_seq : str
            Deoptimized NT sequenced based on a given dictionary of rates

        '''

        if deopt_dict is None:
            deopt_dict = self.codon_dicts.human_codon_frequency_bias_nakamura
        if nt2aa_dict is None:
            nt2aa_dict = self.codon_dicts.aa_table
            aa2nt_dict = self.codon_dicts.aa_table_r

        else:
            aa2nt_dict = {v: k for k, v in nt2aa_dict.items()}

        
        self.__check_sequence_multiple_of_3(nt_seq)
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
        
        aa = []
        for i in range(len(seperated_codons)):
            try:
                aa += [nt2aa_dict[seperated_codons[i]],]
            except KeyError:
                msg = 'Unrecognized codon to amino acid: "%s", at index %i provide an entry for this'\
                    ' codon and its amino acid in '\
                    'the nt2aa dictionary.'%(seperated_codons[i],i)
                raise custom_err.UnrecognizedAAError(msg)
                
            try:
                deopt_dict[seperated_codons[i]]
            except KeyError:
                msg = 'Unrecognized codon to optimize: "%s", at index %i provide an entry for this'\
                    ' codon and its amino acid in '\
                    'the deoptimization dictionary.'%(seperated_codons[i],i) 
                raise custom_err.UnrecognizedCodonError(msg)

        opt_seq = ''
        for i in range(0, len(aa)):
            ind = np.argmin([deopt_dict[x] for x in aa2nt_dict[aa[i]]])
            opt_codon = aa2nt_dict[aa[i]][ind]
            opt_seq = opt_seq + opt_codon
        return opt_seq


    def nt2aa(self, nt_seq):
        '''
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.

        Returns
        -------
        aa : str
            amino acid sequence.

        '''
        self.__check_sequence_multiple_of_3(nt_seq)
        aa = ''
        nt_seq = nt_seq.upper()
        
        for i in range(0, len(nt_seq), 3):
            try:
                aa += self.codon_dicts.aa_table[nt_seq[i:i+3]]
            except:
                aa += 'X'
        return aa


    def get_orfs(self, nt_seq, min_codons=80):
        '''
        convert a nucleotide sequence string into a dictionary of open reading
        frames (orfs). Dictionaries are formatted as follows:
    
            {'0': [(start,stop),...  ],
             '+1': [(start,stop),... ],
             '+2': [(start,stop),... ],
             '-1': [(start,stop),... ]}     
            
        Orf dictionary returns tuple pairs of start and stop indices along a 
        given frame, 0, +1, or +2/-1 (identical).

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence as a string.
        min_codons : int, optional
            minimum number of contiguous nucleotides to consider as
            a valid protein/open reading frame. The default is 80.

        Returns
        -------
        orfs : dictionary
            dictionary of tuple pairs of start and stop indices for sorted by
            keys of '0', '+1', '+2', '-1' frames.

        '''

        if nt_seq == '':
            nt_seq = self.sequence.upper()
        nt_seq = nt_seq.upper()
        allstarts = np.array([m.start() for m in re.finditer(
            '(?=A[TU]G((?:.{3})+?)[TU](?:AG|AA|GA))', nt_seq)])


        #allsegments = re.findall('(?=A[TU]G((?:.{3})+?)[TU](?:AG|AA|GA))',self.sequence_str)
        allstops = np.array(
            [m.start() for m in re.finditer('(?=[TU](?:AG|AA|GA))', nt_seq)])
        start_frames = allstarts%3
        stop_frames = allstops%3
        min_len = min_codons*3
        orf1_starts = allstarts[np.where(start_frames == 0)]
        orf2_starts = allstarts[np.where(start_frames == 1)]
        orf3_starts = allstarts[np.where(start_frames == 2)]

        orf1_stops = allstops[np.where(stop_frames == 0)]
        orf2_stops = allstops[np.where(stop_frames == 1)]
        orf3_stops = allstops[np.where(stop_frames == 2)]

        #self.starts = [orf1_starts, orf2_starts, orf3_starts]
        #self.stops = [orf1_stops, orf2_stops, orf3_stops]
        orfs = {'0':[], '+1':[], '+2':[], '-1':[]}
        orfs = {'0':[], '+1':[], '+2':[], '-1':[]}
        laststop = 0
        for start in orf1_starts:
            nextstop = orf1_stops[np.where(orf1_stops > start)[0][0]]+3
            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    orfs['0'].append((start, nextstop))

                    laststop = nextstop

        laststop = 0
        for start in orf2_starts:
            nextstop = orf2_stops[np.where(orf2_stops > start)[0][0]]+3
            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    orfs['+1'].append((start, nextstop))
                    laststop = nextstop

        laststop = 0
        for start in orf3_starts:
            nextstop = orf3_stops[np.where(orf3_stops > start)[0][0]]+3

            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    orfs['+2'].append((start, nextstop))
                    orfs['-1'].append((start, nextstop))
                    laststop = nextstop



        return orfs
    
    def get_codon_list(self, nt_seq):
        '''
        

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.

        Returns
        -------
        List of codons.

        '''
        #split codons by 3
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)]
        return seperated_codons
    
    def get_codon_count_dict(self, nt_seq):
        '''
        
        
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.

        Returns
        -------
        dict of codon counts.

        '''
        #split codons by 3
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)]
        cdict = {}
        for i in range(len(seperated_codons)):
            if seperated_codons[i] in cdict.keys():
                cdict[seperated_codons[i]] +=1
            else:
                cdict[seperated_codons[i]] = 1
            
        return cdict

    
    def __check_trna_gene_dict(self,trna_gene_dict):
      '''
      check if a dictionary's anticodons are written in 5'-->3' 
      or 3'-->5' order
    
      Returns: index of the 5' anticodon, bool flipped
      '''
      try: 
        tmp_keys = list(trna_gene_dict.keys())  #get the trna dictionary keys
        tmp_keys_lower = [x.lower() for x in tmp_keys] #get the lower strings of them
        ala_key = tmp_keys[tmp_keys_lower.index('ala')] #get the ALA key
        alanines = trna_gene_dict[ala_key]
        if len(set([x[0] for x in alanines])) == 1: # if they are oriented 3' to 5' ie all starting with C's 
            flipped = False
            anti_codon_3rd_nt = 2
        elif len(set([x[0] for x in alanines])) == 3:
            flipped = True
            anti_codon_3rd_nt = 0
        else:
            ValueError('could not recognize what this dictionary is')
      except:
        anti_codon_3rd_nt = 0
        flipped= True
      return anti_codon_3rd_nt, flipped

    def __default_tai_dicts(self):
        '''
        get the default tai_dicts
        '''
        wobble_rules = {'A':['U','C','G','A'],
            'C':['G'],
            'G':['C','U'],
            'U':['A','G','C','U'],
            'I':['A','C','U']}    
        
        bp_rules = {'A':'U', 'G':'C','C':'G','U':'A'}
        
        sij_dict = {'A:U':0, 'A:A':.5, 'A:G':.5, 'A:C':.5,
                      'C:G':0,
                      'G:C':0, 'G:U':.41,
                      'U:A':0, 'U:G':.68, 'U:C':.5, 'U:U':.5,
                      'I:U':0, 'I:C':.28, 'I:A':.99,
                      'L:A':.89, }        
        
        return bp_rules, wobble_rules, sij_dict
      
        
    def get_tai_weights(self, nt_seq,
                        tRNA_gene_copy_dict = None,
                        alt_stopcodons = ['UAG','UGA','UAA'],
                        bp_rules = None,
                        wobble_rules = None, 
                        anticodon_order='5->3',
                        sij_dict = None,
                        ignore_start_codon=True):
        '''
        get the tRNA Adaptation Index (tAI) weights for an arbitrary:
            
        * base pairing dictionary
        * wobble pairng dictionary
        * sij interaction dictionary
        * tRNA gene copy number dictionary

        Parameters
        ----------
        nt_seq : str
            Nucleotide sequence.
        tRNA_gene_copy_dict : dict, optional
            tRNA gene copy count by anticodon key. The default is None.
            If left blank automatically uses hg19 tRNA counts from gTRNAdb,
            stored in aux_cdict.hg19_tRNA_GCN_by_anticodon
        alt_stopcodons : list, optional
            list of stop codons to remove / not consider for TAI. 
            The default is ['UAG','UGA','UAA'].
        bp_rules : dict, optional
            base pair rules. The default is {A:U, G:C, C:G, U:A}.
        wobble_rules : dict, optional
            allowed wobble pairing rules. The default is 
                {'A':['U','C','G','A'],
                'C':['G'],
                'G':['C','U'],
                'U':['A','G','C','U'],
                'I':['A','C','U']}.
        anticodon_order : str, optional
            '5->3' or '3->5' ordering for the keys of
            the tRNA gene copy number dict. The default is '5->3'.
        sij_dict : dict, optional
            wobble base pair efficiencies 0 = 100% efficent (1-sij). The default is
             {'A:U':0, 'A:A':.5, 'A:G':.5, 'A:C':.5,
                'C:G':0,
                'G:C':0, 'G:U':.41,
                'U:A':0, 'U:G':.68, 'U:C':.5, 'U:U':.5,
                'I:U':0, 'I:C':.28, 'I:A':.99,
                'L:A':.89, }
        ignore_start_codon : bool, optional
            ignore the first ATG for calculating tAI. The default is True.

        Returns
        -------
        dict
            codon to tAI weight dictionary.

        '''
        
        bases = ['A','G','C','U']
        all_3mers = [''.join(i) for i in it.product(bases, repeat = 3)]
        
        if tRNA_gene_copy_dict == None:
            tRNA_gene_copy_dict = self.aux_dicts.hg19_tRNA_GCN_by_anticodon
        
        #get the bp/ wobble rules if specified
        if bp_rules == None:
            bp_rules = {'A':'U', 'G':'C','C':'G','U':'A'}
        if wobble_rules == None:
            wobble_rules = {'A':['U','C','G','A'],
            'C':['G'],
            'G':['C','U'],
            'U':['A','G','C','U'],
            'I':['A','C','U']}    
        
        
        #remove stop codons
        for stop in alt_stopcodons:
            all_3mers.remove(stop)
        
        #available trna anticodons from the dictionary given
        available_genes = sorted(list(tRNA_gene_copy_dict.keys()))
        
        #get matches based on the bp/wobble rules and a given codon
        def generate_matches(codon,bp_rules, bp_wobble_rules):
            frnt = bp_rules[codon.upper()[0]] + bp_rules[codon.upper()[1]] 
            return [frnt + x for x in bp_wobble_rules[codon.upper()[2]]] 
        
        #anti_codon_3rd_nt, flipped = check_trna_gene_dict(trna_gene_dict)
        if anticodon_order=='5->3':
            anti_codon_3rd_nt = 0
            flipped=True
          
        codon_to_species = {}
        for i in range(len(all_3mers)): #for every codon
            codon =  all_3mers[i]  
            match_ids = []
            #get all anticodons that match
            for anticodon in available_genes: 
                if flipped: #flip the anticodon if needed
                    matches = generate_matches(anticodon[::-1], bp_rules, wobble_rules)
                    if codon in matches:
                        match_ids =match_ids + [anticodon]
                else:
                    matches =  generate_matches(anticodon, bp_rules, wobble_rules)
                    if codon in matches:
                        match_ids =match_ids + [anticodon]
            
            codon_to_species[codon] = match_ids 

 
        Wi = {}
          
        #interaction efficiency weight dictionary 5' anticodon to 3' codon
        # 0 = 100% efficent 
        if sij_dict==None:
          sij_dict = {'A:U':0, 'A:A':.5, 'A:G':.5, 'A:C':.5,
                        'C:G':0,
                        'G:C':0, 'G:U':.41,
                        'U:A':0, 'U:G':.68, 'U:C':.5, 'U:U':.5,
                        'I:U':0, 'I:C':.28, 'I:A':.99,
                        'L:A':.89, }
          
        codons_given = list(tRNA_gene_copy_dict.keys()) #available anticodons
        codons_to_remove = alt_stopcodons
        codons_given = [x for x in codons_given if x not in codons_to_remove]
          
        #sanity check, check if anticodon is 5'->3' or 3'->5'
        anti_codon_3rd_nt, flipped = self.__check_trna_gene_dict(tRNA_gene_copy_dict) 
        
        for codon in codon_to_species:
            weight = 0
              # Wi = sumj[ (1-sij)*tGCNj]
            for anticodon in codon_to_species[codon]:
                sij = sij_dict[ (str(anticodon[anti_codon_3rd_nt]) + ':' + str(codon[2]))]
                weight += (1-sij)*tRNA_gene_copy_dict[anticodon]
              
            Wi[codon] = weight #store the weights
          
        #get nonzero weights
        Wi_nonzero = np.array(list(Wi.values()))[np.array(list(Wi.values()))!=0] 
        wmean = self.geomean(Wi_nonzero) #geometric mean of nonzero weights
        wmax = np.max(Wi_nonzero) #max weight
        wi = {}
        #if the key = 0 set it to geomean of nonzero weights
        for key in Wi.keys(): 
            if Wi[key] !=0:
                wi[key] = Wi[key]/wmax
            else:
                wi[key] = wmean
            #return the weight dictionary
        return wi

        
    def get_metric_from_weights(self, nt_seq, weight_dict,
                alt_stopcodons = ['UAG','UGA','UAA'],
                ignore_start_codon=True):
        '''
        given an arbitrary dictionary of weights per codon, calculate 
        geometric mean of those weights from a sequence.

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.
        weight_dict : dict
            dict of weights with keys being codons.
        alt_stopcodons : list, optional
            DESCRIPTION. The default is ['UAG','UGA','UAA'].
        ignore_start_codon : bool, optional
            ignore the start codon when calculating the metric. The default is True.

        Returns
        -------
        float
            geometric mean of the weights per codon in the sequence.
        metric_codons : list
            list of weights per codon.

        '''
        
     
        mRNA_sequence_str = nt_seq.upper().replace('T','U') 
        #convert the sequence to uppercase and use U
        #codons_to_remove = ['UGA', 'UAG','UAA'] #remove stop codons
        #get codons from the sequeence and remove stops
        codons = [mRNA_sequence_str[i:i+3] for i in range(0, len(mRNA_sequence_str), 3)]
        codons = [x for x in codons if x not in alt_stopcodons]
        if ignore_start_codon: #ignore the start codon if needed
            codons = codons[1:]
        metric_codons = [] #get the list of tai per codon
        for i in range(len(codons)):
            metric_codons = metric_codons + [weight_dict[codons[i]],]
        #return geometric mean of that list of tai per codons
        
        return self.geomean(metric_codons), metric_codons

    def get_tai(self, nt_seq, tRNA_gene_copy_dict=None,
                sij_dict=None,
                alt_stopcodons = ['UAG','UGA','UAA'],
                bp_rules = None,
                wobble_rules = None, 
                anticodon_order='5->3',
                ignore_start_codon=True):
        
        '''
        get the tRNA Adaptation Index (tAI) weights for an arbitrary:
            
        * base pairing dictionary
        * wobble pairng dictionary
        * sij interaction dictionary
        * tRNA gene copy number dictionary

        Parameters
        ----------
        nt_seq : str
            Nucleotide sequence.
        tRNA_gene_copy_dict : dict, optional
            tRNA gene copy count by anticodon key. The default is None.
            If left blank automatically uses hg19 tRNA counts from gTRNAdb,
            stored in aux_cdict.hg19_tRNA_GCN_by_anticodon
        alt_stopcodons : list, optional
            list of stop codons to remove / not consider for TAI. 
            The default is ['UAG','UGA','UAA'].
        bp_rules : dict, optional
            base pair rules. The default is {A:U, G:C, C:G, U:A}.
        wobble_rules : dict, optional
            allowed wobble pairing rules. The default is 
                {'A':['U','C','G','A'],
                'C':['G'],
                'G':['C','U'],
                'U':['A','G','C','U'],
                'I':['A','C','U']}.
        anticodon_order : str, optional
            '5->3' or '3->5' ordering for the keys of
            the tRNA gene copy number dict. The default is '5->3'.
        sij_dict : dict, optional
            wobble base pair efficiencies 0 = 100% efficent (1-sij). The default is
             {'A:U':0, 'A:A':.5, 'A:G':.5, 'A:C':.5,
                'C:G':0,
                'G:C':0, 'G:U':.41,
                'U:A':0, 'U:G':.68, 'U:C':.5, 'U:U':.5,
                'I:U':0, 'I:C':.28, 'I:A':.99,
                'L:A':.89, }
        ignore_start_codon : bool, optional
            ignore the first ATG for calculating tAI. The default is True.

        Returns
        -------
        float
            geometric mean of the tAI weights per codon in the sequence.
        tai_codons : list
            list of tAI weights per codon.

        '''
        
        
        tai_weight_dict = self.get_tai_weights(nt_seq,
                               tRNA_gene_copy_dict=tRNA_gene_copy_dict,   
                               alt_stopcodons=alt_stopcodons,
                               bp_rules=bp_rules,
                               wobble_rules=wobble_rules,
                               anticodon_order=anticodon_order,
                               ignore_start_codon=ignore_start_codon,
                               sij_dict=sij_dict)
        
        
        mRNA_sequence_str = nt_seq.upper().replace('T','U') 
        #convert the sequence to uppercase and use U
        if alt_stopcodons == None:
            alt_stopcodons = ['UAG','UGA','UAA']
        #codons_to_remove = ['UGA', 'UAG','UAA'] #remove stop codons
        #get codons from the sequeence and remove stops
        codons = [mRNA_sequence_str[i:i+3] for i in range(0, len(mRNA_sequence_str), 3)]
        codons = [x for x in codons if x not in alt_stopcodons]
        if ignore_start_codon: #ignore the start codon if needed
            codons = codons[1:]
        tai_codons = [] #get the list of tai per codon
        for i in range(len(codons)):
            tai_codons = tai_codons + [tai_weight_dict[codons[i]],]
        #return geometric mean of that list of tai per codons
        
        return self.geomean(tai_codons), tai_codons

    

    def get_cai(self, nt_seq, codon_dict=None):
        '''
        get codon adaptation index (CAI)

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.
        codon_dict : dict, optional
            codon frequency dictionary. The default is None.

        Returns
        -------
        float
            codon adaptation index of the sequence.

        '''
        if codon_dict == None:
            codon_dict = self.codon_dicts.human_codon_frequency_bias_nakamura
            
        cai_codons = []
        for i in range(0, len(nt_seq), 3):
            synonmous_codons = self.codon_dicts.aa_table_r[
                self.codon_dicts.aa_table[nt_seq[i:i+3]]]

            max_freq = max([codon_dict[x] for x in synonmous_codons])

            cai_codons.append(codon_dict[nt_seq[i:i+3]] /max_freq)

        return self.geomean(cai_codons)
        
    def codon_usage(self, nt_seq, codon_dict=None):
        '''
        return codon adaptation index (CAI) as well as codon sensitivity and
        cai per codon.
            

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.
        codon_dict : dict, optional
            codon frequency dict. The default is None.

        Returns
        -------
        codon_sensitivity : list
            how codon sensitive a sequence is per codon.
        cai : float
            codon adaptation index.
        cai_codons : list
            codon adaptation index per codon in the sequence.

        '''

        if codon_dict == None:
            codon_dict = self.codon_dicts.human_codon_frequency_bias_nakamura

        codon_usage = np.zeros((1, 21))
        gene_len = len(nt_seq)/3
        aa_seq = self.nt2aa(nt_seq)
        aa_seq = aa_seq.replace('X','') #remove unknown sequences
        for i in range(21):
            if self.codon_dicts.aa_keys[i] != '*':
                codon_usage[0, i] = len(
                    re.findall(self.codon_dicts.aa_keys[i], aa_seq))
            else:
                codon_usage[0, i] = len(re.findall('\*', aa_seq))
        codon_norm = codon_usage/gene_len
        codon_sensitivity = np.round(
            codon_norm*self.codon_dicts.sensitivity_fast_slow, 2)

        cai_codons = []
        for i in range(0, len(nt_seq), 3):
            synonmous_codons = self.codon_dicts.aa_table_r[
                self.codon_dicts.aa_table[nt_seq[i:i+3]]]

            max_freq = max([codon_dict[x] for x in synonmous_codons])

            cai_codons.append(codon_dict[nt_seq[i:i+3]] /max_freq)

        cai = self.geomean(cai_codons)

        return codon_sensitivity, cai, cai_codons


    def clean_seq(self, nt_sequence):
        '''
        Return an mrna sequence of lowercase a,u,c,g from IPUAC substitutions

        .. warning:: 
            this code will replace substitutive nucleotides with
            preferential order a, g , u , c. for example: N (any base) is set
            to A, W (T, U, or A) is set to A, S (C or G) is set to G

        .. warning:: 
            this code will replace any modified nucleotides such as Ψ'
            or m6a to a,u,g,c. The code is not configured for modified nucleotides
            in this version.


        Parameters
        ----------
        seq : str
            nucleotide sequence string to convert to only a, u, g, or c..

        Returns
        -------
        seq : str
            cleaned nucleotide sequence str (only lowercase a,u,c,g).

        '''
        nt_sequence = nt_sequence.lower()

        for key in self.sub_dict.keys():
            nt_sequence = nt_sequence.replace(key, self.sub_dict[key])

        nt_sequence = nt_sequence.replace('t', 'u')
        return nt_sequence
    

    def get_kmer_freq(self, nt_seq, kmer_length, substrings=False):
        '''
        return the K-mer frequency of a nuclotide sequence
    
        seq = 'AACGTACGTAGCTCATG...'
        
        kmer_dict with length 3 = 
        
        {'AAA':1, 'AAC': 2, 'AAG':0 ...}
        

        Parameters
        ----------
        nt_seq : str
            nuclotide sequence to get kmers from.
        kmer_length : int
            the length of the K-mer to generate, for example 3 would generate
            a 3-mer vector of counts of AAA, AAC, AAG... 
        substrings : bool, optional
            generate kmers for all lower length kmers as well, ie
            for k = 3, 2, and 1 instead of just 3. The default is False.

        Returns
        -------
        kmer_freq_vec : 1darray
            kmer frequency per k-mer key.
        kmer_ind : list of str
            list of keys for the k-mer frequency vector.

        '''
          
        
        nt_seq = nt_seq.upper()
        unique_char = list(set(nt_seq))
        total_combos = 0
        for i in range(2,kmer_length+1):
            total_combos += len(unique_char)**(i-1)

        if total_combos > 1e5:
            msg = 'Memory Warning: you have requested %i total combinations of subsequences, this'\
            ' will create large array and may take lots of memory.'%total_combos
            warnings.warn(msg)

        combos =[x for x in it.product(unique_char, repeat=kmer_length)] 
        if substrings:
            for n in range(kmer_length-1, 0, -1):
                combos += [x for x in it.product(unique_char, repeat=n)] 
        kmer_ind = [''.join(y) for y in combos]

        kmer_freq_vec = np.zeros(len(combos)).astype(int)
        for i in range(len(nt_seq)-len(nt_seq)%kmer_length - kmer_length + 1):
            kmer_freq_vec[kmer_ind.index(nt_seq[i:i+kmer_length])] += 1
        if substrings:
            for n in range(kmer_length-1, 0, -1):
                for i in range(len(nt_seq)- len(nt_seq)%n - n + 1):
                    kmer_freq_vec[kmer_ind.index(nt_seq[i:i+n])] += 1            

        return kmer_freq_vec, kmer_ind
    
    
    def clean_nt_seq(self, nt_seq, upper=True, sub_dict=None, t_or_u='u',
                     random_sub=False):
        '''
        Return an mrna sequence of lowercase a,u,c,g from IPUAC substitutions

        .. warning:: this code will replace substitutive nucleotides with
        preferential order a, g , u , c. for example: N (any base) is set
        to A, W (T,U, or A) is set to A, S (C or G) is set to G



        Parameters
        ----------
        seq : str
            sequence string.

        Returns
        -------
        seq : str
            cleaned sequence str (only lowercase a,u,c,g).

        '''
        if sub_dict == None:
            sub_dict = self.codon_dicts.ipuac_nt_t

        seq = nt_seq.lower()

        if random_sub:
            new_str = []
            for i in range(len(seq)):
                if seq[i] not in ['a','t','g','c','u']:
                    idx = np.random.randint( len(sub_dict[seq[i]])   )
                    new_str += [sub_dict[seq[i]][idx] , ]
                else:
                    new_str += [seq[i], ]
            seq = ''.join(new_str)
        else:
            for key in sub_dict.keys():
                seq = seq.replace(key, sub_dict[key][0])
        
        if t_or_u in ['u','U']:
            seq = seq.replace('t', 'u')
        elif t_or_u in ['t','T']:
            seq = seq.replace('u', 't')
        else:
            msg = "Cannot recognize the flag for t_or_u, please use 'u' or 't'"\
                " to indicate which letter to use for the sequence."
            raise custom_err.UnrecognizedFlagError(msg)
            
        if upper:
            seq = seq.upper()
        else:
            seq = seq.lower()
        return seq

    
    def get_kmer_freq_dict(self, nt_seq, kmer_length, substrings=False ):
        '''
        Get K-mer frequency as a dictionary, for example
        
        seq = 'AACGTACGTAGCTCATG...'
        
        kmer_dict with length 3 = 
        
        {'AAA':1, 'AAC': 2, 'AAG':0 ...}
        
        

        Parameters
        ----------
        nt_seq : str
            DESCRIPTION.
        kmer_length : int
            size kmers to generate.
        substrings : bool, optional
            generate kmers for all lower length kmers as well, ie
            for k = 3, 2, and 1 instead of just 3. The default is False.
        
        Returns
        -------
        dict
            dictionary of kmer frequencies.

        '''
        kmer_freq, kmer_key = self.get_kmer_freq(nt_seq, kmer_length, substrings=substrings)
        return  dict(zip(kmer_key, kmer_freq))

    def get_gc_content(self, nt_seq):
        '''
        Get the GC content of a nucleotide sequence

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.

        Returns
        -------
        GC content float (percentage).

        '''
        nt_seq = nt_seq.upper()
        return float((nt_seq.count('G') + nt_seq.count('C')))/len(nt_seq)
        

    def get_tag_loc(self, aa_seq, tag, epitope_loc='front'):
        '''
        

        Parameters
        ----------
        aa_seq : str
            amino acid sequence.
        tag : str
            amino acid epitope sequence.
        epitope_loc : str, optional
            where to consider the binary epitope, at the 'front', 'back' or 'middle'.
            Given the example AAA GGG AAA searching for GGG, the front option 
            will return index 2 as the epitope location, the middle option will return 3,
            and the back option will return 4. The default is 'front'.

        Returns
        -------
        list
            list of binary tag locations per codon.

        '''
        cd = self.codon_dicts
        epitope_loc = epitope_loc.lower()
        if epitope_loc == 'front':
            offset = 0
        if epitope_loc == 'middle':
            offset = int(len(tag)/2)
        if epitope_loc == 'center':
            offset = int(len(tag)/2)
        if epitope_loc == 'back':
            offset = len(tag)

        return [m.start()+1+offset for m in re.finditer(tag, aa_seq)]


    @staticmethod
    def geomean(iterable):
        '''equal weight geometric mean used for codon sensitivity calculations
        '''
        a = np.array(iterable)
        return np.exp(np.sum(np.log(a))/len(a)   )
    
    @staticmethod
    def __check_sequence_multiple_of_3(sequence):
        if len(sequence) %3 !=0:
            msg = 'Given nucleotide/codon Sequence is not a multiple of 3, double check'\
                ' the sequence provided.'
            raise custom_err.InvalidSequenceLengthError(msg)
        else:
            return True
