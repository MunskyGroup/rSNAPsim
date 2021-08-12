# -*- coding: utf-8 -*-
"""
Created on Mon May  3 12:38:15 2021

@author: willi
"""

import numpy as np
import itertools

class tRNAmap():
    def __init__(self):
        self.trna_map = None
        
    
    
    
        # https://www.nature.com/articles/nsmb866
        #5' anticodon base to 3' codon base
        self.bp_rules = {'A':'U', 'G':'C','C':'G','U':'A'}
        
        self.bp_wobble_rules = {'A':['U','C','G','A'],
                                'C':['G'],
                                'G':['C','U'],
                                'U':['A','G','C','U'],
                                'I':['A','C','U']}
        
        self.hg19_trna_gene_number_anticodon = {'AGC': 33, 'CGC':5,
                                             'UGC':11, 'ACG': 8, 'GGC':2,
                                             'CCG':4, 'CCU': 8,
                                             'UCG':6, 'UCU': 6,
                                             'AUU':2, 'GUU':36,
                                             'GUC':19, 'GCA':36,
                                             'CUG':21,'UUG':10,
                                             'CUC':13, 'UUC':14,
                                             'CCC':13, 'GCC':15,
                                             'UCC': 11, 'GUG':10,
                                             'AAU':18, 'GAU':6, 
                                             'UAU':5, 'AAG':13,
                                             'CAA':8, 'CAG':10,
                                             'UAA':6,'UAG':4,
                                             'CUU':26,'UUU':20,
                                             'CAU':11,'GAA':15,
                                             'AGG':11, 'CGG':4,
                                             'UGG':9,'UCA':2,
                                             'AGA':11,'CGA':4,
                                             'GCU':9,'UGA':5,
                                             'CUA':0,'UUA':2,
                                             'AGU':10,'CGU':6,
                                             'UGU':6,'CCA':8,
                                             'AUA':1,'GUA':16,
                                             'AAC':12,'CAC':18, #duplicate for doubled up codons
                                             'UAC':7, }
               
    
    def get_codon_species_dict(self,trna_gene_dict, alt_stopcodons = None):
        bases = ['A','G','C','U']
        all_3mers = [''.join(i) for i in itertools.product(bases, repeat = 3)]
        if alt_stopcodons == None:
            alt_stopcodons = ['UAG','UGA','UAA']
        
        for stop in alt_stopcodons:
            all_3mers.remove(stop)
            
    
        available_genes = sorted(list(trna_gene_dict.keys()))
        
        #get the available tRNA
        #available_genes = available_genes + [ 'IAC', 'IGU', 'ICG', 'IAG','IAU','IGC' ]
        ## check the anticodon orientation                                              
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
            flipped= True
        
        codon_to_species = {}
        for i in range(len(all_3mers)):
            codon =  all_3mers[i]
            match_ids = []
            for anticodon in available_genes:
                if flipped:
                    matches = self.generate_matches(anticodon[::-1])
                    if codon in matches:
                        match_ids =match_ids + [available_genes.index(anticodon)]
                else:
                    matches = self.generate_matches(anticodon)
                    if codon in matches:
                        match_ids =match_ids + [available_genes.index(anticodon)]
                codon_to_species[codon] = match_ids
                
            codon_to_species[codon] = match_ids
        return codon_to_species,available_genes
    
    
    def generate_map_matrix(self,trna_gene_dict, alt_stopcodons = None):
        codon_to_species, available_anticodons  =  self.get_codon_species_dict(trna_gene_dict, alt_stopcodons=alt_stopcodons )
        n_codons = len(codon_to_species)
        n_species = len(available_anticodons)
        codon_to_trna_species_map = np.zeros([n_codons,n_species])
        
        bases = ['A','G','C','U']
        all_3mers = [''.join(i) for i in itertools.product(bases, repeat = 3)]
        
        if alt_stopcodons == None:
            alt_stopcodons = ['UAG','UGA','UAA']       
        for stop in alt_stopcodons:
            all_3mers.remove(stop)
            
        for i,codon in enumerate(all_3mers):
            species = sorted(codon_to_species[codon])
            for idx in species:
                codon_to_trna_species_map[i,idx] = 1
        return codon_to_trna_species_map, all_3mers, available_anticodons
            
        
        
        
        
    def generate_matches(self,codon,rules=None):
        if rules == None:
            rules = self.bp_wobble_rules
        frnt = self.bp_rules[codon.upper()[0]] + self.bp_rules[codon.upper()[1]] 
        return [frnt + x for x in self.bp_wobble_rules[codon.upper()[2]]] 
        