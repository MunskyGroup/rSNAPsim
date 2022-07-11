# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 18:06:56 2020

@author: William Raymond
"""
import numpy as np


'''
This is a dictionary object that handles everything dealing with tag sequences and AA tables
'''


class CodonDictionaries():
    '''
    Attributes
    ----------

    tag_dict: dict
        dictionary of tag name and their epitopes
    tag_colors: dict
        dictionary of tag name and their colors
    tag_full: dict
        dictionary of tag name and their full nt sequence
    aa_keys: list
        keys for amino acids
    aa_table: dict
        codon to aa table
    aa_table_r: dict
        aa to codons table
    human_codon_frequency_bias_nakamura: dict
        codon gene copys for both U|T
    human_codon_frequency_bias_nakamura_single: dict
        a codon gene copy number with only T, for average /
        calculation basis that would be thrown off with including U and T
    trna_ids: list
        tRNA species to ID

    '''


    def __init__(self):
        
        # a, g , u , c. #order of replacement if non random
        self.ipuac_nt_t = {'a':['a'], 't':['t'], 'u':['u'], 'c':['c'],
                         'w':['a','t'], 's':['g','c'], 'm':['a','c'],
                         'k':['g','t'], 'r':['a','g'], 'y':['t','c'],
                         'b':['g','t','c'], 'd':['a','g','t'],
                         'h':['a','t','c'], 'v':['a','g','c'],
                         'n':['a','t','g','c'], 'Ψ':['u'],
            
            }
        
        self.ipuac_nt_u = {'a':['a'], 't':['t'], 'u':['u'], 'c':['c'],
                         'w':['a','u'], 's':['g','c'], 'm':['a','c'],
                         'k':['g','u'], 'r':['a','g'], 'y':['u','c'],
                         'b':['g','u','c'], 'd':['a','g','u'],
                         'h':['a','u','c'], 'v':['a','g','c'],
                         'n':['a','u','g','c'], 'Ψ':['u'],
            
            }
                
        
        #some common epitopes used for protein fluorescent tagging
        self.tag_dict = {'T_SunTag':'EELLSKNYHLENEVARLKK',
                         'T_Flag':'DYKDDDDK',
                         'T_Hemagglutinin':'YPYDVPDYA',
                         'T_Myc':'EQKLISEEDL',
                         'T_Strep':'WSHPQFEK',
                         'T_Hist':'HHHHHH',
                         'T_V5':'GKPIPNPLLGLDST',
                         'T_TC':'CCPGCC'}
        
        #some default colors for those tags
        self.tag_colors = {'T_SunTag':'green',
                           'T_Flag':'blue',
                           'T_Hemagglutinin':'blue',
                           'T_Myc':'green',
                           'T_Strep':'green',
                           'T_Hist':'green',
                           'T_V5':'green',
                           'T_TC':'green'}
        
        #Full NT sequences of tags if known
        self.tag_full = {'T_Flag':('ATGGACTACAAGGACGACGACGACAAAGGTGAC'
                                   'TACAAAGATGATGACGATAAAGGCGACTATA'
                                   'AGGACGATGACGACAAGGGCGGAAACTCACTGA'
                                   'TCAAGGAAAACATGCGGATGAAGGTGGTGAT'
                                   'GGAGGGCTCCGTGAATGGTCACCAGTTCAAGTG'
                                   'CACCGGAGAGGGAGAGGGAAACCCGTACATG'
                                   'GGAACTCAGACCATGCGCATTAAGGTCATCGAA'
                                   'GGAGGTCCGCTGCCGTTCGCTTTCGATATCC'
                                   'TGGCCACTTCGTTCGGAGGAGGGTCGCGCACGTTC'
                                   'ATCAAGTACCCGAAGGGAATCCCGGACTT'
                                   'CTTTAAGCAGTCATTCCCGGAAGGATTCACTTGGG'
                                   'AACGGGTGACCCGGTATGAAGATGGAGGT'
                                   'GTGGTGACTGTCATGCAAGATACTTCGCTGGAGGATGGG'
                                   'TGCCTCGTGTACCACGTCCAAGTCC'
                                   'GCGGAGTGAATTTCCCGTCCAACGGACCAGTGATGCAG'
                                   'AAAAAGACGAAGGGTTGGGAACCTAA'
                                   'TACTGAAATGATGTACCCCGCAGACGGAGGGCTGAGGG'
                                   'GCTACACCCACATGGCGCTGAAGGTC'
                                   'GACGGAGGAGATTACAAGGATGACGACGATAAGCAACAA'
                                   'GATTACAAAGACGATGATGACAAGG'
                                   'GCCAGCAGGGCGACTACAAGGACGACGACGACAAGCAG'
                                   'CAGGACTACAAAGATGACGATGATAA'
                                   'AGGAGGAGGACATCTGTCCTGTTCGTTCGTGACCACCT'
                                   'ACAGATCAAAGAAAACCGTGGGAAAC'
                                   'ATCAAGATGCCGGGCATTCATGCCGTCGACCACCGCCT'
                                   'GGAGCGGCTCGAAGAATCAGACAATG'
                                   'AGATGTTCGTCGTGCAAAGAGAACATGCCGTGGCCAAGTT'
                                   'CGCGGGACTGGGAGGCGGTGGAGG'
                                   'CGATTACAAAGACGATGATGACAAGGGTGACTATAAAGA'
                                   'CGACGATGACAAAGGGGATTACAAG'
                                   'GATGATGATGATAAGGGAGGCGGTGGATCAGGTGGAG'
                                   'GAGGTTCACTGCAG')}

        self.aa_keys = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                        'L', 'K', 'M', 'F',
                        'P', 'S', 'T', 'W', 'Y', 'V', '*', 'X'] #* is stop, X is unspecified

        self.codon_types = dict(zip(
            self.aa_keys, np.ones((1, 22)).flatten().astype(int).tolist()))

        self.aa_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',

            'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
            'ACU':'T',
            'AAU':'N',
            'AGU':'S',
            'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
            'CCU':'P',
            'CAU':'H',
            'CGU':'R',
            'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
            'GCU':'A',
            'GAU':'D',
            'GGU':'G',
            'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
            'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
            'UAC':'Y', 'UAU':'Y', 'UAA':'*', 'UAG':'*',
            'UGC':'C', 'UGU':'C', 'UGA':'*', 'UGG':'W',}

        self.aa_table_r = {
            'A':['GCA', 'GCC', 'GCG', 'GCT', 'GCU'],
            'R':['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA', 'CGU'],
            'N':['AAC', 'AAT', 'AAU'],
            'D':['GAC', 'GAT', 'GAU'],
            'C':['TGC', 'TGT', 'UGC', 'UGU'],
            'Q':['CAA', 'CAG'],
            'E':['GAA', 'GAG'],
            'G':['GGT', 'GGC', 'GGA', 'GGG', 'GGU'],
            'H':['CAC', 'CAT', 'CAU'],
            'I':['ATT', 'ATC', 'ATA', 'AUU', 'AUC', 'AUA'],
            'L':['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG', 'CUA',
                 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'],
            'K':['AAA', 'AAG'],
            'M':['ATG', 'AUG'],
            'F':['TTC', 'TTT', 'UUC', 'UUU'],
            'P':['CCT', 'CCC', 'CCG', 'CCA', 'CCU'],
            'S':['TCA', 'TCC', 'TCG', 'TCT', 'AGC', 'AGT',
                 'UCA', 'UCC', 'UCG', 'UCU', 'AGU'],
            'T':['ACA', 'ACC', 'ACG', 'ACT', 'ACU'],
            'W':['TGG', 'UGG'],
            'Y':['TAT', 'TAC', 'UAC', 'UAU'],
            'V':['GTA', 'GTC', 'GTT', 'GTG', 'GUG', 'GUU',
                 'GUC', 'GUA'],
            '*':['TGA', 'TAG', 'TAA', 'UGA', 'UAG', 'UAA']
            }


        self.human_codon_frequency_bias_nakamura = {
            'TTT': 17.6, 'TCT': 15.2, 'TAT': 12.2, 'TGT': 10.6, 'TTC': 20.3,
            'TCC': 17.7, 'TAC': 15.3, 'TGC': 12.6, 'TTA': 7.7, 'TCA': 12.2,
            'TAA': 1.0, 'TGA': 1.6, 'TTG': 12.9, 'TCG':  4.4, 'TAG': 0.8,
            'TGG': 13.2, 'CTT': 13.2, 'CCT': 17.5, 'CAT': 10.9, 'CGT': 4.5,
            'CTC': 19.6, 'CCC': 19.8, 'CAC': 15.1, 'CGC': 10.4, 'CTA':  7.2,
            'CCA': 16.9, 'CAA': 12.3, 'CGA':  6.2, 'CTG': 39.6, 'CCG':  6.9,
            'CAG': 34.2, 'CGG': 11.4, 'ATT': 16.0, 'ACT': 13.1, 'AAT': 17.0,
            'AGT': 12.1, 'ATC': 20.8, 'ACC': 18.9, 'AAC': 19.1, 'AGC': 19.5,
            'ATA':  7.5, 'ACA': 15.1, 'AAA': 24.4, 'AGA': 12.2, 'ATG': 22.0,
            'ACG': 6.1, 'AAG': 31.9, 'AGG': 12.0, 'GTT': 11.0, 'GCT': 18.4,
            'GAT': 21.8, 'GGT': 10.8, 'GTC': 14.5, 'GCC': 27.7, 'GAC': 25.1,
            'GGC': 22.2, 'GTA':  7.1, 'GCA': 15.8, 'GAA': 29.0, 'GGA': 16.5,
            'GTG': 28.1, 'GCG': 7.4, 'GAG': 39.6, 'GGG': 16.5}

        self.human_codon_frequency_bias_nakamura_single = {
            'TTT': 17.6, 'TCT': 15.2, 'TAT': 12.2, 'TGT': 10.6, 'TTC': 20.3,
            'TCC': 17.7, 'TAC': 15.3, 'TGC': 12.6, 'TTA': 7.7, 'TCA': 12.2,
            'TAA': 1.0, 'TGA': 1.6, 'TTG': 12.9, 'TCG':  4.4, 'TAG': 0.8,
            'TGG': 13.2, 'CTT': 13.2, 'CCT': 17.5, 'CAT': 10.9, 'CGT': 4.5,
            'CTC': 19.6, 'CCC': 19.8, 'CAC': 15.1, 'CGC': 10.4, 'CTA':  7.2,
            'CCA': 16.9, 'CAA': 12.3, 'CGA':  6.2, 'CTG': 39.6, 'CCG':  6.9,
            'CAG': 34.2, 'CGG': 11.4, 'ATT': 16.0, 'ACT': 13.1, 'AAT': 17.0,
            'AGT': 12.1, 'ATC': 20.8, 'ACC': 18.9, 'AAC': 19.1, 'AGC': 19.5,
            'ATA':  7.5, 'ACA': 15.1, 'AAA': 24.4, 'AGA': 12.2, 'ATG': 22.0,
            'ACG': 6.1, 'AAG': 31.9, 'AGG': 12.0, 'GTT': 11.0, 'GCT': 18.4,
            'GAT': 21.8, 'GGT': 10.8, 'GTC': 14.5, 'GCC': 27.7, 'GAC': 25.1,
            'GGC': 22.2, 'GTA':  7.1, 'GCA': 15.8, 'GAA': 29.0, 'GGA': 16.5,
            'GTG': 28.1, 'GCG': 7.4, 'GAG': 39.6, 'GGG': 16.5}



        self.trna_ids = [
            'TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC',
            'TAC', 'TGC', 'TTA', 'TCA', 'TAA',
            'TGA', 'TTG', 'TCG', 'TAG', 'TGG', 'CTT', 'CCT', 'CAT',
            'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CTA',
            'CCA', 'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG', 'ATT',
            'ACT', 'AAT', 'AGT', 'ATC',
            'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ATG',
            'ACG', 'AAG', 'AGG', 'GTT',
            'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA',
            'GCA', 'GAA', 'GGA', 'GTG',
            'GCG', 'GAG', 'GGG']

        self.trna_ids.remove('TAG')
        self.trna_ids.remove('TAA')
        self.trna_ids.remove('TGA')

        trna_temp = [x.replace('T', 'U') for x in self.trna_ids]
        self.trna_ids = self.trna_ids + trna_temp




        self.trna_ids_vals = np.linspace(0, 60, 61).astype(int).tolist() +  np.linspace(0, 60, 61).astype(int).tolist()

        self.trna_dict = dict(zip(self.trna_ids, self.trna_ids_vals))
        self.id_to_trna = dict(map(reversed, self.trna_dict.items()))



        # add the U codons
        for key in list(self.human_codon_frequency_bias_nakamura.keys()):
            if 'T' in key:
                val = self.human_codon_frequency_bias_nakamura[key]
                newkey = key.replace('T', 'U')
                self.human_codon_frequency_bias_nakamura[newkey] = val



        self.human_codon_frequency_bias_nakamura_fast = {
            'GCT': 27.7, 'GCC': 27.7, 'GCA': 27.7, 'GCG': 27.7,  #A
            'CGT': 12.2, 'CGC': 12.2, 'CGA': 12.2, 'CGG': 12.2,
            'AGA': 12.2, 'AGG': 12.2,   # R
            'AAT': 19.1, 'AAC': 19.1,   #N
            'GAT': 25.1, 'GAC': 25.1,   # D
            'TGT': 12.6, 'TGC': 12.6,  # C
            'CAA': 34.2, 'CAG': 34.2,  # Q
            'GAA': 39.6, 'GAG': 39.6,  #E
            'GGT': 22.2, 'GGC': 22.2, 'GGA': 22.2, 'GGG': 22.2,  # G
            'CAT': 15.1, 'CAC': 15.1,  # H
            'ATT': 20.8, 'ATC': 20.8, 'ATA': 20.8,  # I
            'TTA': 39.6, 'TTG': 39.6, 'CTT': 39.6, 'CTC': 39.6,
            'CTA': 39.6, 'CTG': 39.6, # L
            'AAA': 31.9, 'AAG': 31.9,  # K
            'ATG': 22.0,   #M
            'TTT': 20.3, 'TTC': 20.3,    # F
            'CCT': 19.8, 'CCC': 19.8, 'CCA': 19.8, 'CCG': 19.8,  # P
            'TCT': 19.5, 'TCC': 19.5, 'TCA': 19.5, 'TCG': 19.5,
            'AGT': 19.5, 'AGC': 19.5,  # S
            'ACT': 18.9, 'ACC': 18.9, 'ACA': 18.9, 'ACG': 18.9, # T
            'TGG': 13.2,   #W
            'TAT': 15.3, 'TAC': 15.3,  # Y
            'GTT': 28.1, 'GTC': 28.1, 'GTA':28.1, 'GTG': 28.1,  # V
            'TAA': 1.6, 'TAG': 1.6, 'TGA':1.6 #STOP
            }


        for key in list(self.human_codon_frequency_bias_nakamura_fast.keys()):
            if 'T' in key:
                val = self.human_codon_frequency_bias_nakamura_fast[key]
                newkey = key.replace('T', 'U')
                self.human_codon_frequency_bias_nakamura_fast[newkey] = val

        self.human_codon_frequency_bias_nakamura_slow = {
            'GCT': 7.4, 'GCC': 7.4, 'GCA': 7.4, 'GCG': 7.4,  #A
            'CGT': 4.5, 'CGC': 4.5, 'CGA': 4.5, 'CGG': 4.5,
            'AGA':4.5, 'AGG':4.5,   #R
            'AAT': 17.0, 'AAC':17.0,  #%N
            'GAT': 21.8, 'GAC': 21.8,  #D
            'TGT': 10.6, 'TGC':10.6,  #C
            'CAA': 12.3, 'CAG': 12.3,  #Q
            'GAA': 29.0, 'GAG': 29.0,  #E
            'GGT': 10.8, 'GGC': 10.8, 'GGA': 10.8, 'GGG': 10.8,  #G
            'CAT': 10.9, 'CAC':10.9,  #H
            'ATT': 7.5, 'ATC': 7.5, 'ATA': 7.5, #I
            'TTA': 7.2, 'TTG':7.2, 'CTT': 7.2, 'CTC': 7.2,
            'CTA': 7.2, 'CTG': 7.2, #L
            'AAA': 24.4, 'AAG': 24.4, #K
            'ATG': 22.0, #M
            'TTT': 17.6, 'TTC': 17.6, #F
            'CCT': 6.9, 'CCC': 6.9, 'CCA': 6.9, 'CCG': 6.9, #P
            'TCT': 4.4, 'TCC': 4.4, 'TCA': 4.4, 'TCG': 4.4,
            'AGT': 4.4, 'AGC': 4.4, #S
            'ACT': 6.1, 'ACC': 6.1, 'ACA': 6.1, 'ACG': 6.1,#T
            'TGG': 13.2, #W
            'TAT': 12.2, 'TAC': 12.2, #Y
            'GTT': 7.1, 'GTC':7.1, 'GTA': 7.1, 'GTG': 7.1, # V
            'TAA': 0.8, 'TAG': 0.8, 'TGA': 0.8 #STOP CODON}
            }

        for key in list(self.human_codon_frequency_bias_nakamura_slow.keys()):
            if 'T' in key:
                val = self.human_codon_frequency_bias_nakamura_slow[key]
                newkey = key.replace('T', 'U')
                self.human_codon_frequency_bias_nakamura_slow[newkey] = val


        self.fast_codons_value = [27.7, 12.2, 19.1, 25.1, 12.6, 34.2, 39.6, 22.2, 15.1,
                                  20.8, 39.6, 31.9, 22, 20.3, 19.8, 19.5,
                                  18.9, 13.2, 15.3, 28.1, 1.6]


        self.slow_codons_value = [7.4, 4.5, 17, 21.8, 10.6, 12.3, 29, 10.8, 10.9, 7.5, 7.2,
                                  24.4, 22, 17.6, 6.9, 4.4, 6.1, 13.2, 12.2, 7.1, .8]

        self.fullcodonkeys = [
            'GCT', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGT', 'CAT',
            'ATT', 'TTA', 'AAA', 'ATG', 'TTT', 'CCT', 'TCT',
            'ACT', 'TGG', 'TAT', 'GTT', 'TAA',
            'GCU', 'CGU', 'AAU', 'GAU', 'UGU', 'CAA', 'GAA', 'GGU', 'CAU',
            'AUU', 'UUA', 'AAA', 'AUG', 'UUU', 'CCU', 'TCU',
            'ACU', 'UGG', 'UAU', 'GUU', 'UAA',]

        codonkeys = [
            'GCT', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGT', 'CAT',
            'ATT', 'TTA', 'AAA', 'ATG', 'TTT', 'CCT', 'TCT',
            'ACT', 'TGG', 'TAT', 'GTT', 'TAA',]

        self.sensitivity_fast_slow = []
        for i in range(len(codonkeys)):
            self.sensitivity_fast_slow.append(self.human_codon_frequency_bias_nakamura_fast[codonkeys[i]] / self.human_codon_frequency_bias_nakamura_slow[codonkeys[i]])

        self.load_tags()


        self.strCodonFreq = self.human_codon_frequency_bias_nakamura
        self.strCodonFreq_single = self.human_codon_frequency_bias_nakamura_single
        self.strCodonFreq_fast = self.human_codon_frequency_bias_nakamura_fast
        self.strCodonFreq_slow = self.human_codon_frequency_bias_nakamura_slow

    @property
    def mean_genecopynumber(self):
        '''
        Generate the mean copy number from the human_codon_frequency_bias_nakamura dictionary

        Returns
        -------
        float
            Mean gene copy number.

        '''
        keys = ['TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC', 'TTA', 'TCA', 'TAA', 'TGA', 'TTG', 'TCG', 'TAG', 'TGG', 'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CTA', 'CCA', 'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG', 'ATT', 'ACT', 'AAT', 'AGT', 'ATC', 'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ATG', 'ACG', 'AAG', 'AGG', 'GTT', 'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA', 'GCA', 'GAA', 'GGA', 'GTG', 'GCG', 'GAG', 'GGG']

        return sum([self.human_codon_frequency_bias_nakamura[x] for x in keys ])/len(keys)

    def load_tags(self):
        '''
        Opens custom_tags.txt and adds them to the dictionary

        Returns
        -------
        None.

        '''
        try:
            fname = open("custom_tags.txt", "r")
        except:
            return

        raw = fname.readlines()
        previous_tags = []

        for line in raw:
            if line != '\n':
                previous_tags.append(line)

        for line in previous_tags:

            custom_tag = line.strip('\n').split('---')

            if custom_tag[0] not in self.tag_dict.keys():
                self.tag_dict[custom_tag[0]] = custom_tag[2]

                self.tag_full[custom_tag[0]] = custom_tag[1]
        fname.close()


    def add_custom_tag(self, nt_seq, name):
        '''
        add a custom tag to custom_tags.txt and the tag dictionaries

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence of the tag / epitope.
        name : str
            name of the tag / epitope .

        Returns
        -------
        None.

        '''

        fname = open("custom_tags.txt", "r")

        raw = fname.readlines()
        previous_tags = []

        for line in raw:
            if line != '\n':
                previous_tags.append(line)

        if not set(nt_seq.lower()).issubset(set(['a', 't', 'c', 'g', 'u'])):
            print('invalid NT sequence')
            fname.close()
            return


        aa_seq = ''
        for i in range(0, len(nt_seq), 3):
            aa_seq += self.aa_table[nt_seq[i:i+3]]


        newtag = name+'---'+ nt_seq.lower() + '---'+ aa_seq.upper() + '\n'

        if newtag not in previous_tags:
            previous_tags.append(newtag)
        fname.close()

        fname = open("custom_tags.txt","w+")

        for item in previous_tags:
            fname.write('%s' % item)

        fname.close()
