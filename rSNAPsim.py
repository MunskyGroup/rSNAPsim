# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 09:42:24 2018

@author: William
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 06 13:40:30 2018

rSNAPsim CHANGE TOO THIS WHEN PUSHING TO PIP INSTALLLLLLS


@author: William
"""
import re #import regex
import sys
'''
try:
   # sys.path.append('C:\\Users\\wsraymon\\Github\\ssa_cpp\\translation_ssa')

    import ssa_translation
except:
    try:
        sys.path.append('C:\\Users\\william\\Documents\\Github\\ssa_cpp\\translation_ssa')
        import ssa_translation
    except:
        pass

    pass

'''

import os
import time
import json, codecs

from scipy import sparse
from scipy.stats import pearsonr

cwd = os.getcwd()
try:
    os.chdir('../..')
    os.chdir('ssa_cpp')
except:
    os.chdir(cwd)
    os.chdir('ssa_cpp')  

try:
    
    import ssa_translation
    
except:
    pass
os.chdir('..')
os.chdir(cwd)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.patches as mpatches
import matplotlib.animation as animation
from matplotlib.collections import PatchCollection
from matplotlib import cm
from matplotlib import gridspec
from matplotlib.patches import Ellipse
#import scipy.stats.trim_mean as tmean

from scipy.stats import kde

try:
    from Bio import SeqIO
    from Bio import Entrez
except:
    pass




class rSNAPsim():

    """
    The Single Molecule Simulator (SMS) provides a python class for running
    single molecule mRNA translation simulations

    When presented with a valid protein sequence the SMS can find open reading frames
    and simulate intensity trajectories from translation of the protein with given fluorescent tags.

    *model description*

        link to paper here / image


    *main functions*

        -open_seq_file(filepath), opens a txt or .gb file and gets the sequence

        -get_orfs(nt_sequence, min_codons), returns open reading frames of a given
        sequence and a minimum codon length per protein

        -get_temporal_proteins(), gets the proteins after get_orfs

        -analyze_poi(aa_seq,nt_seq), analyzes the proteins of intrest for
        codon sensitivity and elongation rates

        -__.poi(), class to contain proteins of intrest after analyzed

        -run_default(), runs get_orfs, get_temporal proteins, and analyze_poi
        with the first protien found in the sequence



    *attributes*

        **gene_sequence_str** = string of the nucleotide sequence
        **tag_dict** = dictionary with various types of fluorescent tag epitopes
        
        **tag_full** = dictionary of full tag sequences
        
        **aa_keys** = amino acid single letter keys
        
        **codon_types** = flag dictionary of which amino acids are set to Wild-type, fast, or slow
        
        **aa_table** = dictionary of amino acids
        
        **aa_table_r** = reverse dictionary (amino acid letters are the keys)
        
        **strGeneCopy** = dictionary of wild-type tRNA copy numbers
        
        **strGeneCopy_fast** = dictionary of fast tRNA copy numbers
        
        **strGeneCopy_slow** = dictionary of slow tRNA copy numbers
        
        **slow_codons_value** = list of slowest codon tRNA copy numbers
        
        **fast_codons_value** = list of fastest codon tRNA copy numbers
        
        **sensitivity_fast_slow** = list of sensitivity for amino acids
        
        **poi** = Class container for proteins of intrest
        
        **orfs** = dictionary of open reading frames with keys 1,2,3
        
        **seq_str** = sequence string
        
        **proteins** = dictionary of proteins detected in the sequence by ORF
        
        **tagged_proteins** = dictionary of proteins that were detected and tagged


    *POI*

        Protein of intrest has the following attributes:

        **aa_seq** = amino acid sequence
        
        **nt_seq** = nucleotide sequence
        
        **gene_length** = length of the gene
        
        **tag_length** = length of the tags
        
        **total_length** = total length of the full amino acid sequence
        
        **name** = name of the gene
        
        **tag_types** = what types of tags does the protien have
        
        **tag_epitopes** = type of tags and epitope lists per tag
        
        **codon_sensitivity** = how sensitive is the protein per amino acid sequence?
        
        **CAI** = codon activation index
        
        **CAI_codons** = means of the codon activation

    *ssa*

        The ssa container class has the following attributes:


        **no_ribosomes** = number of ribosomes
        
        **n_traj** = number of trajectories
        
        **k** = all kelongation rates (calculated from codon sequence)
        
        **no_rib_per_mrna** = number of ribosomes per mRNA strand on average
        
        **rib_density** = ribosome density
        
        **rib_means** = ribosome means
        
        **rib_vec** = raw ribosome location matrix for each trajectory
        
        **intensity_vec** = fluorescence intensities
        
        **time_vec_fixed** = the time vector
        
        **start_time** = the time the simulation was started

        **evaluating_inhibitor** = was there an inhibitor present?
        
        **evaluating_frap** = was the simulation subjected to a FRAP test
        
        **time_inhibit** = the time of the perturbation

        **autocorr_vec** = autocorrelation vector of intensities
        
        **mean_autocorr** = the average autocorrelations, averaged over trajectories
        
        **error_autocorr** = the standard deviation of the autocorrelation
        
        **dwell_time** = how long do the ribosomes stay on the mRNA strand calculated by the simulation
        
        **ke_sim** = the calculated average elongation rate from the simulations


    """

    def __init__(self):
        self.gene_sequence_str = ''
        

        self.tag_dict = {'T_SunTag':'EELLSKNYHLENEVARLKK',
                         'T_Flag':'DYKDDDDK',
                         'T_Hemagglutinin':'YPYDVPDYA'}
        
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

        self.aa_keys = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
                       'P', 'S', 'T', 'W', 'Y', 'V', '*']

        self.codon_types = dict(zip(self.aa_keys, np.ones((1, 21)).flatten().astype(int).tolist()))

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
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
            'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
            'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
            'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
            'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
            'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
            'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
            'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
            'UAC':'Y', 'UAU':'Y', 'UAA':'*', 'UAG':'*',
            'UGC':'C', 'UGU':'C', 'UGA':'*', 'UGG':'W',}

        self.aa_table_r = {'A':['GCA', 'GCC', 'GCG', 'GCT','GCU'],
                          'R':['CGA', 'CGC', 'CGG', 'CGT','AGG','AGA','CGU'],
                          'N':['AAC', 'AAT','AAU'],
                          'D':['GAC', 'GAT','GAU'],
                          'C':['TGC', 'TGT','UGC','UGU'],
                          'Q':['CAA', 'CAG'],
                          'E':['GAA', 'GAG'],
                          'G':['GGT', 'GGC', 'GGA', 'GGC','GGU'],
                          'H':['CAC', 'CAT','CAU'],
                          'I':['ATT', 'ATC', 'ATA','AUU','AUC','AUA'],
                          'L':['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG','CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'],
                          'K':['AAA', 'AAG'],
                          'M':['ATG','AUG'],
                          'F':['TTC', 'TTT','UUC','UUU'],
                          'P':['CCT', 'CCC', 'CCG', 'CCA','CCU'],
                          'S':['TCA', 'TCC', 'TCG', 'TCT','AGC','AGT','UCA','UCC','UCG'],
                          'T':['ACA', 'ACC', 'ACG', 'ACT','ACU'],
                          'W':['TGG','UGG'],
                          'Y':['TAT', 'TAC','UAC','UAU'],
                          'V':['GTA', 'GTC', 'GTT','GTG','GUG','GUU','GUC','GUA'],
                          '*':['TGA', 'TAG', 'TAA','UGA','UAG','UAA']
                         }


        self.strGeneCopy = {'TTT': 17.6, 'TCT': 15.2, 'TAT': 12.2, 'TGT': 10.6, 'TTC': 20.3,
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

        
        # add the U codons
        for key in list(self.strGeneCopy.keys()):
            if 'T' in key:
                val = self.strGeneCopy[key]
                newkey = key.replace('T','U')
                self.strGeneCopy[newkey] = val



        self.strGeneCopy_fast = {'GCT': 27.7, 'GCC': 27.7, 'GCA': 27.7, 'GCG': 27.7,  #A
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


        for key in list(self.strGeneCopy_fast.keys()):
            if 'T' in key:
                val = self.strGeneCopy_fast[key]
                newkey = key.replace('T','U')
                self.strGeneCopy_fast[newkey] = val

        self.strGeneCopy_slow = {'GCT': 7.4, 'GCC': 7.4, 'GCA': 7.4, 'GCG': 7.4,  #A
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
        
        for key in list(self.strGeneCopy_slow.keys()):
            if 'T' in key:
                val = self.strGeneCopy_slow[key]
                newkey = key.replace('T','U')
                self.strGeneCopy_slow[newkey] = val


        self.fast_codons_value = [27.7, 12.2, 19.1, 25.1, 12.6, 34.2, 39.6, 22.2, 15.1,
                                  20.8, 39.6, 31.9, 22, 20.3, 19.8, 19.5,
                                  18.9, 13.2, 15.3, 28.1, 1.6]


        self.slow_codons_value = [7.4, 4.5, 17, 21.8, 10.6, 12.3, 29, 10.8, 10.9, 7.5, 7.2,
                                  24.4, 22, 17.6, 6.9, 4.4, 6.1, 13.2, 12.2, 7.1, .8]

        fullcodonkeys = ['GCT', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGT', 'CAT',
                     'ATT', 'TTA', 'AAA', 'ATG', 'TTT', 'CCT', 'TCT',
                     'ACT', 'TGG', 'TAT', 'GTT', 'TAA',
                     'GCU', 'CGU', 'AAU', 'GAU', 'UGU', 'CAA', 'GAA', 'GGU', 'CAU',
                     'AUU', 'UUA', 'AAA', 'AUG', 'UUU', 'CCU', 'TCU',
                     'ACU', 'UGG', 'UAU', 'GUU', 'UAA',                     ]

        codonkeys = ['GCT', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGT', 'CAT',
                     'ATT', 'TTA', 'AAA', 'ATG', 'TTT', 'CCT', 'TCT',
                     'ACT', 'TGG', 'TAT', 'GTT', 'TAA']
        self.sensitivity_fast_slow = []
        for i in range(len(codonkeys)):
            self.sensitivity_fast_slow.append(self.strGeneCopy_fast[codonkeys[i]] / self.strGeneCopy_slow[codonkeys[i]])


    def __update_sensitivity(self):

        """
        updates sensitivities for the GUI implementation call
        """

        self.fast_codons_value = []
        for key in self.aa_keys:
            values = []
            codons = self.aa_table_r[key]
            for codon in codons:
                values.append(self.strGeneCopy[codon])

            self.fast_codons_value.append(max(values))

            for codon in codons:
                self.strGeneCopy_fast[codon] = max(values)



        self.slow_codons_value = []
        for key in self.aa_keys:
            values = []
            codons = self.aa_table_r[key]
            for codon in codons:
                values.append(self.strGeneCopy_slow[codon])

            self.slow_codons_value.append(min(values))

            for codon in codons:
                self.strGeneCopy_slow[codon] = min(values)

        codonkeys = ['GCT', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGT', 'CAT', 'ATT',
                     'TTA', 'AAA', 'ATG', 'TTT', 'CCT', 'TCT', 'ACT', 'TGG', 'TAT', 'GTT', 'TAA']

        self.sensitivity_fast_slow = []
        for i in range(len(codonkeys)):
            self.sensitivity_fast_slow.append(self.strGeneCopy_fast[codonkeys[i]] / self.strGeneCopy_slow[codonkeys[i]])


    def nt2aa(self, nt_seq):
        '''
        Translates nucleotides sequences to amino acid sequences

        *args*

            **nt_seq**, nucleotide sequence as a string

        *returns*

            **aa_seq**, amino acid sequence as string
        '''

        aa = ''
        for i in range(0, len(nt_seq), 3):
            aa += self.aa_table[nt_seq[i:i+3]]
        return aa


    def get_orfs(self, nt_seq='', min_codons=80):

        '''
        Returns open reading frames of the nucleotide sequence given

        orfs = {'1':[proteins],
                '2':[proteins],
                '3':[proteins]}

        *keyword args*

            **nt_seq**, nucleotide sequence as a string. If left blank uses
            the self.sequence_str

            **min_codons**, minimum amount of codons to be considered
            a protein in the open reading frame
        '''

        if nt_seq == '':
            nt_seq = self.sequence_str

        allstarts = np.array([m.start() for m in re.finditer('(?=A[TU]G((?:.{3})+?)[TU](?:AG|AA|GA))', nt_seq)])
        
 
        #allsegments = re.findall('(?=A[TU]G((?:.{3})+?)[TU](?:AG|AA|GA))',self.sequence_str)
        allstops = np.array([m.start() for m in re.finditer('(?=[TU](?:AG|AA|GA))', nt_seq)])
        start_frames = allstarts%3
        stop_frames = allstops%3
        min_len = min_codons*3
        orf1_starts = allstarts[np.where(start_frames == 0)]
        orf2_starts = allstarts[np.where(start_frames == 1)]
        orf3_starts = allstarts[np.where(start_frames == 2)]

        orf1_stops = allstops[np.where(stop_frames == 0)]
        orf2_stops = allstops[np.where(stop_frames == 1)]
        orf3_stops = allstops[np.where(stop_frames == 2)]

        self.starts = [orf1_starts, orf2_starts, orf3_starts]
        self.stops = [orf1_stops, orf2_stops, orf3_stops]
        self.orfs = {'1':[], '2':[], '3':[]}


        self.orfs = {'1':[], '2':[], '3':[]}
        laststop = 0
        for start in orf1_starts:
            nextstop = orf1_stops[np.where(orf1_stops > start)[0][0]]
            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    self.orfs['1'].append((start, nextstop))

                    laststop = nextstop

        laststop = 0
        for start in orf2_starts:
            nextstop = orf2_stops[np.where(orf2_stops > start)[0][0]]
            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    self.orfs['2'].append((start, nextstop))
                    laststop = nextstop

        laststop = 0
        for start in orf3_starts:
            nextstop = orf3_stops[np.where(orf3_stops > start)[0][0]]

            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    self.orfs['3'].append((start, nextstop))
                    laststop = nextstop





    def get_k_construct(self, nt_seq, k_init, k_elong_mean, codon_types=None):
        '''
        Returns the k_elongation rates of a given nucleotide sequence under constructed conditions
        given some sort of key describing which amino acids are slow, fast or natural

        *args*

            **nt_seq**, nucleotide sequence to get the propensities of

            **k_init**, initiation rate of starting translation

            **k_elong_mean**, average rate of elongation for the protein translation


        *keyword args*

            **codon_types**, a dictonary or identifier determining which amino acids are slow, fast or natural

                self.codon_types is an example dictionary for the user to change / utilize, if codon_types is left blank
                get_k_construct uses this internal dictonary

                ex: codon_types = 'slow' or 'rare'  all amino acids set to slow
                    codon_types = 'fast' or 'common'  all amino acids set to fast
                    codon_types = 'natural' all amino acids set to fast

                    codon_types = {'A':[0], 'T':[2]}  A set to slow, T set to fast
                    codon_types = {'rare':['A','R'],'common':['L']}  A and R set to slow, L set to fast


        '''



        if codon_types == None:
            codon_types = self.codon_types
        else:
            all_natural = dict(zip(self.aa_keys, np.ones((1, 20)).flatten().astype(int).tolist()))

            if isinstance(codon_types, str):
                if codon_types == 'rare' or codon_types == 'slow':
                    all_natural = dict(zip(self.aa_keys, np.zeros((1, 20)).flatten().astype(int).tolist()))
                if codon_types == 'common' or codon_types == 'fast':
                    all_natural = dict(zip(self.aa_keys, (2*np.ones((1, 20))).flatten().astype(int).tolist()))
            if isinstance(codon_types, dict):
                for key in codon_types.keys():
                    if isinstance(key, str):
                        if key.lower() not in ['rare', 'common', 'natural']:
                            if key.upper() in self.aa_keys:
                                if codon_types[key] in [0, 1, 2]:
                                    all_natural[key] = key
                                if codon_types[key] in ['rare', 'common', 'natural']:
                                    if codon_types[key] == 'rare':
                                        all_natural[key] = 0
                                    if codon_types[key] == 'common':
                                        all_natural[key] = 2
                                    if codon_types[key] == 'natural':
                                        all_natural[key] = 1
                        else:
                            newkeys = codon_types[key]
                            for newkey in newkeys:
                                if newkey.upper() in self.aa_keys:
                                    if key.lower() == 'rare':
                                        all_natural[newkey.upper()] = 0
                                    if key.lower() == 'common':
                                        all_natural[newkey.upper()] = 2
                                    if key.lower() == 'natural':
                                        all_natural[newkey.upper()] = 1


                    if isinstance(key, int):
                        newkeys = codon_types[key]
                        for newkey in newkeys:
                            all_natural[newkey] = key




            codon_types = all_natural


        aa_seq = self.nt2aa(nt_seq)



        tRNA_design = np.zeros((1, len(aa_seq)))
        tRNA_norm = np.zeros((1, len(aa_seq)))

        seperated_codons = [nt_seq[i:i+3] for i in range(0, len(nt_seq), 3)] #split codons by 3


        for i in range(len(seperated_codons)):
            tRNA_norm[0, i] = self.strGeneCopy[seperated_codons[i]]





        for i in range(len(self.aa_keys)-1):

            fs = codon_types[self.aa_keys[i]]
            indexes = [m.start() for m in re.finditer(self.aa_keys[i], aa_seq)]
            for index in indexes:

                if fs == 0:
                    tRNA_design[0, index] = self.slow_codons_value[i]
                if fs == 2:
                    tRNA_design[0, index] = self.fast_codons_value[i]
                if fs == 1:
                    tRNA_design[0, index] = tRNA_norm[0, index]


        tRNA_design[0, -1] = tRNA_norm[0, -1]




        mean_tRNA_copynumber = np.mean(list(self.strGeneCopy.values()))



        k_elongation_design = (tRNA_design / mean_tRNA_copynumber) * k_elong_mean

        all_k_design = [k_init] + k_elongation_design.flatten().tolist() + [k_elong_mean]

        return all_k_design
    
    
    def get_ui(self, nt_seq):
        mean_u = np.mean(self.strGeneCopy.values()) 
        ui = []
        for i in range(0, len(nt_seq), 3):
            ui.append(mean_u/ self.strGeneCopy[nt_seq[i:i+3]])
        return ui
        



    def get_k(self, nt_seq, k_init, k_elong_mean):
        '''
        returns all propensities for a given nucleotide sequence

        *args*

            **nt_seq**,  nucleotide sequence as a string

            **k_initiation**, initiation rate of ribosome binding

            **k_elong_mean**, average rate of elgonation experimentally found


        '''
        codons = nt_seq
        genelength = int(len(codons)/3)
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
        k_elongation = np.zeros((1, genelength))
        tRNA_copynumber = np.zeros((1, genelength))

     
        for i in range(len(seperated_codons)):
            tRNA_copynumber[0, i] = self.strGeneCopy[seperated_codons[i]]

        mean_tRNA_copynumber = np.mean(list(self.strGeneCopy.values()))

        k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean
        all_k = [k_init] + k_elongation.flatten().tolist()[:-1] + [10]

        return all_k



    def get_temporal_proteins(self):
        '''
        gets all the temporal proteins after getting the ORFs

        __.tagged_proteins = dictionary with keys of tag types and a list of proteins
        __.pois = list of proteins of intrest
        __.pois_seq = list of nucleotide sequences of proteins of sequences
        __.proteins = dictonary with keys of 1 2 or 3 orfs

        '''


        self.proteins = {'1':[], '2':[], '3':[]}
        self.tagged_proteins = {a:[] for a in self.tag_dict.keys()}
        self.tagged_protein_seq = {a:[] for a in self.tag_dict.keys()}

        for i in range(len(self.orfs)):
            for j in range(len(self.orfs[str(i+1)])):
                pro = self.nt2aa(self.sequence_str[self.orfs[str(i+1)][j][0]:self.orfs[str(i+1)][j][1]+3])
                nt_seq = self.sequence_str[self.orfs[str(i+1)][j][0]:self.orfs[str(i+1)][j][1]+3]
                self.proteins[str(i+1)].append(pro)
                for tag in self.tag_dict.keys():
                    if self.tag_dict[tag] in pro:
                        self.tagged_protein_seq[tag].append(nt_seq)
                        self.tagged_proteins[tag].append(pro)

        tags = 0
        for key in self.tagged_proteins.keys():
            tags += len(self.tagged_proteins[key])


        self.pois = []
        self.pois_seq = []
        for tag in self.tag_dict.keys():
            for i in range(len(self.tagged_proteins[tag])):
                if self.tagged_proteins[tag][i] not in self.pois:
                    self.pois.append(self.tagged_proteins[tag][i])
                    self.pois_seq.append(self.tagged_protein_seq[tag][i])

        if len(self.pois) == 0:

            POIs = []
            pois_s = []
            pois_nt = []
            for i in range(len(self.gb_obj.features)):

                try:

                    self.gb_obj.features[i].qualifiers['translation']

                    if tags == 0:

                        POIs.append(self.gb_obj.features[i])
                        pois_s.append(self.nt2aa(self.tag_full['T_Flag']) + self.gb_obj.features[i].qualifiers['translation'][0])
                        pois_nt.append(self.tag_full['T_Flag'] + str(self.gb_obj.seq)[int(self.gb_obj.features[i].location.start):int(self.gb_obj.features[i].location.end)])
                    else:

                        POIs.append(self.gb_obj.features[i])
                        pois_s.append(self.gb_obj.features[i].qualifiers['translation'][0])
                        pois_nt.append(str(self.gb_obj.seq)[int(self.gb_obj.features[i].location.start):int(self.gb_obj.features[i].location.end)])

                except:
                    pass


            self.pois = pois_s
            self.pois_seq = pois_nt


    def analyze_poi(self, protein, sequence):
        '''
        Analyzes the protein of intrest and stores it in __.POI

        *args*

            **protein**,  amino acid sequence as a string

            **sequence**, nucleotide sequence that goes with the protein


        '''

        self.POI = poi()
        self.POI.nt_seq = sequence
        self.POI.aa_seq = protein
        self.POI.name = self.sequence_name
        self.POI.total_length = len(protein)

        '''
        for key in self.tagged_proteins:
            if protein in self.tagged_proteins[key]:
                self.POI.tag_types.append(key)
        '''
        self.POI.tag_types = []
        for tag in self.tag_dict.keys():
            if self.tag_dict[tag] in protein:
                self.POI.tag_types.append(tag)

                #''.join(sms.poi[0].split('DYKDDDDK')

        self.POI.tag_epitopes = {a:[] for a in self.POI.tag_types}
        gs = protein


        for i in range(len(self.POI.tag_types)):
            
            try:
                nt_tag = self.tag_full[self.POI.tag_types[i]]
                aa_tag = self.nt2aa(nt_tag)
            except:
                epi = self.tag_dict[self.POI.tag_types[i]]
                firstep = self.POI.aa_seq.find(epi) 
                lastep = len(self.POI.aa_seq) - self.POI.aa_seq[::-1].find(epi[::-1])                
                aa_tag = self.POI.aa_seq[firstep:lastep]
                nt_tag = self.POI.nt_seq[3*firstep:3*lastep]
                
                
                
            self.POI.tag_epitopes[self.POI.tag_types[i]] = [m.start()+1 for m in re.finditer(self.tag_dict[self.POI.tag_types[i]], self.POI.aa_seq)]

            gs = gs.replace(aa_tag, '')





        self.POI.gene_seq = gs
        self.POI.gene_length = len(gs)
        codons = []
        for i in range(0, len(sequence), 3):
            codons.append(sequence[i:i+3])
        self.POI.codons = codons

        self.POI.codon_sensitivity, self.POI.CAI, self.POI.CAI_codons = self.codon_usage(self.POI.nt_seq)


    def open_seq_file(self, seqfile):
        '''
        Reads a sequence file, either a .txt file or a .gb genbank file

        *args*

            **seqfile**, sequence file either in txt, gb, gbk format


        '''
        seq = seqfile
        if '.txt' in seq:
            with open(seq) as f:
                raw = f.readlines()



            raw = ''.join(raw)
           

            
            
            onlychar = re.split(r'[^A-Za-z]', raw)
            validt = ['A', 'G', 'T', 'C']
            validu = ['A', 'G', 'U', 'C']
            namelen = 0
            self.sequence_str = ''
       
            for i in range(len(onlychar)):
                section = onlychar[i]
                if set(section.upper()) == set(validt):
                    self.sequence_str += section.upper()
                    
                
                elif set(section.upper()) == set(validu):
                    self.sequence_str += section.upper()                    
                else:
                    if len(section)>namelen:
                        self.sequence_name = section
                        namelen = len(section)
                    
                    
                    


        if '.gb' in seq:
            gb_record = SeqIO.read(open(seq, "r"), "genbank")
            self.sequence_str = str(gb_record.seq)
            self.sequence_name = gb_record.name
            self.gb_obj = gb_record





    def codon_usage(self, nt_seq):
        '''
        Analyzes codon useage from the nucleotide sequence

        *args*

            **nt_seq**,  nucleotide sequence as a string

        *returns*

            **codon_sensitivity**, a list of codon sensitivity for the nucleotide sequence

            **cai**, cai value

        '''
        codon_usage = np.zeros((1, 21))
        gene_len = len(nt_seq)/3
        aa_seq = self.nt2aa(nt_seq)

        for i in range(len(self.aa_keys)-1):

            codon_usage[0, i] = len(re.findall(self.aa_keys[i], aa_seq))
        codon_usage[0, 20] = len(re.findall('\*', aa_seq))
        codon_norm = codon_usage/gene_len
        codon_sensitivity = np.round(codon_norm*self.sensitivity_fast_slow, 2)

        cai_codons = []
        for i in range(0, len(nt_seq), 3):
            cai_codons.append(self.strGeneCopy[nt_seq[i:i+3]] / self.strGeneCopy_fast[nt_seq[i:i+3]])

        cai = self.geomean(cai_codons)

        return codon_sensitivity, cai, cai_codons


    def get_probvec(self): 
        '''
        returns the probe vectors (epitope positions by codon position) associated with the tagged sequence stored in POI

        *returns*

            **probe_vec**, cumlative probe intensity vector by codon position. Ex: [0,0,0,0,1,1,1,1,2,2,2,3,3,3 etc]

            **probe_loc**, epitope posistion as a binary vector, 1 for epitope pos, 0 for everything else
        '''

        probePositions = []
        
        keylist = list(self.POI.tag_epitopes.keys())
        
        for n in range(len(keylist)):
            probePosition = []
            key = keylist[n]
            
            probePosition = probePosition + self.POI.tag_epitopes[key]
            
            if probePosition != []:
                probePosition = np.unique(probePosition).tolist()
                probePositions.append(probePosition)
                
        
        

        genelength = self.POI.total_length

        pvfull = np.zeros((1, genelength+1)).astype(int).flatten()
        
        if len(probePositions) > 1:
            k = 0
            for n in range(len(keylist)):
                pv = np.zeros((1, genelength+1)).astype(int).flatten()
                key = keylist[n]
                probePosition = probePositions[k] 
                k+=1
                if len(self.POI.tag_epitopes[key]) != 0:
                    for i in range(len(probePosition)):
                        pv[probePosition[i]:] = i+1
                    if n > 0:
                        pvfull = np.vstack((pvfull,pv))
                    else:
                        pvfull = pv
        else:
            probePosition = probePositions[0]
            for n in range(len(keylist)):
                pv = np.zeros((1, genelength+1)).astype(int).flatten()
                key = keylist[n]
                if len(self.POI.tag_epitopes[key]) != 0:
                    for i in range(len(probePosition)):
                        pv[probePosition[i]:] = i+1
                    if n > 0:
                        pvfull = np.vstack((pvfull,pv))
                    else:
                        pvfull = pv
        numtags = 0
        for key in keylist:
            if len(self.POI.tag_epitopes[key]) != 0:
                numtags += 1

        ploc = np.zeros((numtags, self.POI.total_length+1)).astype(int)

        numind = 0
        for n in range(len(keylist)):
            key = keylist[n]
            if len(self.POI.tag_epitopes[key]) != 0:
                ploc[numind][self.POI.tag_epitopes[key]] = 1

                numind += 1

        return pvfull, ploc


    def ssa_solver(self, nt_seq=None, all_k=None, k_elong_mean=10, k_initiation=.03, probePosition=[], n_traj=100, tf=1000, start_time=0, tstep=1000, time_inhibit=0, evaluating_frap=False, evaluating_inhibitor=False,force_python = False):
        '''
        Solve stochastic simulation algorithms (SSA) for the translation simulation.

        *keyword args*

            **nt_seq**, nucleotide sequence to simulate

            **all_k**, the propensity rates for each codon location (obtained via get_k)

            **k_elong_mean**, average elongation rate to normalize by

            **k_initiation**, rate of mRNA translation initiation

            **probePosition**, binary vector of probe positions, i.e. where the tag epitopes start by codon position

            **n_traj**, number of trajectories

            **tf**, final time point

            **tstep**, number of time steps to record from 0 to tf

            **time_inhibit**, inhibition time of translation either, harringtonine assay or FRAP

            **evaluating_frap**, true or false for evaluating frap assay at time_inhibit

            **evaluating_inhibitor**, true or false for evaluating harringtonine at time_inhibit

        *returns*

            **ssa_obj**, a ssa() class containing the raw ribosome posistions simulated and statistics such as intensity vectors from the SSA trajectory group

        '''

        if len(probePosition) == 0:
            '''
            try:
                probePosition = []
                for key in self.POI.tag_epitopes.keys():
                    probePosition = probePosition + self.POI.tag_epitopes[key]
                probePosition = np.unique(probePosition).tolist()
            except:
                print('No POI found')
                #nt_seq = self.tag_full['T_flag'] + nt_seq
            '''
            
            pv,probePosition = self.get_probvec()
        

        if nt_seq == None:
          nt_seq = self.POI.nt_seq
        genelength = int(len(nt_seq)/3)

        if all_k == None:


            codons = nt_seq
            genelength = int(len(codons)/3)
            seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
            k_elongation = np.zeros((1, genelength))
            tRNA_copynumber = np.zeros((1, genelength))

            for i in range(len(seperated_codons)):
                tRNA_copynumber[0, i] = self.strGeneCopy[seperated_codons[i]]

            mean_tRNA_copynumber = np.mean(list(self.strGeneCopy.values()))

            k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean
            all_k = [k_initiation] + k_elongation.flatten().tolist()[:-1] + [10]

        
        if isinstance(probePosition,list):
            probePosition = np.array([probePosition]).astype(int)
            
     

        non_consider_time = start_time
      
        '''
        if probePosition.shape[0] <= 1:
            pv = np.zeros((1, genelength+1)).astype(int).flatten()
            
            for i in range(len(probePosition[0])):
                pv[probePosition[0][i]:] = i+1
        else:
            pv = np.zeros((probePosition.shape[0], genelength+1)).astype(int)
            for j in range(probePosition.shape[0]):
                for i in range(len(probePosition)):
                    pv[j][probePosition[j][i]:] = i+1      
        '''

        npoints = tstep #non_consider_time + tstep
        
        time_vec_fixed = np.linspace(0, npoints-1, npoints, dtype=np.float64)
        truetime = np.linspace(0, tf, tstep, dtype=np.float64)
        
        print(truetime[0],truetime[-1])

        rib_vec = []

        solutions = []
        


        evf = int(evaluating_frap)
        evi = int(evaluating_inhibitor)
        try:
            intime = float(time_inhibit)
        except:
            intime = 0

#        if evaluating_frap == True or evaluating_inhibitor == True:
#            for i in range(nRepetitions):
#
#                soln = self.SSA(all_k,time_vec_fixed,inhibit_time=time_inhibit+non_consider_time,FRAP=evaluating_frap,Inhibitor=evaluating_inhibitor)
#                solutions.append(soln)
#        else:

        solutionssave = []
        
        st = time.time() 
        
        try:
            if force_python == True:
                st[0]
                
            rib_vec = []
    
            solutions = []            
            solutionssave = []
            N_rib = 200
            all_results = np.zeros((n_traj, N_rib*len(time_vec_fixed)), dtype=np.int32)
            all_ribtimes = np.zeros((n_traj,int(1.3*all_k[0]*truetime[-1])),dtype=np.float64)
            result = np.zeros((len(time_vec_fixed)*N_rib), dtype=np.int32)
            nribs = np.array([0],dtype=np.int32)
            k = np.array(all_k)
            seeds = np.random.randint(0, 0x7FFFFFF, n_traj)
            all_frapresults = np.zeros((n_traj,N_rib*len(time_vec_fixed)),dtype=np.int32)
            all_collisions = np.zeros((n_traj,int(1.3*all_k[0]*truetime[-1])),dtype=np.int32)
            all_nribs = np.zeros((n_traj,1))
            all_col_points = []
            for i in range(n_traj):
                result = np.zeros((len(time_vec_fixed)*N_rib), dtype=np.int32)
                ribtimes = np.zeros((int(1.3*k[0]*truetime[-1])),dtype=np.float64)
                frapresult = np.zeros((len(time_vec_fixed)*N_rib),dtype=np.int32)
                coltimes = np.zeros((int(1.3*k[0]*truetime[-1])),dtype=np.int32)
                colpointsx = np.zeros(len(k[1:-1])*(int(1.3*k[0]*truetime[-1])),dtype=np.int32)
                colpointst = np.zeros(len(k[1:-1])*(int(1.3*k[0]*truetime[-1])),dtype=np.float64)
                
                nribs = np.array([0],dtype=np.int32)
                
                ssa_translation.run_SSA(result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
                #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
                all_results[i, :] = result
                all_frapresults[i,:] = frapresult
                all_ribtimes[i,:] = ribtimes
                all_collisions[i,:] = coltimes
                all_nribs[i,:] = nribs
                
                endcolrec = np.where(colpointsx == 0)[0][0]
                
                colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
                all_col_points.append(colpoints.T)
                
                
    
            for i in range(n_traj):
                soln = all_results[i, :].reshape((N_rib, len(time_vec_fixed)))
                validind = np.where(np.sum(soln,axis=1)!=0)[0]
                if np.max(validind) != N_rib-1:
                    validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)
            
                so = soln[(validind,)]
                
                solutionssave.append(so)
                solutions.append(soln)
            
            collisions = np.array([[]])
            watched_ribs = []
            for i in range(n_traj):
                totalrib = all_nribs[i]
            
                if totalrib > all_collisions.shape[1]:
                    collisions = np.append(collisions, all_collisions[i][:])
                    watched_ribs.append(int(all_collisions.shape[1]))
            
                else:
                   
                    collisions = np.append(collisions, all_collisions[i][:int(totalrib[0])])
                    watched_ribs.append(int(totalrib[0]))
            
            sttime = time.time() - st

        except:
            
            print('C++ library failed, Using Python Implementation')
            rib_vec = []
    
            solutions = []            
            solutionssave = []
            N_rib = 200
            collisions = np.array([[]])
            all_results = np.zeros((n_traj, N_rib*len(time_vec_fixed)), dtype=np.int32)
            all_col_points = []
            watched_ribs = []
            for i in range(n_traj):
                
                soln,all_ribtimes,Ncol,col_points = self.SSA(all_k, truetime, inhibit_time=time_inhibit+non_consider_time, FRAP=evaluating_frap, Inhibitor=evaluating_inhibitor)
                #soln = soln.reshape((1, (len(time_vec_fixed)*N_rib)))
                
                collisions = np.append(collisions,Ncol)
                watched_ribs.append(int(len(collisions)))
                validind = np.where(np.sum(soln,axis=1)!=0)[0]
                all_col_points.append(np.array(col_points))
                if np.max(validind) != N_rib-1:
                    validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)
                
                so = soln[(validind,)]
               
                solutionssave.append(so)

                solutions.append(soln)
            
                result = soln.reshape((1, (len(time_vec_fixed)*N_rib)))
                all_results[i, :] = result
            
            sttime = time.time() - st


                #rb = sparse.lil_matrix((len(time_vec_fixed),genelength),dtype=int)
                #for j in range(soln.shape[1]):

                    #if len(np.where(soln[:,j]!=0)[0]) !=0:
                    #print(np.where(soln[:,j]!=0)[0])


                    #rb[j,np.where(soln[:,j]!=0)[0]] = 1


                        #for value in soln[:,j][np.where(soln[:,j]!=0)[0]].astype(int):

                            #rb[j, value-1] = 1

                #rib_vec.append(rb)

        



        no_ribosomes = np.zeros((n_traj, (genelength+1)))
        
        startindex = np.where(truetime >= non_consider_time)[0][0]
        
        #all_results = all_results[:,startindex*N_rib:]

        for i in range(len(solutions)):
            for j in range(len(solutions[0][0][startindex:])):
                rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]
            
                no_ribosomes[i, rib_pos.astype(int)] += 1
        no_ribosomes = no_ribosomes[:, 1:]

        ribosome_means = np.mean(no_ribosomes, axis=0)
        ribosome_density = ribosome_means/npoints

        no_ribosomes_per_mrna = np.mean(no_ribosomes)
        
 

        if probePosition.shape[0] <=1:
            I = np.zeros((n_traj, len(time_vec_fixed[startindex:])))
         
            
        else:
            I = np.zeros((int(probePosition.shape[0]),n_traj, len(time_vec_fixed[startindex:])))
         

        #I = np.zeros((1,tstep+1))
        
        if evaluating_frap == False:
            if probePosition.shape[0] <=1:
                for i in range(n_traj):
        
                    traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
        
                    I[i, :] = np.sum(pv[traj], axis=1)[startindex:].T
            else:
                for j in range(probePosition.shape[0]):
                    for i in range(n_traj):
            
                        traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
            
                        I[j,i, :] = np.sum(pv[j][traj], axis=1)[startindex:].T                
    
    
            intensity_vec = I
        
        else:
            fraptime = time_inhibit
         
            
            inds = np.where(truetime > fraptime)

            inds2 = np.where(truetime  < fraptime+20)
            inds = np.intersect1d(inds,inds2)
            endfrap = inds[-1]-1
            

            
            for i in range(n_traj):
    
                traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
                
                nribs = np.sum(solutionssave[i][:,endfrap]!=0)
             
                #ribloc = solutionssave[i][:,endfrap]
                
                #adj_pv = pv[solutionssave[i][:,inds[-1]][:nribs]]
                frap_app = 20
                print(fraptime+start_time)
                print(fraptime+start_time+frap_app)
                print(len(traj))
                print(genelength)
                print(truetime.shape)
                revI = self.get_negative_intensity(traj,genelength,pv,truetime,fraptime+start_time,fraptime+start_time+frap_app)
                

                I[i, :] = np.sum(pv[traj], axis=1)[startindex:].T
                             
                I[i,inds[0]:inds[0]+20] = 0
                #I[i,endfrap-startindex:] = np.sum(pv[traj],axis=1)[endfrap-startindex:].T

                I[i,inds[0]+frap_app:len(revI)+inds[0]+frap_app] = I[i,inds[0]+frap_app:len(revI)+inds[0]+frap_app] + revI
                
                print(np.min(I))
      
                
                
    
    
            intensity_vec = I




        ssa_obj = ssa()
        ssa_obj.no_ribosomes = no_ribosomes
        ssa_obj.n_traj = n_traj
        ssa_obj.k = all_k
        ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        ssa_obj.rib_density = ribosome_density
        ssa_obj.rib_means = ribosome_means
        ssa_obj.rib_vec = rib_vec
        ssa_obj.intensity_vec = intensity_vec
        ssa_obj.time_vec_fixed = time_vec_fixed
        ssa_obj.time = truetime
        ssa_obj.time_rec = truetime[startindex:]
        ssa_obj.start_time = non_consider_time
        ssa_obj.watched_ribs = watched_ribs
        try:
            ssa_obj.col_points = all_col_points
        except:
            pass


        ssa_obj.evaluating_inhibitor = evaluating_inhibitor
        ssa_obj.evaluating_frap = evaluating_frap
        ssa_obj.time_inhibit = time_inhibit
        ssa_obj.solutions = solutionssave
        ssa_obj.solvetime = sttime
        ssa_obj.collisions = collisions
        
        
        try:
            ssa_obj.ribtimes = all_ribtimes[np.where(all_ribtimes > 0)]
        except:
            pass


        #solt = solutions.T

        fragmented_trajectories = []
        fragtimes = []
        maxlen = 0
    
        fragmentspertraj= []
        for k in range(n_traj):
            ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
            changes = ind[1:] - ind[:-1]
            addindexes = np.where(changes > 0)[0]
            subindexes = np.where(changes < 0)[0]
            
            sub = solutions[k][:,1:] - solutions[k][:,:-1]
            neutralindexes = np.unique(np.where(sub < 0)[1])
            neutralindexes = np.setxor1d(neutralindexes, subindexes)
            
            for index in neutralindexes:
                pre = solutions[k][:,index]
                post = solutions[k][:,index+1]
                changecount = 0
                while len(np.where(post - pre < 0)[0]) > 0:
    
                    post = np.append([genelength],post)
                    pre = np.append(pre,0)
                    
                    changecount+=1
                
                for i in range(changecount):
                    addindexes = np.sort(np.append(addindexes,index))
                    subindexes = np.sort(np.append(subindexes,index))
                    
                changes[index] = -changecount
                ind[index] += changecount
             
                
            for index in np.where(np.abs(changes)>1)[0]:
                if changes[index] < 0:
                    for i in range(np.abs(changes[index])-1):
                        subindexes = np.sort(np.append(subindexes,index))
                else:
                    for i in range(np.abs(changes[index])-1):
                        addindexes = np.sort(np.append(addindexes,index))   
                
            truefrags = len(subindexes)
     
                
        
           
            if len(subindexes) < len(addindexes):
                subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(len(truetime)-1)).astype(int))
                
            
            fragmentspertraj.append(len(subindexes))
            
            for m in range(min(len(subindexes),len(addindexes))):
                traj = solutions[k][:, addindexes[m]:subindexes[m]+1]
                traj_ind = changes[addindexes[m]:subindexes[m]+1]
                
                startind = ind[addindexes[m]]
                minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
                fragment = np.array([])
            
                    
                
                iterind = startind
                
                if subindexes[m]-addindexes[m] > 0:
                    if len(minusloc) > 1:
                        if m <= truefrags:
                            for n in range(len(minusloc)-1):
                                iterind = iterind + min(0,traj_ind[minusloc[n]])
                                fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                                
                                
                                
                  
                
                      
                            
                            fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                            
                        else:
                            for n in range(len(minusloc)-1):

                                iterind = iterind + min(0,traj_ind[minusloc[n]])
                                
                                fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                  
                                
                            fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
          
                        
                    
                    else:

                        fragment = solutions[k][startind][addindexes[m]:subindexes[m]+1].flatten()
                   
                
                    
                    fragtimes.append(addindexes[m]+1)
                       
                    
                    fragmented_trajectories.append(fragment)
                    #if m <= truefrags:
                        #kes.append(genelength/truetime[len(fragment)])
            
                    if len(fragment) > maxlen:
                        maxlen = len(fragment)
                    
    
            fragarray = np.zeros((len(fragmented_trajectories), maxlen))
            for i in range(len(fragmented_trajectories)):
                fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
            
        ssa_obj.fragments = fragarray
        ssa_obj.fragtimes = fragtimes
        ssa_obj.frag_per_traj = fragmentspertraj
        ssa_obj.full_frags = truefrags
        ssa_obj.all_results = all_results
        
        if probePosition.shape[0] > 1:
            for i in range(probePosition.shape[0]):
                if i > 0:
                    autocorr_vec2, mean_autocorr2, error_autocorr2, dwelltime2, ke_sim2  = self.get_autocorr(intensity_vec[i], truetime, 0, genelength)
                    autocorr_vec = np.vstack((autocorr_vec,autocorr_vec2))
                    mean_autocorr = np.vstack((mean_autocorr,mean_autocorr2))
                    error_autocorr = np.vstack((error_autocorr,error_autocorr2))
                    dwelltime.append(dwelltime2)
                    ke_sim.append(ke_sim2)
                else:
                    autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.get_autocorr(intensity_vec[i], truetime, 0, genelength)
                    dwelltime = [dwelltime]
                    ke_sim = [ke_sim]
            
        else:
            autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.get_autocorr(intensity_vec, truetime, 0, genelength)
        
        
        ssa_obj.autocorr_vec = autocorr_vec
        ssa_obj.mean_autocorr = mean_autocorr
        ssa_obj.error_autocorr = error_autocorr
        ssa_obj.dwelltime = dwelltime
        ssa_obj.ke_sim = ke_sim
        ssa_obj.ke_true = float(genelength)/np.mean(ssa_obj.ribtimes)
        ssa_obj.probe = probePosition

        return ssa_obj


    def get_negative_intensity(self,solution,gene_length,pv,tvec,ti,stop_frap):
        
        startindex = np.where(tvec >= ti)[0][0]
        stop_frap = np.where(tvec >= stop_frap)[0][0]
      
            
        solution = solution.T
        fragmented_trajectories = []
        fragtimes = []
        endfragtimes = []
        maxlen = 0
        
        fragmentspertraj= []
  
        ind = np.array([next(j for j in range(0,solution.shape[0]) if int(solution[j, i]) == 0 or int(solution[j, i]) == -1) for i in range(0, solution.shape[1])])
        changes = ind[1:] - ind[:-1]
        addindexes = np.where(changes > 0)[0]
        subindexes = np.where(changes < 0)[0]
        
        sub = solution[:,1:] - solution[:,:-1]
        neutralindexes = np.unique(np.where(sub < 0)[1])
        neutralindexes = np.setxor1d(neutralindexes, subindexes)
        
        for index in neutralindexes:
            pre = solution[:,index]
            post = solution[:,index+1]
            changecount = 0
            while len(np.where(post - pre < 0)[0]) > 0:
    
                post = np.append([gene_length],post)
                pre = np.append(pre,0)
                
                changecount+=1
            
            for i in range(changecount):
                addindexes = np.sort(np.append(addindexes,index))
                subindexes = np.sort(np.append(subindexes,index))
                
            changes[index] = -changecount
            ind[index] += changecount
         
            
        for index in np.where(np.abs(changes)>1)[0]:
            if changes[index] < 0:
                for i in range(np.abs(changes[index])-1):
                    subindexes = np.sort(np.append(subindexes,index))
            else:
                for i in range(np.abs(changes[index])-1):
                    addindexes = np.sort(np.append(addindexes,index))   
            
        truefrags = len(subindexes)
     
            
    
       
        if len(subindexes) < len(addindexes):
            subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(len(tvec)-1)).astype(int))
            
        
        fragmentspertraj.append(len(subindexes))
        
        for m in range(min(len(subindexes),len(addindexes))):
            traj = solution[:, addindexes[m]:subindexes[m]+1]
            traj_ind = changes[addindexes[m]:subindexes[m]+1]
            
            startind = ind[addindexes[m]]
            minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
            fragment = np.array([])
        
                
            
            iterind = startind
            
            if subindexes[m]-addindexes[m] > 0:
                if len(minusloc) > 1:
                    if m <= truefrags:
                        for n in range(len(minusloc)-1):
                            iterind = iterind + min(0,traj_ind[minusloc[n]])
                            fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                            
                            
                            
              
            
                  
                        
                        fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                        
                    else:
                        for n in range(len(minusloc)-1):
    
                            iterind = iterind + min(0,traj_ind[minusloc[n]])
                            
                            fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
              
                            
                        fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
      
                    
                
                else:
    
                    fragment = solution[startind][addindexes[m]:subindexes[m]+1].flatten()
               
            
                
                fragtimes.append(addindexes[m]+1)
                if addindexes[m]+1  + len(fragment) > len(tvec):
                    endfragtimes.append(len(tvec))
                else:
                    endfragtimes.append(addindexes[m]+1  + len(fragment))
                   
                
                fragmented_trajectories.append(fragment)
                #if m <= truefrags:
                    #kes.append(genelength/truetime[len(fragment)])
        
                if len(fragment) > maxlen:
                    maxlen = len(fragment)
                
    
        fragarray = np.zeros((len(fragmented_trajectories), maxlen))
        for i in range(len(fragmented_trajectories)):
            fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
            
            
        affected_frags = []
        fragindexes = []
        
        
        for i in range(len(fragtimes)):
        
           if  np.sum([fragtimes[i]> np.array([startindex, stop_frap]), endfragtimes[i] > np.array([startindex, stop_frap])]) in [1,2,3]:
               affected_frags.append(i)
               fragindexes.append([fragtimes[i],endfragtimes[i]])
            
        
        #affected_frags = np.intersect1d(np.where(np.array(fragtimes) >=  startindex), np.where(np.array(fragtimes)<= stop_frap))
        if len(fragindexes)> 0:
            findexes = np.array(fragindexes)
  
            frange = findexes[:,1]-stop_frap
            afterfrapribs = findexes[np.where(frange > 0 )]
            
            
            relevantfrags = np.array(affected_frags)[np.where(frange > 0 )]
            if len(relevantfrags) > 0:
                cooked_ribs = 0#(len(affected_frags) - len(relevantfrags))*max(pv)
      
                stopfrapindex = stop_frap - afterfrapribs[:,0]
                
                rfrags = fragarray[relevantfrags]
                np.diag(rfrags[:,stopfrapindex])
                laglen = afterfrapribs[:,1] - stop_frap
                posistions_at_end_of_FRAP = np.diag(rfrags[:,stopfrapindex])
             
                offset = pv[posistions_at_end_of_FRAP.astype(int)]
             
                trailing_intensity = np.zeros((max(laglen)))
                
                for i in range(len(laglen)):
                    trailing_intensity[:laglen[i]] -= offset[i] 
                    
                trailing_intensity= trailing_intensity-cooked_ribs
            else:
                trailing_intensity = np.array([0])
        else:
            trailing_intensity = np.array([0])
        
        return trailing_intensity


    def ssa_solver_append(self, ssa_obj, n=100):

        nRepetitions = ssa_obj.n_traj
        all_k = ssa_obj.k
        no_ribosomes_per_mrna = ssa_obj.no_rib_per_mrna
        ribosome_density = ssa_obj.rib_density
        ribosome_means = ssa_obj.rib_means
        rib_vec = ssa_obj.rib_vec
        intensity_vec = ssa_obj.intensity_vec
        time_vec_fixed = ssa_obj.time_vec_fixed
        non_consider_time = ssa_obj.start_time

        evaluating_inhibitor = ssa_obj.evaluating_inhibitor
        evaluating_frap = ssa_obj.evaluating_frap
        time_inhibit = ssa_obj.time_inhibit



        try:
            probePosition = []
            for key in self.POI.tag_epitopes.keys():
                probePosition = probePosition + self.POI.tag_epitopes[key]
            probePosition = np.unique(probePosition).tolist()
        except:
            print('No POI found')
                #nt_seq = self.tag_full['T_flag'] + nt_seq


        nt_seq = self.POI.nt_seq
        genelength = int(len(nt_seq)/3)



        pv = np.zeros((1, genelength)).astype(int).flatten()

        for i in range(len(probePosition)):
            pv[probePosition[i]:] = i





        npoints = len(time_vec_fixed)
        tstep = npoints-non_consider_time
        for i in range(n):

            soln = self.SSA(all_k, time_vec_fixed, inhibit_time=time_inhibit+non_consider_time, FRAP=evaluating_frap, Inhibitor=evaluating_inhibitor)

            rb = sparse.lil_matrix((len(time_vec_fixed), genelength), dtype=int)
            for j in range(soln.shape[1]):

                #if len(np.where(soln[:,j]!=0)[0]) !=0:
                #print(np.where(soln[:,j]!=0)[0])


                #rb[j,np.where(soln[:,j]!=0)[0]] = 1


                    for value in soln[:, j][np.where(soln[:, j] != 0 )[0]].astype(int):

                        rb[j, value-1] = 1

            rib_vec.append(rb)


        no_ribosomes = np.zeros((len(rib_vec), genelength))



        for i in range(len(rib_vec)):
            no_ribosomes[i] = np.sum(rib_vec[i].todense()[non_consider_time:], axis=0).flatten()

        ribosome_means = np.mean(no_ribosomes, axis=0)
        ribosome_density = ribosome_means/npoints

        no_ribosomes_per_mrna = np.mean(no_ribosomes)

        intensity_vec = np.zeros((len(rib_vec), tstep+1))

        I = np.zeros((1, tstep+1))
        for i in range(len(rib_vec)):
            for j in range(tstep):
                temp_output = rib_vec[i][non_consider_time + j, :].todense()

                I[0, j] = np.sum(pv * temp_output.flatten().T)
            intensity_vec[i] = I



        ssa_obj = ssa()

        ssa_obj.n_traj = nRepetitions + n
        ssa_obj.k = all_k
        ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        ssa_obj.rib_density = ribosome_density
        ssa_obj.rib_means = ribosome_means
        ssa_obj.rib_vec = rib_vec
        ssa_obj.intensity_vec = intensity_vec
        ssa_obj.time_vec_fixed = time_vec_fixed
        ssa_obj.start_time = non_consider_time
        ssa_obj.probe = probePosition
        ssa_obj.evaluating_inhibitor = evaluating_inhibitor
        ssa_obj.evaluating_frap = evaluating_frap
        ssa_obj.time_inhibit = time_inhibit



        if evaluating_inhibitor == False:
            autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.get_autocorr(intensity_vec, time_vec_fixed, 0, genelength)
            ssa_obj.autocorr_vec = autocorr_vec
            ssa_obj.mean_autocorr = mean_autocorr
            ssa_obj.error_autocorr = error_autocorr
            ssa_obj.dwelltime = dwelltime
            ssa_obj.ke_sim = ke_sim

        return ssa_obj



    def multitau_acc(self, ivec, n, sampling_rate, sample_rate_seconds):
        '''
        Multi-tau acc
        '''
        sigmas = 3
        acc = np.array([[]])
        for i in range(0, n):
            tempdata = ivec[i, :].flatten()
            tempdata[np.where(tempdata > tmean(tempdata, 10)) + sigmas*np.std(tempdata)] = 0
            tempdata[np.where(tempdata < tmean(tempdata, 10)) - sigmas*np.std(tempdata)] = 0

            if np.isnan(tempdata[0]):
                tempdata = tempdata[1:]
            if np.isnan(tempdata[-1]):
                tempdata = tempdata[:-1]

            outliers = np.where(tempdata == 0)[0]
            if outliers[-1] == len(tempdata)-1:
                outliers = outliers[:-1]
            if outliers[0] == 0:
                outliers = outliers[1:]

            tempdata[outliers] = 1/2*(tempdata[outliers-1] + tempdata[outliers+1])
            tempdata = tempdata-np.mean(tempdata)

            preacc = self.get_acc2(tempdata)
            if i == 0:
                acc = preacc
            else:
                acc = np.hstack((acc, preacc))
        for i in range(0, n):
            data = acc[i]
            data[0:sample_rate_seconds] = []

            binnedData_1 = data







    def geomean(self, iterable):
        '''geometric mean used for codon sensitivity calculations
        '''
        a = np.array(iterable)
        return a.prod()**(1.0/len(a))


    def SSA(self, k, t_array, inhibit_time=0, FRAP=False, Inhibitor=False):
        '''
        mRNA Translation simulation python implementation

        given a propensity vector k, time array to record, and inhibitory conditions, run a single trajectory of translation simulation

        The simulation is described here: [PUT LINK HERE TO PAPER]

        *args*

            **k**, propensity vector of size gene length + 2, [initiation rate,  Codon dependent rates,  completion rate / unbinding rate]
            for reference the codon dependent rates are refering to the time rate of a ribosome to move on to the next codon

            **t_array**, time points to record the ribosome posistions at

        *keyword args*

            **inhibit_time**, the time to start inhibition assays if FRAP or Inhibitor (harringtonine) == True

            **FRAP**, True or false to apply Fluorescence Recovery After Photobleaching (FRAP) https://en.wikipedia.org/wiki/Fluorescence_recovery_after_photobleaching

            **Inhibitor**, True or false to apply harringtonine at inhibit_time. Harringtonine acts as a protien translation initiation inhibitor

        '''

        #SSA params and propensities
        R = 10 #exclusion volume (ribosome footprint), ribosomes cant be less than 10 codons apart because of their physical size
        kelong = np.array([k[1:-1]]).T  #rates for ribosomes moving to the next codon, based on tRNA concentrations

        N = len(kelong)  #Number of codons in the mRNA
        kbind = k[0]   #rate for a ribosome to bind and start translation
        kcompl = k[-1]     #rate for a ribosome at the end of the mRNA to unbind
        X = np.array([0, 0], dtype=int)   #the updating ribosome posistion vector that is changed in the simulation


        Ncol = np.zeros((1,0))
        
        #example X arrays and how its formatted:
        # X = [423 30 10 0 ]  read from left to right theres a ribosome in position 423 30 and 10, with a 0 kept as a buffer for simulation

        t = t_array[0]  #time point
        Nt = len(t_array)  #number of time points to record over
        tf = t_array[-1]  #final time point
        N_rib = 200  #Maximum number of ribosomes on a single mRNA (hard limit for the simulation not a physical constant)
        col = np.zeros((1,N_rib))
        X_array = np.zeros((N_rib, Nt))  #recording array that records the ribosome posistions over time array points
        NR = 0  #number of ribosomes bound
        it = 1  #number of iterations
        Sn_p = np.eye(max(NR+1, 2), dtype=int) #stoichiometry for the SSA
        wn_p = np.zeros((X.shape[0], 1)) # propensities for the SSA
        
        T = np.array([0, 0], dtype=float)
        ribtimes = np.array([[0,0]],dtype=float)
        col_points = []
        #wn_p = np.zeros((1,X.shape[0])).flatten()
        wshape = len(wn_p)
        Inhibit_condition = 1  #set up inhibitor flags
        while t < tf:


            if Inhibitor == True:
                if t >= inhibit_time:

                    Inhibit_condition = 0
                else:

                    Inhibit_condition = 1
            else:
                Inhibit_condition = 1
            if FRAP == True :   #if the Photobleaching is happening, "remove" ribosomes
                if t >= inhibit_time and t < inhibit_time + 20:
                    #X = np.array([0, 0])
                    a=1
                    #T = np.array([0,0])



            oldNR = NR

            #other options for NR calc
            #NR = len(np.where(X>0)[0])
            #NR = len(np.where(X!=0)[0])
            #NR = len(np.argwhere(X))
            #NR = np.nonzero(X)[0].shape[0]
            #NR = max(0,len(X)-1)
            #NR = np.sum(X!=0)
            #NR = np.where(X!=0)[0][-1]+1
            #NR = np.flatnonzero(X).shape[0]

            NR = len(np.flatnonzero(X)) #each iteration get the number of ribosomes on the mRNA



            if X.shape[0] < NR+1:  #if the last reaction added a ribosome put a 0 on the end of X vec

                X = np.append(X, [0])
                T = np.append(T, [0])
                T[-2] = t


            X[-1] = 0
            T[-1] = 0

            X = X[0:max(NR, 1)+1]  #clear any additional 0's on the end
            T = T[0:max(NR, 1)+1]

            if oldNR != NR:     #if the number of ribosomes has changed reallocate the sizes of stoich and propensities
                Sn_p = np.eye(max(NR+1, 2), dtype=int)
                wn_p = np.zeros((X.shape[0], 1))

                wshape = len(wn_p)
                

            Sn = Sn_p
            wn = wn_p


            #get indices of where X vecs are > 0 ie where the ribosome values are
            inds = X > 0


            wn[inds] = kelong[X[inds]-1]  #update propensities



            if X[0] == N:  #if the ribosome in the 0 position is at the end of the mRNA set propensities to the reaction for completion

                Sn[:, 0] = (np.append(X[1:], np.array([0]))-X[0:])
                

                wn[0] = kcompl


            #if there are no ribosomes or when there is enough room for a new ribosome to bind add the propensity for binding
            if NR == 0:

                wn[NR] = kbind*Inhibit_condition
                
                
                

            if NR > 0 and X[NR-1] > R:
                wn[NR] = kbind*Inhibit_condition

            REST = np.less(X[1:]+10, X[0:-1])  #apply the footprint condition ie set any propensities where it violates the > 10 codons apart rule to 0


            wn[1:] = (wn[1:].T*REST).T  #apply that logical^ to propensities

            w0 = sum(wn.flat)  #get the sum of propensities
            randnum = np.random.random_sample(2) #update time to point of next reaction (exponential waiting time distb)
            t = (t-np.log(randnum[0])/w0)

            while it < Nt and t > t_array[it]:  #record state if past timepoint
                X_array[0:len(X), it] = X
                it += 1

            if t < tf:  #if still running simulation pick which reaction happened via random number and propensity sum
                r2 = w0*randnum[1]
                tmp = 0

                for i in range(wshape):
                    tmp = tmp + wn[i]
                    if tmp >= r2:
                        event = i
                        break

            X = (X + Sn[:, event].T)  #update X vector for new ribosome state
            if np.sum(Sn[:,event]) < 0 :
                
                ribtimes = np.vstack((ribtimes,[T[0],t]))
                T[:-1] = T[1:]
                Ncol = np.append(Ncol,col[0][0] )
                col = np.atleast_2d(np.append(col[:,1:],[0]))
                
            else:
                if X[event-1] == X[event] + R:
                    col[0][event] +=1
                    col_points.append( (X[event],t) )
                    
                
            
        return X_array,ribtimes[1:,:],Ncol,col_points  #return the completed simulation


    def get_acc2(self, data, trunc=False):
        '''
        Get autocorrelation function

        *NOT* multi-tau
        '''
        N = len(data)
        fvi = np.fft.fft(data, n=2*N)
        acf = fvi*np.conjugate(fvi)
        acf = np.fft.ifft(acf)
        acf = np.real(acf[:N])/float(N)
        if trunc:
            acf[acf < 0]=0
            for i in range(1, len(acf)):
                if acf[i] > acf[i-1]:
                    acf[i] = acf[i-1]
        return acf


    def elongation_animation(self, ti=0, tf=1000, tstep=1000, cell_radius=50, imagesize=5, dpi=90, filename='simulated_cell', ssa_obj=None, fcolor='#00FF00' ,rnacolor='#FF0000', xkcd=False):
        '''
        function that creates a mrna translation animation
        '''

        custom_cmap = ['#69dd42', '#e5361b', '#db11c7']
        def rpts(x, y, angle):
            nx = np.cos(angle)*x - np.sin(angle)*y
            ny = np.sin(angle)*x + np.cos(angle)*y
            return nx, ny


        def update_line(num, xpos, ypos, line):  #function for the FuncAnimation
            if num != 0:
                ax.get_lines()[-1].remove()

                for child in ax.get_children():  #remove the previous patch collection (green spots)

                    if isinstance(child, PatchCollection):
                        child.remove()


            patches = []
            gp = []
            ep = []
            radi = np.ones(xpos[:, inds[num]].shape)*4 #create a max radius of 3 for intensity vecs
            ypos = np.ones(xpos[:, inds[num]].shape)*(ytop+3)
            x = xpos[:, inds[num]]
            x[np.where(x == 0)] = x[np.where(x == 0)] - 300

            for x1, y1, r in zip(xpos[:, inds[num]], ypos, radi):   #make circle objects of radius based on ivec
                circle = mpatches.Circle((x1, y1), r, facecolor='#FF0000', edgecolor='k')
                patches.append(circle)

            pcolor = custom_cmap[0]

            for i in range(len(x.flatten())):

                if x[i] > 0:
                    xpts = np.linspace(0, int(x[i])-1, int(x[i]))
                    ypts = 5*np.sin(1/10*np.linspace(0, int(x[i])-1, int(x[i])))
                    xpts, ypts = rpts(ypts, xpts, 1)
                    ypts = ypts+ytop+3
                    xpts = xpts+x[i]
                    radi = np.ones(xpts.shape)*1
                    k = 0
                    ypts = np.fliplr(np.atleast_2d(ypts))
                    ypts = ypts.flatten()
                    xpts = np.fliplr(np.atleast_2d(xpts))
                    xpts = xpts.flatten()

                    for x2, y2, r2 in zip(xpts, ypts, radi):
                        probloc = False
                        j = 0
                        for key in epitopes.keys():
                            if k in epitopes[key]:
                                probloc = True
                                pcolor = custom_cmap[j]
                                j += 1
                        rx = np.random.rand()*2
                        ry = np.random.rand()*2
                        if probloc == False:

                            circle = mpatches.Circle((x2+rx, y2+ry), r2, facecolor='#0000FF', edgecolor='#FFFFFF', lw=2, ls='solid')
                            gp.append(circle)
                        else:
                            circle = mpatches.Circle((x2+rx, y2+ry), r2*3, facecolor='#00FF00', edgecolor='#000000', lw=2, ls='solid')
                            ep.append(circle)

                        k += 1


                #fig.gca().add_artist(circle)
            '''
            xs = np.flip(np.sort(xpos[:,inds[num]][0].flatten()),axis=0)
            for i in range(max_ribs):
                line.set_data(xpos[:,inds[num]],ypos[inds[num]])
                line.set_linewidth(0)
                line.set_marker('o')
                line.set_markersize(3)
            '''
            p = PatchCollection(patches, facecolors=('#FF0000',), zorder=5)  #create a patch collection to add to axis

            m = PatchCollection(gp, facecolors=('#0000FF',), lw=2, zorder=3)  #create a patch collection to add to axis
            e = PatchCollection(ep, facecolors=(pcolor,), zorder=4)



            n = num

            ax.plot(np.linspace(0, tag_length, len(ssa_obj.time_vec_fixed[ssa_obj.start_time:]))[:n], 3*ssa_obj.intensity_vec.flatten()[:n]+total_length, color=pcolor)

            fldot = mpatches.Ellipse((total_length-30, total_length+40), width=ssa_obj.intensity_vec.flatten()[n], height=ssa_obj.intensity_vec.flatten()[n]*1.0, color=pcolor)
            f = [fldot]
            fe = PatchCollection(f, facecolors=(pcolor,), zorder=4)
            ax.add_collection(p)  #adds the circles to axis
            ax.add_collection(m)  #adds the circles to axis
            ax.add_collection(e)
            ax.add_collection(fe)
            plt.xlabel(str(inds[num]))  #update time label
            return line,

        if ssa_obj == None:
            ssa_obj = self.ssa_solver(n_traj=1, tf=tf, tstep=tstep)
        if xkcd == True:
            plt.xkcd()



        fig1 = plt.figure(figsize=(imagesize+5, imagesize), dpi=dpi)  #make figure
        fig1.tight_layout()

        ax = fig1.add_subplot('111')
        ax.set_aspect(1)

        tag_length = self.POI.tag_length
        total_length = self.POI.total_length

        epitopes = self.POI.tag_epitopes

        tag_length = total_length - self.POI.gene_length


        ax.cla()
        ybot = 90
        ytop = 110
        ax.plot([0, total_length], [ybot, ybot], color='white', zorder=3)
        ax.plot([0, total_length], [ytop, ytop], color='white', zorder=3)
        ax.plot([0, 0], [ybot, ytop], color='white', zorder=3)
        ax.plot([total_length, total_length], [ybot, ytop], color='white', zorder=3)
        ax.axis([-10, total_length+10, 80, total_length+np.max(ssa_obj.intensity_vec)*3+20])
        ax.plot([tag_length, tag_length], [ybot, ytop], color='white', linewidth=1, zorder=3)
        k = 0
        for key in epitopes.keys():
            for i in range(len(epitopes[key])):
                ax.plot([epitopes[key][i], epitopes[key][i]], [ybot, ytop], color=custom_cmap[k], linewidth=2, zorder=3)
            rect = mpatches.Rectangle(xy=(tag_length, ybot), width=total_length-tag_length, height=ytop-ybot, color='#0000FF')
            #ax.fill_between([tag_length,tag_length,total_length,total_length],[ybot,ytop,ytop,ybot],color='#00FF00')
            ax.add_patch(rect)
            k += 1

        ticks = np.linspace(0, total_length, 10).astype(int)
        ax.set_xticks(ticks)
        ax.set_xlabel('Codon Position')
        ax.get_yaxis().set_visible(False)
        ax.set_facecolor('k')
        filename = 'elong.gif'

        Writer = animation.writers['pillow']

        print('making movie...')
        max_ribs = np.max(np.nonzero(ssa_obj.solutions[0])[0])

        l, = plt.plot([], [], 'r-')
        t = ssa_obj.time_vec_fixed[ssa_obj.start_time:]
        inds = np.linspace(0, len(t)-1, len(t)).astype(int)
        xpos = np.zeros((max_ribs, len(ssa_obj.time_vec_fixed[ssa_obj.start_time:])))

        ypos = np.ones((1, len(ssa_obj.time_vec_fixed[ssa_obj.start_time:]))).flatten()

        xpos[:, :] = ssa_obj.solutions[0][:max_ribs, ssa_obj.start_time:len(ssa_obj.time_vec_fixed)]

        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
        line_ani = animation.FuncAnimation(fig1, update_line, len(ssa_obj.time_vec_fixed[ssa_obj.start_time:]), fargs=(xpos, ypos, l),
                                           interval=50, blit=True)

        line_ani.save((filename), writer=writer)  #save the animation




    def simulate_cell(self, diffusion_constant, kon, koff, kRNA, kdecay, ti=0, tf=1000, tstep=1000, cell_radius=50, imagesize=5, dpi=90, filename='simulated_cell', ssa_obj=None, fcolor='#00FF00', rnacolor='#FF0000'):
        '''
        [DNA] ==kRNA==> [RNA] <==koff== [RNA*] ==translation simulation==> [Protein]===> null
                        // ||             /\
                        || `'=====kon====='`
                        ||
                        \/
                        null

        '''
        print('simulating RNA creation....')
        t = np.linspace(ti, tf, tstep)


        dna_s = np.array([[ 0,  0],
                          [ 1, -1]])

        dna_w1 = np.array([[kRNA, 0],
                           [0, 0]],dtype=float)


        dna_w0 = np.array([[0], [0]])


        dna_si = GenericSSA(type='linear' )
        dna_si.W1 = dna_w1
        dna_si.W0 = dna_w0
        dna_si.S = dna_s

        dna_si.ti = t[0]
        dna_si.tf = t[-1]
        dna_si.n = 1
        xi = np.zeros((2, 1))
        xi[0] = 1
        dna_si.xi = xi
        dna_si.ptimes = len(t)

        dna_si.time_variant = False
        dna_si._solve(1)
        rna_creation_data = dna_si.data




        stoich = np.array([[  0,    0,  1],
                           [  -1,  1, -1],
                           [  1, -1, 0]])

        propensity = np.array([[0, kon, 0],
                              [0, 0,koff],
                              [0,kdecay, 0]], dtype=float)

        w0 = np.array([[0],[0],[0]])

        solver_instance = GenericSSA(type='linear' )
        solver_instance.W1 = propensity
        solver_instance.W0 = w0
        solver_instance.S = stoich

        solver_instance.ti = t[0]
        solver_instance.tf = t[-1]
        solver_instance.n = 1
        xi = np.zeros((3,1))
        xi[1] = 1
        solver_instance.xi = xi
        solver_instance.ptimes = len(t)

        solver_instance.time_variant = False




        print('simulating RNA activation....')



        R = cell_radius
        squarelen = float(R/np.sqrt(2))

        n_RNA_t = np.zeros((len(t),int(np.max(rna_creation_data[1]))))

        nRNA = 0
        nparticles = (int(np.max(rna_creation_data[1])))
        for i in range(len(t)):

            while nRNA != rna_creation_data[1][i]:
                data = solver_instance._solve(1)


                rnaonoff = data[2] + 1 - data[0]



                n_RNA_t[i:, nRNA] = rnaonoff[:-i].flatten()
                nRNA += 1


        rna_particles = n_RNA_t.T
        rna_exist = np.where(rna_particles >0,1,0)
        rnaex = data



        print('simulating RNA motion....')
        rna_locations = np.empty((nparticles, len(t), 2))

        dt = t[-1]/len(t)

        delta = diffusion_constant



        def linecirc(m, b, xc, yc, r):

            if np.isinf(m) == False:
                a = 1+m**2
                e = 2*(m*(b-yc)-xc)
                c = yc**2+xc**2 + b**2-2*yc*b-r**2
                x = np.roots([a, e, c])

                if np.isreal(x).all() == False:
                    x = [np.nan, np.nan]
                    y = [np.nan, np.nan]
                else:
                    y = [b + m*x[0], b+m*x[1]]

            elif abs(xc-b) > r:
                x = [np.nan, np.nan]
            else:
                x = [b, b]
                step = np.sqrt(r**2-(b-xc)**2)
                y = [yc + step, yc-step]

            return [x[0], y[0]], [x[1], y[1]]


        def dist(x1, y1, x2, y2):
            return np.sqrt((x1-x2)**2+(y1-y2)**2)


        for i in range(nparticles):
            x = np.empty((2,len(t) - np.where(rna_exist[i] != 0 )[0][0]  ))
            centers = np.zeros(x.shape)
            x[:,0] = np.random.random()*squarelen
            x0 = [  ((R+squarelen/4) - (R-squarelen/4))*np.random.random() + (R-squarelen/4),((R+squarelen/4) - (R-squarelen/4))*np.random.random() + (R-squarelen/4) ]


            x0 = x0 - np.array([R, R])
            x[:,0] =x0
            r = norm.rvs(size=np.array(x0).shape + (len(t) - np.where(rna_exist[i] !=0 )[0][0],), scale=delta*np.sqrt(dt))



            out = np.empty(r.shape)

            np.cumsum(r, axis=-1, out=out)
            out += np.expand_dims(np.array(x0), axis=-1)

            #out = np.array([[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36],
                            #[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36]])

            centers = np.zeros(out.shape)
            dists = np.zeros((x.shape[1], 1)).flatten()

            incirc = np.hypot(out.T[:, 0]-centers.T[:, 0], out.T[:, 1]-centers.T[:, 1])
            dists[np.where(out[0] != 0)] = incirc[np.where(out[0] != 0)]

            while len(np.where(dists>R)[0]) != 0:   #trajectory left the cell
                out = out.T
                left_cell = np.where(dists > R)[0][0]



                pts = [[out[left_cell][0], out[left_cell][1]], [out[left_cell-1][0], out[left_cell-1][1]]]

                p = np.polyfit([out[left_cell][0], out[left_cell-1][0]], [out[left_cell][1], out[left_cell-1][1]], 1)
                m = p[0]
                b = p[1]

                intercepts = linecirc(m, b, 0, 0, R)
                if dist(*tuple(intercepts[0])+tuple(pts[0])) > dist(*tuple(intercepts[1])+tuple(pts[0])):
                    inter = np.array(intercepts[1])
                else:
                    inter = np.array(intercepts[0])

                a = out[left_cell] - inter


                out[left_cell-1:] = out[left_cell-1:] - 2*(np.dot(inter, a)/np.linalg.norm(inter)**2)*inter




                dists = np.zeros((x.shape[1], 1)).flatten()
                out = out.T
                incirc = np.hypot(out.T[:, 0]-centers.T[:, 0], out.T[:, 1]-centers.T[:, 1])
                dists[np.where(out[0] != 0)] = incirc[np.where(out[0] != 0)]



            data = ((out.T).T*rna_exist[i][np.where(rna_exist[i] != 0)[0][0]:].T).T
            data[np.where(rna_exist[i] != 0)[0][-1]- np.where(rna_exist[i] != 0)[0][0]+1 :] = -R


            rna_locations[i, np.where(rna_exist[i] != 0)[0][0]:, :] =  data
            rna_locations[i, :np.where(rna_exist[i] != 0)[0][0], :] =  -R


        print(nparticles)
        rna_loc_compressed = rna_locations[np.where(np.sum(np.sum(rna_locations+R, axis=1), axis=1) > 0)]

        if ssa_obj == None:
            print('no ssa data given')
            print('simulating translation....')
            print(int(rna_loc_compressed.shape[0]))
            ssa_obj = self.ssa_solver(n_traj=int(rna_loc_compressed.shape[0]),tf=tf,tstep=tstep)

            ivec = ssa_obj.intensity_vec/np.max(ssa_obj.intensity_vec)
            ivec = ivec.T  #get the intensity vec for the "fluorescence"


        else:
            print('Translation data given')
            print('Given ' + str(ssa_obj.n_traj) + ' Needed '+str(int(rna_loc_compressed.shape[0])) )
            if ssa_obj.n_traj  < int(rna_loc_compressed.shape[0]):
                print('simulating ' + str(int(rna_loc_compressed.shape[0]) - ssa_obj.n_traj) + ' additional trajectories....')
                ssa_obj = self.ssa_solver_append(ssa_obj, n=int(rna_loc_compressed.shape[0]) - ssa_obj.n_traj)
                ivec = ssa_obj.intensity_vec/np.max(ssa_obj.intensity_vec)
                ivec = ivec.T  #get the intensity vec for the "fluorescence"

            else:
                ivec = ssa_obj.intensity_vec[0:int(rna_loc_compressed.shape[0])]/np.max(ssa_obj.intensity_vec[0:int(rna_loc_compressed.shape[0])])
                ivec = ivec[0:int(rna_loc_compressed.shape[0])].T  #get the intensity vec for the "fluorescence"







        print('making movie...')
        #simulate brownian motion
        def update_line(num, xpos,ypos, line):  #function for the FuncAnimation
            if num !=0:

                for child in ax.get_children():  #remove the previous patch collection (green spots)

                    if isinstance(child, PatchCollection):
                        child.remove()
                    if isinstance(child, mpatches.Ellipse):
                        child.remove()

            patches = []
            radi = 3*ivec[inds[num]]   #create a max radius of 3 for intensity vecs


            for x1, y1, r in zip(xpos[inds[num]],ypos[inds[num]], radi):   #make circle objects of radius based on ivec
                circle = mpatches.Circle((x1, y1), r, color=fcolor)
                patches.append(circle)
                #fig.gca().add_artist(circle)

            line.set_data(xpos[inds[num]], ypos[inds[num]])
            line.set_linewidth(0)
            line.set_marker('o')
            line.set_markersize(3)
            line.set_color(rnacolor)
            line.set
            p = PatchCollection(patches, zorder=3, facecolors=(fcolor,))  #create a patch collection to add to axis
            ax.add_collection(p)  #adds the circles to axis


            p = mpatches.Circle((0,0), radius=R, color='black')  #add the black circle
            ax.add_patch(p)


            whitep = mpatches.Ellipse((-R, -R), width=7, height=7, color='white', zorder=5)  #add the black circle
            ax.add_patch(whitep)

            plt.xlabel(str(inds[num]))  #update time label


            return line,
        xpos = rna_loc_compressed.T[0]
        ypos = rna_loc_compressed.T[1]


        filetype='.mov'

        if filetype == '.gif':

            Writer = animation.writers['pillow']
        if filetype == '.html':
            Writer = animation.writers['html']
        if filetype == '.gif':

            Writer = animation.writers['FFMpeg']

        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

        fig1 = plt.figure(figsize=(imagesize, imagesize),dpi=dpi)  #make figure
        fig1.tight_layout()

        ax= fig1.add_subplot('111')
        plt.yticks([])
        plt.xticks([])
        p = mpatches.Circle((0, 0), radius=R, color='black')  #add the black circle
        ax.add_patch(p)
        plt.gca().set_aspect('equal', adjustable='box')

        l, = plt.plot([], [], 'r-')
        plt.xlim(-R-10, R+10)
        plt.ylim(-R-10, R+10)
        plt.xlabel('0')
        plt.title('Simulated Cell')

        inds = np.linspace(0, len(t)-1, len(t)).astype(int)
        #creates the animation
        line_ani = animation.FuncAnimation(fig1, update_line, tstep, fargs=(xpos,ypos, l),
                                           interval=50, blit=True)

        line_ani.save((filename + filetype), writer=writer)  #save the animation



        #return solver_instance,n_RNA_t,rna_creation_data,data,rna_locations
        return rna_locations, rna_loc_compressed, rna_particles, rna_creation_data, rna_exist, rnaonoff, rnaex


    def get_simulated_mov(self, ssa_obj, filename, filetype):
        '''
        Create a gif or html file of the simulated circlular cell from any sms ssa object

        '''

        R = 50   #set up the random circle and points
        num = ssa_obj.n_traj
        r1 = np.zeros((1, num)).flatten()
        theta1 = np.zeros((1, num)).flatten()
        x1 = np.zeros((1, num)).flatten()
        y1 = np.zeros((1, num)).flatten()


        r2 = np.zeros((1, num)).flatten()
        theta2 = np.zeros((1, num)).flatten()
        x2 = np.zeros((1, num)).flatten()
        y2 = np.zeros((1, num)).flatten()

        for n in range(0, num):  #for  all trajectories make initial points
            r1[n] = R*np.sqrt(np.random.random(1))
            r2[n] = R*np.sqrt(np.random.random(1))
            theta1[n] = 2*np.pi*np.random.random(1)
            theta2[n] = 2*np.pi*np.random.random(1)
            x1[n] = np.cos(theta1[n])*r1[n]
            x2[n] = np.cos(theta2[n])*r2[n]
            y1[n] = np.sin(theta1[n])*r1[n]
            y2[n] = np.sin(theta1[n])*r2[n]

        movement = .7
        xpos = np.zeros((len(ssa_obj.time_vec_fixed[ssa_obj.start_time:]), num))
        ypos = np.zeros((len(ssa_obj.time_vec_fixed[ssa_obj.start_time:]), num))

        #over time perterub the simulated ribosome posistions and save to xpos , ypos
        for j in range(len(ssa_obj.time_vec_fixed[ssa_obj.start_time:])):
            if j == 0:
                xpos[0] = x1
                ypos[0] = y1
            else:
                for i in range(0,num):
                    xpos[j, i] = xpos[j-1,i]-movement + 2*movement*np.random.random(1)
                    if xpos[j, i] > 52:
                        xpos[j, i] = 51
                    if xpos[j, i] < -52:
                        xpos[j, i] = -51
                    ypos[j, i] = ypos[j-1,i]-movement + 2*movement*np.random.random(1)
                    if ypos[j, i] > 52:
                        ypos[j, i] = 51
                    if ypos[j, i] < -52:
                        ypos[j, i] = -51

        ivec = ssa_obj.intensity_vec/np.max(ssa_obj.intensity_vec)
        ivec = ivec.T  #get the intensity vec for the "fluorescence"
        k = 0
        def update_line(num, xpos, ypos, line):  #function for the FuncAnimation
            if num !=0:
                for child in ax.get_children():  #remove the previous patch collection (green spots)

                    if isinstance(child, PatchCollection):
                        child.remove()

            patches = []
            radi = 3*ivec[inds[num]]   #create a max radius of 3 for intensity vecs


            for x1, y1, r in zip(xpos[inds[num]],ypos[inds[num]], radi):   #make circle objects of radius based on ivec
                circle = mpatches.Circle((x1, y1), r, color='#00FF00')
                patches.append(circle)
                #fig.gca().add_artist(circle)

            line.set_data(xpos[inds[num]],ypos[inds[num]])
            line.set_linewidth(0)
            line.set_marker('o')
            line.set_markersize(.5)
            p = PatchCollection(patches, zorder=2, facecolors=('#00FF00',))  #create a patch collection to add to axis





            ax.add_collection(p)  #adds the circles to axis
            plt.xlabel(str(inds[num]))  #update time label





            return line,
        if filetype == '.gif':

            Writer = animation.writers['pillow']
        if filetype == '.html':
            Writer = animation.writers['html']
        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

        fig1 = plt.figure()  #make figure
        fig1.tight_layout()

        ax= fig1.add_subplot('111')
        plt.yticks([])
        plt.xticks([])
        p = mpatches.Circle((0, 0), radius=65, color='black')  #add the black circle
        ax.add_patch(p)
        plt.gca().set_aspect('equal', adjustable='box')

        l, = plt.plot([], [], 'r-')
        plt.xlim(-70, 70)
        plt.ylim(-70, 70)
        plt.xlabel('0')
        plt.title('Simulated Cell')
        inds = np.linspace(0, xpos.shape[0]-1, 120).astype(int)
        #creates the animation
        line_ani = animation.FuncAnimation(fig1, update_line, 120, fargs=(xpos,ypos, l),
                                           interval=50, blit=True)

        line_ani.save((filename + filetype), writer=writer)  #save the animation


    def analyze_seq_file(self, filename):
        '''
        General catch all to run all functions necessary before a SSA and store the first POI found from any given sequence

        *args*

            **filename** a txt or gb file to be read and analyzed

        '''


        self.open_seq_file(filename)
        self.get_orfs(self.sequence_str, min_codons=80)
        self.get_temporal_proteins()
        self.analyze_poi(self.pois[0], self.pois_seq[0])
        self.POI.k = self.get_k(self.POI.nt_seq, .03, 10)
        probe_vec,probe_loc = self.get_probvec()
        self.POI.probe_vec = probe_vec
        self.POI.probe_loc = probe_loc

    def run_default(self):

        self.get_orfs(self.sequence_str, min_codons=80)
        self.get_temporal_proteins()
        self.analyze_poi(self.pois[0], self.pois_seq[0])



    def get_gb_file(self, accession_number, savetofile=False):
        '''
        A function to poll genbank given an accession number and pull the relevant gb file

        *args*

            **accession_number**, the accession number of the sequence to find.
            http://www.nslc.wustl.edu/elgin/genomics/bio4342/1archives/2006/AccReference.pdf

        *keyword args*

            **savetofile**, true or false to save the gb file in the same directory as sms for future use



        '''



        Entrez.email = "wsraymon@rams.colostate.edu"
        Entrez.tool = 'SingleMoleculeSimulator'

        er = False

        try:
            handle =  Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=accession_number)
            gb_record = SeqIO.read(handle, "genbank") #using "gb" as an alias for "genbank"
            handle.close()
        except:
            er = True


        time.sleep(2)

        if er == True:
            print('HTTP Error: Could not find specified ascession ID')

            return


        self.gb_rec = gb_record
        self.gb_obj = gb_record

        self.sequence_str = str(gb_record.seq)
        self.sequence_name = gb_record.name

        if savetofile:
            filename = self.sequence_name
            f = open(filename, 'w')


            f.write(self.gb_rec.format('gb'))

            f.close()



    def tau_plot(self,ssa_obj,t,tau,plot_type='contour', plot_all = False):

        
        
        stime = ssa_obj.time_rec-ssa_obj.start_time
        idx_t = (np.abs(stime - t)).argmin()
        idx_tau = (np.abs(stime - tau)).argmin()
        
        diff = idx_tau - idx_t
        difftime = t-tau
        
        if plot_type == 'Average':

            fig,ax= plt.subplots()
            for i in range(len(stime)-idx_tau,0,-4):
                idx_tau = (np.abs(stime- (stime[i]+difftime ))).argmin() 
                Itau = ssa_obj.intensity_vec[:,idx_tau]
                x,y = np.mean(ssa_obj.intensity_vec[:,idx_tau]/np.sum(ssa_obj.probe)),np.mean(ssa_obj.intensity_vec[:,idx_tau+diff]/np.sum(ssa_obj.probe))

        if plot_type == 'window':
            
            minx = 10000000
            maxx = 0
            
            miny = 10000000
            maxy = 0
            
            fig,ax= plt.subplots()
            for i in range(len(stime)-idx_tau,0,-10):
                idx_tau = (np.abs(stime - (idx_t+i))).argmin()  
                Itau = ssa_obj.intensity_vec[:,idx_tau]
                x,y = np.mean(ssa_obj.intensity_vec[:,idx_tau]/np.sum(ssa_obj.probe)),np.mean(ssa_obj.intensity_vec[:,idx_tau+diff]/np.sum(ssa_obj.probe))
                minx = min(np.min(x),minx)
                miny = min(np.min(y),miny)
                maxx = max(np.max(x),maxx)
                maxy = max(np.max(y),maxy)
                
                ax.scatter(x, y,zorder=3,color= cm.viridis_r(1.*i/len(stime)))


            c_map_ax = fig.add_axes([.95, 0.1, 0.1, 0.8])
          

            c_map_ax.axes.get_xaxis().set_visible(False)

            cbar = mpl.colorbar.ColorbarBase(c_map_ax, cmap=cm.viridis_r, orientation = 'vertical')
            
            
            cbar.ax.set_yticklabels(np.linspace(idx_t,stime[-1],6).astype(int) )
            cbar.ax.set_title('t')
            
            ax.plot([min(minx,miny),max(maxx,maxy)],[min(minx,miny),max(maxx,maxy)], color='red',ls='--')
            
            ax.set_ylabel(('<I(t=' + 't + tau'+')>'))
            ax.set_xlabel(('<I(t=' +'t'+')>'))
            ax.set_title(( 'Average I(t) vs Average I(t+tau) for tau = ' + str(diff) ) )
            
            
        if plot_type == 'density':
            fig,ax= plt.subplots()
            nbins = int(np.max(ssa_obj.intensity_vec/np.sum(ssa_obj.probe)))+2
            x, y = ssa_obj.intensity_vec[:,idx_t]/np.sum(ssa_obj.probe),ssa_obj.intensity_vec[:,idx_tau]/np.sum(ssa_obj.probe)
            k = kde.gaussian_kde([x,y])
            xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(), yi.flatten()])) 
            
            R = pearsonr(x,y)[0]
            ax.set_title(('Density Plot' + ' R = ' + str(np.round(R,3))))
            ax.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.viridis)
            ax.contour(xi, yi, zi.reshape(xi.shape) )   
            ax.set_ylabel(('I(t=' + str(tau)+')'))
            ax.set_xlabel(('I(t=' + str(t)+')'))
            fig.show()    
            
            
        if plot_type == 'set_tau':
            fig,ax= plt.subplots()
            for i in range(len(stime)-diff-idx_t):
                    idx_tau = (np.abs(stime - (idx_t+i))).argmin()                
                    plt.scatter(ssa_obj.intensity_vec[:,i]/np.sum(ssa_obj.probe), ssa_obj.intensity_vec[:,i+diff]/np.sum(ssa_obj.probe),c= cm.viridis(1.*i/len(stime)),alpha=.5  )
            plt.ylabel('I(t + s)')
            plt.xlabel(('I(t)'))
            plt.title(('Set tau, all times s = ' + str(diff) ))            

            c_map_ax = fig.add_axes([.95, 0.1, 0.1, 0.8])
          

            c_map_ax.axes.get_xaxis().set_visible(False)

            cbar = mpl.colorbar.ColorbarBase(c_map_ax, cmap=cm.viridis, orientation = 'vertical')
            
            cbar.ax.set_yticklabels(np.linspace(idx_t,stime[-1],6).astype(int) )
        
        if plot_type == 'scatter':
            if not plot_all:
                
                plt.scatter(ssa_obj.intensity_vec[:,idx_t]/np.sum(ssa_obj.probe), ssa_obj.intensity_vec[:,idx_tau]/np.sum(ssa_obj.probe) )
                plt.ylabel(('I(t=' + str(tau)+')'))
                
                
            else:
               
                for i in range(idx_t,len(stime)):
                    idx_tau = (np.abs(stime - (idx_t+i))).argmin()                
                    plt.scatter(ssa_obj.intensity_vec[:,idx_t]/np.sum(ssa_obj.probe), ssa_obj.intensity_vec[:,idx_tau]/np.sum(ssa_obj.probe),c= cm.viridis(1.*i/len(stime)),alpha=.1  )
                    plt.ylabel('I(tau)')
            plt.xlabel(('I(t=' + str(t)+')'))
            
        if plot_type == 'contour':
            fig,ax= plt.subplots()
            if not plot_all:
                It = ssa_obj.intensity_vec[:,idx_t]
                Itau = ssa_obj.intensity_vec[:,idx_tau]
                
                cov = np.cov(It,Itau)
                
                eigs, v = np.linalg.eig(cov)
                eigs = np.sqrt(eigs)
                plt.ylabel(('I(t=' + str(tau)+')'))
                colors = [cm.viridis(1.0),cm.viridis(.5),cm.viridis(0.0),cm.viridis(0.0)]
      
                
                for j in xrange(3, 0,-1):
                   
                    ell_artist = Ellipse(xy=(np.mean(It), np.mean(Itau)),
                                  width=eigs[0]*j*2, height=eigs[1]*j*2,
                                  angle=np.rad2deg(np.arccos(v[0, 0])))
                    
                    ell_artist.set_linewidth(2)
                    ell_artist.set_edgecolor(colors[j-1])
                    ell_artist.set_color(colors[j-1])
                    ax.add_patch(ell_artist)
                    
                ax.autoscale()      
                ax.set_xlim(0)
                ax.set_ylim(0)
                ax.scatter(It, Itau,zorder=3,alpha=0.3,color='red',marker='.')
                fig.show()
            else:
                plt.ylabel('I(tau)')
                It = ssa_obj.intensity_vec[:,idx_t]
                for i in range(len(stime)-idx_t,0,-10):
                    idx_tau = (np.abs(stime - (idx_t+i))).argmin()  
                    Itau = ssa_obj.intensity_vec[:,idx_tau]
                   
                    cov = np.cov(It,Itau)
                    
                    eigs, v = np.linalg.eig(cov)
                    eigs = np.sqrt(eigs)
                    
                    
                    j = 3
                    ell_artist = Ellipse(xy=(np.mean(It), np.mean(Itau)),
                                  width=eigs[0]*j*2, height=eigs[1]*j*2,
                                  angle=np.rad2deg(np.arccos(v[0, 0])))
                    
                    ell_artist.set_linewidth(2)
                    ell_artist.set_edgecolor( cm.viridis_r(1.*i/len(stime)))
                    ell_artist.set_color( cm.viridis_r(1.*i/len(stime)))
                    ax.autoscale()    
                    ax.add_patch(ell_artist)
                    ax.figure.canvas.draw()
            
                
                    
                plt.xlabel(('I(t=' + str(t)+')'))
                ax.set_xlim(0)
                ax.set_ylim(0)
                c_map_ax = fig.add_axes([.95, 0.1, 0.1, 0.8])
              

                c_map_ax.axes.get_xaxis().set_visible(False)

                cbar = mpl.colorbar.ColorbarBase(c_map_ax, cmap=cm.viridis_r, orientation = 'vertical')
                
                cbar.ax.set_yticklabels(np.linspace(idx_t,stime[-1],6).astype(int) )
                
                fig.show()       
            
    

    def kymograph(self,ssa_obj,n_traj,bg_intense=True,show_intense = True,tag = 0, show_col=True,col_size = 1.5, custom_fig = None, facecolor='white', *args,**kwargs):
        '''
        Constructs a kymograph of ribosome locations
        '''
        startfrags = 0
        for i in range(n_traj):
            startfrags += ssa_obj.frag_per_traj[i]
            
        endfrags = startfrags + ssa_obj.frag_per_traj[n_traj]
        fragments = ssa_obj.fragments[startfrags:endfrags]

        

        time = ssa_obj.time#[0:len(ssa_obj.time_rec)-1]
        

        
        
        
        
        if len(ssa_obj.intensity_vec.shape) ==3:
            ivec = ssa_obj.intensity_vec[tag][n_traj]
        else:        
            ivec = ssa_obj.intensity_vec[n_traj]
        ftimes = ssa_obj.fragtimes[startfrags:startfrags+endfrags]


        nfrag = fragments.shape[0]
        maxlen= fragments.shape[1]

        
        #plt.figure(figsize=(5,10))
        if show_intense == True:
            gs = gridspec.GridSpec(1, 2, custom_fig, width_ratios=[3, 1]) 
        else:
            gs = gridspec.GridSpec(1, 1)
        
        plt.subplot(gs[0])
        lenplot = np.max(fragments)
        maxin = np.max(ivec)
        ax = plt.gca()
        ax.set_facecolor(facecolor)

        if bg_intense == True:
            for i in range(len(time)):
                plt.plot([0,lenplot],[time[i],time[i]],color = cm.summer(1.*ivec[i]/maxin),lw=1)
            
        for i in range(nfrag):
            
            
            
            if maxlen <= np.where(fragments[i] > 0 )[0][-1]:       
                timeseg = time[ftimes[i]:ftimes[i]+maxlen]
                
                plt.plot(fragments[i][0:len(timeseg)] ,timeseg[::-1] )
                
            else:
                timeseg = time[ftimes[i]:]
                
                stop = np.where(fragments[i] > 0 )[0][-1]
                timelen = len(fragments[i][0:stop]) 

                plt.plot(fragments[i][0:stop]   ,timeseg[0:timelen],**kwargs )

        plt.xlabel('Ribosome position')
        plt.ylabel('Time (sec)')
        segtime = ssa_obj.time[0:len(ssa_obj.time_rec)]
        plt.ylim(ssa_obj.time_rec[-1], ssa_obj.time_rec[0])
                
        if show_col == True:
            try:
                col = ssa_obj.col_points[n_traj]

                plt.plot(col[:,0],col[:,1],color='#00ff00',markersize=col_size,linestyle='none',marker='o')
            except:
                pass

        if show_intense == True:
            plt.subplot(gs[1])
            ax = plt.gca()
            ax.set_facecolor(facecolor)
            print((np.sum(ssa_obj.probe)))
            plt.plot(ivec.T/ np.sum(ssa_obj.probe),segtime,**kwargs)
            plt.xlabel('Intensity (ump)')
            plt.xlim(0,30)
            plt.ylim(segtime[-1], segtime[0])
            
        
            
            plt.tight_layout()
                

            


    def get_autocorr(self, intensity_vec, time_vec, totalSimulationTime, geneLength):
        '''
        returns the autocorrelations
        '''

        autocorr_vec = np.zeros((intensity_vec.shape))
        
        for i in range(intensity_vec.shape[0]):
            autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-np.mean(intensity_vec[i]))

        normalized_autocorr = autocorr_vec.T/ autocorr_vec[:,0]
        mean_autocorr = np.mean(normalized_autocorr, axis=1)
        error_autocorr = np.std(normalized_autocorr, axis=1)/np.sqrt(intensity_vec.shape[0])
        
 
        dwelltime = None

        try:
            dwelltime = time_vec[np.where(mean_autocorr < .01)[0][0]]

        except:
            try:
                dwelltime = time_vec[np.where(mean_autocorr < .05)[0][0]]
            except:
                dwelltime = 1




        ke_exp = np.round(geneLength/dwelltime ,1)

        return autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_exp


    def get_crosscorr(self, iv1,iv2):
        '''
        returns the autocorrelations
        '''

        i = 0
        slen = np.correlate(iv1[i]-np.mean(iv1[i]),iv2[i]-np.mean(iv2[i]),'full').shape[0]
        crosscorr_vec = np.zeros((iv1.shape[0],slen))
        
        for i in range(iv1.shape[0]):
            crosscorr_vec[i,:] = np.correlate(iv1[i]-np.mean(iv1[i]),iv2[i]-np.mean(iv2[i]),'full')

        normalized_autocorr = crosscorr_vec.T/ crosscorr_vec[:,0]
        mean_autocorr = np.mean(normalized_autocorr, axis=1)

        return crosscorr_vec, mean_autocorr





class poi():
    '''
    Protein of Intrest class

    Holds all the information for analyzed proteins
    '''
    def __init__(self):
        self.aa_seq = ''    #amino sequence
        self.nt_seq = ''    #nucleotide sequence
        self.gene_length = 0   #length of the gene
        self.tag_length = 0   #length of the tags
        self.total_length = 0  #total length of the full amino acid sequence
        self.name = ''         #name of the gene
        self.tag_types = []
        self.tag_epitopes = {}  #type of tags and epitope lists per tag

class ssa():
    '''
    SSA container class

    holds intensity / ribosome data as well as the propensities used

    __.n_traj = number of trajectories
    __.k = propensities used for the simulation
    __.rib_density = ribosome density per mRNA strand
    __.ribosome_means


    '''
    def __init__(self):
        self.n_traj = 0   #number trajectories
        self.k = []       #propensities
        self.no_rib_per_mrna = 0    #number of ribosomes per mrna strand
        self.rib_density = 0      #ribosome density
        self.ribosome_means = 0  #mean ribosomes
        self.rib_vec = 0          #actual vector of ribosome locations
        self.intensity_vec = []   #intensity vectors per SSA trajectory
        self.time_vec_fixed = []   #time array




    def save_txt(self,filename):

        if '.txt' in filename:
            f = open(filename, 'a')
            for key in self.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    f.write((key + '\r\n'))
                    np.savetxt(f, np.atleast_2d(self.__dict__[key]), delimiter=',', fmt='%s')
                    f.write(('\r\n'))

        else:
            filename = filename + '.txt'
            f = open(filename,'a')
            for key in self.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    f.write((key + '\r\n'))
                    np.savetxt(f, np.atleast_2d(self.__dict__[key]), delimiter=',', fmt='%s')
                    f.write(('\r\n'))
        f.close()

    def load_from_txt(self, filename):
        if '.txt' in filename:
            ssa_obj = np.loadtxt(filename, dtype=str,delimiter='\n')
            solutions = []
            for i in range(0,len(ssa_obj)-1):
                label = ssa_obj[i]
                
                
                if label in ['rib_means',
                             'rib_vec',
                             'n_traj',
                             'start_time',
                             'k',
                             'time_vec_fixed',
                             'dwelltime',
                             'mean_autocorr',
                             'no_rib_per_mrna',
                             'ke_sim',
                             'autocorr_vec',
                             'ribosome_means',
                             'error_autocorr',
                             'rib_density',
                             'time',
                             'ke_true',
                             'evaluating_inhibitor',
                             'time_inhibit',
                             'evaluating_frap']:

                    if label in ['start_time', 'no_rib_per_mrna', 'ke_sim', 'dwelltime','ke_true','time_inhibit']:
                        
                        array = np.fromstring(ssa_obj[i+1], dtype=float, sep=',')[0]
                        exec(('self.'+label+ '=array'))
                    elif label in ['n_traj']:
                        array = int(np.fromstring(ssa_obj[i+1], dtype=float, sep=',')[0])
                        exec(('self.'+label+ '=array'))
                    else:
                        array = np.fromstring(ssa_obj[i+1], dtype=float, sep=',')
                        exec(('self.'+label+ '=array'))

                if label in ['evaluating_inhibitor','evaluating_frap']:
                    if 'False' in ssa_obj[i+1]:                         
                        exec(('self.'+label+ '=False'))
                    if 'True' in ssa_obj[i+1]:                         
                        exec(('self.'+label+ '=True'))


            for i in range(0,len(ssa_obj)-1):
                label = ssa_obj[i]                    
                    
                if label == 'intensity_vec':

                    tvec = self.time_vec_fixed[np.where(self.time_vec_fixed >= self.start_time)]
                    i_vec = np.zeros((self.n_traj, len(self.time)))

                    for j in range(self.n_traj):
                        array = np.fromstring(ssa_obj[i+j+1], dtype=float,sep=',')
                        i_vec[j] = array
                        
                    exec(('self.'+label+ '=i_vec'))
                     
                if label == 'solutions':    
                    for j in range(self.n_traj):
                        array = np.fromstring(ssa_obj[i+j+1], dtype=float,sep=',')
                        solutions.append(array)
                        
                    exec(('self.'+label+ '=solutions'))
                        
                  
                    
                    
                     
                     
                     

    def save_from_json(self, filename):

        if '.json' in filename:

            ssadict = {}
            for key in self.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa_harr.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa_harr.__dict__[key]

            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

        else:
            filename =  filename + '.json'

            ssadict = {}
            for key in self.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa_harr.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa_harr.__dict__[key]

            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)


    def load_json(self,filename):
        if '.json' in filename:

            obj_text = codecs.open(filename, 'r', encoding='utf-8').read()
            ssadict = json.loads(obj_text)

            for key in ssadict.keys():
                if key in ['rib_means','time_vec_fixed','intensity_vec','mean_autocorr','autocorr_vec','error_autocorr','rib_density']:

                    self.__dict__[key] = np.array(ssadict[key])
                else:
                    self.__dict__[key] = ssadict[key]


class GenericSSA():

    '''
    Generic SSA solver - used for the simulated cell animations
    '''

    def __init__(self,type='linear'):

        self.time_variant = False
        self.xi=np.array([])
        self.ti= None
        self.tf=None
        self.S=np.array([])
        self.type=type
        self.ptimes=100
        self.params={}
        if type=='linear':
            #self.fast_rxn = 0.5
            self.W0=np.array([])
            self.W1=np.array([])
        if type == 'nonlinear':
            #self.fast_rxn = 0.5
            self.P=lambda x,t:None



    def gettvec(self):
        return np.linspace(self.ti,self.tf,self.ptimes)

    def updateparams(self):
        for k, v in self.params.iteritems():
            setattr(self, k, v)

    def _run_trajectory(self):#!!!!!!!!!!!!!!!!!renamed run to solve(big deal)
        x=self.xi
        t=self.ti
        __n=len(x)
        self.time=self.gettvec()
        data=np.zeros((len(self.xi),self.ptimes))
        ip=0

        if self.type=='linear':
            if self.time_variant == False:
                while t<self.tf:
                    rate=np.atleast_2d(np.dot(self.W1,x))+self.W0
                    rate=np.cumsum(rate)
                    with np.errstate(divide='ignore', invalid='ignore'):
                        t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()
                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            data[:,ip]=x.reshape(__n)
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i
                            break
                    x=x+np.atleast_2d(self.S[:,event]).T

            else:


                x = np.concatenate((np.array([0]),self.xi.flatten()))

                #x = self.xi
                __n=len(x)
                self.time=self.gettvec()
                data=np.zeros((len(self.xi),self.ptimes))
                a,b = self.S.shape
                S = np.vstack((np.zeros(b),self.S))
                S = np.hstack((np.zeros((a+1,1)),S))
                while t<self.tf:
                    __n=len(x)
                    self.time=self.gettvec()
                    data=np.zeros((len(self.xi),self.ptimes))
                    a,b = self.S.shape
                    S = np.vstack((np.zeros(b),self.S))
                    S = np.hstack((np.zeros((a+1,1)),S))
                    while t<self.tf:
                        trate=self.get_P(x[1:],t)

                        rate = np.concatenate((np.array([self.fast_rxn]),trate))
                        rate=np.cumsum(rate)


                        t=(t-np.log(np.random.rand(1))/rate[-1])
                        ro=rate[-1]*np.random.rand()

                        while t>self.time[ip]:
                            if t>self.tf:
                                b = len(self.time[ip:])
                                fill = np.repeat(x[1:],b)
                                data[:,ip:]=fill.reshape(__n-1,b)
                                return data
                            else:
                                #data[:,ip]=x.reshape(__n)
                                data[:,ip]=x[1:]
                                ip=ip+1
                        for i in range(len(rate)):
                            if rate[i]>=ro:
                                event=i

                                break

                        x=x+S[:,event].ravel()
                    '''
                    rate=np.atleast_2d(np.dot(self.W1(t),x))+self.W0(t)

                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    print(t)
                    ro=rate[-1]*np.random.rand()
                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            data[:,ip]=x.reshape(__n)
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i
                            break
                    x=x+np.atleast_2d(self.S[:,event]).T

                    '''



        elif self.type=='nonlinear':
            if self.time_variant == True:  #if time variant use fast reaction
                x = np.concatenate((np.array([0]),self.xi.flatten()))

                #x = self.xi
                __n=len(x)
                self.time = self.gettvec()
                data = np.zeros((len(self.xi), self.ptimes))
                a, b = self.S.shape
                S = np.vstack((np.zeros(b), self.S))
                S = np.hstack((np.zeros((a+1, 1)), S))
                while t < self.tf:
                    trate=self.get_P(x[1:],t)
                    rate = np.concatenate((np.array([self.fast_rxn]),trate))
                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()

                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x[1:],b)
                            data[:,ip:]=fill.reshape(__n-1,b)
                            return data
                        else:
                            #data[:,ip]=x.reshape(__n)
                            data[:,ip]=x[1:]
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i

                            break

                    x=x+S[:,event].ravel()

            else:   #if not time variant ignore fast reaction

                x = self.xi.flatten()

                #x = self.xi
                __n=len(x)
                self.time=self.gettvec()

                while t<self.tf:
                    rate=self.get_P(x,t)

                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()

                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            #data[:,ip]=x.reshape(__n)
                            data[:,ip]=x
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i

                            break

                    x=x+self.S[:,event].ravel()



        else:
            'Error'
        self.data=data
        return data

    def _solve(self,n):
        __data=np.zeros((len(self.xi),self.ptimes,n))

        for i in range(n):
            __d=self._run_trajectory()
            __data[:,:,i]=__d
        self.data = __data

        return __data

    def setpar(self,key,val):
        self.params[key]=val

    def get_dist(self,specID=0):
        '''
        build distribution (non-normalized and pdf)
        of rna for the model)
        '''
        n_specs, n_times, n_traj = self.data.shape
        n_specs = int(n_specs)
        n_times = int(n_times)
        n_traj = int(n_traj)
        specID = int(specID)
        max_rna = int(np.max(self.data[specID,:,:]))
        self.pdf = np.zeros((n_times,max_rna+1))
        self.fdist = np.zeros((n_times,max_rna+1))
        for i in range(n_times):
            ind = int(i)
            for j in range(n_traj):
                jnd = int(j)
                self.fdist[ind,int(self.data[int(specID),ind,jnd])] +=1
            self.pdf[ind,:] = self.fdist[ind,:] / np.sum(self.fdist[ind,:])

    def get_traj(self,specID=0,ntraj='all'):
        n_specs, n_times, n_traj = self.data.shape
        n_specs = int(n_specs)
        n_times = int(n_times)




        if isinstance(specID,int):
            if ntraj == 'all':
                return self.data[specID],ntraj
            else:
                try:
                    ntraj = ntraj.flatten().astype(int).tolist()
                except:
                    ntraj = int(ntraj)
                    pass

                return self.data[specID][:,ntraj],ntraj
        else:
            if specID == 'all':
                if ntraj == 'all':
                    return self.data,ntraj
                else:

                    try:
                        ntraj = ntraj.flatten().astype(int).tolist()
                    except:
                        pass

                    return self.data,ntraj

    def get_means(self,specID=0):
        '''
        get the first moment.
        '''
        n_specs, n_times, n_traj = self.data.shape
        max_rna = int(np.max(self.data[int(specID),:,:])+1)
        self.means = np.zeros(n_times)
        for i in range(n_times):
            self.means[i] = np.sum(np.arange(max_rna)*self.pdf[i,:])

    def get_variances(self,specID=0):
        '''
        get the second moment.
        '''
        self.get_dist()
        n_specs, n_times, n_traj = self.data.shape
        max_rna = int(np.max(self.data[specID,:,:])+1)
        self.variances = np.zeros(n_times)
        self.covariances = np.zeros((n_specs,n_specs,n_times))
        for i in range(n_times):
            self.variances[i] = np.sum((np.arange(max_rna)**2)*self.pdf[i,:])-(np.sum(np.arange(max_rna)*self.pdf[i,:])**2)
            self.covariances[:,:,i] = np.cov(self.data[:,i,:])

    def return_all_var(self):
        all_members = self.__dict__.keys()

        return [ (item, self.__dict__[item]) for item in all_members if not item.startswith("_")]

    def return_names(self):
        all_members = self.__dict__.keys()
        return [ item for item in all_members if not item.startswith("_")]







'''
if __name__ == '__main__': #do the main
    sms = sms()
    sms.open_seq_file("KDM5B.txt")
    sms.get_orfs(sms.sequence_str, min_codons = 80)
    sms.get_temporal_proteins()
    sms.analyze_poi(sms.pois[0],sms.pois_seq[0])
'''



'''



Unimplemented codes


def SSA_c(self):

    a no argument version of SSA for multithreading

    *****Depreciated****

    k = self.k
    t_array = self.t
    #SSA params and propensities
    R=10 #exclusion volume
    kelong = np.array([k[1:-1]]).T

    N = len(kelong)
    kbind = k[0]
    kcompl = k[-1]
    X = np.array([0,0])

    t = t_array[0]
    Nt = len(t_array)
    tf = t_array[-1]
    N_rib = 200
    X_array = np.zeros((N_rib,Nt))

    it = 1
    Inhibitor = False
    FRAP = False
    inhibit_time = 0

    while t < tf:

        if Inhibitor == True:
            if t >= inhibit_time:
                inhibit_off = False
                Inhibit_condition = 0
            else:
                inhibit_off = True
                Inhibit_condition = 1
        else:
            Inhibit_condition = 1
        if FRAP == True :
            if t >= inhibit_time and t < inhibit_time + 20:
                X = np.array([0,0])


        NR = len(np.where(X>0)[0])

        if len(X)< NR+1:

            X = np.append(X,[0])

        X[max(NR,1)] = 0
        X = X[0:max(NR,1)+1]
        Sn = np.eye(max(NR+1,2))
        wn = np.zeros((X.shape[0],1))
        inds = np.where(X>0)
        try:
            inds[0]
            elongpass = False
        except:
            elongpass = True
        if elongpass == False:
            if len(X[inds[0]]) != 0:

                wn[inds[0]] = kelong[inds[0]]

        if X[0] == N:

            Sn[:,0] = (np.append(X[1:],np.array([0]))  -X[0:])

            wn[0] = kcompl


        if NR == 0 or X[NR-1] > R:

            wn[NR] = kbind*Inhibit_condition


        x = X
        REST = np.less(x[1:]+10,x[0:-1])


        wn[1:] = (wn[1:].T*REST).T

        w = wn
        S = Sn
        w0=np.cumsum(w)
        #print(w0)
        t=(t-np.log(np.random.rand(1))/w0[-1])
        #print(t)

        while it < Nt and t > t_array[it]:
            X_array[0:len(X),it] = X
            it +=1

        if t < tf:
            r2 = w0[-1]*np.random.rand(1)[0]

            for i in range(len(w0)):
            #event = np.where(w0>=r2)[0][0]

                if w0[i]>=r2:
                    event=i
                    break


        X = (X.T + S[:,event]).T

    return X_array










def SSA_2(self,k,t_array,inhibit_time=0,FRAP=False,Inhibitor=False):

    Translation simulation

    *****Depreciated, was used for testing 2 versions of ssa against eachother for efficency testing****


    #SSA params and propensities
    R=10 #exclusion volume
    kelong = np.array([k[1:-1]]).T

    N = len(kelong)
    kbind = k[0]
    kcompl = k[-1]
    X = np.array([0,0],dtype=int)

    t = t_array[0]
    Nt = len(t_array)
    tf = t_array[-1]
    N_rib = 200
    X_array = np.zeros((N_rib,Nt))
    NR = 0
    it = 1
    Sn_p = np.eye(max(NR+1,2),dtype=int)
    wn_p = np.zeros((X.shape[0],1))
    Inhibit_condition = 1
    while t < tf:


        if Inhibitor == True:
            if t >= inhibit_time:

                Inhibit_condition = 0
            else:

                Inhibit_condition = 1
        else:
            Inhibit_condition = 1
        if FRAP == True :
            if t >= inhibit_time and t < inhibit_time + 20:
                X = np.array([0,0])



        oldNR = NR
        NR = len(np.where(X>0)[0])


        if len(X)< NR+1:

            X = np.append(X,[0])


        X[max(NR,1)] = 0
        X = X[0:max(NR,1)+1]



        #Sn = np.eye(max(NR+1,2))
        #wn = np.zeros((X.shape[0],1))

        if oldNR !=NR:
            Sn_p = np.eye(max(NR+1,2),dtype=int)
            wn_p = np.zeros((X.shape[0],1))

        Sn = Sn_p
        wn = wn_p

        inds = X>0

        """
        try:
            inds[0]
            elongpass = False
        except:
            elongpass = True
        if elongpass == False:
        """

        X_inds = X[inds]

        #if X_inds.shape[0] > 0 :

        wn[inds] = kelong[X_inds-1]



        if X[0] == N:

            Sn[:,0] = (np.append(X[1:],np.array([0]))  -X[0:])

            wn[0] = kcompl


        if NR == 0 or X[NR-1] > R:

            wn[NR] = kbind*Inhibit_condition



        REST = np.less(X[1:]+10,X[0:-1])


        wn[1:] = (wn[1:].T*REST).T

        w = wn

        S = Sn
        w0=np.cumsum(w)

        #w0= np.sum(w)

        #print(w0)
        randnum = np.random.random(2)

        t=(t-np.log(randnum[0])/w0[-1])
        #t=(t-np.log(randnum[0])/w0)
        #print(t)

        while it < Nt and t > t_array[it]:
            X_array[0:len(X),it] = X
            it +=1

        if t < tf:
           # r2 = w0*np.random.rand(1)

            r2 = w0[-1]*randnum[1]





            event = np.where(w0>=r2)[0][0]


        X = (X.T + S[:,event]).T



    return X_array


'''




'''
def SSA_f(self,k,t_array,inhibit_time=0,k_on = .3,k_off = .1, FSloc = 200,FRAP=False,Inhibitor=False):
    """
    Translation simulation

    **unimplemented frameshift in python code**

    """

    #SSA params and propensities
    R=10 #exclusion volume
    kelong = np.array([k[1:-1]]).T

    N = len(kelong)/2
    kbind = k[0]
    kcompl = k[-1]
    X = np.array([0,0],dtype=int)

    t = t_array[0]
    Nt = len(t_array)
    tf = t_array[-1]
    N_rib = 200
    X_array = np.zeros((N_rib,Nt))
    NR = 0
    it = 1

    z = np.random.binomial(1,  k_on/(k_on+k_off)  )

    Sn_p = np.eye(max(NR+1,2),dtype=int)
    wn_p = np.zeros((X.shape[0],1))

    #wn_p = np.zeros((1,X.shape[0])).flatten()
    wshape = len(wn_p)
    Inhibit_condition = 1
    while t < tf:


        if Inhibitor == True:
            if t >= inhibit_time:

                Inhibit_condition = 0
            else:

                Inhibit_condition = 1
        else:
            Inhibit_condition = 1
        if FRAP == True :
            if t >= inhibit_time and t < inhibit_time + 20:
                X = np.array([0,0])



        oldNR = NR
        #NR = len(np.where(X>0)[0])
        #NR = len(np.where(X!=0)[0])



        #NR = len(np.argwhere(X))
        #NR = np.nonzero(X)[0].shape[0]
        #NR = max(0,len(X)-1)
        #NR = np.sum(X!=0)
        #NR = np.where(X!=0)[0][-1]+1
        #NR = np.flatnonzero(X).shape[0]
        NR = len(np.flatnonzero(X))



        if X.shape[0]< NR+1:

            X = np.append(X,[0])


        #X[max(NR,1)] = 0
        X[-1] = 0

        X = X[0:max(NR,1)+1]




        #Sn = np.eye(max(NR+1,2))
        #wn = np.zeros((X.shape[0],1))

        if oldNR !=NR:
            Sn_p = np.eye(max(NR+1,2),dtype=int)
            wn_p = np.zeros((X.shape[0],1))
            #wn_p = np.zeros((1,X.shape[0])).flatten()
            wshape = len(wn_p)

        Sn = Sn_p
        wn = wn_p

        inds = X>0

        """
        try:
            inds[0]
            elongpass = False
        except:
            elongpass = True
        if elongpass == False:
        """

        #X_inds = X[inds]


        #if X_inds.shape[0] > 0 :

        wn[inds] = kelong[X[inds]-1]



        if X[0] == N:

            Sn[:,0] = (np.append(X[1:],np.array([0]))  -X[0:])

            wn[0] = kcompl


        if NR == 0 or X[NR-1] > R:
            wn[NR] = kbind*Inhibit_condition

        REST = np.less(X[1:]+10,X[0:-1])


        wn[1:] = (wn[1:].T*REST).T

        #w = wn#.flatten()
        #S = Sn
        w0= sum(wn.flat)
        randnum = np.random.random_sample(2)
        t=(t-np.log(randnum[0])/w0)

        while it < Nt and t > t_array[it]:
            X_array[0:len(X),it] = X
            it +=1

        if t < tf:
            r2 = w0*randnum[1]
            tmp = 0

            for i in range(wshape):
                tmp = tmp + wn[i]
                if tmp >=r2:
                    event=i
                    break

        X = (X + Sn[:,event].T)



    return X_array

'''

'''
vertical stacking ssa, slower but correct
def SSA_stacking(self,k,t_array,inhibit_time=0,FRAP=False,Inhibitor=False):

    #SSA params and propensities
    R=10 #exclusion volume
    kelong = np.array([k[1:-1]]).T

    N = len(kelong)
    kbind = k[0]
    kcompl = k[-1]
    X = np.array([[0],[0]])

    t = t_array[0]
    Nt = len(t_array)
    tf = t_array[-1]
    N_rib = 200
    X_array = np.zeros((N_rib,Nt))

    it = 1

    while t < tf:
        if Inhibitor == True:
            if t >= inhibit_time:
                inhibit_off = False
                Inhibit_condition = 0
            else:
                inhibit_off = True
                Inhibit_condition = 1
        else:
            Inhibit_condition = 1
        if FRAP == True :
            if t >= inhibit_time and t < inhibit_time + 20:
                X = np.array([[0]])

        NR = len(np.where(X>0)[0])

        if X.shape[0]< NR+1:

            X = np.vstack((X,np.array([[0]])))

        X[max(NR,1),0] = 0
        X = X[0:max(NR,1)+1]
        Sn = np.eye(max(NR+1,2))
        wn = np.zeros((X.shape[0],1))
        inds = np.where(X>0)
        try:
            inds[0]
            elongpass = False
        except:
            elongpass = True
        if elongpass == False:
            if len(X[inds[0]]) != 0:

                wn[inds[0]] = kelong[inds[0]]

        if X[0] == N:

            Sn[:,0] = (np.vstack((X[1:],np.array([[0]])))  -X[0:]).flatten()

            wn[0] = kcompl


        if NR == 0 or X[NR-1] > R:

            wn[NR] = kbind*Inhibit_condition


        x = X.flatten()
        REST = np.less(x[1:]+10,x[0:-1])


        wn[1:] = (wn[1:].T*REST).T

        w = copy.deepcopy(wn)
        S = copy.deepcopy(Sn)
        w0=np.cumsum(w)
        #print(w0)
        t=(t-np.log(np.random.rand(1))/w0[-1])
        #print(t)

        while it < Nt and t > t_array[it]:
            X_array[0:len(X),it] = X.flatten()
            it +=1

        if t < tf:
            r2 = w0[-1]*np.random.rand(1)[0]

            for i in range(len(w0)):

                if w0[i]>=r2:
                    event=i
                    break

        X = (X.T + S[:,event]).T

    return X_array

'''


