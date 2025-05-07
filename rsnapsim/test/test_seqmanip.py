# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 12:54:52 2021

@author: willi
"""

##############################################################################
# Testing file for the sequence manipulation class
# Current functions tested:
#   optimization - converting a codon sequence to highest value 
#                  given a codon to value and codon to aa dictionary
#   deoptimization - reverse of above
#
#
#
##############################################################################


import os
cwd = os.getcwd()
os.chdir('../../')
print(os.getcwd())
import rsnapsim as rsnp
from rsnapsim import seqmanip

import numpy as np
import time


os.chdir(cwd)


import unittest

example_mRNA = '''ATGGCGAACCTTGGCTGCTGGATGCTGGTTCTCTTTGTGGCCACATGGAGTGACCTGGGC
                    CTCTGCAAGAAGCGCCCGAAGCCTGGAGGATGGAACACTGGGGGCAGCCGATACCCGGGG
                    CAGGGCAGCCCTGGAGGCAACCGCTACCCACCTCAGGGCGGTGGTGGCTGGGGGCAGCCT
                    CATGGTGGTGGCTGGGGGCAGCCTCATGGTGGTGGCTGGGGGCAGCCCCATGGTGGTGGC
                    TGGGGACAGCCTCATGGTGGTGGCTGGGGTCAAGGAGGTGGCACCCACAGTCAGTGGAAC
                    AAGCCGAGTAAGCCAAAAACCAACATGAAGCACATGGCTGGTGCTGCAGCAGCTGGGGCA
                    GTGGTGGGGGGCCTTGGCGGCTACATGCTGGGAAGTGCCATGAGCAGGCCCATCATACAT
                    TTCGGCAGTGACTATGAGGACCGTTACTATCGTGAAAACATGCACCGTTACCCCAACCAA
                    GTGTACTACAGGCCCATGGATGAGTACAGCAACCAGAACAACTTTGTGCACGACTGCGTC
                    AATATCACAATCAAGCAGCACACGGTCACCACAACCACCAAGGGGGAGAACTTCACCGAG
                    ACCGACGTTAAGATGATGGAGCGCGTGGTTGAGCAGATGTGTATCACCCAGTACGAGAGG
                    GAATCTCAGGCCTATTACCAGAGAGGATCGAGCATGGTCCTCTTCTCCTCTCCACCTGTG
                    ATCCTCCTGATCTCTTTCCTCATCTTCCTGATAGTGGGATGA'''.replace('\n','').replace(' ','')    

example_protein = 'MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG*'




class TestSeqmanip(unittest.TestCase):


##############################################################################            
# Sequence Opening Tests
    
   
    def test_parsing_multiline_fasta(self):
        aa = 'MKGPILLGTCTSYPGASILSTSTGWASTTASAPAGLFTCLVEASISFRSLRKGILMVRCMRPRKTALPSWSSSTCGRSTPVRCWRAPGSSMSCPTTEAGSTCWTRRSTGSPSTGVQLPQLSSLSAALWSDDTDAAKRWLALSSK*'
        example_file_paths = './test_gene_files/'
        files = os.listdir(example_file_paths)
        for f in files:
            if f == 'multiline_fasta.fasta':
                a,c,b,d = rsnp.seqmanip.open_seq_file(example_file_paths + f, add_tag=True)    
        self.assertEqual(a['0'][0], aa) #check that both sequences are right from each multiline
        self.assertEqual(a['0'][1], aa)
            
            
        
##############################################################################            
# Sequence Optimization Tests
    
    def test_unknown_codon_to_aa_optimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon is provided and cant be decoded
        '''
        test_sequence = 'aaacccggguuuaax'
        with self.assertRaises(rsnp.custom_errors.UnrecognizedAAError): 
            rsnp.seqmanip.optimize_ntseq(test_sequence)
            
    def test_invalid_length_sequence_nt_optimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon for optimization
        '''
        test_sequence = 'aaacccggguuuaaau'
        with self.assertRaises(rsnp.custom_errors.InvalidSequenceLengthError): 
            rsnp.seqmanip.optimize_ntseq(test_sequence)
        
    def test_unknown_codon_to_optimize_nt_optimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon for optimization, codon missing from the optimization dictionary
        '''
        test_sequence = 'aaacccggguuuaau'
        opt_dict = {'AAA':1,'CCC':2,'GGG':3,'UUU':3,}
        with self.assertRaises(rsnp.custom_errors.UnrecognizedCodonError): 
            rsnp.seqmanip.optimize_ntseq(test_sequence, opt_dict=opt_dict)
            
##############################################################################            
# Sequence DeOptimization Tests
    
    def test_unknown_codon_to_aa_deoptimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon is provided and cant be decoded
        '''
        test_sequence = 'aaacccggguuuaax'
        with self.assertRaises(rsnp.custom_errors.UnrecognizedAAError): 
            rsnp.seqmanip.deoptimize_ntseq(test_sequence)
            
    def test_invalid_length_sequence_nt_deoptimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon for optimization
        '''
        test_sequence = 'aaacccggguuuaaau'
        with self.assertRaises(rsnp.custom_errors.InvalidSequenceLengthError): 
            rsnp.seqmanip.deoptimize_ntseq(test_sequence)
        
    def test_unknown_codon_to_optimize_nt_deoptimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon for optimization, codon missing from the optimization dictionary
        '''
        test_sequence = 'aaacccggguuuaau'
        opt_dict = {'AAA':1,'CCC':2,'GGG':3,'UUU':3,}
        with self.assertRaises(rsnp.custom_errors.UnrecognizedCodonError): 
            rsnp.seqmanip.deoptimize_ntseq(test_sequence, deopt_dict=opt_dict)
        

##############################################################################            
# kmer Tests    

    def test_kmer_freq_wo_substrings(self):
        test_sequence = 'aaacccuuugggaaa'
        
        kmer_vec, kmer_inds = seqmanip.get_kmer_freq(test_sequence, 3)
        self.assertEqual(np.sum(kmer_vec), 13)
        self.assertEqual(kmer_vec[kmer_inds.index('UUU')], 1) #uuu = 1
        self.assertEqual(kmer_vec[kmer_inds.index('UUG')], 1) #uug = 1
        self.assertEqual(kmer_vec[kmer_inds.index('AAA')], 2) #aaa = 2 for first and last
        self.assertEqual(len(kmer_vec), 64) # total possible str combos
    
    def test_kmer_freq_w_substrings(self):
        test_sequence = 'aaacccuuugggaaa'
        
        kmer_vec, kmer_inds = seqmanip.get_kmer_freq(test_sequence, 3, substrings=True)
        self.assertEqual(np.sum(kmer_vec), 41)
        self.assertEqual(kmer_vec[kmer_inds.index('UUU')], 1) #uuu = 1
        self.assertEqual(kmer_vec[kmer_inds.index('UUG')], 1) #uug = 1
        self.assertEqual(kmer_vec[kmer_inds.index('A')], 6) #a = 6 a's total
        self.assertEqual(len(kmer_vec), 84) # total possible str combos
    
    '''
    def test_kmer_freq_large_warning(self):
        test_sequence = 'aaacccuuugggaaannneee'
        with self.assertWarns(Warning):
            kmer_vec, kmer_inds = seqmanip.get_kmer_freq(test_sequence, 8, substrings=True)

    '''

##############################################################################            
# nt2aa Tests        

    def test_nt2aa(self):
        self.assertEqual(seqmanip.nt2aa(example_mRNA),example_protein)

    def test_nt2aa_invalid_length(self):
        with self.assertRaises(rsnp.custom_errors.InvalidSequenceLengthError): 
            seqmanip.nt2aa(example_mRNA + 'A')

if __name__ == '__main__':
    unittest.main()