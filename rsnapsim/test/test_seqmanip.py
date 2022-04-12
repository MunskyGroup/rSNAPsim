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
os.chdir('../../..')

import rsnapsim as rss
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
# Sequence Optimization Tests
    
    def test_unknown_codon_to_aa_optimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon is provided and cant be decoded
        '''
        test_sequence = 'aaacccggguuuaax'
        with self.assertRaises(rss.custom_errors.UnrecognizedAAError): 
            rss.seqmanip.optimize_ntseq(test_sequence)
            
    def test_invalid_length_sequence_nt_optimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon for optimization
        '''
        test_sequence = 'aaacccggguuuaaau'
        with self.assertRaises(rss.custom_errors.InvalidSequenceLengthError): 
            rss.seqmanip.optimize_ntseq(test_sequence)
        
    def test_unknown_codon_to_optimize_nt_optimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon for optimization, codon missing from the optimization dictionary
        '''
        test_sequence = 'aaacccggguuuaau'
        opt_dict = {'AAA':1,'CCC':2,'GGG':3,'UUU':3,}
        with self.assertRaises(rss.custom_errors.UnrecognizedCodonError): 
            rss.seqmanip.optimize_ntseq(test_sequence, opt_dict=opt_dict)
            
##############################################################################            
# Sequence DeOptimization Tests
    
    def test_unknown_codon_to_aa_deoptimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon is provided and cant be decoded
        '''
        test_sequence = 'aaacccggguuuaax'
        with self.assertRaises(rss.custom_errors.UnrecognizedAAError): 
            rss.seqmanip.deoptimize_ntseq(test_sequence)
            
    def test_invalid_length_sequence_nt_deoptimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon for optimization
        '''
        test_sequence = 'aaacccggguuuaaau'
        with self.assertRaises(rss.custom_errors.InvalidSequenceLengthError): 
            rss.seqmanip.deoptimize_ntseq(test_sequence)
        
    def test_unknown_codon_to_optimize_nt_deoptimization(self):
        '''
        Test that the unrecognized codon Error is raised when given an unknown 
        codon for optimization, codon missing from the optimization dictionary
        '''
        test_sequence = 'aaacccggguuuaau'
        opt_dict = {'AAA':1,'CCC':2,'GGG':3,'UUU':3,}
        with self.assertRaises(rss.custom_errors.UnrecognizedCodonError): 
            rss.seqmanip.deoptimize_ntseq(test_sequence, deopt_dict=opt_dict)
        

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
        with self.assertRaises(rss.custom_errors.InvalidSequenceLengthError): 
            seqmanip.nt2aa(example_mRNA + 'A')
        
    '''
    def test_pull_file(self):
        if os.path.isfile(os.path.join('.', 'content', 'MN908947.gb')):
            os.remove(os.path.join('.', 'content', 'MN908947.gb'))
                               
                               
        rss.seqmanip.get_gb_file(Ascession_Number, '.\content')
        
        self.assertTrue(os.path.isfile(os.path.join('.',
                                                         'content',
                                                         'MN908947.gb')))

    def test_wrong_path(self):
        with self.assertRaises(rss.SequenceManipMethods.PathDoesNotExistError) as context:
            rss.seqmanip.get_gb_file('aaa', '.\content')
            msg = 'Specified save path does not exist, double check the path'\
            ' specified.'
            self.assertEqual(
                context.exception.msg,
                msg)
            
    def test_wrong_asc(self):
        with self.assertRaises(rss.SequenceManipMethods.AscNumDoesNotExistError) as context:
            rss.seqmanip.get_gb_file('113', '.\contents')
            msg = 'Cannot find given ascession number for genbank, file re'\
                'quest failed.'
            
            print('')
            print(context.exception)
            print('')
            self.assertEqual(
                context.exception.msg,
                msg)
    '''

if __name__ == '__main__':
    unittest.main()