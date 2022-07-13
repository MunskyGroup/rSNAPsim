# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:26:32 2022

@author: willi
"""

##############################################################################
# Testing file for the protein class
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

import rsnapsim as rss
from rsnapsim import inta

import numpy as np
import time
import collections
compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

os.chdir(cwd)

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
                    ATCCTCCTGATCTCTTTCCTCATCTTCCTGATAGTGGGATGA'''.replace('\n','').replace(' ','') .replace('T','U')   


tag = '''ATGGACTACAAGGACGACGACGACAAAGGTGACTACAAAGATGATGACGATAAAGGCGACTATAAGGACGATGACGACAAGGGCGGAAACTCACTGATCAAGGAAAACATGCGGATGAAGGTGGTGATGGAGGGCTCCGTGAATGGTCACCAGTTCAAGTGCACCGGAGAGGGAGAGGGAAACCCGTACATGGGAACTCAGACCATGCGCATTAAGGTCATCGAAGGAGGTCCGCTGCCGTTCGCTTTCGATATCCTGGCCACTTCGTTCGGAGGAGGGTCGCGCACGTTCATCAAGTACCCGAAGGGAATCCCGGACTTCTTTAAGCAGTCATTCCCGGAAGGATTCACTTGGGAACGGGTGACCCGGTATGAAGATGGAGGTGTGGTGACTGTCATGCAAGATACTTCGCTGGAGGATGGGTGCCTCGTGTACCACGTCCAAGTCCGCGGAGTGAATTTCCCGTCCAACGGACCAGTGATGCAGAAAAAGACGAAGGGTTGGGAACCTAATACTGAAATGATGTACCCCGCAGACGGAGGGCTGAGGGGCTACACCCACATGGCGCTGAAGGTCGACGGAGGAGATTACAAGGATGACGACGATAAGCAACAAGATTACAAAGACGATGATGACAAGGGCCAGCAGGGCGACTACAAGGACGACGACGACAAGCAGCAGGACTACAAAGATGACGATGATAAAGGAGGAGGACATCTGTCCTGTTCGTTCGTGACCACCTACAGATCAAAGAAAACCGTGGGAAACATCAAGATGCCGGGCATTCATGCCGTCGACCACCGCCTGGAGCGGCTCGAAGAATCAGACAATGAGATGTTCGTCGTGCAAAGAGAACATGCCGTGGCCAAGTTCGCGGGACTGGGAGGCGGTGGAGGCGATTACAAAGACGATGATGACAAGGGTGACTATAAAGACGACGATGACAAAGGGGATTACAAGGATGATGATGATAAGGGAGGCGGTGGATCAGGTGGAGGAGGTTCACTGCAG'''

example_protein = 'MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG*'
poi = rss.seqmanip.seq_to_protein_obj(example_mRNA, add_tag=False)['0'][0]

tagged_poi = rss.seqmanip.seq_to_protein_obj(example_mRNA, add_tag=True)['0'][0]



import unittest

class TestPOI(unittest.TestCase):


        
    def test_get_aa_seq(self):
        self.assertEqual(example_protein, poi.aa_seq)
        
    def test_get_nt_seq(self):
        self.assertEqual(example_mRNA, poi.nt_seq)
        
    def test_get_orf(self):
        self.assertEqual('0', poi.orf)
        
    def test_get_5p(self):
        self.assertEqual('', poi.UTR_5p)
        
    def test_get_3p(self):
        self.assertEqual('', poi.UTR_3p)
        
    def test_get_all_k(self):
        self.assertAlmostEqual(3311.515214162576 , sum(poi.all_k))
        
    def test_get_codon_sensitivity(self):
        self.assertAlmostEqual(2.42, np.sum(poi.codon_sensitivity))
        
    def test_get_codons(self):
        codons = [example_mRNA[i:i+3] for i in range(0, len(example_mRNA), 3)]
        self.assertTrue(compare(codons, poi.codons))
        
    def test_get_gc_content(self):
        gc_content =(example_mRNA.count('C') + example_mRNA.count('G'))/len(example_mRNA)
        self.assertAlmostEqual(gc_content, poi.gc_content)        
        
    def test_get_protein_length(self):
        self.assertEqual(254, poi.gene_length)     
        
    def test_get_ke(self):
        self.assertEqual(10,poi.ke)   
        
    def test_get_ke_mu(self):
        self.assertEqual(10,poi.ke_mu)   
        
    def test_get_loc(self):
        self.assertEqual(0, poi.loc[0]) 
        self.assertEqual(762, poi.loc[1])   
        
    def test_get_name(self):
        self.assertEqual('', poi.name)   
        
    def test_get_source_seq(self):
        self.assertEqual(example_mRNA, poi.source_seq)   
        

    def test_get_tag_added(self):
        self.assertEqual(False, poi.tag_added)
        
    def test_get_tag_length(self):
        self.assertEqual(0, poi.tag_length)
        
    def test_get_tag_epitopes(self):
        self.assertEqual(0, len(poi.tag_epitopes))
        
    def test_get_total_length(self):
        self.assertEqual(254, poi.total_length)


class TestPOI_tagged(unittest.TestCase):

        
    def test_get_nt_seq(self):
        self.assertEqual(tag + example_mRNA, tagged_poi.nt_seq)
        
    def test_get_orf(self):
        self.assertEqual('0', tagged_poi.orf)
        
    def test_get_5p(self):
        self.assertEqual('', tagged_poi.UTR_5p)
        
    def test_get_3p(self):
        self.assertEqual('', tagged_poi.UTR_3p)
        
    def test_get_all_k(self):
        self.assertAlmostEqual(8068.131710648465, sum(tagged_poi.all_k))
        
    def test_get_codon_sensitivity(self):
        self.assertAlmostEqual(2.25, np.sum(tagged_poi.codon_sensitivity))
        
    def test_get_codons(self):
        codons = [(tag + example_mRNA)[i:i+3] for i in range(0, len((tag + example_mRNA)), 3)]
        self.assertTrue(compare(codons, tagged_poi.codons))
        
    def test_get_gc_content(self):
        gc_content =((tag + example_mRNA).count('C') + (tag + example_mRNA).count('G'))/len((tag + example_mRNA))
        self.assertAlmostEqual(gc_content, tagged_poi.gc_content)        
        
    def test_get_protein_length(self):
        self.assertEqual(254, tagged_poi.gene_length)     
        
    def test_get_ke(self):
        self.assertEqual(10,tagged_poi.ke)   
        
    def test_get_ke_mu(self):
        self.assertEqual(10,tagged_poi.ke_mu)   
        
    def test_get_loc(self):
        self.assertEqual(0, tagged_poi.loc[0]) 
        self.assertEqual(762, tagged_poi.loc[1])   
        
    def test_get_name(self):
        self.assertEqual('', tagged_poi.name)   
        
    def test_get_source_seq(self):
        self.assertEqual(example_mRNA, tagged_poi.source_seq)   
        
    def test_get_tag_added(self):
        self.assertEqual(True, tagged_poi.tag_added)
        
    def test_get_tag_length(self):
        # length between first and last epitope
        self.assertEqual(316, tagged_poi.tag_length)
        
    def test_get_tag_epitopes(self):
        self.assertTrue(compare([2, 11, 20, 196, 206, 218, 228, 300, 309, 318], tagged_poi.tag_epitopes['T_Flag']))
        
    def test_get_total_length(self):
        self.assertEqual(254+337, tagged_poi.total_length)


if __name__ == '__main__':
    unittest.main()