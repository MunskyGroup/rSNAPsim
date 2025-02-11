# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 15:08:34 2025

@author: willi
"""

##############################################################################
# Testing file for the saving custom solutions
#
#
#
##############################################################################


import os
cwd = os.getcwd()
os.chdir('../../')
print(os.getcwd())
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





class TestSSAsoln_saving(unittest.TestCase):

    aa = 'MKGPILLGTCTSYPGASILSTSTGWASTTASAPAGLFTCLVEASISFRSLRKGILMVRCMRPRKTALPSWSSSTCGRSTPVRCWRAPGSSMSCPTTEAGSTCWTRRSTGSPSTGVQLPQLSSLSAALWSDDTDAAKRWLALSSK*'
    example_file_paths = './test_gene_files/'
    files = os.listdir(example_file_paths)
    for f in files:
        if f == 'H2B_withTags.txt':
            a,c,b,d = rss.seqmanip.open_seq_file(example_file_paths + f, add_tag=False)        

    poi = b['0'][0]
    t = np.linspace(0,500,501)
    
    example_kelong = np.ones(100)
    example_init = .03
    
    def test_save_default_model_npy(self):
        st = time.time()
        soln = rss.solver.solve_ssa(self.poi,self.t, n_traj=1, seed=1 )
        soln.save('./test_solution_save', fmt='.npy')
        return     

    def test_load_poi_model_npy(self):
        st = time.time()
        soln = rss.solver.solve_ssa(self.poi,self.t, n_traj=1, seed=1 )
        soln.save('./test_solution_save2', fmt='.npy')

        soln2 = rss.solver.load_soln('./test_solution_save2.npy')
        
        error_sum = [np.sum(soln.ribosome_array.flatten() != soln2.ribosome_array.flatten()),
        np.sum(soln.state_array.flatten() != soln2.state_array.flatten()),
        np.sum(soln.resource_array.flatten() != soln2.resource_array.flatten()),
        np.sum(soln.kelong_mat.flatten() != soln2.kelong_mat.flatten()),
        np.sum(soln.t.flatten() != soln2.t.flatten())]
        self.assertAlmostEqual(np.sum(error_sum), 0)
    
    def test_load_poi_model_npz(self):
        st = time.time()
        soln = rss.solver.solve_ssa(self.poi,self.t, n_traj=1, seed=1 )
        soln.save('./test_solution_save2', fmt='.npz')

        soln2 = rss.solver.load_soln('./test_solution_save2.npz')
        
        error_sum = [np.sum(soln.ribosome_array.flatten() != soln2.ribosome_array.flatten()),
        np.sum(soln.state_array.flatten() != soln2.state_array.flatten()),
        np.sum(soln.resource_array.flatten() != soln2.resource_array.flatten()),
        np.sum(soln.kelong_mat.flatten() != soln2.kelong_mat.flatten()),
        np.sum(soln.t.flatten() != soln2.t.flatten())]
        self.assertAlmostEqual(np.sum(error_sum), 0)


if __name__ == '__main__':
    unittest.main()