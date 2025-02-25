# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:06:11 2024

@author: willi
"""

##############################################################################
# Testing file for the default models
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





class TestDefaultSSAModel(unittest.TestCase):

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
    

##############################################################################
# Default SSA models

    def test_default_model_run_py_single(self):
        st = time.time()
        rss.solver.solve_ssa(self.poi,self.t, n_traj=1, seed=1 )
        print('1 python trajectory of 500s took: %s s'%(time.time() - st))
        return 


    def test_default_model_set_seed_py(self):
        st = time.time()
        soln = rss.solver.solve_ssa(self.poi,self.t, n_traj=1, seed=1 )
        soln2 = rss.solver.solve_ssa(self.poi,self.t, n_traj=1, seed=1 )
        self.assertTrue(soln.ribosome_array[0,25].sum() == soln2.ribosome_array[0,25].sum())
        return 


    def test_default_model_run_py_multiple(self):
        st = time.time()
        rss.solver.solve_ssa(self.poi,self.t, n_traj=50, seed=1 )
        print('50 python trajectories of 500s took: %s s'%(time.time() - st))
        return 
    

    # def test_default_model_run_py_multiprocessing(self):
    #     st = time.time()
    #     rss.solver.solve_ssa(self.poi,self.t, n_traj=50, seed=1, parallel=True, cores=16)
    #     print('50 multiproccessing python trajectories of 500s took: %s s'%(time.time() - st))
    #     return 
    
    def test_single_model_run_cpp(self):
        return 
    
    
    def test_multiple_model_run_cpp(self):
        return 


    def test_multithreaded_model_run_cpp(self):
        return 


if __name__ == '__main__':
    unittest.main()