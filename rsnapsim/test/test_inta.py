# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 12:54:52 2021

@author: willi
"""

##############################################################################
# Testing file for the intensity analyses class
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
from rsnapsim import inta

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


#Generate some example signals to test correlation calculations with 

freq1 = .01
freq2 = .03
np.random.seed(0)
t = np.linspace(0,2999,3000)
example_signal = np.zeros([2,3000,100])
example_signal[0,:,:] = np.sin(freq1*np.vstack([t]*100).T + np.random.randint(0,1000,100))
example_signal[1,:,:] = np.sin(freq2*np.vstack([t]*100).T + np.random.randint(0,1000,100))



example_signal_ragged = []
for j in range(2):
    freq = [freq1,freq2][j]
    signal_array = []
    for i in range(100):
        signal_array.append(np.sin(freq*t + np.random.randint(0,1000))[:np.random.randint(2500,3000)])

    example_signal_ragged.append(signal_array)




freq1 = .01
freq2 = .03
np.random.seed(0)
t = np.linspace(0,2999,3000)
example_signal_cc = np.zeros([2,3000,100])
example_signal_cc[0,:,:] = np.sin(freq1*np.vstack([t]*100).T + 0)
example_signal_cc[1,:,:] = np.sin(freq2*np.vstack([t]*100).T + 15)

example_signal_ragged_cc = []
for j in range(2):
    freq = [freq1,freq2][j]
    offset = [0,15][j]
    signal_array = []
    np.random.seed(0)
    for i in range(100):
        signal_array.append(np.sin(freq*t + offset)[:np.random.randint(2500,3000)])

    example_signal_ragged_cc.append(signal_array)


class TestIntA(unittest.TestCase):


##############################################################################            
# autocorrelation Tests
    
    def test_autocov_ind(self):
        
        #Check individual normalization
        acov, acov_error = inta.get_autocov(example_signal,norm='ind')
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)
        
    def test_autocov_global(self):
        
        #Check Global normalization
        acov, acov_error = inta.get_autocov(example_signal,norm='global')
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)
        
    def test_autocov_ind_scalefix(self):
        
        #Check Global normalization
        acov, acov_error = inta.get_autocov(example_signal,norm='ind',scale_fix=True)
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)


    def test_autocov_jagged_ind(self):
        #Check individual normalization
        acov, acov_error = inta.get_autocov(example_signal_ragged,norm='ind')
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)

    def test_autocov_jagged_global(self):
        #Check Global normalization
        acov, acov_error = inta.get_autocov(example_signal_ragged,norm='global')
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)
        
    def test_autocov_jagged_ind_scalefix(self):
        #Check Global normalization
        acov, acov_error = inta.get_autocov(example_signal_ragged,norm='ind',scale_fix=True)
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)

##############################################################################            
# cross correlation Tests

    def test_autocov_ind(self):
        
        #Check individual normalization
        acov, acov_error = inta.get_autocov(example_signal,norm='ind')
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)
        
    def test_autocov_global(self):
        
        #Check Global normalization
        acov, acov_error = inta.get_autocov(example_signal,norm='global')
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)
        
    def test_autocov_ind_scalefix(self):
        
        #Check Global normalization
        acov, acov_error = inta.get_autocov(example_signal,norm='ind',scale_fix=True)
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)


    def test_autocov_jagged_ind(self):
        #Check individual normalization
        acov, acov_error = inta.get_autocov(example_signal_ragged,norm='ind')
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)

    def test_autocov_jagged_global(self):
        #Check Global normalization
        acov, acov_error = inta.get_autocov(example_signal_ragged,norm='global')
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)
        
    def test_autocov_jagged_ind_scalefix(self):
        #Check Global normalization
        acov, acov_error = inta.get_autocov(example_signal_ragged,norm='ind',scale_fix=True)
        mean_acov = np.mean(acov[0],axis=1)
        
        #cycle values where the autocorrelation function for color 1 should be decorrelated
        x = [int(6.28/freq1)*x + int(6.28/freq1/4) for x in [0,1,2,3]]
        sum_cycles = np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        mean_acov = np.mean(acov[1],axis=1)
        
        #cycle values where the autocorrelation function for color 2 should be decorrelated
        x = [int(6.28/freq2)*x + int(6.28/freq2/4) for x in [0,1,2,3]]
        sum_cycles += np.sum(np.abs(mean_acov[x])) #these should be close to 0 when added up
        
        
        self.assertAlmostEqual(sum_cycles, 0, delta=.3)


if __name__ == '__main__':
    unittest.main()