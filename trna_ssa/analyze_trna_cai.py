# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 10:07:09 2020

@author: willi
"""

import os

import sys
sys.path.append('C:/Users/willi/Documents/GitHub/rSNAPsim')

os.chdir('..')
from rss import SequenceManipMethods as smm
from rss import TranslationSolvers as tss
from rss import rSNAPsim as rss
from rss import CodonOptimizer as copt
from rss import IntensityAnalyses as ia

os.chdir('./trna_ssa')


import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
import time

cais = np.load('cais.npy')
wt_taus = np.load('wt_taus.npy')
opt_taus = np.load('opt_taus.npy')

with open('sequences_cai.txt', 'r') as f:
    sequences = f.readlines()





strGeneCopy = {'TTT': 17.6, 'TCT': 15.2, 'TAT': 12.2, 'TGT': 10.6, 'TTC': 20.3,
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

st = time.time()
strGeneCopy.pop('TAG')
strGeneCopy.pop('TAA')
strGeneCopy.pop('TGA')

k_trna = np.array(list(strGeneCopy.values()))

t = np.linspace(0,1400,2801)
k_diffusion = 1
k_bind = .033
elong_scale = 1
optimizer = copt()
rss = rss()
solver = tss()

seq = sequences[1392][:-1]

optimized_seq = optimizer.optimize_ntseq(seq)

poi = rss.seq_to_protein_obj(seq)['1'][0]
print(seq[:40])


predicted_tau = int(np.ceil(len(seq)/10))


t = np.linspace(0,predicted_tau*5 ,2*(predicted_tau*5)+1)

optimized_poi = rss.seq_to_protein_obj(optimized_seq)['1'][0]

solver.protein  = poi

trna_ssa_soln_wt = solver.solve_ssa_trna( np.array(poi.ktrna_id) , k_diffusion,k_bind,elong_scale ,t,n_traj=40,k_trna= k_trna)

solver.protein = optimized_poi
trna_ssa_soln_opt = solver.solve_ssa_trna(np.array(optimized_poi.ktrna_id), k_diffusion,k_bind,elong_scale ,t,n_traj=40,k_trna= k_trna)

print('time for simulation: ')
print(time.time()-st)

acov,err_acov = ia().get_autocov(trna_ssa_soln_wt.intensity_vec[:,predicted_tau*2:,:],norm='global')
acc_wt,err_acorr_wt = ia().get_autocorr(acov)

acov,err_acov = ia().get_autocov(trna_ssa_soln_opt.intensity_vec[:,predicted_tau*2:,:],norm='global')
acc_opt,err_acorr_opt = ia().get_autocorr(acov)

tau_acc_opt = t[np.where(np.mean(acc_opt[0],axis=1) < 0.05)[0][0]]
tau_acc_wt = t[np.where(np.mean(acc_wt[0],axis=1) < 0.05)[0][0]]