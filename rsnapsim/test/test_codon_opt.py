# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 11:58:56 2021

@author: willi
"""

import os
cwd = os.getcwd()
os.chdir('../../..')

import rsnapsim as rss
from rsnapsim import seqmanip

import numpy as np
import time
import matplotlib.pyplot as plt

os.chdir(cwd)

aa_seqs, poi_objs, tagged_proteins, raw_seq = rss.seqmanip.open_seq_file('./test_gene_files/MN908947.gb', add_tag=False)
wt_spike = poi_objs['2'][1]

cai, sensitivity, cai_codons = rss.seqmanip.codon_usage(wt_spike.nt_seq)

ccount = rss.seqmanip.get_codon_count_dict(wt_spike.nt_seq)

t = {'AAA': 0.44680851063829785,
 'AAC': 1.0,
 'AAG': 0.7106382978723403,
 'AAU': 0.7072340425531916,
 'ACA': 0.23404255319148937,
 'ACC': 0.1702127659574468,
 'ACG': 0.27489361702127657,
 'ACU': 0.2765957446808511,
 'AGA': 0.1276595744680851,
 'AGC': 0.2553191489361702,
 'AGG': 0.21106382978723404,
 'AGU': 0.1768085106382979,
 'AUA': 0.2978723404255319,
 'AUC': 0.3723404255319149,
 'AUG': 0.45957446808510644,
 'AUU': 0.5114893617021277,
 'CAA': 0.2127659574468085,
 'CAC': 0.3191489361702128,
 'CAG': 0.5148936170212766,
 'CAU': 0.23191489361702128,
 'CCA': 0.30851063829787234,
 'CCC': 0.2127659574468085,
 'CCG': 0.2634042553191489,
 'CCU': 0.32978723404255317,
 'CGA': 0.2127659574468085,
 'CGC': 0.14893617021276595,
 'CGG': 0.21106382978723404,
 'CGU': 0.23404255319148937,
 'CUA': 0.22340425531914893,
 'CUC': 0.18085106382978725,
 'CUG': 0.37829787234042556,
 'CUU': 0.3191489361702128,
 'GAA': 0.2978723404255319,
 'GAC': 0.5531914893617021,
 'GAG': 0.3719148936170213,
 'GAU': 0.38744680851063834,
 'GCA': 0.5851063829787234,
 'GCC': 0.5106382978723404,
 'GCG': 0.5323404255319149,
 'GCU': 0.8442553191489361,
 'GGA': 0.23404255319148937,
 'GGC': 0.43617021276595747,
 'GGG': 0.35148936170212763,
 'GGU': 0.30531914893617024,
 'GUA': 0.2765957446808511,
 'GUC': 0.20212765957446807,
 'GUG': 0.5582978723404255,
 'GUU': 0.32978723404255317,
 'UAC': 0.3723404255319149,
 'UAU': 0.24340425531914897,
 'UCA': 0.22340425531914893,
 'UCC': 0.1702127659574468,
 'UCG': 0.23617021276595743,
 'UCU': 0.2872340425531915,
 'UGC': 0.7872340425531915,
 'UGG': 0.18382978723404256,
 'UGU': 0.4731914893617022,
 'UUA': 0.1276595744680851,
 'UUC': 0.3829787234042553,
 'UUG': 0.21106382978723404,
 'UUU': 0.2521276595744681}

opt_spike = rss.seqmanip.optimize_ntseq(wt_spike.nt_seq,opt_dict=t)
deopt_spike = rss.seqmanip.deoptimize_ntseq(wt_spike.nt_seq)


t = np.linspace(0,2000,2001)
ki = .033
rss.solver.protein= wt_spike
blank_probe = np.zeros([1,len(wt_spike.kelong)],dtype=int)
#blank_probe[0,1] = 1
blank_vec = np.cumsum(blank_probe,axis=1)

wt_traj = rss.solver.solve_ssa(wt_spike.kelong,t,ki=ki,n_traj=1, probe_loc=blank_probe, probe_vec=blank_vec,low_memory=False )

