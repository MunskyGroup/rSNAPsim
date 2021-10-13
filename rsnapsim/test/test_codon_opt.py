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

opt_spike = rss.seqmanip.optimize_ntseq(wt_spike.nt_seq)
deopt_spike = rss.seqmanip.deoptimize_ntseq(wt_spike.nt_seq)


t = np.linspace(0,2000,2001)
ki = .033
rss.solver.protein= wt_spike
blank_probe = np.zeros([1,len(wt_spike.kelong)],dtype=int)
#blank_probe[0,1] = 1
blank_vec = np.cumsum(blank_probe,axis=1)

wt_traj = rss.solver.solve_ssa(wt_spike.kelong,t,ki=ki,n_traj=1, probe_loc=blank_probe, probe_vec=blank_vec,low_memory=False )

