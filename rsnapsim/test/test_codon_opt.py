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