# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 11:53:44 2021

@author: wsraymon
"""

import os
import rsnapsim
import numpy as np
os.chdir('../..')
print(os.getcwd())

proteins, protein_strs, poi_objs, seqs = rsnapsim.seqmanip.open_seq_file('./rsnapsim/gene_files/Bactin_withTags.txt' )
rsnapsim.solver.protein = poi_objs['1'][0]
t = np.linspace(0,500,501)
bactin = poi_objs['1'][0]
soln = rsnapsim.solver.solve_ssa(pro.kelong, t,n_traj=1, low_memory=False)

X = soln.ribosome_locations

bactin_mrna_blank_mw = rsnapsim.diffcalc.calculate_rna_strand_base_mw(bactin.nt_seq, n_loops=24 )

mw_vec = rsnapsim.diffcalc.convert_rib_pos_tensor(X, bactin.nt_seq, bactin_mrna_blank_mw)