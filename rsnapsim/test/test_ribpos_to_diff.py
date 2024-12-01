# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 11:53:44 2021

@author: wsraymon
"""

import os
print(os.getcwd())
os.chdir('../..')
import rsnapsim
import numpy as np
import matplotlib.pyplot as plt
print(os.getcwd())

proteins, protein_strs, poi_objs, seqs = rsnapsim.seqmanip.open_seq_file('./rsnapsim/gene_files/Bactin_withTags.txt' )
rsnapsim.solver.protein = poi_objs['1'][0]
pro = poi_objs['1'][0]
t = np.linspace(0,500,501)
bactin = poi_objs['1'][0]
soln = rsnapsim.solver.solve_ssa(pro.kelong, t,n_traj=1,ki=.033, low_memory=False)

X = soln.ribosome_locations
bactin_mrna_blank_mw = rsnapsim.diffcalc.calculate_rna_strand_base_mw(bactin.nt_seq, n_loops=24, fluorophore=0, coat_protein=0 )
aa_mw_vec = rsnapsim.diffcalc.calculate_single_rib_mw(bactin.aa_seq,bactin.probe_loc,)


pl = np.vstack([bactin.probe_loc,bactin.probe_loc])
mw_vec = rsnapsim.diffcalc.mw_over_time(X, bactin.nt_seq, pl, bactin_mrna_blank_mw)
print(np.mean(mw_vec))
df_vec = rsnapsim.diffcalc.calculate_diffusion_constant(mw_vec)
print(np.mean(df_vec))
plt.plot(mw_vec.T)
plt.figure()
plt.plot(df_vec.T)
