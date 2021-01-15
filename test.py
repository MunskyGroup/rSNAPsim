# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 15:16:37 2020

@author: willi
"""

#test_everything

import rsnapsim
import numpy as np

#### test a full example

## 1 color sim
proteins, protein_strs, poi_objs = rsnapsim.seqmanip.open_seq_file('./rsnapsim/gene_files/Bactin_withTags.txt' )
rsnapsim.solver.protein = poi_objs['1'][0]
t = np.linspace(0,500,501)
pro = poi_objs['1'][0]
soln = rsnapsim.solver.solve_ssa(pro.kelong, t,n_traj=10)

# 2 color sim
proteins, protein_strs, poi_objs = rsnapsim.seqmanip.open_seq_file('./rsnapsim/gene_files/H2B_2tags.txt' )
pro = poi_objs['1'][0]
pro.tag_epitopes['T_Flag'] = [10,20,30,40,50,60,70]
pro.tag_epitopes['T_Hemagglutinin'] = [300,330,340,350]
rsnapsim.solver.protein = pro
t = np.linspace(0,500,501)

soln2 = rsnapsim.solver.solve_ssa(pro.kelong, t,n_traj=10)
pro.visualize_probe()

acov,err_acov = rsnapsim.inta.get_autocov(soln2.intensity_vec,norm='ind')
acc,acc_err = rsnapsim.inta.get_autocorr(acov)
cross_corr,inds = rsnapsim.inta.get_crosscorr(soln2.intensity_vec)