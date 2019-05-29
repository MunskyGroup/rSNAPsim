# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 20:36:30 2019

@author: willi
"""
import numpy as np

import matplotlib.pyplot as plt
import rSNAPsim

sms2 = rSNAPsim.rSNAPsim()
#sms.strGeneCopy['ACC'] = 3
sms2.open_seq_file("gene_files/Bactin_withtags.txt")
sms2.get_orfs(sms2.sequence_str, min_codons = 80)
sms2.get_temporal_proteins()
sms2.analyze_poi(sms2.pois[0],sms2.pois_seq[0])
rates2 = sms2.get_k(sms2.POI.nt_seq,.03,10)

true_kes = []
kis = np.linspace(0.01,10,35)
1/0
for ki in kis:
    testssa2 = sms2.ssa_solver(n_traj=1,all_k=rates2,tf=1000,tstep=3001)