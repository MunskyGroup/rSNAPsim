# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 14:07:33 2019

@author: willi
"""
import rSNAPsim
import matplotlib.pyplot as plt
sms = rSNAPsim.rSNAPsim()
#sms.strGeneCopy['ACC'] = 3
sms.open_seq_file("gene_files/sunKif8b.txt")
sms.get_orfs(sms.sequence_str, min_codons = 80)
sms.get_temporal_proteins()
sms.analyze_poi(sms.pois[0],sms.pois_seq[0])
rates = sms.get_k(sms.POI.nt_seq,.08,3.1)

testssa = sms.ssa_solver(n_traj=1,all_k=rates,tf=2000,tstep=3001)

sms.kymograph(testssa,0,bg_intense=False,show_intense=False,color='white',lw=1)


sms2 = rSNAPsim.rSNAPsim()
#sms.strGeneCopy['ACC'] = 3
sms2.open_seq_file("gene_files/KDM5B_withtags.txt")
sms2.get_orfs(sms2.sequence_str, min_codons = 80)
sms2.get_temporal_proteins()
sms2.analyze_poi(sms2.pois[0],sms2.pois_seq[0])
rates2 = sms2.get_k(sms2.POI.nt_seq,.03,10)


testssa2 = sms2.ssa_solver(n_traj=1,all_k=rates2,tf=2000,tstep=3001)
plt.figure()
sms2.kymograph(testssa2,0,bg_intense=False,show_intense=False,color='white',lw=1)