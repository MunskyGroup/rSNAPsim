# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 14:07:33 2019

@author: willi
"""
import numpy as np
import rSNAPsim
import matplotlib.pyplot as plt
sms = rSNAPsim.rSNAPsim()
#sms.strGeneCopy['ACC'] = 3
sms.open_seq_file("gene_files/sunKif8b.txt")
sms.get_orfs(sms.sequence_str, min_codons = 80)
sms.get_temporal_proteins()
sms.analyze_poi(sms.pois[0],sms.pois_seq[0])
rates = sms.get_k(sms.POI.nt_seq,.08,3.1)

testssa = sms.ssa_solver(n_traj=1,all_k=rates,tf=1000,tstep=3001)

sms.kymograph(testssa,0,bg_intense=False,show_intense=False,color='white',lw=1)


sms2 = rSNAPsim.rSNAPsim()
#sms.strGeneCopy['ACC'] = 3
sms2.open_seq_file("gene_files/KDM5B_withtags.txt")
sms2.get_orfs(sms2.sequence_str, min_codons = 80)
sms2.get_temporal_proteins()
sms2.analyze_poi(sms2.pois[0],sms2.pois_seq[0])
rates2 = sms2.get_k(sms2.POI.nt_seq,.03,10)


testssa2 = sms2.ssa_solver(n_traj=1,all_k=rates2,tf=1000,tstep=3001,force_python=True)
plt.figure()
sms2.kymograph(testssa2,0,bg_intense=False,show_intense=False,color='white',lw=1)





testki  = np.linspace(.02,.06,12)
testrates = np.linspace(2,14,12)

cases = np.meshgrid(testrates,testki)
colcount = np.zeros((len(testrates),len(testrates)))
colcountzeros = np.zeros((len(testrates),len(testrates)))
colcountnonzero = np.zeros((len(testrates),len(testrates)))

for i in range(len(testrates)):
    print(i)
    for j in range(len(testki)):
        rate = testrates[i]
        ki = testki[j]
        rates = sms2.get_k(sms2.POI.nt_seq,ki,rate)
        ssa = sms2.ssa_solver(n_traj=10,all_k=rates,tf=2000,tstep=2000)
        colcount[i,j] = np.mean(ssa.collisions)
        colcountzeros[i,j] = np.sum(ssa.collisions == 0)
        colcountnonzero[i,j] = np.mean(ssa.collisions[np.where(ssa.collisions !=0)])
        print(np.mean(ssa.collisions))
        
