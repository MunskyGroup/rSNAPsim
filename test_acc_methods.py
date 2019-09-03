# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 10:12:25 2019

@author: wsraymon
"""

#test all autocorrelation methods


import rSNAPsim as rss
import numpy as np
import matplotlib.pyplot as plt

rss =rSNAPsim.rSNAPsim()

rss.open_seq_file('gene_files/KDM5B.txt')
rss.run_default()
ssa_traj = rss.ssa_solver(n_traj=10)

ivec = ssa_traj.intensity_vec
ug = np.mean(ivec)
stdg = np.std(ivec)
varg = np.var(ivec)
ntraj = ivec.shape[0]

autocorr_ui = np.zeros((ivec.shape))
for i in range(ivec.shape[0]):
    autocorr_ui[i,:] = rss.get_acc2(ivec[i]-np.mean(ivec[i]))
    
autocorr_ug = np.zeros((ivec.shape))
for i in range(ivec.shape[0]):
    autocorr_ug[i,:] = rss.get_acc2(ivec[i]-ug)
    
    

mean_autocorr_ug = np.mean(autocorr_ug.T, axis=1)
mean_autocorr_ui = np.mean(autocorr_ui.T, axis=1)

mean_autocorr_ug_norm = np.mean(autocorr_ug.T/varg, axis=1)


autocorr_ui_norm = np.zeros((ivec.shape))
for i in range(ivec.shape[0]):
    autocorr_ui_norm[i,:] = rss.get_acc2(ivec[i]-np.mean(ivec[i])) 
    autocorr_ui_norm[i,:] = autocorr_ui_norm[i,:]/np.var(ivec[i,:])
    
mean_autocorr_ui_norm = np.mean(autocorr_ui_norm.T, axis=1)

sem_autocorr_ui_norm = 1.0/np.sqrt(ntraj)*np.std(autocorr_ui_norm.T,ddof=1,axis=1)
sem_autocorr_ug_norm = 1.0/np.sqrt(ntraj)*np.std(autocorr_ug.T/varg,ddof=1,axis=1)

sem_autocorr_ui = 1.0/np.sqrt(ntraj)*np.std(autocorr_ui.T,ddof=1,axis=1)
sem_autocorr_ug = 1.0/np.sqrt(ntraj)*np.std(autocorr_ug.T,ddof=1,axis=1)


plt.plot(mean_autocorr_ug,'g')
plt.plot(mean_autocorr_ui,'m')

plt.plot(mean_autocorr_ug-sem_autocorr_ug,'g--')
plt.plot(mean_autocorr_ug+sem_autocorr_ug,'g--')

plt.plot(mean_autocorr_ui-sem_autocorr_ui,'m--')
plt.plot(mean_autocorr_ui+sem_autocorr_ui,'m--')

plt.plot([0,1000],[0,0],'r--')
plt.legend(['ug','ui'])
    
plt.figure()
plt.plot(mean_autocorr_ug_norm,'g')
plt.plot(autocorr_ui_norm.T,'m')
plt.plot(mean_autocorr_ug_norm-sem_autocorr_ug_norm,'g--')
plt.plot(mean_autocorr_ug_norm+sem_autocorr_ug_norm,'g--')

plt.plot(mean_autocorr_ui_norm-sem_autocorr_ui_norm,'m--')
plt.plot(mean_autocorr_ui_norm+sem_autocorr_ui_norm,'m--')
plt.plot([0,1000],[0,0],'r--')
plt.legend(['varg','vari'])