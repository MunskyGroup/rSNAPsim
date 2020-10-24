# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 11:19:42 2020

@author: willi
"""

import numpy as np
import os



from rss import rSNAPsim
from rss import ProbeVectorFactory as pvf
from rss import PropensityFactory as pff
from rss import TranslationSolvers as tss
from rss import FragmentSeperator as fss
from rss import CodonDictionaries as cdict
from rss import SSA_Soln as solution_obj
from rss import CodonOptimizer as copt
from rss import SequenceManipMethods as smm
import matplotlib.pyplot as plt
import time
from rss import IntensityAnalyses as ia
from matplotlib import cm

from cycler import cycler

#############################
# Global setup

colors = ['#5e1150', '#d90000', '#167f39','#ff8c00','#04756f']

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

plt.rcParams.update({'font.size': 12, 'font.weight':'bold','font.family':'normal'  }   )
plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})

plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})
plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})


plt.rcParams.update({'xtick.major.width'   : 2.8 })
plt.rcParams.update({'xtick.labelsize'   : 14 })



plt.rcParams.update({'ytick.major.width'   : 2.8 })
plt.rcParams.update({'ytick.labelsize'   : 14})


plt.rcParams.update({'axes.linewidth':2.8})
plt.rcParams.update({'axes.labelpad':8})
plt.rcParams.update({'axes.titlepad':10})
plt.rcParams.update({'figure.dpi':300})


#def figure_1():
'''
Comparison of tRNA and old model
'''

rsim = rSNAPsim()
rsim.open_seq_file('./gene_files/Bactin_withTags.txt')
poi = rsim.proteins['1'][0]  #protein object
solver = tss()  #solver class
solver.protein=poi


Bactin_opt = copt().optimize_ntseq(poi.nt_seq)
Bactin_deopt = copt().deoptimize_ntseq(poi.nt_seq)

ktrna_id_opt = pff().get_trna_ids(Bactin_opt)
ktrna_id_deopt = pff().get_trna_ids(Bactin_deopt)

colors = cm.viridis(np.linspace(0,1,60))

k_index = np.array(poi.ktrna_id)

(tid,counts) = np.unique(poi.ktrna_id, return_counts=True)
fig,ax = plt.subplots(1,2)
ax[0].pie(counts,colors = colors,explode=.05*np.ones(len(counts)))
#ax[0].set_title('WildType trna ID distribution')

ax[1].hist(poi.ktrna_id,density=True,bins=61)
ax[1].set_ylim([0,.13])
fig.tight_layout()


(tid,counts) = np.unique(ktrna_id_opt, return_counts=True)
fig,ax = plt.subplots(1,2)
ax[0].pie(counts,colors = colors[tid],explode=.05*np.ones(len(counts)))
#ax[0].set_title('Optimized trna ID distribution')

ax[1].hist(ktrna_id_opt,density=True,bins=61)
ax[1].set_ylim([0,.13])
fig.tight_layout()


(tid,counts) = np.unique(ktrna_id_deopt, return_counts=True)
fig,ax = plt.subplots(1,2)
ax[0].pie(counts,colors = colors[tid],explode=.05*np.ones(len(counts)))
#ax[0].set_title('Deoptimized trna ID distribution')

ax[1].hist(ktrna_id_deopt,density=True,bins=61)
ax[1].set_ylim([0,.13])
fig.tight_layout()



bactin = rsim.proteins['1'][0]  #protein object

rsim2 = rSNAPsim()
rsim2.open_seq_file('./gene_files/Bactin_withTags.txt')
bactin_deopt = rsim2.proteins['1'][0]  #protein object

rsim3 = rSNAPsim()
rsim3.open_seq_file('./gene_files/Bactin_withTags.txt')
bactin_opt = rsim3.proteins['1'][0]  #protein object

bactin_opt.nt_seq = Bactin_opt
bactin_opt.aa_seq = smm('').nt2aa(Bactin_opt)

bactin_deopt.nt_seq = Bactin_deopt
bactin_deopt.aa_seq = smm('').nt2aa(Bactin_deopt)





wt_solver = tss()  #solver class
wt_solver.protein=bactin
wt_solver.default_conditions['burnin'] = 1000

opt_solver = tss()  #solver class
opt_solver.protein=bactin_opt
wt_solver.default_conditions['burnin'] = 1000

deopt_solver = tss()  #solver class
deopt_solver.protein=bactin_deopt
wt_solver.default_conditions['burnin'] = 1000

bactin_mean = 14.0

t = np.linspace(0,2000,2001)
bactin.ke_mu = bactin_mean
bactin_opt.ke_mu = bactin_mean
bactin_deopt.ke_mu = bactin_mean
ssa_soln_wt = wt_solver.solve_ssa([.033] + bactin.kelong + [10],t,n_traj=50, low_memory=False, record_stats=True)
ssa_soln_opt = opt_solver.solve_ssa([.033] + bactin_opt.kelong + [10],t,n_traj=50)
ssa_soln_deopt = deopt_solver.solve_ssa([.033] + bactin_deopt.kelong + [10],t,n_traj=50)





acov,err_acov = ia().get_autocov(ssa_soln_wt.intensity_vec[:,1000:,:],norm='global')
acc,err_acorr = ia().get_autocorr(acov)

plt.figure()

plt.plot(np.mean(acc[0],axis=1),color='seagreen');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='seagreen', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='seagreen', label='_nolegend_',alpha=.3)
plt.plot([0,500],[0,0],'r--', label='_nolegend_')
plt.xlim([0,500])


acov,err_acov = ia().get_autocov(ssa_soln_opt.intensity_vec[:,1000:,:],norm='global')
acc,err_acorr = ia().get_autocorr(acov)
plt.plot(np.mean(acc[0],axis=1),color='indigo');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='indigo', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='indigo', label='_nolegend_',alpha=.3)

plt.xlim([0,500])



acov,err_acov = ia().get_autocov(ssa_soln_deopt.intensity_vec[:,1000:,:],norm='global')
acc,err_acorr = ia().get_autocorr(acov)

plt.plot(np.mean(acc[0],axis=1),color='violet');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='violet', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='violet', label='_nolegend_',alpha=.3)

plt.xlim([0,200])

plt.legend(['wt','opt','deopt',])



strGeneCopy = {'TTT': 17.6, 'TCT': 15.2, 'TAT': 12.2, 'TGT': 10.6, 'TTC': 20.3,
                            'TCC': 17.7, 'TAC': 15.3, 'TGC': 12.6, 'TTA': 7.7, 'TCA': 12.2,
                            'TAA': 1.0, 'TGA': 1.6, 'TTG': 12.9, 'TCG':  4.4, 'TAG': 0.8,
                            'TGG': 13.2, 'CTT': 13.2, 'CCT': 17.5, 'CAT': 10.9, 'CGT': 4.5,
                            'CTC': 19.6, 'CCC': 19.8, 'CAC': 15.1, 'CGC': 10.4, 'CTA':  7.2,
                            'CCA': 16.9, 'CAA': 12.3, 'CGA':  6.2, 'CTG': 39.6, 'CCG':  6.9,
                            'CAG': 34.2, 'CGG': 11.4, 'ATT': 16.0, 'ACT': 13.1, 'AAT': 17.0,
                            'AGT': 12.1, 'ATC': 20.8, 'ACC': 18.9, 'AAC': 19.1, 'AGC': 19.5,
                            'ATA':  7.5, 'ACA': 15.1, 'AAA': 24.4, 'AGA': 12.2, 'ATG': 22.0,
                            'ACG': 6.1, 'AAG': 31.9, 'AGG': 12.0, 'GTT': 11.0, 'GCT': 18.4,
                            'GAT': 21.8, 'GGT': 10.8, 'GTC': 14.5, 'GCC': 27.7, 'GAC': 25.1,
                            'GGC': 22.2, 'GTA':  7.1, 'GCA': 15.8, 'GAA': 29.0, 'GGA': 16.5,
                            'GTG': 28.1, 'GCG': 7.4, 'GAG': 39.6, 'GGG': 16.5}


strGeneCopy.pop('TAG')
strGeneCopy.pop('TAA')
strGeneCopy.pop('TGA')

k_trna = np.array([strGeneCopy[cdict().trna_ids[x]] for x in range(0,61)]) 
k_trna = np.array(list(strGeneCopy.values()))*10

t = np.linspace(0,3000,3001)
k_diffusion = .1
k_trna = np.array([strGeneCopy[cdict().trna_ids[x]] for x in range(0,61)])*k_diffusion
k_bind = .033
elong_scale = 1
kcompl = 10


trna_ssa_soln_wt = wt_solver.solve_ssa_trna( np.array(bactin.ktrna_id) , k_diffusion,k_bind,elong_scale,kcompl ,t,n_traj=10,k_trna= k_trna)

trna_ssa_soln_opt = opt_solver.solve_ssa_trna(np.array(bactin_opt.ktrna_id), k_diffusion,k_bind,elong_scale,kcompl ,t,n_traj=10,k_trna= k_trna)

trna_ssa_soln_deopt = deopt_solver.solve_ssa_trna(np.array(bactin_deopt.ktrna_id), k_diffusion,k_bind,elong_scale,kcompl ,t,n_traj=10,k_trna= k_trna)




plt.gca()
acov,err_acov = ia().get_autocov(trna_ssa_soln_wt.intensity_vec[:,1000:,:],norm='global')
acc,err_acorr = ia().get_autocorr(acov)



plt.plot(np.mean(acc[0],axis=1),color='royalblue');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='royalblue', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='royalblue', label='_nolegend_',alpha=.3)
plt.plot([0,500],[0,0],'r--', label='_nolegend_')
plt.xlim([0,500])


acov,err_acov = ia().get_autocov(trna_ssa_soln_opt.intensity_vec[:,1000:,:],norm='global')
acc,err_acorr = ia().get_autocorr(acov)
plt.plot(np.mean(acc[0],axis=1),color='crimson');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='crimson', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='crimson', label='_nolegend_',alpha=.3)

plt.xlim([0,500])



acov,err_acov = ia().get_autocov(trna_ssa_soln_deopt.intensity_vec[:,1000:,:],norm='global')
acc,err_acorr = ia().get_autocorr(acov)

plt.plot(np.mean(acc[0],axis=1),color='orange');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='orange', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='orange', label='_nolegend_',alpha=.3)

plt.xlim([0,200])

#plt.legend(['wt','opt','deopt',])
plt.title('tRNA solutions')
