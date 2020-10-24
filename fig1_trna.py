# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 13:33:03 2020

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
import matplotlib.pyplot as plt
import time
from rss import IntensityAnalyses as ia





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


k_index = np.array(poi.ktrna_id)

codons = [cdict().trna_ids[x] for x in poi.ktrna_id]
bactin_mean_trna_elongation = poi.total_length/np.sum(1/ np.array([cdict().strGeneCopy[y] for y in codons]))

poi.ke_mu = 14.00


t = np.linspace(0,2000,2001)
wt_solver = tss()  #solver class
wt_solver.protein=poi
wt_solver.default_conditions['burnin'] = 1000

ssa_soln_wt = wt_solver.solve_ssa([.033] + poi.kelong + [10],t,n_traj=50, low_memory=False, record_stats=True)

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
k_diff = 3
k_trna = np.array([strGeneCopy[cdict().trna_ids[x]] for x in range(0,61)])*k_diff
k_bind = .033
kelong = 1
k_compl = 10


t = np.linspace(0,3000,3001)
#ssa_soln = solution_obj()
#ssa_soln.load('trna_kdiff3_50traj')

ssa_soln3 = solver.solve_ssa_trna(k_index,  k_diff,k_bind,kelong,k_compl ,t,n_traj=50, k_trna=k_trna)


k_diff = 1
k_trna = np.array([strGeneCopy[cdict().trna_ids[x]] for x in range(0,61)])*k_diff
k_bind = .033
kelong = 1
k_compl = 10


t = np.linspace(0,3000,3001)
#ssa_soln = solution_obj()
#ssa_soln.load('trna_kdiff3_50traj')
ssa_soln1 = solver.solve_ssa_trna(k_index,  k_diff,k_bind,kelong,k_compl ,t,n_traj=50, k_trna=k_trna)

k_diff = .5
k_trna = np.array([strGeneCopy[cdict().trna_ids[x]] for x in range(0,61)])*k_diff
k_bind = .033
kelong = 1
k_compl = 10


t = np.linspace(0,3000,3001)
#ssa_soln = solution_obj()
#ssa_soln.load('trna_kdiff3_50traj')
ssa_solnpt5 = solver.solve_ssa_trna(k_index,  k_diff,k_bind,kelong,k_compl ,t,n_traj=50, k_trna=k_trna)



k_diff = .1
k_trna = np.array([strGeneCopy[cdict().trna_ids[x]] for x in range(0,61)])*k_diff
k_bind = .033
kelong = 1
k_compl = 10


t = np.linspace(0,3000,3001)
#ssa_soln = solution_obj()
#ssa_soln.load('trna_kdiff3_50traj')
ssa_solnpt1 = solver.solve_ssa_trna(k_index,  k_diff,k_bind,kelong,k_compl ,t,n_traj=50, k_trna=k_trna)


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
##########################






plt.plot(ssa_soln3.solutions[0].T,'#1cfff7')


plt.figure(dpi=300)
ft, fa = fss().get_fragments(ssa_soln3.solutions[0])
for i in range(len(fa)): 
    
    frag = fa[i][fa[i] > 0]
    
    timeseg = t[ft[i]: ft[i] + len(frag) ]
    
    
    plt.plot(frag ,timeseg[::-1],'seagreen',alpha=.8  )
    


ft, fa = fss().get_fragments(ssa_soln1.solutions[0])
for i in range(len(fa)): 
    
    frag = fa[i][fa[i] > 0]
    
    timeseg = t[ft[i]: ft[i] + len(frag) ]
    
    
    plt.plot(frag ,timeseg[::-1],'violet',alpha=.8  )
    
plt.xlabel('position')
plt.ylabel('time')




ft, fa = fss().get_fragments(ssa_solnpt5.solutions[0])
for i in range(len(fa)): 
    
    frag = fa[i][fa[i] > 0]
    
    timeseg = t[ft[i]: ft[i] + len(frag) ]
    
    
    plt.plot(frag ,timeseg[::-1],'royalblue',alpha=.8  )
    
plt.xlabel('position')
plt.ylabel('time')



ft, fa = fss().get_fragments(ssa_solnpt1.solutions[0])
for i in range(len(fa)): 
    
    frag = fa[i][fa[i] > 0]
    
    timeseg = t[ft[i]: ft[i] + len(frag) ]
    
    
    plt.plot(frag ,timeseg[::-1],'crimson',alpha=.8 )
    
plt.xlabel('position')
plt.ylabel('time')
plt.ylim([2500,3000])






#######################################

plt.figure()
plt.hist(ssa_soln3.intensity_vec[0].flatten(),alpha=1,histtype='step',density=True, lw=3,color='seagreen')
plt.hist(ssa_soln3.intensity_vec[0].flatten(),alpha=.1,density=True, color='seagreen')

plt.hist(ssa_soln_wt.intensity_vec[0].flatten(),alpha=1,histtype='step',density=True, lw=3,color='violet')
plt.hist(ssa_soln_wt.intensity_vec[0].flatten(),alpha=.3,density=True, color='violet')

plt.hist(ssa_soln1.intensity_vec[0].flatten(),alpha=1,histtype='step',density=True, lw=3,color='royalblue')
plt.hist(ssa_soln1.intensity_vec[0].flatten(),alpha=.3,density=True, color='royalblue')

plt.hist(ssa_solnpt5.intensity_vec[0].flatten(),alpha=1,histtype='step',density=True, lw=3,color='crimson')
plt.hist(ssa_solnpt5.intensity_vec[0].flatten(),alpha=.3,density=True, color='crimson')

plt.hist(ssa_solnpt1.intensity_vec[0].flatten(),alpha=1,histtype='step',density=True, lw=3,color='indigo')
plt.hist(ssa_solnpt1.intensity_vec[0].flatten(), alpha=.3,density=True, color='indigo')


#######################################





vmaxv = 45

a = ssa_soln3.all_trna_results[0,:].reshape((61,len(t)))
plt.rcParams.update({'font.size': 8, 'font.weight':'bold','font.family':'normal'  }   )

trna_counts = np.mean(np.mean(ssa_soln3.all_trna_results.reshape((50,61,len(t))),axis=2),axis=0)
t2 = np.append(trna_counts,[0,0,0]).reshape(8,8)
fig = plt.figure()
ax = fig.add_subplot(111)

size = 8

x_start = 0.0
x_end = 8.0
y_start = 0.0
y_end = 8.0

values =  np.round(np.fliplr(np.flip(t2)),1)
extent = [x_start, x_end, y_start, y_end]
im = ax.imshow(t2,extent=extent,vmax=vmaxv)

jump_x = (x_end - x_start) / (2.0 * size)
jump_y = (y_end - y_start) / (2.0 * size)
x_positions = np.linspace(start=x_start, stop=x_end, num=size, endpoint=False)
y_positions = np.linspace(start=y_start, stop=y_end, num=size, endpoint=False)

for y_index, y in enumerate(y_positions):
    for x_index, x in enumerate(x_positions):
        label = values[y_index, x_index]
        if label == 0.0:
            label = '-'
        text_x = x + jump_x
        text_y = y + jump_y
        ax.text(text_x, text_y, label, color='black', ha='center', va='center')




fig.colorbar(im)
plt.show()




a = ssa_soln3.all_trna_results[0,:].reshape((61,len(t)))
plt.rcParams.update({'font.size': 8, 'font.weight':'bold','font.family':'normal'  }   )

trna_counts = np.mean(np.mean(ssa_soln1.all_trna_results.reshape((50,61,len(t))),axis=2),axis=0)
t2 = np.append(trna_counts,[0,0,0]).reshape(8,8)
fig = plt.figure()
ax = fig.add_subplot(111)

size = 8

x_start = 0.0
x_end = 8.0
y_start = 0.0
y_end = 8.0

values =  np.round(np.fliplr(np.flip(t2)),1)
extent = [x_start, x_end, y_start, y_end]
im = ax.imshow(t2,extent=extent,vmax=vmaxv)

jump_x = (x_end - x_start) / (2.0 * size)
jump_y = (y_end - y_start) / (2.0 * size)
x_positions = np.linspace(start=x_start, stop=x_end, num=size, endpoint=False)
y_positions = np.linspace(start=y_start, stop=y_end, num=size, endpoint=False)

for y_index, y in enumerate(y_positions):
    for x_index, x in enumerate(x_positions):
        label = values[y_index, x_index]
        if label == 0.0:
            label = '-'
        text_x = x + jump_x
        text_y = y + jump_y
        ax.text(text_x, text_y, label, color='black', ha='center', va='center')




fig.colorbar(im)
plt.show()


a = ssa_soln3.all_trna_results[0,:].reshape((61,len(t)))
plt.rcParams.update({'font.size': 8, 'font.weight':'bold','font.family':'normal'  }   )

trna_counts = np.mean(np.mean(ssa_solnpt5.all_trna_results.reshape((50,61,len(t))),axis=2),axis=0)
t2 = np.append(trna_counts,[0,0,0]).reshape(8,8)
fig = plt.figure()
ax = fig.add_subplot(111)

size = 8

x_start = 0.0
x_end = 8.0
y_start = 0.0
y_end = 8.0

values =  np.round(np.fliplr(np.flip(t2)),1)
extent = [x_start, x_end, y_start, y_end]
im = ax.imshow(t2,extent=extent,vmax=vmaxv)

jump_x = (x_end - x_start) / (2.0 * size)
jump_y = (y_end - y_start) / (2.0 * size)
x_positions = np.linspace(start=x_start, stop=x_end, num=size, endpoint=False)
y_positions = np.linspace(start=y_start, stop=y_end, num=size, endpoint=False)

for y_index, y in enumerate(y_positions):
    for x_index, x in enumerate(x_positions):
        label = values[y_index, x_index]
        if label == 0.0:
            label = '-'
        text_x = x + jump_x
        text_y = y + jump_y
        ax.text(text_x, text_y, label, color='black', ha='center', va='center')




fig.colorbar(im)
plt.show()


a = ssa_soln3.all_trna_results[0,:].reshape((61,len(t)))
plt.rcParams.update({'font.size': 8, 'font.weight':'bold','font.family':'normal'  }   )

trna_counts = np.mean(np.mean(ssa_solnpt1.all_trna_results.reshape((50,61,len(t))),axis=2),axis=0)
t2 = np.append(trna_counts,[0,0,0]).reshape(8,8)
fig = plt.figure()
ax = fig.add_subplot(111)

size = 8

x_start = 0.0
x_end = 8.0
y_start = 0.0
y_end = 8.0

values =  np.round(np.fliplr(np.flip(t2)),1)
extent = [x_start, x_end, y_start, y_end]
im = ax.imshow(t2,extent=extent,vmax=vmaxv)

jump_x = (x_end - x_start) / (2.0 * size)
jump_y = (y_end - y_start) / (2.0 * size)
x_positions = np.linspace(start=x_start, stop=x_end, num=size, endpoint=False)
y_positions = np.linspace(start=y_start, stop=y_end, num=size, endpoint=False)

for y_index, y in enumerate(y_positions):
    for x_index, x in enumerate(x_positions):
        label = values[y_index, x_index]
        if label == 0.0:
            label = '-'
        text_x = x + jump_x
        text_y = y + jump_y
        ax.text(text_x, text_y, label, color='black', ha='center', va='center')




fig.colorbar(im)
plt.show()



#################################################


plt.figure()
acov,err_acov = ia().get_autocov(ssa_soln3.intensity_vec[:,1000:,:],norm='global')
acc,err_acorr = ia().get_autocorr(acov)



plt.plot(np.mean(acc[0],axis=1),color='seagreen');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='seagreen', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='seagreen', label='_nolegend_',alpha=.3)
plt.plot([0,500],[0,0],'r--', label='_nolegend_')
plt.xlim([0,100])

acov,err_acov = ia().get_autocov(ssa_soln_wt.intensity_vec,norm='global')
acc,err_acorr = ia().get_autocorr(acov)



plt.plot(np.mean(acc[0],axis=1),color='violet');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='violet', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='violet', label='_nolegend_',alpha=.3)
plt.plot([0,500],[0,0],'r--', label='_nolegend_')
plt.xlim([0,100])


acov,err_acov = ia().get_autocov(ssa_soln1.intensity_vec[:,1000:,:],norm='global')
acc,err_acorr = ia().get_autocorr(acov)



plt.plot(np.mean(acc[0],axis=1),color='royalblue');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='royalblue', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='royalblue', label='_nolegend_',alpha=.3)
plt.plot([0,500],[0,0],'r--', label='_nolegend_')
plt.xlim([0,100])



acov,err_acov = ia().get_autocov(ssa_solnpt5.intensity_vec[:,1000:,:],norm='global')
acc,err_acorr = ia().get_autocorr(acov)



plt.plot(np.mean(acc[0],axis=1),color='crimson');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='crimson', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='crimson', label='_nolegend_',alpha=.3)
plt.plot([0,500],[0,0],'r--', label='_nolegend_')
plt.xlim([0,100])

acov,err_acov = ia().get_autocov(ssa_solnpt1.intensity_vec[:,1000:,:],norm='global')
acc,err_acorr = ia().get_autocorr(acov)



plt.plot(np.mean(acc[0],axis=1),color='indigo');
plt.plot(np.mean(acc[0],axis=1) - err_acorr[0],'--',color='indigo', label='_nolegend_',alpha=.3)
plt.plot(np.mean(acc[0],axis=1)+ err_acorr[0],'--',color='indigo', label='_nolegend_',alpha=.3)
plt.plot([0,500],[0,0],'r--', label='_nolegend_')
plt.xlim([0,100])

plt.ylabel('Autocorrelation')
plt.xlabel('$tau$')
plt.legend(['kdiff 3','Auguilera2019','kdiff 1','kdiff .5','kdiff .1'])