# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:09:50 2020

@author: willi
"""

#get all CAIs

import os

import sys
sys.path.append('C:/Users/willi/Documents/GitHub/rSNAPsim')

os.chdir('..')
from rss import SequenceManipMethods as smm
from rss import TranslationSolvers as tss
from rss import rSNAPsim as rss
from rss import CodonOptimizer as copt
from rss import IntensityAnalyses as ia

os.chdir('./trna_ssa')


import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

sequences = []
with open('CCDS_nucleotide.current.fna', 'r') as f:
    for record in SeqIO.parse(f, "fasta"):
        sequences.append(record)
        
seq_manip = smm('')

print('getting CAI....')
CAI_ccds = []
cai_sequences = []
for i in range(len(sequences)):
    try:
        CAI = seq_manip.codon_usage(str(sequences[i].seq))[1]
        cai_sequences.append(str(sequences[i].seq))
        CAI_ccds.append( CAI)
    except:
        pass
    
CAI_ccds = np.array(CAI_ccds)
#CAI_ccds = [ seq_manip.codon_usage(str(x.seq))[1] for x in sequences]


sorted_order = np.argsort(CAI_ccds)

valid_cai = np.where( CAI_ccds[np.argsort(CAI_ccds)] >0)[0] #CAI above 0

tRNA_ai_bp_table = {'IU':0, 'IC': .28, 'IA': .9999, 'GU':.41,'GC':0,'UG':.68,'UA':0,'CG':0,'LA':.89 }






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

k_trna = np.array(list(strGeneCopy.values()))

t = np.linspace(0,1400,2801)
k_diffusion = 1
k_bind = .033
elong_scale = 1
optimizer = copt()
rss = rss()
solver = tss()


print('running simulations...')

import time

st = time.time()

wt_taus = []
opt_taus = []
cais = []
sliceval = 10
checked_seqs = []
total_seq = len(sorted_order[valid_cai][::sliceval] )


for i in range(0,total_seq):
    
    print('-----------------')
    print('current seq: ')
    print('%i out of %i' %(i,total_seq ))
    
    try:
        seq = cai_sequences[   sorted_order[valid_cai][::sliceval][i]  ]
        
        optimized_seq = optimizer.optimize_ntseq(seq)
        
        poi = rss.seq_to_protein_obj(seq)['1'][0]
        print(seq[:40])
        
        
        predicted_tau = int(np.ceil(len(seq)/10))
        
        
        t = np.linspace(0,predicted_tau*5 ,2*(predicted_tau*5)+1)
        
        optimized_poi = rss.seq_to_protein_obj(optimized_seq)['1'][0]
        
        solver.protein  = poi
    
        trna_ssa_soln_wt = solver.solve_ssa_trna( np.array(poi.ktrna_id) , k_diffusion,k_bind,elong_scale ,t,n_traj=10,k_trna= k_trna)
        
        solver.protein = optimized_poi
        trna_ssa_soln_opt = solver.solve_ssa_trna(np.array(optimized_poi.ktrna_id), k_diffusion,k_bind,elong_scale ,t,n_traj=10,k_trna= k_trna)
        
        print('time for simulation: ')
        print(time.time()-st)
        
        acov,err_acov = ia().get_autocov(trna_ssa_soln_wt.intensity_vec[:,predicted_tau*2:,:],norm='global')
        acc_wt,err_acorr_wt = ia().get_autocorr(acov)
        
        acov,err_acov = ia().get_autocov(trna_ssa_soln_opt.intensity_vec[:,predicted_tau*2:,:],norm='global')
        acc_opt,err_acorr_opt = ia().get_autocorr(acov)
        
        tau_acc_opt = t[np.where(np.mean(acc_opt[0],axis=1) < 0.05)[0][0]]
        tau_acc_wt = t[np.where(np.mean(acc_wt[0],axis=1) < 0.05)[0][0]]
        cais.append(CAI_ccds[   sorted_order[valid_cai][::sliceval][i]  ])
        
        wt_taus.append(tau_acc_wt)
        opt_taus.append(tau_acc_opt)
        
        checked_seqs.append(seq)
    except:
        pass
    
    
    
    
np.save('cais',np.array(cais))
np.save('wt_taus',np.array(wt_taus))
np.save('opt_taus',np.array(wt_taus))


with open('sequences_cai.txt', 'w') as f:
    for item in checked_seqs:
        f.write("%s\n" % item)


fig,ax = plt.subplots(1,dpi=300) #.figure()
ax.scatter(cais,(np.array(wt_taus) - np.array(opt_taus) )/np.array(wt_taus)   ,s=1)
ax.plot([.6,1],[0,0],'r--')
ax.set_xlabel('CAI')
ax.set_ylabel('(wt_tau - opt_tau) / wt_tau')
fig.show()