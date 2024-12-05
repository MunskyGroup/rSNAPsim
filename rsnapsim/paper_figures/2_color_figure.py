# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 20:01:27 2021

@author: willi
"""
import os
cwd = os.getcwd()
os.chdir('../../..')

import rsnapsim as rss
from rsnapsim import seqmanip

import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


os.chdir(cwd)

cb_friendly_8 = ['#000000', '#2271B2', '#3DB7E9', '#F748A5', '#359B73', '#D55E00', '#E69F00', 
                 '#F0E442']

cb_friendly_15 = ['#68023F', '#008169', '#EF0096', '#00DCB5', '#FFCFE2', '#003C86', '#9400E6',
                  '#009FFA', '#FF71FD', '#7CFFFA', '#6A0213', '#008607', '#F60239', '#00E307','#FFDC3D']

cb_15_order_fluorescent = [3,2,9,4,13,14]
cb_8_order_fluorescent = [2,3,6,7]


import matplotlib.pyplot as plt
from cycler import cycler
########################################
dark = True
if not dark:
    colors = ['#ef476f', '#073b4c','#06d6a0','#7400b8','#073b4c', '#118ab2', '#00A090',
              "00DCB5"]
    
    colors = [cb_friendly_15[x] for x in cb_15_order_fluorescent]
    
else:
    plt.style.use('dark_background')
    #colors = ['#118ab2','#57ffcd', '#ff479d', '#ffe869','#ff8c00','#04756f']
    colors = ['#ef476f', '#073b4c','#06d6a0','#7400b8','#073b4c', '#118ab2', '#00A090',
              "00DCB5"]
    
    colors = [cb_friendly_15[x] for x in cb_15_order_fluorescent]
    
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

save = True

plt.rcParams.update({'font.size': 12, 'font.weight':'bold','font.family':'normal'  }   )
plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})

plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})
plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})


plt.rcParams.update({'xtick.major.width'   : 2.8 })
plt.rcParams.update({'xtick.labelsize'   : 12 })



plt.rcParams.update({'ytick.major.width'   : 2.8 })
plt.rcParams.update({'ytick.labelsize'   : 12})

plt.rcParams.update({'axes.titleweight'   : 'bold'})
plt.rcParams.update({'axes.titlesize'   : 10})
plt.rcParams.update({'axes.labelweight'   : 'bold'})
plt.rcParams.update({'axes.labelsize'   : 12})

plt.rcParams.update({'axes.linewidth':2.8})
plt.rcParams.update({'axes.labelpad':8})
plt.rcParams.update({'axes.titlepad':10})
plt.rcParams.update({'figure.dpi':300})

##########################################


cmap = plt.get_cmap('Spectral')

#my_cmap = cmap(np.arange(cmap.N))
#my_cmap[:, -1] = .9
#my_cmap = ListedColormap(my_cmap)

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

new_cmap = truncate_colormap(cmap, 0.1, 1)
my_cmap = new_cmap(np.arange(cmap.N))
my_cmap[:, -1] = .9
my_cmap = ListedColormap(my_cmap)


suntag_kif = '''atgataataGCTAGCTCTGGAGGAGAAGAACTTTTGAGCAAGAATTATCATCTTGAGAACGAAGTGGCTCGTCTTAAGAAAGGTTCTGGCAGTGGAGAAGAACTGCTTTCAAAGAATTACCACCTGGAAAATGAGGTAGCTAGACTGAAAAAGGGGAGCGGAAGTGGGGAGGAGTTGCTGAGCAAAAATTATCATTTGGAGAACGAAGTAGCACGACTAAAGAAAGGGTCCGGATCGGGTGAGGAGTTACTCTCGAAAAATTATCATCTCGAAAACGAAGTGGCTCGGCTAAAAAAGGGCAGTGGTTCTGGAGAAGAGCTATTATCTAAAAACTACCACCTCGAAAATGAGGTGGCACGCTTAAAAAAGGGAAGTGGCAGTGGTGAAGAGCTACTATCCAAGAATTATCATCTTGAGAACGAGGTAGCGCGTTTGAAGAAGGGTTCCGGCTCAGGAGAGGAACTGCTCTCGAAGAACTATCATCTTGAAAATGAGGTCGCTCGATTAAAAAAGGGATCGGGCAGTGGTGAGGAACTACTTTCAAAGAATTACCACCTCGAAAACGAAGTAGCTCGATTAAAGAAAGGTTCAGGGTCGGGTGAAGAATTACTGAGTAAAAATTATCATCTGGAAAATGAGGTAGCGAGACTAAAAAAGGGGAGTGGTTCTGGCGAAGAGTTGCTATCGAAAAATTATCATCTTGAGAACGAAGTTGCTAGGCTCAAAAAGGGCTCAGGCTCAGGCGAGGAGTTGCTCTCGAAAAACTACCACTTGGAAAATGAGGTCGCGAGGTTGAAAAAGGGGAGCGGGTCGGGCGAGGAGTTATTGAGCAAAAACTATCATTTAGAGAACGAAGTCGCGCGCTTAAAGAAAGGCTCGGGCTCGGGCGAAGAACTCTTATCGAAGAACTACCACCTCGAAAATGAGGTCGCCAGGTTGAAAAAGGGCAGTGGCAGCGGGGAGGAACTCTTGAGCAAGAACTACCACTTGGAGAATGAGGTCGCGAGATTGAAGAAAGGGTCGGGGAGCGGCGAGGAATTGCTCAGCAAGAATTATCATTTGGAGAACGAAGTCGCCAGGCTCAAGAAAGGCTCGGGGTCGGGGGAGGAATTGTTGAGTAAAAACTACCACTTGGAAAATGAAGTCGCCAGGCTCAAAAAAGGGAGTGGGAGCGGCGAAGAGTTATTGAGCAAAAATTACCACTTGGAGAACGAAGTGGCAAGGCTCAAGAAAGGGAGCGGCAGCGGGGAGGAGCTCTTATCGAAGAACTACCACTTAGAGAATGAAGTCGCCCGCTTGAAGAAAGGCTCGGGGAGCGGGGAAGAGCTCTTGAGCAAGAACTACCACTTGGAAAATGAGGTGGCGCGCTTGAAGAAAGGGAGCGGGAGCGGGGAAGAGTTACTATCTAAGAATTATCATCTCGAGAACGAGGTGGCTCGACTAAAGAAGGGCTCCGGCAGTGGGGAGGAACTCCTGTCGAAGAACTATCATCTTGAAAATGAGGTTGCAAGACTTAAAAAGGGGTCCGGATCAGGTGAGGAACTACTCAGTAAGAATTACCACCTGGAAAACGAAGTTGCACGTTTGAAGAAAGGATCAGGATCAGGCGAAGAACTGCTCTCAAAAGATTATCATTTGGAAAATGAGGTTGCACGTTTAAAAAAGGGAAGTGGCAGTGGTGAGGAACTTCTGTCGAAAAATTATCATCTCGAGAATGAAGTAGCCCGACTTAAAAAGGGTTCTGGCTCGGGTCAGCGGCCGCACCGGTCAGCAGTGGAGGACAGCACGCTGCAAGTAGTGGTACGGGTGCGGCCCCCCACCCCTCGGGAGCTGGACAGTCAGCGGCGGCCAGTGGTTCAGGTGGTGGACGAGCGGGTGCTGGTGTTTAACCCTGAGGAGCCCGATGGAGGGTTCCCTGGCCTGAAATGGGGTGGCACCCATGATGGCCCCAAGAAGAAGGGCAAAGACCTGACGTTTGTCTTTGACCGGGTCTTTGGCGAGGCGGCCACCCAACAGGACGTGTTCCAGCACACCACGCACAGCGTCCTGGACAGCTTCCTCCAGGGCTACAACTGCTCAGTGTTTGCCTACGGGGCCACCGGGGCTGGGAAGACACACACCATGCTGGGAAGGGAGGGGGACCCCGGCATCATGTACCTGACCACCGTGGAACTGTACAGGCGCCTGGAGGCCCGCCAGCAGGAGAAGCACTTCGAGGTGCTCATCAGCTACCAGGAGGTGTATAATGAACAGATCCATGACCTCCTGGAGCCCAAGGGGCCCCTTGCCATCCGCGAGGACCCCGACAAGGGGGTGGTGGTGCAAGGACTTTCTTTCCACCAGCCAGCCTCAGCCGAGCAGCTGCTGGAGATACTGACCAGGGGGAACCGTAACCGCACGCAGCACCCCACTGATGCCAACGCGACTTCCTCCCGCTCCCATGCCATCTTCCAGATCTTTGTGAAGCAGCAGGACCGGGTTCCAGGACTGACCCAGGCTGTCCAGGTGGCCAAGATGAGCCTGATTGACCTGGCTGGCTCAGAGCGGGCATCCAGCACCCATGCGAAGGGGGAGCGGCTGCGGGAGGGGGCCAACATCAACCGCTCTCTGCTGGCGCTCATCAACGTCCTCAATGCCTTGGCCGATGCAAAGGGCCGCAAGACCGCTGTGCCCTACGCGGACAGCGCACTGACCCGCCTGCTCAAAGACTCCCTCGGGGGCAACTGCCGCACAGTGATGATCGCTGCCATCAGCCCCTCCAGCCTGACCTACGAGGACACGTATAATACCCTCAAATATGCCGACCGGGCCAAGGAGATCAGGCTCTCGCTGAAGAGCAATGTGACCAGCCTGGACTGTCACATCAGCCAGTATGCTACCATCTGCCAACAGCTCCAGGCTGAGGTAGCCGCTCTGAGGAAGAAGCTCCAAGTGTATGAGGGGGGAGGCCAGCCCCCACCACAGGACCTCCCAGGATCTCCCAAGTCGGGACCACCACCAGAACACCTTCCCAGCTCCCCCTTGCCACCCCACCCTCCCAGCCAGCCCTGCACCCCAGAGCTCCCTGCAGGGCCTAGAGCCCTTCAAGAGGAGAGTCTGGGGATGGAGGCCCAGGTGGAGAGGGCCATGGAAGGGAACTCTTCAGACCAGGAGCAGTCCCCAGAGGATGAGGATGAAGGCCCAGCTGAGGAGGTTCCAACCCAGATGCCAGAGCAGAACCCCACACATGCACTGCCAGAGTCCCCTCGCCTGACCCTGCAGCCCAAGCCAGTCGTGGGCCACTTCTCAGCACGGGAACTGGATGGGGACCGTTCTAAGCGGTTGGCCCTAAAGGTGCTGTGCGTTGCCCAGCGGCAGTACTCCCTGCTCCAAGCAGCCAACCTCCTGACGCCCGACATGATCACAGAGTTTGAGACCCTACAGCAGCTGGTGCAAGAGGAAAAAATTGAGCCTGGGGCAGAGGCCTTGAGGACTTCAGGCCTGGCCAGGGGGGCACCTCTGGCTCAGGAGCTGTGTTCAGAGTCAAAGCCTCCAGGATACACTGGCCCTGTGACCCGGACTATGGCGAGGCGACTGAGTGGCCCCCTGCACACCCTGGGAATCCCGCCTGGACCCAACTGCACCCCAGCCCAGGGGTCCCGATGGCCCATGGAGAAGAAGAGGAGGAGACCAAGCGCCTTGGAGGCAGACAGTCCCATGGCCCCAAAGCGGGGCACCAAGCGCCAGCGCCAGTCCTTCCTGCCCTGCCTAAGGAGAGGGTCTCTGCCTGACACCCAACCTTCACAGGGGCCCAGCACCCCCAAAGGAGAAAGGGCCTCCTCCCCCTGCCATTCCCCTCGCGTTTGCCCAGCCACAGTCATCAAAAGCCGGGTGCCCCTGGGCCCTTCCGCCATGCAGAACTGCTCCACCCCGCTGGCTCTGCCCACTCGAGACCTCAATGCCACCTTTGATCTCTCTGAGGAGCCTCCCTCAAAGCCCAGTTTCCATGAATGCATTGGCTGGGACAAAATACCCCAGGAGCTGAGCAGGCTGGACCAGCCCTTCATCCCCAGGGCACCTGTGCCCCTGTTCACCATGAAGGGCCCCAAGCCAACATCTTCCCTCCCTGGGACCTCTGCCTGCAAGAAGAAGCGCGTTGCGAGTTCCTCAGTCTCCCATGGCCGCAGCCGCATCGCCCGCCTCCCCAGCAGCACTTTGAAGAGGCCAGCTGGGCCCCTTGTACTCCCAGAGCTGCCCTTGAGTCCCCTGTGCCCTAGCAACCGGAGGAATGGAAAGGACCTCATCAGGGTGGGGAGAGCACTCTCAGCAGGGAACGGCGTCACCAAGGTGTCCGATAAGGACCTAGGCGGACTGTTACTGAGCTGCGTTTTACACCCTTTCTTTGACAAAACCTAA'''
bactin_seq = '''ATGGACTACAAGGACGACGACGACAAAGGTGACTACAAAGATGATGACGATAAAGGCGACTATAAGGACGATGACGACAAGGGCGGAAACTCACTGATCAAGGAAAACATGCGGATGAAGGTGGTGATGGAGGGCTCCGTGAATGGTCACCAGTTCAAGTGCACCGGAGAGGGAGAGGGAAACCCGTACATGGGAACTCAGACCATGCGCATTAAGGTCATCGAAGGAGGTCCGCTGCCGTTCGCTTTCGATATCCTGGCCACTTCGTTCGGAGGAGGGTCGCGCACGTTCATCAAGTACCCGAAGGGAATCCCGGACTTCTTTAAGCAGTCATTCCCGGAAGGATTCACTTGGGAACGGGTGACCCGGTATGAAGATGGAGGTGTGGTGACTGTCATGCAAGATACTTCGCTGGAGGATGGGTGCCTCGTGTACCACGTCCAAGTCCGCGGAGTGAATTTCCCGTCCAACGGACCAGTGATGCAGAAAAAGACGAAGGGTTGGGAACCTAATACTGAAATGATGTACCCCGCAGACGGAGGGCTGAGGGGCTACACCCACATGGCGCTGAAGGTCGACGGAGGAGATTACAAGGATGACGACGATAAGCAACAAGATTACAAAGACGATGATGACAAGGGCCAGCAGGGCGACTACAAGGACGACGACGACAAGCAGCAGGACTACAAAGATGACGATGATAAAGGAGGAGGACATCTGTCCTGTTCGTTCGTGACCACCTACAGATCAAAGAAAACCGTGGGAAACATCAAGATGCCGGGCATTCATGCCGTCGACCACCGCCTGGAGCGGCTCGAAGAATCAGACAATGAGATGTTCGTCGTGCAAAGAGAACATGCCGTGGCCAAGTTCGCGGGACTGGGAGGCGGTGGAGGCGATTACAAAGACGATGATGACAAGGGTGACTATAAAGACGACGATGACAAAGGGGATTACAAGGATGATGATGATAAGGGAGGCGGTGGATCAGGTGGAGGAGGTTCACTGCAGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCCCCAGGCACCAGGGCGTGATGGTGGGCATGGGTCAGAAGGATTCCTATGTGGGCGACGAGGCCCAGAGCAAGAGAGGCATCCTCACCCTGAAGTACCCCATCGAGCACGGCATCGTCACCAACTGGGACGACATGGAGAAAATCTGGCACCACACCTTCTACAATGAGCTGCGTGTGGCTCCCGAGGAGCACCCCGTGCTGCTGACCGAGGCCCCCCTGAACCCCAAGGCCAACCGCGAGAAGATGACCCAGATCATGTTTGAGACCTTCAACACCCCAGCCATGTACGTTGCTATCCAGGCTGTGCTATCCCTGTACGCCTCTGGCCGTACCACTGGCATCGTGATGGACTCCGGTGACGGGGTCACCCACACTGTGCCCATCTACGAGGGGTATGCCCTCCCCCATGCCATCCTGCGTCTGGACCTGGCTGGCCGGGACCTGACTGACTACCTCATGAAGATCCTCACCGAGCGCGGCTACAGCTTCACCACCACGGCCGAGCGGGAAATCGTGCGTGACATTAAGGAGAAGCTGTGCTACGTCGCCCTGGACTTCGAGCAAGAGATGGCCACGGCTGCTTCCAGCTCCTCCCTGGAGAAGAGCTACGAGCTGCCTGACGGCCAGGTCATCACCATTGGCAATGAGCGGTTCCGCTGCCCTGAGGCACTCTTCCAGCCTTCCTTCCTGGGCATGGAGTCCTGTGGCATCCACGAAACTACCTTCAACTCCATCATGAAGTGTGACGTGGACATCCGCAAAGACCTGTACGCCAACACAGTGCTGTCTGGCGGCACCACCATGTACCCTGGCATTGCCGACAGGATGCAGAAGGAGATCACTGCCCTGGCACCCAGCACAATGAAGATCAAGATCATTGCTCCTCCTGAGCGCAAGTACTCCGTGTGGATCGGCGGCTCCATCCTGGCCTCGCTGTCCACCTTCCAGCAGATGTGGATCAGCAAGCAGGAGTATGACGAGTCCGGCCCCTCCATCGTCCACCGCAAATGCTTCTAG'''
#the original suntag_kif has 24x Epitopes, so we are going to manually delete 14 of them with this sequence:
suntags_todelete = 'GSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKDYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGQ'

hairpin = 'ACGUGCCACGAUUCAACGUGGCACAGC'
bactin_only = bactin_seq[337*3:] #getting just B-actin
suntag = suntag_kif[:591*3] #deleting KIF
suntag = suntag[:729] + suntag[1737:] #deleting 14x Suntag epitopes
flagtag = bactin_seq[:337*3] #getting just the Flag tag 

bactin_codons = [bactin_only[i:i+3] for i in range(0, len(bactin_only), 3)]  #convert b-actin to pairs of codons
suntag_codons = [suntag[i:i+3] for i in range(0, len(suntag), 3)] #suntag codons
flagtag_codons = [flagtag[i:i+3] for i in range(0, len(flagtag), 3)] #flagtag codons
hairpin_codons = [hairpin[i:i+3] for i in range(0, len(hairpin), 3)] #hairpin codons

construct1 = flagtag_codons + suntag_codons + bactin_codons 
construct_sequence_str1 = ''.join(construct1)

construct2 = flagtag_codons + bactin_codons[:-1] + suntag_codons + ['TGA',]
construct_sequence_str2 = ''.join(construct2)


bactin_1 = rss.seqmanip.seq_to_protein_obj(construct_sequence_str1,add_tag=False)['0'][0]
bactin_2 = rss.seqmanip.seq_to_protein_obj(construct_sequence_str2,add_tag=False)['0'][0]



def movmean(a, n=3): #rolling average
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

f,ax = plt.subplots(1,1,tight_layout=True)
ax.plot(movmean(bactin_1.kelong,n=10))
ax.set_xlabel('codon')
ax.set_ylabel('stepping rate')
if save:
    plt.savefig('kelong_bactin1.svg')

bin_means = [np.mean(bactin_1.kelong[x*10:x*10+10]) for x in range(0,int(967/10))]
f,ax = plt.subplots(1,1,tight_layout=True)
a = plt.imshow(np.array([bin_means]), aspect=10, cmap=my_cmap)
ax.set_xlabel('codons (binned by 10)')
ax.set_title('Stepping rate calculation')
plt.colorbar(a, ax=ax, orientation='horizontal')
plt.savefig('kelong_bactin1.svg')



f,ax = plt.subplots(1,1,tight_layout=True)
a = plt.imshow(np.array([movmean(bactin_1.kelong,10)]), aspect=100, cmap=my_cmap)
ax.set_xlabel('codons')
ax.set_title('Stepping rate calculation')
plt.colorbar(a, ax=ax, orientation='horizontal')
plt.savefig('kelong_bactin2.svg')




viridis = cm.get_cmap('viridis', 100)
newcolors = viridis(np.linspace(0, 1, 100))
c1 = np.array([239/256, 0/256, 150/256, 1])
c2 = np.array([0/256, 220/256, 181/256, 1])
c3 = np.array([255/256, 207/256, 226/256, 1])
c4 = np.array([255/256, 220/256, 61/256, 1])
newcolors[:25, :] = c1
newcolors[25:50, :] = c2
newcolors[50:75, :] = c3
newcolors[75:, :] = c4
newcmp = ListedColormap(newcolors)

f,ax = plt.subplots(1,1,tight_layout=True, dpi=300)

w = 25
h = 25
seqconv = [{'A':1,'U':2,'C':3,'G':4}[i] for i in list(bactin_1.source_seq)]
plt.imshow(np.array(seqconv[:w*h]).reshape(w,h), cmap=newcmp)
lmat = np.array(seqconv[:w*h]).reshape(w,h)
for i in range(h):
    for j in range(w):
        text = ax.text(j, i, {1:'A',2:'U',3:'C',4:'G'}[lmat[i, j]],
                       ha="center", va="center", color="black",fontsize=4)
plt.savefig('seq.svg')

t = np.linspace(0,2000,2000)
rss.solver.colors=2
bactin_soln1 = rss.solver.solve_ssa(bactin_1.kelong,t, ki=.033, n_traj=1000, probe_vec= bactin_1.probe_vec, probe_loc = bactin_1.probe_loc )
bactin_soln2 = rss.solver.solve_ssa(bactin_2.kelong,t, ki=.033, n_traj=1000, probe_vec= bactin_2.probe_vec, probe_loc = bactin_2.probe_loc )

f,ax = plt.subplots(1,2,tight_layout=True)

x,bins = np.histogram(bactin_soln1.intensity_vec[0,200::200,:].flatten()/10,bins=np.linspace(0,10,11), density=True)
ax[0].hist(x,bins, color=colors[0],lw=2, histtype = 'step') 
x,bins = np.histogram(bactin_soln2.intensity_vec[0,200::200,:].flatten()/10,bins=np.linspace(0,10,11), density=True)
ax[0].hist(x,bins, color=colors[1],lw=2, histtype='step') 

ax[0].hist(bactin_soln1.intensity_vec[0,200::200,:].flatten()/10, bins=np.linspace(0,10,11),alpha=.3, density=True)
ax[0].hist(bactin_soln2.intensity_vec[0,200::200,:].flatten()/10, bins=bins,alpha=.3, density=True)


x,bins = np.histogram(bactin_soln1.intensity_vec[1,200::200,:].flatten()/10,bins=np.linspace(0,10,11), density=True)
ax[1].hist(x,bins, color=colors[0], lw=2, histtype='step') 
x,bins = np.histogram(bactin_soln2.intensity_vec[1,200::200,:].flatten()/10,bins=np.linspace(0,10,11), density=True)
ax[1].hist(x,bins, color=colors[1], lw=2, histtype='step') 
ax[1].hist(bactin_soln1.intensity_vec[1,200::200,:].flatten()/10,bins=np.linspace(0,10,11),alpha=.3, density=True)
ax[1].hist(bactin_soln2.intensity_vec[1,200::200,:].flatten()/10,bins=np.linspace(0,10,11),alpha=.3, density=True)



ax[0].set_xlabel('Intensity (UMP)')
ax[1].set_xlabel('Intensity (UMP)')

ax[0].set_ylabel('Probability')
ax[1].set_ylabel('Probability')

ax[0].set_title('SunTag')
ax[1].set_title('FlagTag')

plt.legend(['C1','C2'])
if save:
    plt.savefig('2_color_intensity_hist_fig2.svg')


plt.figure()
bactin_1.visualize_probe(colors=[colors[0],colors[1]])
plt.plot([len(flagtag_codons), len(flagtag_codons)], [0,1] ,'#00ff00')
plt.plot([-5, -5], [0,1] ,'#00ff00')
plt.plot([-5, len(flagtag_codons)], [1,1] ,'#00ff00')
plt.plot([-5, len(flagtag_codons)], [0,0] ,'#00ff00')


plt.plot([len(suntag_codons + flagtag_codons), len(suntag_codons + flagtag_codons)], [0,1] ,'#00ffff')
plt.plot([len( flagtag_codons)+1, len(flagtag_codons)+1], [0,1] ,'#00ffff')
plt.plot([len( flagtag_codons)+1, len(suntag_codons + flagtag_codons)], [1,1] ,'#00ffff')
plt.plot([len( flagtag_codons)+1,len(suntag_codons + flagtag_codons)], [0,0] ,'#00ffff')
plt.text(0,1.5,'10X FLAG-Tag')
plt.text(len( flagtag_codons),1.5,'10X Sun-Tag')

plt.text(len( flagtag_codons + suntag_codons)+100,.3,r'$\beta$-Actin')
if save:
    plt.savefig('construct1_fig2.svg')


plt.figure()
bactin_2.visualize_probe(colors=[colors[0],colors[1]])
plt.plot([len(flagtag_codons), len(flagtag_codons)], [0,1] ,'#00ff00')
plt.plot([-5, -5], [0,1] ,'#00ff00')
plt.plot([-5, len(flagtag_codons)], [1,1] ,'#00ff00')
plt.plot([-5, len(flagtag_codons)], [0,0] ,'#00ff00')


plt.plot([len(construct2), len(construct2)], [0,1] ,'#00ffff')
plt.plot([len(bactin_codons + flagtag_codons), len(bactin_codons + flagtag_codons)], [0,1] ,'#00ffff')
plt.plot([len(bactin_codons + flagtag_codons), len(construct2)], [1,1] ,'#00ffff')
plt.plot([len(bactin_codons + flagtag_codons), len(construct2)], [0,0] ,'#00ffff')
plt.text(0,1.5,'10X FLAG-Tag')
plt.text(len(bactin_codons + flagtag_codons),1.5,'10X Sun-Tag')

plt.text(len( flagtag_codons)+100,.3,r'$\beta$-Actin')
if save:
    plt.savefig('construct2_fig2.svg')


f,ax = plt.subplots(1,2,tight_layout=True)

a,e = rss.inta.get_autocov(bactin_soln1.intensity_vec[:,1000:,:],norm='global')
acc, acc_error = rss.inta.get_autocorr(a,shot_noise_type='G0')
ax[0].plot([0,200],[0,0],'g--', label='_nolegend_')
ax[1].plot([0,200],[0,0],'g--', label='_nolegend_')

ax[0].plot(np.mean(acc[0],axis=1)[:200],color=colors[0] )

ax[0].plot(np.mean(acc[0],axis=1)[:200] -acc_error[0][:200],color=colors[0],ls='--',  label='_nolegend_' )
ax[0].plot(np.mean(acc[0],axis=1)[:200] + acc_error[0][:200],color=colors[0],ls='--',  label='_nolegend_' )

ax[1].plot(np.mean(acc[1],axis=1)[:200],color=colors[0] )
ax[1].plot(np.mean(acc[1],axis=1)[:200] - acc_error[1][:200],color=colors[0],ls='--', label='_nolegend_' )
ax[1].plot(np.mean(acc[1],axis=1)[:200] + acc_error[1][:200],color=colors[0],ls='--', label='_nolegend_' )



a,e = rss.inta.get_autocov(bactin_soln2.intensity_vec[:,1000:,:],norm='global')
acc, acc_error = rss.inta.get_autocorr(a,shot_noise_type='G0')
ax[0].plot(np.mean(acc[0],axis=1)[:200],color=colors[1] )
ax[0].plot(np.mean(acc[0],axis=1)[:200] - acc_error[0][:200],color=colors[1],ls='--' )
ax[0].plot(np.mean(acc[0],axis=1)[:200] + acc_error[0][:200],color=colors[1],ls='--' )



ax[1].plot(np.mean(acc[1],axis=1)[:200],color=colors[1], )
ax[1].plot(np.mean(acc[1],axis=1)[:200] - acc_error[1][:200],color=colors[1],ls='--' )
ax[1].plot(np.mean(acc[1],axis=1)[:200] + acc_error[1][:200],color=colors[1],ls='--' )

ax[0].set_xlabel('Tau (s)')
ax[1].set_xlabel('Tau (s)')

ax[0].set_ylabel('Autocorrelation')
ax[1].set_ylabel('Autocorrelation')

ax[0].set_title('SunTag')
ax[1].set_title('FlagTag')
plt.legend(['C1','C2'])
if save:
    plt.savefig('2_color_acc_fig2.svg')



##################### Figure 3


t = np.linspace(0,2000,2000)
bactin_soln1 = rss.solver.solve_ssa(bactin_1.kelong,t, ki=.033, n_traj=1000, probe_vec= bactin_1.probe_vec, probe_loc = bactin_1.probe_loc, low_memory=False, record_stats=True )
bactin_soln2 = rss.solver.solve_ssa(bactin_2.kelong,t, ki=.033, n_traj=1000, probe_vec= bactin_2.probe_vec, probe_loc = bactin_2.probe_loc, low_memory=False, record_stats=True  )

############################


f,ax = plt.subplots(1,1,tight_layout=True)

ax.hist(bactin_soln1.ribtimes[bactin_soln1.ribtimes !=0].flatten(),bins=30,density=True,align='mid')
ax.set_xlabel('Ribosomal Dwell Time (s)')
ax.set_ylabel('Probability')
ax.set_title('')
if save:
    plt.savefig('dwell_fig3.svg')


f,ax = plt.subplots(1,1,tight_layout=True)
ax.hist(bactin_soln1.collisions,bins=40,density=True,align='mid')
ax.set_xlabel('N Collisions per Ribosome')
ax.set_ylabel('Probability')
ax.set_title('')
plt.savefig('collision_fig3.svg')

f,ax = plt.subplots(1,1,tight_layout=True)
ax.hist(bactin_soln1.collisions[bactin_soln1.collisions!=0],bins=40,density=True,align='mid')
ax.set_xlabel('N Collisions per Ribosome')
ax.set_ylabel('Probability')
ax.set_title('')
if save:
    plt.savefig('collision_fig3_sub.svg')


def movmean(a, n=3): #rolling average
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

f,ax = plt.subplots(1,1,tight_layout=True)
ax.plot(movmean(bactin_soln1.rib_density,n=10))
ax.set_xlabel('Codon (10 rolling average)')
ax.set_ylabel('Occupation Probability')
ax.set_title('')
if save:
    plt.savefig('chip_fig3_sub.svg')


f,ax = plt.subplots(1,1,tight_layout=True)
ax.plot(bactin_soln1.ribosome_locations[0][500:1000],t[500:1000],'.',color=colors[0],markersize=2)
ax.set_xlabel('Codon')
ax.set_ylabel('Time')
ax.set_title('')
ax.invert_yaxis()
ax.set_xlim([1,1000])
if save:
    plt.savefig('kymograph_fig3_sub.svg')


cols = np.hstack([bactin_soln1.col_points[x][:,0] for x in range( len(bactin_soln1.col_points))])
x,bins = np.histogram(cols,bins = np.linspace(0,len(construct1), int((len(construct1)+1)/10) ), density=True )
f,ax = plt.subplots(1,1,tight_layout=True)

import matplotlib.patches as patches
bactin_1.visualize_probe(colors=[colors[0],colors[1]])

plt.text(0,1.5,'10X FLAG-Tag')
plt.text(len( flagtag_codons),1.5,'10X Sun-Tag')

#plt.text(len( flagtag_codons + suntag_codons)+100,.3,r'$\beta$-Actin')

plt.xlabel('Collision Probability per Location')

plt.title('')
ax =plt.gca()
nx = (x - np.min(x)) / (np.max(x) - np.min(x))
import matplotlib.cm as cm
for i in range(len(x)):
    c = cm.viridis(nx[i])
    n = 10.17894737
    rect = patches.Rectangle((n*i, 0), n,1, linewidth=0, edgecolor='none', facecolor=c,zorder=10)
    ax.add_patch(rect)

plt.plot([len(suntag_codons + flagtag_codons), len(suntag_codons + flagtag_codons)], [0,1] ,'#ff00ff',zorder=22)
plt.plot([len( flagtag_codons)+1, len(flagtag_codons)+1], [0,1] ,'#ff00ff',zorder=22)
plt.plot([len( flagtag_codons)+1, len(suntag_codons + flagtag_codons)], [1,1] ,'#ff00ff',zorder=22)
plt.plot([len( flagtag_codons)+1,len(suntag_codons + flagtag_codons)], [0,0] ,'#ff00ff',zorder=22)

plt.plot([len(flagtag_codons), len(flagtag_codons)], [0,1] ,'#ff00ff',zorder=22)
plt.plot([-5, -5], [0,1] ,'#ff00ff',zorder=22)
plt.plot([-5, len(flagtag_codons)], [1,1] ,'#ff00ff',zorder=22)
plt.plot([-5, len(flagtag_codons)], [0,0] ,'#ff00ff',zorder=22)

plt.plot([len(construct1), len(construct1)], [0,1] ,'#000000',zorder=21)
plt.plot([-4, -4], [0,1] ,'#000000',zorder=21)
plt.plot([-4, len(construct1)], [1,1] ,'#000000',zorder=21)
plt.plot([-4, len(construct1)], [0,0] ,'#000000',zorder=21)







############################ Figure 4 
plt.figure()
I = bactin_soln2.intensity_vec

plt.plot(I[0,500:,2]/10, color=colors[0])
plt.plot(I[1,500:,2]/10, color=colors[1])

plt.ylabel('Time (s)')
plt.xlabel('Intensity (UMP)')
plt.legend(['color 1', 'color 2'])
if save:
    plt.savefig('intensity_fig4.svg')


plt.figure()
I = bactin_soln2.intensity_vec

plt.plot([0,200],[0,0],'g--', label='_nolegend_')
a,e = rss.inta.get_autocov(bactin_soln2.intensity_vec[:,1000:,:],norm='global')
acc, acc_error = rss.inta.get_autocorr(a,g0='G0')
plt.plot(np.mean(acc[0],axis=1)[:200],color=colors[0] )
plt.plot(np.mean(acc[0],axis=1)[:200] - acc_error[0][:200],color=colors[0],ls='--', label='_nolegend_' )
plt.plot(np.mean(acc[0],axis=1)[:200] + acc_error[0][:200],color=colors[0],ls='--', label='_nolegend_' )


a,e = rss.inta.get_autocov(bactin_soln2.intensity_vec[:,1000:,:],norm='global')
acc, acc_error = rss.inta.get_autocorr(a,g0='G0')
plt.plot(np.mean(acc[1],axis=1)[:200],color=colors[1] )
plt.plot(np.mean(acc[1],axis=1)[:200] - acc_error[1][:200],color=colors[1],ls='--', label='_nolegend_' )
plt.plot(np.mean(acc[1],axis=1)[:200] + acc_error[1][:200],color=colors[1],ls='--', label='_nolegend_' )


plt.ylabel('Autocorrelation')
plt.xlabel('Tau (s)')
plt.legend(['color 1', 'color 2'])
if save:
    plt.savefig('acc_fig4.svg')


plt.figure()
plt.plot([-100,200],[0,0],'g--', label='_nolegend_')
plt.plot([0,0],[0,.7],'g--', label='_nolegend_')
a =rss.inta.get_crosscorr(bactin_soln2.intensity_vec[:,1000:,:],norm='global',g0='indiv_max')


plt.plot(np.linspace(-100,199,300), np.mean(a[0][1],axis=1)[900:1200],  color=colors[0])

plt.plot( np.linspace(-100,199,300),  np.mean(a[0][1],axis=1)[900:1200] - a[1][1,900:1200]   ,color=colors[0],ls='--', label='_nolegend_' )
plt.plot( np.linspace(-100,199,300),  np.mean(a[0][1],axis=1)[900:1200] + a[1][1,900:1200]   ,color=colors[0],ls='--', label='_nolegend_' )



plt.plot(np.linspace(-100,199,300), np.mean(a[0][0],axis=1)[900:1200],  color=colors[1])

plt.plot( np.linspace(-100,199,300),  np.mean(a[0][0],axis=1)[900:1200] - a[1][1,900:1200]   ,color=colors[1],ls='--', label='_nolegend_' )
plt.plot( np.linspace(-100,199,300),  np.mean(a[0][0],axis=1)[900:1200] + a[1][1,900:1200]   ,color=colors[1],ls='--', label='_nolegend_' )



plt.plot(np.linspace(-100,199,300), np.mean(a[0][3],axis=1)[900:1200],  color=colors[2])

plt.plot( np.linspace(-100,199,300),  np.mean(a[0][3],axis=1)[900:1200] - a[1][1,900:1200]   ,color=colors[2],ls='--', label='_nolegend_' )
plt.plot( np.linspace(-100,199,300),  np.mean(a[0][3],axis=1)[900:1200] + a[1][1,900:1200]   ,color=colors[2],ls='--', label='_nolegend_' )




plt.ylabel('Cross-Correlation')
plt.xlabel('Tau (s)')
plt.legend(['Flag-Sun cc', 'Sun acc', 'Flag acc'])
if save:
    plt.savefig('cc_fig4.svg')


plt.figure()
x,bins = np.histogram(bactin_soln2.intensity_vec[0,200::200,:].flatten()/10,bins=np.linspace(0,10,11)-.5, density=True)
plt.hist(x,bins, color=colors[0], lw=2) 
x,bins = np.histogram(bactin_soln2.intensity_vec[1,200::200,:].flatten()/10,bins=np.linspace(0,10,11)-.5, density=True)
plt.hist(x,bins, color=colors[1], lw=2) 
plt.hist(bactin_soln2.intensity_vec[0,200::200,:].flatten()/10,bins=np.linspace(0,10,11)-.5,alpha=.3, density=True, align='mid')
plt.hist(bactin_soln2.intensity_vec[1,200::200,:].flatten()/10,bins=np.linspace(0,10,11)-.5,alpha=.3, density=True, align='mid')

plt.xlabel('Intensity (UMP)')

plt.ylabel('Probability')


plt.legend(['SunTag','FlagTag'])


if save:
    plt.savefig('hist_fig4.svg')

n = 25
f,ax = plt.subplots(1,3,tight_layout=True)
    
tau = 3

ax[0].scatter(bactin_soln2.intensity_vec[1,200:-tau, :n]/10, bactin_soln2.intensity_vec[1,200+tau:, :n]/10,alpha=.1,marker='.', facecolor=colors[1] )
ax[0].scatter(bactin_soln2.intensity_vec[0,200:-tau, :n]/10, bactin_soln2.intensity_vec[0,200+tau:, :n]/10,alpha=.1,marker='.' , facecolor=colors[0])
ax[0].plot([0,8],[0,8],'r--')
ax[0].set_xlabel('I(t) (UMP)')
ax[0].set_ylabel('I(t + tau) (UMP)')
ax[0].set_title('tau = 3s')

tau = 50
ax[1].scatter(bactin_soln2.intensity_vec[1,200:-tau, :n]/10, bactin_soln2.intensity_vec[1,200+tau:, :n]/10,alpha=.1,marker='.', facecolor=colors[1] )
ax[1].scatter(bactin_soln2.intensity_vec[0,200:-tau, :n]/10, bactin_soln2.intensity_vec[0,200+tau:, :n]/10,alpha=.1,marker='.', facecolor=colors[0] )
ax[1].plot([0,7],[0,8],'r--')
ax[1].set_xlabel('I(t) (UMP)')
#ax[1].set_ylabel('Intensity(t + tau) (UMP)')
ax[1].set_title('tau = 50s')

    
tau = 300
ax[2].scatter(bactin_soln2.intensity_vec[1,200:-tau, :n]/10, bactin_soln2.intensity_vec[1,200+tau:, :n]/10,alpha=.1,marker='.', facecolor=colors[1] )
ax[2].scatter(bactin_soln2.intensity_vec[0,200:-tau, :n]/10, bactin_soln2.intensity_vec[0,200+tau:, :n]/10,alpha=.1,marker='.', facecolor=colors[0] )
ax[2].plot([0,8],[0,8],'r--')
ax[2].set_xlabel('I(t) (UMP)')
#ax[2].set_ylabel('Intensity(t + tau) (UMP)')
ax[2].set_title('tau = 300s')


ax[0].set_aspect(aspect='equal')
ax[1].set_aspect(aspect='equal')
ax[2].set_aspect(aspect='equal')


if save:
    plt.savefig('icomp1_fig4.svg')

f,ax = plt.subplots(1,3,tight_layout=True,)
    
tau = 3

ax[0].scatter(bactin_soln2.intensity_vec[1,200:-tau, :n]/10, bactin_soln2.intensity_vec[0,200+tau:, :n]/10,alpha=.1,marker='.', facecolor=colors[2] )
ax[0].plot([0,8],[0,8],'r--')
ax[0].set_xlabel('I1(t) (UMP)')
ax[0].set_ylabel('I2(t + tau) (UMP)')
ax[0].set_title('tau = 3s')


tau = 50
ax[1].scatter(bactin_soln2.intensity_vec[1,200:-tau, :n]/10, bactin_soln2.intensity_vec[0,200+tau:, :n]/10,alpha=.1,marker='.', facecolor=colors[2] )
ax[1].plot([0,8],[0,8],'r--')
ax[1].set_xlabel('I1(t) (UMP)')
#ax[1].set_ylabel('Intensity(t + tau) (UMP)')
ax[1].set_title('tau = 50s')

tau = 500
ax[2].scatter(bactin_soln2.intensity_vec[1,200:-tau, :n]/10, bactin_soln2.intensity_vec[0,200+tau:, :n]/10,alpha=.1,marker='.', facecolor=colors[2] )
ax[2].plot([0,8],[0,8],'r--')
ax[2].set_xlabel('I1(t) (UMP)')
#ax[2].set_ylabel('Intensity(t + tau) (UMP)')
ax[2].set_title('tau = 300s')


ax[0].set_aspect(aspect='equal')
ax[1].set_aspect(aspect='equal')
ax[2].set_aspect(aspect='equal')

    
if save:
    plt.savefig('icomp2_fig4.svg')




def norm(sig):
   return(sig - np.min(sig))/ (np.max(sig) - np.min(sig))


m0 = np.mean(bactin_soln2.intensity_vec[0])/10
m1 = np.mean(bactin_soln2.intensity_vec[1])/10
n = 100
f,ax = plt.subplots(2,3,tight_layout=True)
nbins = 12
tau = 3
x,y = norm(bactin_soln2.intensity_vec[1,200:-tau, :n]/10) , norm(bactin_soln2.intensity_vec[1,200+tau:, :n]/10)
heatmap,xedges,yedges  = np.histogram2d(x.flatten(), y.flatten(), bins=nbins)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

ax[0,0].imshow(heatmap.T, extent=extent, origin='lower' )
ax[0,0].plot([0,5],[0,5],'r--')
x,y = norm(bactin_soln2.intensity_vec[1,200:-tau, :n]/10), norm(bactin_soln2.intensity_vec[0,200+tau:, :n]/10)
heatmap,xedges,yedges  = np.histogram2d(x.flatten(), y.flatten(), bins=nbins)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

ax[1,0].imshow(heatmap.T, extent=extent, origin='lower' )


#ax[0].scatter(bactin_soln2.intensity_vec[0,200:-tau, :n]/10, bactin_soln2.intensity_vec[0,200+tau:, :n]/10,alpha=.1,marker='.' , facecolor=colors[0])
ax[1,0].plot([0,5],[0,5],'r--')
ax[1,0].set_xlabel('I1(t) (UMP)')
ax[0,0].set_ylabel('I1(t + tau) (UMP)')
ax[1,0].set_ylabel('I2(t + tau) (UMP)')
ax[0,0].set_title('tau = 3s')





tau = 50

x,y = norm(bactin_soln2.intensity_vec[1,200:-tau, :n]/10), norm(bactin_soln2.intensity_vec[1,200+tau:, :n]/10)
heatmap,xedges,yedges  = np.histogram2d(x.flatten(), y.flatten(), bins=nbins)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

ax[0,1].imshow(heatmap.T, extent=extent, origin='lower' )
ax[0,1].plot([0,5],[0,5],'r--')
x,y = norm(bactin_soln2.intensity_vec[1,200:-tau, :n]/10), norm(bactin_soln2.intensity_vec[0,200+tau:, :n]/10)
heatmap,xedges,yedges  = np.histogram2d(x.flatten(), y.flatten(), bins=nbins)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

ax[1,1].imshow(heatmap.T, extent=extent, origin='lower' )


#ax[0].scatter(bactin_soln2.intensity_vec[0,200:-tau, :n]/10, bactin_soln2.intensity_vec[0,200+tau:, :n]/10,alpha=.1,marker='.' , facecolor=colors[0])
ax[1,1].plot([0,5],[0,5],'r--')
ax[1,1].set_xlabel('I1(t) (UMP)')

ax[0,1].set_title('tau = 50s')


tau = 500
x,y = norm(bactin_soln2.intensity_vec[1,200:-tau, :n]/10), norm(bactin_soln2.intensity_vec[1,200+tau:, :n]/10)
heatmap,xedges,yedges  = np.histogram2d(x.flatten(), y.flatten(), bins=nbins)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

ax[0,2].imshow(heatmap.T, extent=extent, origin='lower' )
ax[0,2].plot([0,5],[0,5],'r--')
x,y = norm(bactin_soln2.intensity_vec[1,200:-tau, :n]/10), norm(bactin_soln2.intensity_vec[0,200+tau:, :n]/10)
heatmap,xedges,yedges  = np.histogram2d(x.flatten(), y.flatten(), bins=nbins)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

ax[1,2].imshow(heatmap.T, extent=extent, origin='lower' )


#ax[0].scatter(bactin_soln2.intensity_vec[0,200:-tau, :n]/10, bactin_soln2.intensity_vec[0,200+tau:, :n]/10,alpha=.1,marker='.' , facecolor=colors[0])
ax[1,2].plot([0,5],[0,5],'r--')
ax[1,2].set_xlabel('I1(t) (UMP)')
ax[0,2].set_title('tau = 500s')



ax[0,0].set_aspect(aspect='equal')
ax[0,1].set_aspect(aspect='equal')
ax[0,2].set_aspect(aspect='equal')
ax[0,0].set_xlim([0,1])
ax[0,1].set_xlim([0,1])
ax[0,2].set_xlim([0,1])
ax[0,0].set_ylim([0,1])
ax[0,1].set_ylim([0,1])
ax[0,2].set_ylim([0,1])

ax[1,0].set_aspect(aspect='equal')
ax[1,1].set_aspect(aspect='equal')
ax[1,2].set_aspect(aspect='equal')

ax[1,0].set_xlim([0,1])
ax[1,1].set_xlim([0,1])
ax[1,2].set_xlim([0,1])
ax[1,0].set_ylim([0,1])
ax[1,1].set_ylim([0,1])
ax[1,2].set_ylim([0,1])


if save:
    plt.savefig('icomp1_fig4.svg')
####################


plt.figure()

base_mw = rss.diffcalc.calculate_rna_strand_base_mw(bactin_1.nt_seq + 'A'*100, )
mw_trajectory = rss.diffcalc.mw_over_time(bactin_soln1.ribosome_locations[0], bactin_1.nt_seq, bactin_1.probe_loc, base_mw, fluorophore=['Cy3','GFP' ] )
plt.plot(mw_trajectory)

