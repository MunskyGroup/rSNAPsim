# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 20:38:17 2021

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

os.chdir(cwd)

os.chdir('..')
os.chdir('fjc')
import fjc_translation_cpp
os.chdir(cwd)

import matplotlib.pyplot as plt
from cycler import cycler
########################################
dark = False
if not dark:
    colors = ['#ef476f', '#073b4c','#06d6a0','#7400b8','#073b4c', '#118ab2',]
else:
    plt.style.use('dark_background')
    colors = ['#57ffcd', '#118ab2', '#ff479d', '#ffe869','#ff8c00','#04756f']

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

save = False

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


bactin_1 = rss.seqmanip.seq_to_protein_obj(construct_sequence_str1,add_tag=False)['1'][0]
bactin_2 = rss.seqmanip.seq_to_protein_obj(construct_sequence_str2,add_tag=False)['1'][0]

t = np.linspace(0,2000,4000)
bactin_soln1 = rss.solver.solve_ssa(bactin_1.kelong,t, ki=.02, low_memory=False, n_traj=1, probe_vec= bactin_1.probe_vec, probe_loc = bactin_1.probe_loc )
bactin_soln2 = rss.solver.solve_ssa(bactin_2.kelong,t, ki=.02, low_memory=False, n_traj=1, probe_vec= bactin_2.probe_vec, probe_loc = bactin_2.probe_loc )



construct_L = bactin_1.probe_vec.shape[1]*3 + 100
footprint = 10

pL =.33
sL = 5
tscale = .01
seed = np.random.randint(0x7fffff)
seed = 0

nn = 1000
pos_tensor = bactin_soln1.ribosome_locations[:,:nn,:]

pt = np.ascontiguousarray(pos_tensor[0]*3, dtype=np.int32)

xi = np.zeros((construct_L, len(t)) ,dtype=np.float64)    
yi = np.zeros(( construct_L, len(t)),dtype=np.float64)    
st = time.time()
fjc_translation_cpp.generate_wlc_2d(xi,yi,pt,pL,sL,footprint,construct_L,0,t,tscale)
print(time.time()-st)

xi = np.array([xi.T])
yi = np.array([yi.T])

import time
st = time.time()
#xi,yi, = rss.riba.WLC_positions_2D(bactin_soln1.ribosome_locations[:,:nn,:], bactin_1.probe_vec,
                                   #.08, compressed_length=.3,)
print(time.time() - st)

_, intensity_tensor, utr_intensity = rss.riba.linear_probe_positions(bactin_soln1.ribosome_locations[:,:nn,:], bactin_1.probe_vec)


for i in range(nn):
    com_x = np.sum(xi[0,i,:]) / len(xi[0,i,:])
    com_y = np.sum(yi[0,i,:]) / len(yi[0,i,:])
    xi[0,i,:] = xi[0,i,:] - com_x
    yi[0,i,:] = yi[0,i,:] - com_y  



from matplotlib.animation import FuncAnimation

max_rib = intensity_tensor.shape[-1]

fig, ax = plt.subplots()
ax.set(facecolor = '#555555')
xdata, ydata = [], []
ln, = plt.plot([], [], 'b-',lw=.5)
ln2, = plt.plot([], [], 'ko', markersize=5)
ln3, = plt.plot([], [], 'ro', markersize=5)

color1 = plt.scatter([0,]*max_rib , [0,]*max_rib, s=0, c='#00FF00', zorder=3, alpha=.5)
color2 = plt.scatter([0,]*max_rib, [0,]*max_rib, s=0, c='#00FFFF',zorder=3, alpha=.5)
bbox = [np.min(xi), np.min(yi), np.max(xi), np.max(yi) ]

def init():
    ax.set_xlim(bbox[0], bbox[2])
    ax.set_ylim(bbox[1], bbox[3])
    ax.set_xlabel('nm')
    ax.set_ylabel('nm')
    return ln,

def update(frame):
    ln.set_data(xi[0,frame], yi[0,frame])
    occupied = pos_tensor[0,frame,:][pos_tensor[0,frame,:] !=0]
    if len(occupied) != 0:
        ln2.set_data( xi[0,frame, occupied.astype(int)*3 ], yi[0,frame, occupied.astype(int)*3]  )
    else:
        ln2.set_data([],[])
    ln3.set_data(xi[0,frame, -10], yi[0,frame, -10])
    
    if len(occupied) != 0:
        cint1  = intensity_tensor[0,0,frame,:][pos_tensor[0,frame,:] !=0]*4
        cint2 = intensity_tensor[1,0,frame,:][pos_tensor[0,frame,:] !=0]*4
        color1.set_offsets(np.array([xi[0,frame, occupied.astype(int)*3 ]-.3, yi[0,frame, occupied.astype(int)*3]-.3]).T )
        color2.set_offsets(np.array([xi[0,frame, occupied.astype(int)*3 ]+.3, yi[0,frame, occupied.astype(int)*3]+.3]).T )

        color1.set_sizes(cint1)
        color2.set_sizes(cint2)
    else:
        color1.set_offsets( np.array([[0,]*max_rib,[0,]*max_rib]).T)
        color2.set_offsets(np.array([[0,]*max_rib,[0,]*max_rib]).T)
        color1.set_sizes([0,]*max_rib)
        color2.set_sizes([0,]*max_rib)
        
    return ln,ln2,ln3,color1,color2

ani = FuncAnimation(fig, update, frames=[x for x in range(nn)],
                    init_func=init, blit=True)
ani.save('test_wlc_long.mp4', writer='ffmpeg', fps=12)


     
'''        

x,y,I, utr_I = rss.riba.circular_probe_positions(bactin_soln1.ribosome_locations, bactin_1.probe_vec)
plt.figure()
plt.scatter(x[0].flatten(), y[0].flatten(), s=I[1][0].flatten()**2.5,color='#00FF00' );
plt.scatter(x[0].flatten(), y[0].flatten(), s=I[0][0].flatten()**2.5,color='#00FFFF' );
plt.scatter(utr_I[0],utr_I[1],s=utr_I[2]*10,color='r')
plt.xlabel('nm'); plt.ylabel('nm')


x,y,I, utr_I = rss.riba.circular_probe_positions(bactin_soln2.ribosome_locations, bactin_2.probe_vec)
plt.figure()
plt.scatter(x[0].flatten(), y[0].flatten(), s=I[1][0].flatten()**2.5,color='#00FF00' ,alpha=.5);
plt.scatter(x[0].flatten(), y[0].flatten(), s=I[0][0].flatten()**2.5,color='#00FFFF',alpha=.5 );
plt.scatter(utr_I[0],utr_I[1],s=utr_I[2]*10,color='r')
plt.xlabel('nm'); plt.ylabel('nm')

###################################

x,I, utr_I = rss.riba.linear_probe_positions(bactin_soln1.ribosome_locations, bactin_1.probe_vec)
plt.figure()
plt.scatter(x[0].flatten(), np.ones(x[0].flatten().shape), s=I[1][0].flatten()**2.5,color='#00FF00' ,alpha=.5);
plt.scatter(x[0].flatten(), np.ones(x[0].flatten().shape), s=I[0][0].flatten()**2.5,color='#00FFFF',alpha=.5 );
plt.scatter(utr_I[0],1,s=utr_I[1]*10,color='r')
plt.xlabel('nm'); plt.ylabel('nm')




x,I, utr_I = rss.riba.linear_probe_positions(bactin_soln2.ribosome_locations, bactin_2.probe_vec)
plt.figure()
plt.scatter(x[0].flatten(), np.ones(x[0].flatten().shape), s=I[1][0].flatten()**2.5,color='#00FF00' ,alpha=.5);
plt.scatter(x[0].flatten(), np.ones(x[0].flatten().shape), s=I[0][0].flatten()**2.5,color='#00FFFF',alpha=.5 );
plt.scatter(utr_I[0],1,s=utr_I[1]*10,color='r')
plt.xlabel('nm'); plt.ylabel('nm')



####################################

x,y,I, utr_I = rss.riba.spring_probe_positions(bactin_soln1.ribosome_locations, bactin_1.probe_vec, stretched_length=100)
plt.figure()
plt.scatter(x[0].flatten(), y[0].flatten(), s=I[1][0].flatten()**2.5,color='#00FF00',alpha=.05 );
plt.scatter(x[0].flatten(), y[0].flatten(), s=I[0][0].flatten()**2.5,color='#00FFFF',alpha=.05 );
plt.scatter(utr_I[0],utr_I[1],s=utr_I[2]*10,color='r')
plt.xlabel('nm'); plt.ylabel('nm')


x,y,I, utr_I = rss.riba.spring_probe_positions(bactin_soln2.ribosome_locations, bactin_2.probe_vec, stretched_length=100)
plt.figure()
plt.scatter(x[0].flatten(), y[0].flatten(), s=I[1][0].flatten()**2.5,color='#00FF00' ,alpha=.05);
plt.scatter(x[0].flatten(), y[0].flatten(), s=I[0][0].flatten()**2.5,color='#00FFFF',alpha=.05 );
plt.scatter(utr_I[0],utr_I[1],s=utr_I[2]*10,color='r')
plt.xlabel('nm'); plt.ylabel('nm')

'''