# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 12:33:25 2021

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
        'size'   : 8}

save = False

plt.rcParams.update({'font.size': 12, 'font.weight':'bold','font.family':'normal'  }   )
plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})

plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})
plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})


plt.rcParams.update({'xtick.major.width'   : 2.8 })
plt.rcParams.update({'xtick.labelsize'   : 8 })



plt.rcParams.update({'ytick.major.width'   : 2.8 })
plt.rcParams.update({'ytick.labelsize'   : 8})

plt.rcParams.update({'axes.titleweight'   : 'bold'})
plt.rcParams.update({'axes.titlesize'   : 10})
plt.rcParams.update({'axes.labelweight'   : 'bold'})
plt.rcParams.update({'axes.labelsize'   : 8})

plt.rcParams.update({'axes.linewidth':2.8})
plt.rcParams.update({'axes.labelpad':8})
plt.rcParams.update({'axes.titlepad':10})
plt.rcParams.update({'figure.dpi':300})




suntag_kif = '''atgataataGCTAGCTCTGGAGGAGAAGAACTTTTGAGCAAGAATTATCATCTTGAGAACGAAGTGGCTCGTCTTAAGAAAGGTTCTGGCAGTGGAGAAGAACTGCTTTCAAAGAATTACCACCTGGAAAATGAGGTAGCTAGACTGAAAAAGGGGAGCGGAAGTGGGGAGGAGTTGCTGAGCAAAAATTATCATTTGGAGAACGAAGTAGCACGACTAAAGAAAGGGTCCGGATCGGGTGAGGAGTTACTCTCGAAAAATTATCATCTCGAAAACGAAGTGGCTCGGCTAAAAAAGGGCAGTGGTTCTGGAGAAGAGCTATTATCTAAAAACTACCACCTCGAAAATGAGGTGGCACGCTTAAAAAAGGGAAGTGGCAGTGGTGAAGAGCTACTATCCAAGAATTATCATCTTGAGAACGAGGTAGCGCGTTTGAAGAAGGGTTCCGGCTCAGGAGAGGAACTGCTCTCGAAGAACTATCATCTTGAAAATGAGGTCGCTCGATTAAAAAAGGGATCGGGCAGTGGTGAGGAACTACTTTCAAAGAATTACCACCTCGAAAACGAAGTAGCTCGATTAAAGAAAGGTTCAGGGTCGGGTGAAGAATTACTGAGTAAAAATTATCATCTGGAAAATGAGGTAGCGAGACTAAAAAAGGGGAGTGGTTCTGGCGAAGAGTTGCTATCGAAAAATTATCATCTTGAGAACGAAGTTGCTAGGCTCAAAAAGGGCTCAGGCTCAGGCGAGGAGTTGCTCTCGAAAAACTACCACTTGGAAAATGAGGTCGCGAGGTTGAAAAAGGGGAGCGGGTCGGGCGAGGAGTTATTGAGCAAAAACTATCATTTAGAGAACGAAGTCGCGCGCTTAAAGAAAGGCTCGGGCTCGGGCGAAGAACTCTTATCGAAGAACTACCACCTCGAAAATGAGGTCGCCAGGTTGAAAAAGGGCAGTGGCAGCGGGGAGGAACTCTTGAGCAAGAACTACCACTTGGAGAATGAGGTCGCGAGATTGAAGAAAGGGTCGGGGAGCGGCGAGGAATTGCTCAGCAAGAATTATCATTTGGAGAACGAAGTCGCCAGGCTCAAGAAAGGCTCGGGGTCGGGGGAGGAATTGTTGAGTAAAAACTACCACTTGGAAAATGAAGTCGCCAGGCTCAAAAAAGGGAGTGGGAGCGGCGAAGAGTTATTGAGCAAAAATTACCACTTGGAGAACGAAGTGGCAAGGCTCAAGAAAGGGAGCGGCAGCGGGGAGGAGCTCTTATCGAAGAACTACCACTTAGAGAATGAAGTCGCCCGCTTGAAGAAAGGCTCGGGGAGCGGGGAAGAGCTCTTGAGCAAGAACTACCACTTGGAAAATGAGGTGGCGCGCTTGAAGAAAGGGAGCGGGAGCGGGGAAGAGTTACTATCTAAGAATTATCATCTCGAGAACGAGGTGGCTCGACTAAAGAAGGGCTCCGGCAGTGGGGAGGAACTCCTGTCGAAGAACTATCATCTTGAAAATGAGGTTGCAAGACTTAAAAAGGGGTCCGGATCAGGTGAGGAACTACTCAGTAAGAATTACCACCTGGAAAACGAAGTTGCACGTTTGAAGAAAGGATCAGGATCAGGCGAAGAACTGCTCTCAAAAGATTATCATTTGGAAAATGAGGTTGCACGTTTAAAAAAGGGAAGTGGCAGTGGTGAGGAACTTCTGTCGAAAAATTATCATCTCGAGAATGAAGTAGCCCGACTTAAAAAGGGTTCTGGCTCGGGTCAGCGGCCGCACCGGTCAGCAGTGGAGGACAGCACGCTGCAAGTAGTGGTACGGGTGCGGCCCCCCACCCCTCGGGAGCTGGACAGTCAGCGGCGGCCAGTGGTTCAGGTGGTGGACGAGCGGGTGCTGGTGTTTAACCCTGAGGAGCCCGATGGAGGGTTCCCTGGCCTGAAATGGGGTGGCACCCATGATGGCCCCAAGAAGAAGGGCAAAGACCTGACGTTTGTCTTTGACCGGGTCTTTGGCGAGGCGGCCACCCAACAGGACGTGTTCCAGCACACCACGCACAGCGTCCTGGACAGCTTCCTCCAGGGCTACAACTGCTCAGTGTTTGCCTACGGGGCCACCGGGGCTGGGAAGACACACACCATGCTGGGAAGGGAGGGGGACCCCGGCATCATGTACCTGACCACCGTGGAACTGTACAGGCGCCTGGAGGCCCGCCAGCAGGAGAAGCACTTCGAGGTGCTCATCAGCTACCAGGAGGTGTATAATGAACAGATCCATGACCTCCTGGAGCCCAAGGGGCCCCTTGCCATCCGCGAGGACCCCGACAAGGGGGTGGTGGTGCAAGGACTTTCTTTCCACCAGCCAGCCTCAGCCGAGCAGCTGCTGGAGATACTGACCAGGGGGAACCGTAACCGCACGCAGCACCCCACTGATGCCAACGCGACTTCCTCCCGCTCCCATGCCATCTTCCAGATCTTTGTGAAGCAGCAGGACCGGGTTCCAGGACTGACCCAGGCTGTCCAGGTGGCCAAGATGAGCCTGATTGACCTGGCTGGCTCAGAGCGGGCATCCAGCACCCATGCGAAGGGGGAGCGGCTGCGGGAGGGGGCCAACATCAACCGCTCTCTGCTGGCGCTCATCAACGTCCTCAATGCCTTGGCCGATGCAAAGGGCCGCAAGACCGCTGTGCCCTACGCGGACAGCGCACTGACCCGCCTGCTCAAAGACTCCCTCGGGGGCAACTGCCGCACAGTGATGATCGCTGCCATCAGCCCCTCCAGCCTGACCTACGAGGACACGTATAATACCCTCAAATATGCCGACCGGGCCAAGGAGATCAGGCTCTCGCTGAAGAGCAATGTGACCAGCCTGGACTGTCACATCAGCCAGTATGCTACCATCTGCCAACAGCTCCAGGCTGAGGTAGCCGCTCTGAGGAAGAAGCTCCAAGTGTATGAGGGGGGAGGCCAGCCCCCACCACAGGACCTCCCAGGATCTCCCAAGTCGGGACCACCACCAGAACACCTTCCCAGCTCCCCCTTGCCACCCCACCCTCCCAGCCAGCCCTGCACCCCAGAGCTCCCTGCAGGGCCTAGAGCCCTTCAAGAGGAGAGTCTGGGGATGGAGGCCCAGGTGGAGAGGGCCATGGAAGGGAACTCTTCAGACCAGGAGCAGTCCCCAGAGGATGAGGATGAAGGCCCAGCTGAGGAGGTTCCAACCCAGATGCCAGAGCAGAACCCCACACATGCACTGCCAGAGTCCCCTCGCCTGACCCTGCAGCCCAAGCCAGTCGTGGGCCACTTCTCAGCACGGGAACTGGATGGGGACCGTTCTAAGCGGTTGGCCCTAAAGGTGCTGTGCGTTGCCCAGCGGCAGTACTCCCTGCTCCAAGCAGCCAACCTCCTGACGCCCGACATGATCACAGAGTTTGAGACCCTACAGCAGCTGGTGCAAGAGGAAAAAATTGAGCCTGGGGCAGAGGCCTTGAGGACTTCAGGCCTGGCCAGGGGGGCACCTCTGGCTCAGGAGCTGTGTTCAGAGTCAAAGCCTCCAGGATACACTGGCCCTGTGACCCGGACTATGGCGAGGCGACTGAGTGGCCCCCTGCACACCCTGGGAATCCCGCCTGGACCCAACTGCACCCCAGCCCAGGGGTCCCGATGGCCCATGGAGAAGAAGAGGAGGAGACCAAGCGCCTTGGAGGCAGACAGTCCCATGGCCCCAAAGCGGGGCACCAAGCGCCAGCGCCAGTCCTTCCTGCCCTGCCTAAGGAGAGGGTCTCTGCCTGACACCCAACCTTCACAGGGGCCCAGCACCCCCAAAGGAGAAAGGGCCTCCTCCCCCTGCCATTCCCCTCGCGTTTGCCCAGCCACAGTCATCAAAAGCCGGGTGCCCCTGGGCCCTTCCGCCATGCAGAACTGCTCCACCCCGCTGGCTCTGCCCACTCGAGACCTCAATGCCACCTTTGATCTCTCTGAGGAGCCTCCCTCAAAGCCCAGTTTCCATGAATGCATTGGCTGGGACAAAATACCCCAGGAGCTGAGCAGGCTGGACCAGCCCTTCATCCCCAGGGCACCTGTGCCCCTGTTCACCATGAAGGGCCCCAAGCCAACATCTTCCCTCCCTGGGACCTCTGCCTGCAAGAAGAAGCGCGTTGCGAGTTCCTCAGTCTCCCATGGCCGCAGCCGCATCGCCCGCCTCCCCAGCAGCACTTTGAAGAGGCCAGCTGGGCCCCTTGTACTCCCAGAGCTGCCCTTGAGTCCCCTGTGCCCTAGCAACCGGAGGAATGGAAAGGACCTCATCAGGGTGGGGAGAGCACTCTCAGCAGGGAACGGCGTCACCAAGGTGTCCGATAAGGACCTAGGCGGACTGTTACTGAGCTGCGTTTTACACCCTTTCTTTGACAAAACCTAA'''
bactin_seq = '''ATGGACTACAAGGACGACGACGACAAAGGTGACTACAAAGATGATGACGATAAAGGCGACTATAAGGACGATGACGACAAGGGCGGAAACTCACTGATCAAGGAAAACATGCGGATGAAGGTGGTGATGGAGGGCTCCGTGAATGGTCACCAGTTCAAGTGCACCGGAGAGGGAGAGGGAAACCCGTACATGGGAACTCAGACCATGCGCATTAAGGTCATCGAAGGAGGTCCGCTGCCGTTCGCTTTCGATATCCTGGCCACTTCGTTCGGAGGAGGGTCGCGCACGTTCATCAAGTACCCGAAGGGAATCCCGGACTTCTTTAAGCAGTCATTCCCGGAAGGATTCACTTGGGAACGGGTGACCCGGTATGAAGATGGAGGTGTGGTGACTGTCATGCAAGATACTTCGCTGGAGGATGGGTGCCTCGTGTACCACGTCCAAGTCCGCGGAGTGAATTTCCCGTCCAACGGACCAGTGATGCAGAAAAAGACGAAGGGTTGGGAACCTAATACTGAAATGATGTACCCCGCAGACGGAGGGCTGAGGGGCTACACCCACATGGCGCTGAAGGTCGACGGAGGAGATTACAAGGATGACGACGATAAGCAACAAGATTACAAAGACGATGATGACAAGGGCCAGCAGGGCGACTACAAGGACGACGACGACAAGCAGCAGGACTACAAAGATGACGATGATAAAGGAGGAGGACATCTGTCCTGTTCGTTCGTGACCACCTACAGATCAAAGAAAACCGTGGGAAACATCAAGATGCCGGGCATTCATGCCGTCGACCACCGCCTGGAGCGGCTCGAAGAATCAGACAATGAGATGTTCGTCGTGCAAAGAGAACATGCCGTGGCCAAGTTCGCGGGACTGGGAGGCGGTGGAGGCGATTACAAAGACGATGATGACAAGGGTGACTATAAAGACGACGATGACAAAGGGGATTACAAGGATGATGATGATAAGGGAGGCGGTGGATCAGGTGGAGGAGGTTCACTGCAGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCCCCAGGCACCAGGGCGTGATGGTGGGCATGGGTCAGAAGGATTCCTATGTGGGCGACGAGGCCCAGAGCAAGAGAGGCATCCTCACCCTGAAGTACCCCATCGAGCACGGCATCGTCACCAACTGGGACGACATGGAGAAAATCTGGCACCACACCTTCTACAATGAGCTGCGTGTGGCTCCCGAGGAGCACCCCGTGCTGCTGACCGAGGCCCCCCTGAACCCCAAGGCCAACCGCGAGAAGATGACCCAGATCATGTTTGAGACCTTCAACACCCCAGCCATGTACGTTGCTATCCAGGCTGTGCTATCCCTGTACGCCTCTGGCCGTACCACTGGCATCGTGATGGACTCCGGTGACGGGGTCACCCACACTGTGCCCATCTACGAGGGGTATGCCCTCCCCCATGCCATCCTGCGTCTGGACCTGGCTGGCCGGGACCTGACTGACTACCTCATGAAGATCCTCACCGAGCGCGGCTACAGCTTCACCACCACGGCCGAGCGGGAAATCGTGCGTGACATTAAGGAGAAGCTGTGCTACGTCGCCCTGGACTTCGAGCAAGAGATGGCCACGGCTGCTTCCAGCTCCTCCCTGGAGAAGAGCTACGAGCTGCCTGACGGCCAGGTCATCACCATTGGCAATGAGCGGTTCCGCTGCCCTGAGGCACTCTTCCAGCCTTCCTTCCTGGGCATGGAGTCCTGTGGCATCCACGAAACTACCTTCAACTCCATCATGAAGTGTGACGTGGACATCCGCAAAGACCTGTACGCCAACACAGTGCTGTCTGGCGGCACCACCATGTACCCTGGCATTGCCGACAGGATGCAGAAGGAGATCACTGCCCTGGCACCCAGCACAATGAAGATCAAGATCATTGCTCCTCCTGAGCGCAAGTACTCCGTGTGGATCGGCGGCTCCATCCTGGCCTCGCTGTCCACCTTCCAGCAGATGTGGATCAGCAAGCAGGAGTATGACGAGTCCGGCCCCTCCATCGTCCACCGCAAATGCTTCTAG'''
#the original suntag_kif has 24x Epitopes, so we are going to manually delete 14 of them with this sequence:
suntags_todelete = 'GSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGEELLSKDYHLENEVARLKKGSGSGEELLSKNYHLENEVARLKKGSGSGQ'

hairpin = 'ACGUGCCACGAUUCAACGUGGCACAGC'

bactin_only = bactin_seq[337*3:] #getting just B-actin
suntag = suntag_kif[:1012*3] #deleting KIF
suntag = suntag[:729] + suntag[1737:] #deleting 14x Suntag epitopes
flagtag = bactin_seq[:337*3] #getting just the Flag tag 


bactin_codons = [bactin_only[i:i+3] for i in range(0, len(bactin_only), 3)]  #convert b-actin to pairs of codons
suntag_codons = [suntag[i:i+3] for i in range(0, len(suntag), 3)] #suntag codons
flagtag_codons = [flagtag[i:i+3] for i in range(0, len(flagtag), 3)] #flagtag codons
hairpin_codons = [hairpin[i:i+3] for i in range(0, len(hairpin), 3)] #hairpin codons
construct = flagtag_codons + hairpin_codons + suntag_codons + bactin_codons 

construct_sequence_str = ''.join(construct)

bactin_only = bactin_seq[337*3:] #getting just B-actin
suntag = suntag_kif[:1012*3] #deleting KIF
suntag = suntag[:729] + suntag[1737:] #deleting 14x Suntag epitopes
flagtag = bactin_seq[:337*3] #getting just the Flag tag 

poi = rss.seqmanip.seq_to_protein_obj(construct_sequence_str) 
bactin = poi['1'][0] #getting the poi object
bactin.tag_length = len(flagtag_codons) + len(hairpin_codons) + len(suntag_codons) 



koff = .04
kon = .005
kin = .033
kout = 10
hair_pin_location = len(flagtag_codons)
construct_end = len(bactin.kelong)
print(construct_end)
parameters = np.array([kon, koff, kin, kout, hair_pin_location, construct_end], dtype=np.float)


forward_rates = bactin.kelong + [0,] #forward rates of the construct (with an extra 0 at the end)

#stoichiometry for the states  
stoich_states = np.array([[-1,  1,  ],
                          [ 1,  -1,  ],], dtype = np.int32)

#stoichiometry for the lattice ractions
stoich_lattice = np.zeros([2,len(forward_rates)], dtype=np.int32)
stoich_lattice[0,0] = 1
stoich_lattice[1,-1] = -1

#inital conditions
xi_states = np.array([[1,0]], dtype= np.int32)
xi_lattice = np.zeros([1,len(forward_rates)], dtype=np.int32)


additional_rules = '''
#               0      1    2     3           4             5
#parameters = [kon,  koff, kin, kout, hairpin_location, end_location]

#propensity function for hairpin model
int hairpin_location = cast_to_int(parameters[4]) # location of the hairpin
int end_loc = cast_to_int(parameters[5]) #end of the construct

if state[0] == 1:
    #if the hairpin formed, state = off
    wn[0] = parameters[0]
    step[ hairpin_location-1  ] = 0
    
if state[1] == 1:
    #if the hairpin formed, state = on
    if sum(X[hairpin_location:hairpin_location+10  ]) < 1: 
        wn[1] = parameters[1]

if free[0]:
    #if the front location is free, allow kin
    wn[2] = parameters[2]

if X[end_loc] == 1:
    # if its in the final location, add kout
    step[end_loc] = 0
    wn[3] = parameters[3]

'''

# parsed_rules = rss.model_builder.rss.().parse_rules(additional_rules)
# #mf = ModelFactory()
# #mf.compile_model('cap_ires_model', overwrite=True, rules = additional_rules)
# from models import cap_ires_model
# model = cap_ires_model.cap_ires_model

#rss.model_builder.compile_model('hairpin_model',overwrite=False, rules = additional_rules)
model = rss.model_builder.get_model('hairpin_model')


t_array = np.linspace(0,10000,10000)
seeds = np.random.randint(0,0x7fff, 1000)
print('seed:')
print(seeds[0])


t_array = np.linspace(0,4000,4001)
n_traj = 10

import time
st = time.time()


koff = .04
kon = .001
kin = .033
kout = 10
hair_pin_location = len(flagtag_codons)
construct_end = len(bactin.kelong)
print(construct_end)
parameters = np.array([kon, koff, kin, kout, hair_pin_location, construct_end], dtype=np.float)


ssa_soln = rss.solver.solve_custom_model(model, parameters, forward_rates, t_array,
                           stoich_lattice, stoich_states,
                           xi_lattice, xi_states, n_traj=1,
                           probe_loc=bactin.probe_loc, poi=None)



Ncolors = 2
# import time
# import matplotlib.pyplot as plt
# tf = t_array[-1]
# Nt = len(t_array)
# max_rib = 100 #int(len(kelong)/8 + 5)
# n_traj = 1
# st = time.time()
# all_results = np.zeros([n_traj,  Nt, max_rib,],dtype=np.intc)
# all_intensities = np.zeros([n_traj,  Nt, Ncolors,],dtype=np.intc)
# all_states = np.zeros([n_traj,  Nt, max(xi_states.shape),],dtype=np.intc)
# for i in range(0,n_traj):
#     result = np.zeros([Nt, max_rib,],dtype=np.intc)
#     intensity = np.zeros([Nt, Ncolors,],dtype=np.intc)
#     states = np.zeros([ Nt, max(xi_states.shape),],dtype=np.intc)

 
#     a = model.run_ssa_cpp(result, intensity, states, stoich_states,
#                           stoich_lattice, parameters,
#                           np.array(kelong,dtype=np.float),  t_array,
#                           xi_lattice, xi_states,probe_locations,
#                           length, seeds[i], n_total_reactions)
    
#     all_results[i] = result.reshape(max_rib,Nt).T
#     all_intensities[i] = intensity.reshape(Ncolors,Nt).T
#     all_states[i] = states.reshape(max_rib,Nt).T
    
# print('time for 100 - 1000t trajectories (c++): %f'% (time.time() - st))


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


ribosome_position_tensor = ssa_soln.ribosome_locations
intensity_tensor = ssa_soln.intensity_vec
state_tensor = ssa_soln.states
max_rib = ssa_soln.ribosome_locations.shape[-1]
Nt = ssa_soln.ribosome_locations.shape[1]
Ncolors = 2

b = ribosome_position_tensor[0]
c = state_tensor[0]
d = intensity_tensor[0]
tasep_vis = np.zeros([Nt, np.max(b)+1])
for i in range(Nt):
    tasep_vis[i, b[i,:]] = 1
    
tasep_vis[:,0] = 0
plt.rcParams['figure.dpi'] = 120
plt.figure()
#plt.matshow(tasep_vis, aspect = (np.max(b)+1)/Nt,  cmap='gray_r')
x,y = np.where(tasep_vis == 1)
plt.scatter(y.astype(int)[::1], x.astype(int)[::1],marker='.', color='gray', s=.05)

plt.gca().invert_yaxis()
plt.gca().set_aspect( (np.max(b)+1)/Nt)


plt.plot(  c[:,0]*hair_pin_location, t_array ,'r.' )

for i in (np.array([bactin.tag_epitopes['T_SunTag']])).flatten():
    plt.plot(np.ones(Nt)*i, t_array, color=colors[1],lw=.5 )
    
for i in (np.array([bactin.tag_epitopes['T_Flag']])).flatten():
    plt.plot(np.ones(Nt)*i, t_array, color=colors[1],lw=.5  )

plt.savefig('hairpin_kym.png',transparent=True)

plt.xlabel('Location')
plt.ylabel('Time')

plt.figure()

plt.plot(d[:,1]/10, )
plt.plot(d[:,0]/10, )
plt.xlabel('Time')
plt.ylabel('Intensity (UMP)')
plt.legend(['FlagTag', 'SunTag'])
plt.savefig('hairpin_int.svg')


plt.figure()
state_plot = np.zeros((1,Nt))
for i in range(Nt):
    if np.sum(c[i,:]) == 0:
    
        state_plot[0,i] = 1
    if np.sum(c[i,:]) == 1:
        if c[i,:][0] ==1:
            state_plot[0,i] = 1
        else:
            state_plot[0,i] = 0


plt.plot(state_plot[0],marker='.',lw=1)

plt.xlabel('Time')
plt.yticks([0,1], labels=['Hairpin Off','Hairpin On'])
plt.savefig('hairpin_states.svg')


'''
tf = t_array[-1]
Nt = len(t_array)
max_rib = 100 #int(len(kelong)/8 + 5)
# n_traj = 1
# st = time.time()

b = ssa_soln.ribosome_locations[0].T
c = ssa_soln.states[0].T
d = ssa_soln.intensity_vec[0].T
tasep_vis = np.zeros([Nt, np.max(b)+1])
for i in range(Nt):
    tasep_vis[i, b[:,i]] = 1
    
tasep_vis[:,0] = 0
plt.rcParams['figure.dpi'] = 300
plt.figure()
plt.matshow(tasep_vis, aspect = (np.max(b)+1)/Nt, cmap='gray_r')
x,y = np.where(tasep_vis == 1)
plt.scatter(y.astype(int)[::1], x.astype(int)[::1],marker='.', color='gray', s=.05)


#plt.plot(  c[0]*50, t_array ,'r' )
for i in (np.array([suntags['T_SunTag']])+IRES_start).flatten():
    plt.plot(np.ones(Nt)*i, t_array, color=colors[0],lw=.5 )
    
for i in (np.array([flagtags['T_Flag']])).flatten():
    plt.plot(np.ones(Nt)*i, t_array, color=colors[1],lw=.5  )

plt.xlabel('location')
plt.ylabel('time')

plt.savefig('cap_ires_kym.svg')

plt.figure()

plt.plot(d.T[:,1]/10,)
plt.plot(d.T[:,0]/10,)
plt.xlabel('time')
plt.ylabel('intensity (UMP)')
plt.legend(['IRES (SunTag)', 'CAP (FlagTag)'])
plt.savefig('cap_ires_ints.svg')

plt.figure()
state_plot = np.zeros((1,Nt))
for i in range(Nt):
    if np.sum(c[:,i]) == 0:
    
        state_plot[0,i] =0
    if np.sum(c[:,i]) == 1:
        if c[:,i][0] ==1:
            state_plot[0,i] = 1
        else:
            state_plot[0,i] = 2
    if np.sum(c[:,i]) == 2:
        state_plot[0,i] = 3

plt.plot(state_plot[0],marker='o',lw=1)

plt.xlabel('time')
plt.yticks([0,1,2,3,], labels=['Off', 'CAP on', 'IRES on', 'CAP+IRES on'])

plt.savefig('cap_ires_states.svg')


'''
