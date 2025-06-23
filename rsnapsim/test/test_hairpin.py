# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 13:46:20 2025

@author: willi
"""



import sys, os
#print(os.listdir('.'))
print(os.getcwd())
sys.path.append('.')
cwd = os.getcwd()
os.chdir('../..')
import rsnapsim as rsnp
os.chdir(cwd)
import numpy as np

import multiprocessing
import time
import inspect


example_mRNA = '''ATGGCGAACCTTGGCTGCTGGATGCTGGTTCTCTTTGTGGCCACATGGAGTGACCTGGGC
                    CTCTGCAAGAAGCGCCCGAAGCCTGGAGGATGGAACACTGGGGGCAGCCGATACCCGGGG
                    CAGGGCAGCCCTGGAGGCAACCGCTACCCACCTCAGGGCGGTGGTGGCTGGGGGCAGCCT
                    CATGGTGGTGGCTGGGGGCAGCCTCATGGTGGTGGCTGGGGGCAGCCCCATGGTGGTGGC
                    TGGGGACAGCCTCATGGTGGTGGCTGGGGTCAAGGAGGTGGCACCCACAGTCAGTGGAAC
                    AAGCCGAGTAAGCCAAAAACCAACATGAAGCACATGGCTGGTGCTGCAGCAGCTGGGGCA
                    GTGGTGGGGGGCCTTGGCGGCTACATGCTGGGAAGTGCCATGAGCAGGCCCATCATACAT
                    TTCGGCAGTGACTATGAGGACCGTTACTATCGTGAAAACATGCACCGTTACCCCAACCAA
                    GTGTACTACAGGCCCATGGATGAGTACAGCAACCAGAACAACTTTGTGCACGACTGCGTC
                    AATATCACAATCAAGCAGCACACGGTCACCACAACCACCAAGGGGGAGAACTTCACCGAG
                    ACCGACGTTAAGATGATGGAGCGCGTGGTTGAGCAGATGTGTATCACCCAGTACGAGAGG
                    GAATCTCAGGCCTATTACCAGAGAGGATCGAGCATGGTCCTCTTCTCCTCTCCACCTGTG
                    ATCCTCCTGATCTCTTTCCTCATCTTCCTGATAGTGGGATGA'''.replace('\n','').replace(' ','')

# Make the mRNA object
poi = rsnp.seqmanip.seq_to_CDS_obj(example_mRNA) # convert a given sequence to a protein of interest object
mRNA = poi['0'][0]                              # pull out the main open reading frame
mRNA_length = len(mRNA.kelong)                  # get the length of the mRNA
mRNA.generate_3frame_tags()                     # call to generate all open reading frames

# manually adding epitopes for two colors
mRNA.multiframe_epitopes[0] = {'T_Flag': [1, 10, 19, 195, 205, 217, 227, 299, 308, 317], 'T_HA':[400,410,420,430]}



hairpin_model = rsnp.tasep_model(mRNA,'hairpin') # model object

hairpin_location = 200

# Make the kelong mat (manually adding an extra location that is equal to zero, so particles dont run over the simulation)
kelong_mat = np.zeros([3, mRNA_length+1])
kelong_mat[0, :-1] = rsnp.propf.get_k(mRNA.nt_seq, .033, 10, 10)[1:-1]
kelong_mat[1, :-2] = rsnp.propf.get_k(mRNA.nt_seq[1:-2], .1, 10, 10)[1:-1]
kelong_mat[2, :-2] = rsnp.propf.get_k(mRNA.nt_seq[2:-1], .1, 10, 10)[1:-1]
kelong_mat[0, -2] = 0

hairpin_model._kelong_mat = kelong_mat #override the current kelongation mat

hairpin_model.add_states(2, state0=[1,0], names=['off','on'])

parameters = [0.03, 10, 0.05, 0.05, 0, 10000, hairpin_location]

# ribosomal initiation
footprint = 9
# first add the reaction, in this case, we want a lattice reaction at the first
# location for a ribosome to bind (excluded)
init = lambda k,t,p,ke,o,l,pr,s,r,nr: ~np.any(l[0:0+footprint])*k[0]*(t<k[5])
hairpin_model.add_lattice_reaction(init, parameters, rxn_name='init', frame=0, loc=0, dexist=1,)


# Now we need a reaction for ribosomes to leave the lattice at the end (location 590)
leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[590]*k[1] #(lattice location 590 = 1) * parameter
hairpin_model.add_lattice_reaction(leave, parameters, rxn_name='termination', frame=0, loc=590, dexist=-1,)


hairpin_on = lambda k,t,p,ke,o,l,pr,s,r,nr: (sum(l[200:(200+50)]) == 0)*s[0]*k[2]
hairpin_model.add_state_reaction(hairpin_on, parameters, rxn_name='hairpin_on', inds = [0,1], dstates=[-1,1])

hairpin_off = lambda k,t,p,ke,o,l,pr,s,r,nr: s[1]*k[3]
hairpin_model.add_state_reaction(hairpin_off, parameters, rxn_name='hairpin_off', inds=[0,1], dstates=[1,-1])

# DEFAULT STEPPING OF ELONGATION USING THE ELONGATION MATRIX
elongation = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+footprint])) for i in range(nr)]

hairpin_model.add_ribosome_reaction(elongation, parameters, rxn_name='elongation', dloc=1) #default stepping

# finally specify which reactions are ribosome specific
hairpin_model._ribosome_reactions = [4,]
hairpin_model._constant_reactions = [0,1,2,3]
hairpin_model._lattice_arr0 = np.zeros([hairpin_model._length+1], dtype=int)
hairpin_model._state_arr0[0] = 1

hairpin_model.compile_model_c()

#hairpin_model.load_model_c('hairpin')

n_model_runs = 1
t = np.linspace(0,1000,1001)
hairpin_model_soln = rsnp.solver.solve_ssa(hairpin_model, t, n_traj=n_model_runs, burnin=0, verbose=True, cplus=True)

import matplotlib.pyplot as plt
plt.plot(hairpin_model_soln.state_array[0])