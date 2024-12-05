# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 15:10:10 2024

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
poi = rsnp.seqmanip.seq_to_protein_obj(example_mRNA) # convert a given sequence to a protein of interest object
mRNA = poi['0'][0]                              # pull out the main open reading frame
mRNA_length = len(mRNA.kelong)                  # get the length of the mRNA
mRNA.generate_3frame_tags()                     # call to generate all open reading frames

# manually adding epitopes for two colors
mRNA.multiframe_epitopes[0] = {'T_Flag': [1, 10, 19, 195, 205, 217, 227, 299, 308, 317], 'T_HA':[400,410,420,430]}

#mRNA.visualize_mrna_strand() # plot open reading frame 1

model = rsnp.tasep_model(mRNA, 'base') # model object




# Make the kelong mat (manually adding an extra location that is equal to zero, so particles dont run over the simulation)
kelong_mat = np.zeros([3, mRNA_length+1])
kelong_mat[0, :-1] = rsnp.propf.get_k(mRNA.nt_seq, .033, 10, 10)[1:-1]
kelong_mat[1, :-2] = rsnp.propf.get_k(mRNA.nt_seq[1:-2], .1, 10, 10)[1:-1]
kelong_mat[2, :-2] = rsnp.propf.get_k(mRNA.nt_seq[2:-1], .1, 10, 10)[1:-1]
kelong_mat[0, -2] = 0

model._kelong_mat = kelong_mat #override the current kelongation mat

footprint = 9
parameters = [0.03, 10, .4, 0]

init = lambda k,t,p,ke,o,l,pr,s,r,nr: (1-np.any(l[0:0+footprint]))*k[0]
model.add_lattice_reaction(init, parameters, rxn_name = 'initiation', exclusion=1, frame=0, loc=0, dexist=1,)


leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[L-1]*k[1] #(lattice location 590 = 1) * parameter
model.add_lattice_reaction(leave, parameters, rxn_name='termination', exclusion=0, frame=0, loc=590, dexist=-1,)

elong_step = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+footprint])) for i in range(nr)]
model.add_ribosome_reaction(elong_step, parameters, rxn_name='elongation', exclusion=1, dloc=1) 

model._ribosome_reactions = [2]
model._constant_reactions = [0,1,]
model._lattice_arr0 = np.zeros([model._length+1], dtype=int)


model.compile_model_c()

solver = rsnp.Solver.Solver()
t = np.linspace(0,500,5001)
for i in range(1):
    soln = solver.solve_ssa(model, t, cplus=True)



