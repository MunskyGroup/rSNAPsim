# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 12:54:08 2025

@author: willi
"""



# FILE TO GENERATE ALL FIGURES FOR RSNAPSIM PAPER


import numpy as np
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
os.chdir('../../')
import rsnapsim as rsnp
os.chdir(cwd)


import numpy as np
import time

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch
import multiprocessing
from joblib import Parallel, delayed
import inspect
from itertools import product

import tqdm

##############################################################################
# Plotting style
import matplotlib.pyplot as plt
from cycler import cycler
import matplotlib.cm as cm
from matplotlib.lines import Line2D

colors = [ '#06d6a0','#ef476f', '#7400b8','#073b4c', '#118ab2',]
#colors = ['#fa8174', '#b3de69', '#bc82bd','#ccebc4','#ffed6f','#81b1d2']
font = {'family' : 'sans-serif',
        'weight' : 'bold',
        'size'   : 12}
plt.rcParams.update({'font.size': 12, 'font.weight':'bold', 'font.family':'sans-serif'  }   )
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
##############################################################################

#### SET THE RANDOM SEED
np.random.seed(42)
seed = 42

##############################################################################
#### script options and flags

global_dpi = 300

regenerate = False   # remake/rerun the simulations
resave = False       # resave the simulations and overwrite
data_save_folder = './figuredata/' # where to save the simulations

#figure 3 figsize
model_figsize = (7,5)

figure_folder = '.' #where to save the figures
figure_format = '.svg' #what format for the figures

n_model_runs = 500 #for figure 3 how many times to run the models
use_cplus = False # use c++ models vs python for generation of the data (will change the outcomes due to RNG algorithm differences)

recompile = True
##############################################################################



###############################################################################
# Experiment Design example.
###############################################################################

#######################################
#Model one, base mRNA, no changes, codon dependent tasep with stepping rates from human
# codon frequency.

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



#@title model setup
base_model = rsnp.tasep_model(mRNA,'base1') # model object


# Make the kelong mat (manually adding an extra location that is equal to zero, so particles dont run over the simulation)
kelong_mat = np.zeros([3, mRNA_length+1])
kelong_mat[0, :-1] = rsnp.propf.get_k(mRNA.nt_seq, .033, 10, 10)[1:-1]
kelong_mat[1, :-2] = rsnp.propf.get_k(mRNA.nt_seq[1:-2], .1, 10, 10)[1:-1]
kelong_mat[2, :-2] = rsnp.propf.get_k(mRNA.nt_seq[2:-1], .1, 10, 10)[1:-1]
kelong_mat[0, -2] = 0

base_model._kelong_mat = kelong_mat #override the current kelongation mat

parameters = [0.03, 10, 1000]

# ribosomal initiation
footprint = 9
# first add the reaction, in this case, we want a lattice reaction at the first
# location for a ribosome to bind (excluded)
init = lambda k,t,p,ke,o,l,pr,s,r,nr: ~np.any(l[0:0+footprint])*k[0]*(t <= 1000)
base_model.add_lattice_reaction(init, parameters, rxn_name='init', frame=0, loc=0, dexist=1,)


# Now we need a reaction for ribosomes to leave the lattice at the end (location 590)
leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[590]*k[1] #(lattice location 590 = 1) * parameter
base_model.add_lattice_reaction(leave, parameters, rxn_name='termination', frame=0, loc=590, dexist=-1,)

# DEFAULT STEPPING OF ELONGATION USING THE ELONGATION MATRIX
elongation = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+footprint]))  for i in range(nr)]

base_model.add_ribosome_reaction(elongation, 0, rxn_name='elongation', dloc=1) #default stepping

# finally specify which reactions are ribosome specific
base_model._ribosome_reactions = [2,]
base_model._constant_reactions = [0,1,]
base_model._lattice_arr0 = np.zeros([base_model._length+1], dtype=int)

#base_model._parameters[0] = .03

#######################################
#Model two, hairpin model with exclusion
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

parameters = [0.03, 10, 0.05, 0.05, 0, 10000, 200]

# ribosomal initiation
footprint = 9
# first add the reaction, in this case, we want a lattice reaction at the first
# location for a ribosome to bind (excluded)
init = lambda k,t,p,ke,o,l,pr,s,r,nr: ~np.any(l[0:0+footprint])*k[0]*(t<k[5])
hairpin_model.add_lattice_reaction(init, parameters, rxn_name='init', frame=0, loc=0, dexist=1,)


# Now we need a reaction for ribosomes to leave the lattice at the end (location 590)
leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[590]*k[1] #(lattice location 590 = 1) * parameter
hairpin_model.add_lattice_reaction(leave, parameters, rxn_name='termination', frame=0, loc=590, dexist=-1,)


hairpin_on = lambda k,t,p,ke,o,l,pr,s,r,nr: (np.sum(l[hairpin_location:hairpin_location+50]) > 0)*s[0]*k[2]
hairpin_model.add_state_reaction(hairpin_on, parameters, rxn_name='hairpin_on', inds = [0,1], dstates=[-1,1])

hairpin_off = lambda k,t,p,ke,o,l,pr,s,r,nr: s[1]*k[3]
hairpin_model.add_state_reaction(hairpin_off, parameters, rxn_name='hairpin_off', inds=[0,1], dstates=[1,-1])

# DEFAULT STEPPING OF ELONGATION USING THE ELONGATION MATRIX
elongation = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+footprint]))  for i in range(nr)]

hairpin_model.add_ribosome_reaction(elongation, parameters, rxn_name='elongation', dloc=1) #default stepping

# finally specify which reactions are ribosome specific
hairpin_model._ribosome_reactions = [4,]
hairpin_model._constant_reactions = [0,1,2,3]
hairpin_model._lattice_arr0 = np.zeros([hairpin_model._length+1], dtype=int)
hairpin_model._state_arr0 = np.array([1,0])

#######################################
#Model three, constant ribosomal drop off rate

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


dropoff_model = rsnp.tasep_model(mRNA,'dropoff') # model object
parmeters = [0.03, 10, 0.02, 10000]

# Make the kelong mat (manually adding an extra location that is equal to zero, so particles dont run over the simulation)
kelong_mat = np.zeros([3, mRNA_length+1])
kelong_mat[0, :-1] = rsnp.propf.get_k(mRNA.nt_seq, .033, 10, 10)[1:-1]
kelong_mat[1, :-2] = rsnp.propf.get_k(mRNA.nt_seq[1:-2], .1, 10, 10)[1:-1]
kelong_mat[2, :-2] = rsnp.propf.get_k(mRNA.nt_seq[2:-1], .1, 10, 10)[1:-1]
kelong_mat[0, -2] = 0

dropoff_model._kelong_mat = kelong_mat #override the current kelongation mat

# ribosomal initiation
footprint = 9
# first add the reaction, in this case, we want a lattice reaction at the first
# location for a ribosome to bind (excluded)
init = lambda k,t,p,ke,o,l,pr,s,r,nr: ~np.any(l[0:0+footprint])*k[0]*(t<k[3])
dropoff_model.add_lattice_reaction(init, parameters, rxn_name='init', frame=0, loc=0, dexist=1,)


# Now we need a reaction for ribosomes to leave the lattice at the end (location 590)
leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[590]*k[1] #(lattice location 590 = 1) * parameter
dropoff_model.add_lattice_reaction(leave, parameters, rxn_name='termination', frame=0, loc=590, dexist=-1,)

# DEFAULT STEPPING OF ELONGATION USING THE ELONGATION MATRIX
elongation = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+footprint]))  for i in range(nr)]
dropoff_model.add_ribosome_reaction(elongation, parameters, rxn_name='elongation', dloc=1) #default stepping

# Dropping off
drop_off = lambda k,t,p,ke,o,l,pr,s,r,nr: [k[2] for i in range(nr)]
dropoff_model.add_ribosome_reaction(drop_off, parameters, rxn_name='drop_off', dexist=-1)

# finally specify which reactions are ribosome specific
dropoff_model._ribosome_reactions = [2,3]
dropoff_model._constant_reactions = [0,1,]
dropoff_model._lattice_arr0 = np.zeros([dropoff_model._length+1], dtype=int)


if recompile:
    base_model.compile_model_c()
    hairpin_model.compile_model_c()
    dropoff_model.compile_model_c()
else:
    base_model.load_model_c('base1')
    hairpin_model.load_model_c('hairpin')
    dropoff_model.load_model_c('dropoff')

n_model_runs = 100
t = np.linspace(0,6000,6001)

if regenerate: 
    #for i in tqdm.tqdm(range(n_model_runs)):
    base_model_soln = rsnp.solver.solve_ssa(base_model, t, n_traj=n_model_runs, burnin=0, verbose=True)
    hairpin_model_soln = rsnp.solver.solve_ssa(hairpin_model, t, n_traj=n_model_runs, burnin=100, verbose=True)
    dropoff_model_soln = rsnp.solver.solve_ssa(dropoff_model, t, n_traj=n_model_runs, burnin=0, verbose=True)
    if resave:
        base_model_soln.save(data_save_folder + 'base_model',fmt='.npz')
        hairpin_model_soln.save(data_save_folder + 'hairpin_model',fmt='.npz')
        dropoff_model_soln.save(data_save_folder + 'dropoff_model',fmt='.npz')
else:
    
    base_model_soln = rsnp.solver.load_soln(data_save_folder + 'base_model.npz')
    hairpin_model_soln = rsnp.solver.load_soln(data_save_folder + 'hairpin_model.npz')
    dropoff_model_soln = rsnp.solver.load_soln(data_save_folder + 'dropoff_model.npz')


## Generate the heatmaps

def get_acc(intensity, FR, start, color=0):
  acc, acc_err = rsnp.inta.get_autocov(np.swapaxes(intensity[:,start::FR,:], -1,0), norm='global')
  acc, acc_err = rsnp.inta.get_autocorr(acc)
  return np.mean(acc[color,:,:],axis=-1), acc_err[color], acc

def get_mean_int(intensity, n, FR, start, stop, color=0):
  return np.mean(intensity[:n,start:stop:FR,color], axis=0), np.std(intensity[:n,start:stop:FR,color], axis=0)/np.sqrt(n)

def movmean(a, w=3):
    ret = np.cumsum(a, dtype=float)
    ret[w:] = ret[w:] - ret[:-w]
    return ret[w - 1:] / w

def get_ribosomal_occupancy(lattices, binning = 10):
  return movmean(np.mean(lattices,axis=0),w=binning)

def get_LL(acc1, acc2, acc1_err, acc2_err, pts):
  return -1/len(pts) * np.sum((acc1[pts] - acc2[pts])**2 / (np.sqrt(acc1_err[pts])*np.sqrt(acc2_err[pts])) )


# 5 experiments by 3 models

distance_metrics = np.zeros([5,3])

##########################
# Experiment 1 and 2, FCS


framerate_1 = 5      # 5 second framerate
framerate_2 = 2        # 2 second framerate
start_t = 1000
n_spots = 50

Is = [base_model_soln.I, hairpin_model_soln.I, dropoff_model_soln.I]
mean_accs_5 = []
mean_accs_2 = []
mean_accs_5b = []
mean_accs_2b = []

# get 50 ACCs from the intensities of each model for A-B comparison
for i in range(3):
    mean_accs_5.append((get_acc(Is[i][:n_spots,:,], framerate_1, start_t)))
    mean_accs_2.append((get_acc(Is[i][:n_spots,:,], framerate_2, start_t)))
    
# get 50 more ACCs from the intensities of each model for A-A comparison to normalize
for i in range(3):
    mean_accs_5b.append((get_acc(Is[i][n_spots:2*n_spots,:,], framerate_1, start_t)))
    mean_accs_2b.append((get_acc(Is[i][n_spots:2*n_spots,:,], framerate_2, start_t)))
    

# Get LL normalized metric for both frame rates
combos = [(0,1),(1,2),(0,2), (0,0), (1,1), (2,2)]
LLs = []
for i in range(len(combos)):
    LLs.append(get_LL(mean_accs_5[combos[i][0]][0], mean_accs_5b[combos[i][1]][0], mean_accs_5[combos[i][0]][1], mean_accs_5b[combos[i][1]][1], [x for x in range(1,15)] ))
    
distance_metrics[0,0] = -LLs[0] / (max(-LLs[3],-LLs[4]))
distance_metrics[0,1] = -LLs[1] / (max(-LLs[4],-LLs[5]))
distance_metrics[0,2] = -LLs[2] / (max(-LLs[3],-LLs[5]))


plt.figure()
plt.plot(mean_accs_5[0][0][:50]);
plt.plot(mean_accs_5[1][0][:50]);
plt.plot(mean_accs_5[2][0][:50]);

LLs = []
for i in range(len(combos)):
    LLs.append(get_LL(mean_accs_2[combos[i][0]][0], mean_accs_2b[combos[i][1]][0], mean_accs_2[combos[i][0]][1], mean_accs_2b[combos[i][1]][1], [x for x in range(1,15)] ))
    
distance_metrics[1,0] = -LLs[0] / (max(-LLs[3],-LLs[4]))
distance_metrics[1,1] = -LLs[1] / (max(-LLs[4],-LLs[5]))
distance_metrics[1,2] = -LLs[2] / (max(-LLs[3],-LLs[5]))


####################
# Experiment 3, mean intensity differences

mean_ints = []
mean_intsb = []
for i in range(3):
    mean_ints.append((get_mean_int(Is[i][:n_spots], 100, 5, 1400,1800)))
    mean_intsb.append((get_mean_int(Is[i][n_spots:2*n_spots], 100, 5, 1400,1800)))
LLs = []
for i in range(len(combos)):
    LLs.append(get_LL(mean_ints[combos[i][0]][0], mean_intsb[combos[i][1]][0], mean_ints[combos[i][0]][1], mean_intsb[combos[i][1]][1], np.arange(10,40) ))
    

distance_metrics[2,0] = LLs[0] / (max(LLs[3],LLs[4]))
distance_metrics[2,1] = LLs[1] / (max(LLs[4],LLs[5]))
distance_metrics[2,2] = LLs[2] / (max(LLs[3],LLs[5]))
