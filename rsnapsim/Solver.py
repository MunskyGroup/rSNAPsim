# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:52:28 2020

@author: willi
"""
import numpy as np
from . import CodonDictionaries
from . import translation_models as models
from . import ODE_Soln
from . import SSA_Soln
from . import custom_errors as custom_err

from . import seqmanip
from . import propf
from . import TranslationModel


SSA_Soln = SSA_Soln.SSA_Soln
ODE_Soln = ODE_Soln.ODE_Soln



import tqdm as tqdm
import time
import warnings
import os
import multiprocessing
from joblib import Parallel, delayed
import inspect

## Search locally for ssa_cpp
path_to_cpp = ''
path_to_gen = ''
path_to_trna = ''

file_path = os.path.dirname(os.path.realpath(__file__))

for root, dirs, files in os.walk(file_path, topdown=False):
   for branch in dirs:

       if 'ssa_cpp' in branch:
           path_to_cpp = os.path.join(root, branch)

       if 'generalized_cpp' in branch:
           path_to_gen = os.path.join(root, branch)
       if 'trna_ssa' in branch:
           path_to_trna = os.path.join(root, branch)
if path_to_cpp != '':
    try:
        cwd = os.getcwd()
        os.chdir(path_to_cpp)

        #import ssa_translation
        import ssa_translation_lowmem
        os.chdir(cwd)
    except:
        os.chdir(cwd)


if path_to_gen != '':
    try:
        cwd = os.getcwd()

        os.chdir(path_to_gen)
        import ssa_translation_lowmem
        os.chdir(cwd)
    except:
        os.chdir(cwd)

if path_to_trna !='':
    try:
        cwd = os.getcwd()

        os.chdir(path_to_trna)
        import ssa_trna
        os.chdir(cwd)
    except:
        os.chdir(cwd)

try:
    import ssa_translation_lowmem
except:
    pass


class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class CythonMissingError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
class Solver():
    def __init__(self, NUMBER_OF_CORES=None):
        if NUMBER_OF_CORES == None:
            self.NUMBER_OF_CORES = int(multiprocessing.cpu_count()/2) #claim half the cores
        else:
            self.NUMBER_OF_CORES = 1
        self.default_footprint_size = 9


    def load_soln(self, filename):
        
        if filename[-4:] == '.npz':
            
            loaded = np.load(filename)
            solve_time, n_traj, L, n_colors, burnin = loaded['constants']
            
            n_traj = int(n_traj)
            L = int(L)
            n_colors = int(n_colors)

            ribosome_array = loaded['ribosome_array']
            state_array = loaded['state_array']
            resource_array = loaded['resource_array']
            kelong_mat = loaded['kelong_mat']
            t = loaded['t']
            
            soln =  CustomSSASoln([0], n_colors, ribosome_array, state_array,
                                  resource_array, t, burnin, n_traj, solve_time) 
            
            #overwrite these entries since we initialized the object with a blank list.
            soln.L = L
            soln.kelong_mat = kelong_mat
            
            return soln
            
        
        if filename[-4:] == '.npy':
            linear_array = np.load(filename)
            barcode_size = int(linear_array[0]) #read the size of every compressed array
            
            # decompress the constant values
            solve_time = linear_array[1+barcode_size]
            n_traj = int(linear_array[1+barcode_size+1])
            L = int(linear_array[1+barcode_size+2])
            n_colors = int(linear_array[1+barcode_size+3])
            burnin = linear_array[1+barcode_size+4]
            
            start = 1+barcode_size+4+1
            barcode = linear_array[1:barcode_size+1].astype(int)
            ribosome_array_shape = barcode[1:5]
            state_array_shape = barcode[5:8]
            resource_array_shape = barcode[8:11]
            kelong_mat_shape = barcode[11:13]
            time_shape = barcode[13]
            
            ribosome_array = linear_array[start:ribosome_array_shape[0]*ribosome_array_shape[1]*ribosome_array_shape[2]*ribosome_array_shape[3]+start].reshape(ribosome_array_shape)
            
            start = ribosome_array_shape[0]*ribosome_array_shape[1]*ribosome_array_shape[2]*ribosome_array_shape[3]+start
            state_array = linear_array[start:state_array_shape[0]*state_array_shape[1]*state_array_shape[2]+start].reshape(state_array_shape)
            
            
            if 0 not in state_array_shape:
                start = 1+state_array_shape[0]*state_array_shape[1]*state_array_shape[2]+start
            else:
                start = state_array_shape[0]*state_array_shape[1]*state_array_shape[2]+start
            resource_array = linear_array[start:resource_array_shape[0]*resource_array_shape[1]*resource_array_shape[2]+start].reshape(resource_array_shape)

            if 0 not in resource_array_shape:
                start = 1+resource_array_shape[0]*resource_array_shape[1]*resource_array_shape[2]+start
            else:
                start = resource_array_shape[0]*resource_array_shape[1]*resource_array_shape[2]+start
            kelong_mat = linear_array[start:kelong_mat_shape[0]*kelong_mat_shape[1]+start].reshape(kelong_mat_shape)
            
            if 0 not in kelong_mat_shape:
                start = kelong_mat_shape[0]*kelong_mat_shape[1]+start
            else:
                start = kelong_mat_shape[0]*kelong_mat_shape[1]+start
            t = linear_array[start:time_shape+start].reshape(time_shape)
            
                        
            soln =  CustomSSASoln([0], n_colors, ribosome_array, state_array,
                                  resource_array, t, burnin, n_traj, solve_time) 
            
            #overwrite these entries since we initialized the object with a blank list.
            soln.L = L
            soln.kelong_mat = kelong_mat
            
        return soln

    def solve_ssa(self, mRNA_model, t, n_traj=1, burnin=0, seed=None, parallel=False, cplus=False, cores=4, 
                  probe_list=None, ki=None, kt=None, verbose=False):
        
        # check if the user passed a model object or an mRNA object
        
        
        
        st = time.time()
        if seed == None:
            seeds = np.random.randint(0x7FFFFF, size=n_traj)
        else:
            np.random.seed(seed)
            seeds = np.random.randint(0x7FFFFF, size=n_traj)
        
        if not cplus:
            
            constants = self.__setup_python_sim(mRNA_model, probe_list, ki, kt)
            if parallel:
                NUMBER_OF_CORES = cores
                solns = Parallel(n_jobs = NUMBER_OF_CORES, prefer="threads")(delayed(self.__run)(*args, **kwargs)
                                                                             for args, kwargs in ([[(*constants, t), {'burnin':burnin, 'seed':seeds[i]}] for i in range(0,n_traj)]))
    
            else:
                solns = []
                if verbose:
                    for i in tqdm.tqdm(range(n_traj)):
                        solns.append(self.__run(*constants, t, burnin, seeds[i]))
                else:
                    for i in range(n_traj):
                        solns.append(self.__run(*constants, t, burnin, seeds[i]))
                
            rib_array = np.array([solns[i][0] for i in range(len(solns))])
            resource_array = np.array([solns[i][2] for i in range(len(solns))])
            state_array = np.array([solns[i][1] for i in range(len(solns))])
            solve_time = time.time() - st
            n_colors = constants[4]
            soln = CustomSSASoln(mRNA_model, n_colors, rib_array, state_array, resource_array, t, burnin, n_traj, solve_time) 
            
            
        if cplus:
            
            # see if we were given a poi object, a list, or a custom model object
            try:
                mRNA_model.cmodel
                is_custom_model = True
            except:
                is_custom_model = False
            
            if is_custom_model:
                # C++ solver for custom TASEP model objects
                # most of their information is stored in the model object.
                parameters = mRNA_model._parameters
                pars = []
                npars = []
                for i in range(len(parameters)):
                    if isinstance(parameters[i], list):
                        npars.append(len(parameters[i]))
                        for j in range(len(parameters[i])):
                            pars.append(parameters[i][j])
                    else:
                        npars.append(1)
                        pars.append(parameters[i])
                    
                
                if len(mRNA_model._state_arr0) == 0:
                    temp_state_arr0 = np.zeros([1], dtype=np.int32)
                else:
                    temp_state_arr0 = mRNA_model._state_arr0
                    
                if len(mRNA_model._resource_arr0) == 0:
                    temp_resource_arr0 = np.zeros([1], dtype=np.int32)
                else:
                    temp_resource_arr0 = mRNA_model._resource_arr0
                n_states = 0
                n_resources = 0 
                
                st = time.time()
                def convert_c_pa(pa):
                    return pa.swapaxes(1,0).reshape([1,len(t), mRNA_model._rib_arr0.shape[1], mRNA_model._rib_arr0.shape[0]]).swapaxes(-1,-2)
                
                
                ribosome_array = np.zeros([n_traj, len(t), mRNA_model._rib_arr0.shape[0], mRNA_model._rib_arr0.shape[1] ])
                state_array = np.zeros([n_traj, len(t), max(len(mRNA_model._state_arr0),1)]) # minimum of shape one for C++, wont allow 0 shaped arrays
                resource_array = np.zeros([n_traj, len(t),  max(len(mRNA_model._resource_arr0),1)])
                
                for i in range(n_traj):
                    seed = seeds[i]
    
                    pa = np.zeros([ mRNA_model._rib_arr0.shape[0]*mRNA_model._rib_arr0.shape[1], len(t)], dtype=np.int32, order='C')
                    sa = np.zeros([len(t), max(len(mRNA_model._state_arr0),1) ], dtype=np.int32)
                    ra = np.zeros([len(t), max(len(mRNA_model._resource_arr0),1)], dtype=np.int32)
                    
                    mRNA_model.cmodel.run_ssa_cpp(pa, sa, ra, mRNA_model._rib_arr0, temp_state_arr0, temp_resource_arr0,
                                               mRNA_model._rxn_mat.astype(np.int32), np.array(mRNA_model._constant_reactions + mRNA_model._ribosome_reactions),
                                               mRNA_model._kelong_mat,
                                               mRNA_model._probe_mat.astype(np.int32), np.array(pars),
                                               t,
                                               mRNA_model._rib_arr0.shape[0], mRNA_model._n_states, mRNA_model._n_resources, len(mRNA_model._constant_reactions),
                                               len(mRNA_model._ribosome_reactions),
                                               burnin, seed, )
                    
                    ribosome_array[i] = convert_c_pa(pa)
                    state_array[i] = np.copy(sa)
                    resource_array[i] = np.copy(ra)
                    
                    #print(pa)
                solve_time = time.time() - st
                soln = CustomSSASoln(mRNA_model, 0, ribosome_array, state_array, resource_array, t, burnin, n_traj, solve_time)         
            
            else:
                # C++ SOLVING WHEN GIVEN A LIST OR CDS OBJECT
                # if the use gives only a list, or a cds object to solve taseps, make a new blank model
                # and load the default compiled c++ model to use.
                
                
                # make a new default_model and use it to load the precompiled c++ one
                mRNA_model = self.__cds_or_list_to_mRNA_model(mRNA_model, probe_list, ki, kt)
                mRNA_model.load_model_c('default') #MUST BE NAMED DEFAULT FOR NOW
                L = mRNA_model._length
                parameters = [ki, kt, L]
                
                pars = []
                npars = []
                for i in range(len(parameters)):
                    if isinstance(parameters[i], list):
                        npars.append(len(parameters[i]))
                        for j in range(len(parameters[i])):
                            pars.append(parameters[i][j])
                    else:
                        npars.append(1)
                        pars.append(parameters[i])
                    
                
                if len(mRNA_model._state_arr0) == 0:
                    temp_state_arr0 = np.zeros([1], dtype=np.int32)
                else:
                    temp_state_arr0 = mRNA_model._state_arr0
                    
                if len(mRNA_model._resource_arr0) == 0:
                    temp_resource_arr0 = np.zeros([1], dtype=np.int32)
                else:
                    temp_resource_arr0 = mRNA_model._resource_arr0
                n_states = 0
                n_resources = 0 
                
                st = time.time()
                def convert_c_pa(pa):
                    return pa.swapaxes(1,0).reshape([1,len(t), mRNA_model._rib_arr0.shape[1], mRNA_model._rib_arr0.shape[0]]).swapaxes(-1,-2)
                
                
                ribosome_array = np.zeros([n_traj, len(t), mRNA_model._rib_arr0.shape[0], mRNA_model._rib_arr0.shape[1] ])
                state_array = np.zeros([n_traj, len(t), max(len(mRNA_model._state_arr0),1)]) # minimum of shape one for C++, wont allow 0 shaped arrays
                resource_array = np.zeros([n_traj, len(t),  max(len(mRNA_model._resource_arr0),1)])
                
                for i in range(n_traj):
                    seed = seeds[i]
    
                    pa = np.zeros([ mRNA_model._rib_arr0.shape[0]*mRNA_model._rib_arr0.shape[1], len(t)], dtype=np.int32, order='C')
                    sa = np.zeros([len(t), max(len(mRNA_model._state_arr0),1) ], dtype=np.int32)
                    ra = np.zeros([len(t), max(len(mRNA_model._resource_arr0),1)], dtype=np.int32)
                    
                    mRNA_model.cmodel.run_ssa_cpp(pa, sa, ra, mRNA_model._rib_arr0, temp_state_arr0, temp_resource_arr0,
                                               mRNA_model._rxn_mat.astype(np.int32), np.array(mRNA_model._constant_reactions + mRNA_model._ribosome_reactions),
                                               mRNA_model._kelong_mat,
                                               mRNA_model._probe_mat.astype(np.int32), np.array(pars),
                                               t,
                                               mRNA_model._rib_arr0.shape[0], mRNA_model._n_states, mRNA_model._n_resources, len(mRNA_model._constant_reactions),
                                               len(mRNA_model._ribosome_reactions),
                                               burnin, seed, )
                    
                    ribosome_array[i] = convert_c_pa(pa)
                    state_array[i] = np.copy(sa)
                    resource_array[i] = np.copy(ra)
                    
                    #print(pa)
                solve_time = time.time() - st
                soln = CustomSSASoln(mRNA_model, 0, ribosome_array, state_array, resource_array, t, burnin, n_traj, solve_time)         

        return soln
    

        
        
    def __initalize_trajs(self, mRNA_model, n_colors):

        try:
            particle_size = mRNA_model.particle_size
            is_model_obj = True
        except:
            is_model_obj = False
        
        if is_model_obj:
            rib_arr, lattice_arr, state_arr, resource_arr = mRNA_model._x0()
            NR = np.sum(lattice_arr)
            occupied = rib_arr[:,3]
        else:
            if isinstance(mRNA_model, list):
                particle_size = self.default_footprint_size 
                rib_arr = np.zeros([int(len(mRNA_model)/particle_size)+5, 4 + 3 + n_colors], dtype=np.int32)
                lattice_arr = np.zeros(len(mRNA_model), dtype=np.int32)
                state_arr = np.array([],dtype=np.int32)
                resource_arr = np.array([],dtype=np.int32)
                
                NR = np.sum(lattice_arr)
                occupied = rib_arr[:,3]
            else:
                particle_size = self.default_footprint_size 
                n_colors = int(len(mRNA_model.tag_epitopes))
                rib_arr = np.zeros([int(mRNA_model.total_length/particle_size)+5, 4 + 3 + n_colors], dtype=np.int32)
                lattice_arr = np.zeros(len(mRNA_model.kelong), dtype=np.int32)
                state_arr = np.array([],dtype=np.int32)
                resource_arr = np.array([],dtype=np.int32)
                
                NR = np.sum(lattice_arr)
                occupied = rib_arr[:,3]
            
        return np.copy(rib_arr), np.copy(lattice_arr), np.copy(state_arr), np.copy(resource_arr), np.copy(occupied), NR
        
    def __constants(self, mRNA_model, probe_list, ki, kt):
        #if the user passes an object that is not a model, use the default model
        #otherwise, read the constants stored in the model.
        try:
            particle_size = mRNA_model.particle_size
            is_model_obj = True
        except:
            is_model_obj = False
            
        if is_model_obj:
            particle_size = mRNA_model.particle_size    
            rxn_mat = mRNA_model._rxn_mat.astype(int)
            n_colors = int(np.max(mRNA_model._probe_mat))
            n_states = mRNA_model._n_states
            n_resources = mRNA_model._n_resources
            n_rxns = int(rxn_mat.shape[0])
            n_ribosome_reactions = np.sum(rxn_mat[:,0] == 0)
            max_rib = mRNA_model._max_particles
            n_constant_reactions = len(rxn_mat) - n_ribosome_reactions
            L = mRNA_model._kelong_mat.shape[1]
            kelong_mat = mRNA_model._kelong_mat
            kelong_mat[0,-1] = 0
            probe_mat = mRNA_model._probe_mat.astype(int)
            
            constant_props = [mRNA_model._propensities[i] for i in mRNA_model._constant_reactions]
            constant_parameters = [mRNA_model._parameters[i] for i in mRNA_model._constant_reactions]
    
            ribosome_props = [mRNA_model._propensities[i] for i in mRNA_model._ribosome_reactions]
            ribosome_parameters = [mRNA_model._parameters[i] for i in mRNA_model._ribosome_reactions]
            reaction_ids = mRNA_model._constant_reactions + mRNA_model._ribosome_reactions 
            probe_function = mRNA_model._probe_function
            probe_parameters = mRNA_model._probe_parameters
        else:
            if isinstance(mRNA_model, list):
                #overwrite to a default model if just a stepping list was used
                particle_size = self.default_footprint_size 
                n_colors = int(len(probe_list)) # number of color tags
                n_states = 0
                n_resources = 0
    
                max_rib = int(len(mRNA_model)/particle_size)+5 #total number of ribosomes possible
                
    
                L = len(mRNA_model)
                #build the default reaction matrix
                rxn_mat = np.array([[  2.,   1.,   0.,   0.,   1.,   0.,   0.,   0.,   0.],  #initiation
                                    [  2.,   0.,   0., L-1,  -1.,   0.,   0.,   0.,   0.],     #termination
                                    [  0.,   1.,   0.,   0.,   0.,   0.,   1.,   0.,   0.]]).astype(np.int32) #stepping
    
                n_rxns = int(rxn_mat.shape[0])
                n_ribosome_reactions = np.sum(rxn_mat[:,0] == 0)
                n_constant_reactions = len(rxn_mat) - n_ribosome_reactions
                #kelong_mat = mRNA_model.kelong_mat
                kelong_mat = np.zeros([3,L])
                kelong_mat[0,:] = mRNA_model
                kelong_mat[0,-1] = 0
                #probe_mat = mRNA_model.probe_mat
                probe_mat = np.zeros([3,L])
                k = 1
                for i in range(len(probe_list)):
                    probe_mat[0,probe_list[i]] = k
                    k+=1
                
                
                init = lambda k,t,p,ke,o,l,pr,s,r,nr: (1-np.any(l[0:0+int(particle_size)]))*k[0]
                leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[L-1]*k[1]
                elong_step = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+particle_size])) for i in range(nr)]
                probe_function = lambda k,t,p,ke,o,l,pr,s,r,nr: 1
                constant_props = [init, leave]
                constant_parameters = [[ki], [0,kt]]
                probe_parameters = 0
                ribosome_props = [elong_step]
                ribosome_parameters = [[1]]
                reaction_ids =[0,1,2]                
            else:
                #overwrite to a default model if just a protein obj was passed
                particle_size = self.default_footprint_size 
                n_colors = int(len(mRNA_model.tag_epitopes)) # number of color tags
                n_states = 0
                n_resources = 0
    
                max_rib = int(mRNA_model.total_length/particle_size)+5 #total number of ribosomes possible
                
    
                L = len(mRNA_model.kelong)
                #build the default reaction matrix
                rxn_mat = np.array([[  2.,   1.,   0.,   0.,   1.,   0.,   0.,   0.,   0.],  #initiation
                                    [  2.,   0.,   0., L-1,  -1.,   0.,   0.,   0.,   0.],     #termination
                                    [  0.,   1.,   0.,   0.,   0.,   0.,   1.,   0.,   0.]]).astype(np.int32) #stepping
    
                n_rxns = int(rxn_mat.shape[0])
                n_ribosome_reactions = np.sum(rxn_mat[:,0] == 0)
                n_constant_reactions = len(rxn_mat) - n_ribosome_reactions
                kelong_mat = mRNA_model.kelong_mat
                kelong_mat[0,-1] = 0
                probe_mat = mRNA_model.probe_mat
                
                
                init = lambda k,t,p,ke,o,l,pr,s,r,nr: (1-np.any(l[0:0+int(particle_size)]))*k[0]
                leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[L-1]*k[1]
                elong_step = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+particle_size])) for i in range(nr)]
                probe_function = lambda k,t,p,ke,o,l,pr,s,r,nr: 1
                constant_props = [init, leave]
                constant_parameters = [[mRNA_model.ki], [0,mRNA_model.kt]]
                probe_parameters = 0
                ribosome_props = [elong_step]
                ribosome_parameters = [[mRNA_model.ke_mu]]
                reaction_ids =[0,1,2]
            
        return particle_size, max_rib, rxn_mat, n_colors, n_states, n_rxns, n_resources, n_constant_reactions, n_ribosome_reactions, L, kelong_mat, probe_mat, constant_props, constant_parameters, ribosome_props, ribosome_parameters, probe_function, probe_parameters, reaction_ids

    def __cds_or_list_to_mRNA_model(self, mRNA_model, probe_list, ki, kt):
        
                                 
        if isinstance(mRNA_model, list):
            
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
            poi = seqmanip.seq_to_CDS_obj(example_mRNA) # convert a given sequence to a protein of interest object
            mRNA = poi['0'][0]                              # pull out the main open reading frame
            mRNA_length = len(mRNA.kelong)                  # get the length of the mRNA
            mRNA.generate_3frame_tags()                     # call to generate all open reading frames
            
            # manually adding epitopes for two colors
            mRNA.multiframe_epitopes[0] = {'T_Flag': [1, 10, 19, 195, 205, 217, 227, 299, 308, 317], 'T_HA':[400,410,420,430]}
            
                        # Make the mRNA object
            #poi = seqmanip.seq_to_CDS_obj(example_mRNA) # convert a given sequence to a protein of interest object
            #mRNA = poi['0'][0]                              # pull out the main open reading frame
            mRNA_length = len(mRNA_model)                  # get the length of the mRNA
            #mRNA.generate_3frame_tags()                     # call to generate all open reading frames
            
            # manually adding epitopes for two colors
            #mRNA.multiframe_epitopes[0] = {'T_Flag': [1, 10, 19, 195, 205, 217, 227, 299, 308, 317], 'T_HA':[400,410,420,430]}
            
            
            model = TranslationModel.TranslationModel(mRNA, 'default') # model object
            
            
            # Make the kelong mat (manually adding an extra location that is equal to zero, so particles dont run over the simulation)
            kelong_mat = np.zeros([3, mRNA_length+1])
            kelong_mat[0, :-1] = mRNA_model
            #kelong_mat[1, :-2] = propf.get_k(mRNA.nt_seq[1:-2], .1, 10, 10)[1:-1]
            #kelong_mat[2, :-2] = propf.get_k(mRNA.nt_seq[2:-1], .1, 10, 10)[1:-1]
            kelong_mat[0, -2] = 0
            
            model._kelong_mat = kelong_mat #override the current kelongation mat
            
            footprint = 9
            parameters = [0.03, 10, 591]
            
            # DEFAULT STEPPING INITIATION
            init = lambda k,t,p,ke,o,l,pr,s,r,nr: (1-np.any(l[0:0+footprint]))*k[0]
            model.add_lattice_reaction(init, parameters, rxn_name = 'initiation', exclusion=1, frame=0, loc=0, dexist=1,)
            
            
            # DEFAULT TERMINATION 
            leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[k[2]-1]*k[1] #(lattice location 590 = 1) * parameter
            model.add_lattice_reaction(leave, parameters, rxn_name='termination', exclusion=0, frame=0, loc=590, dexist=-1,)
            
            
            # DEFAULT STEPPING OF ELONGATION USING THE ELONGATION MATRIX
            #default stepping
            elong_step = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+footprint])) for i in range(nr)]
            model.add_ribosome_reaction(elong_step, parameters, rxn_name='elongation', exclusion=1, dloc=1) 
            
            
            # finally specify which reactions are ribosome specific
            model._ribosome_reactions = [2]
            model._constant_reactions = [0,1,]
            model._lattice_arr0 = np.zeros([model._length+1], dtype=int)

            return model

            # #overwrite to a default model if just a stepping list was used
            # particle_size = self.default_footprint_size 
            # n_colors = int(len(probe_list)) # number of color tags
            # n_states = 0
            # n_resources = 0

            # max_rib = int(len(mRNA_model)/particle_size)+5 #total number of ribosomes possible
            

            # L = len(mRNA_model)
            # #build the default reaction matrix
            # rxn_mat = np.array([[  2.,   1.,   0.,   0.,   1.,   0.,   0.,   0.,   0.],  #initiation
            #                     [  2.,   0.,   0., L-1,  -1.,   0.,   0.,   0.,   0.],     #termination
            #                     [  0.,   1.,   0.,   0.,   0.,   0.,   1.,   0.,   0.]]).astype(np.int32) #stepping

            # n_rxns = int(rxn_mat.shape[0])
            # n_ribosome_reactions = np.sum(rxn_mat[:,0] == 0)
            # n_constant_reactions = len(rxn_mat) - n_ribosome_reactions
            # #kelong_mat = mRNA_model.kelong_mat
            # kelong_mat = np.zeros([3,L])
            # kelong_mat[0,:] = mRNA_model
            # kelong_mat[0,-1] = 0
            # #probe_mat = mRNA_model.probe_mat
            # probe_mat = np.zeros([3,L])
            # k = 1
            # for i in range(len(probe_list)):
            #     probe_mat[0,probe_list[i]] = k
            #     k+=1
            
            
            # init = lambda k,t,p,ke,o,l,pr,s,r,nr: (1-np.any(l[0:0+int(particle_size)]))*k[0]
            # leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[L-1]*k[1]
            # elong_step = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+particle_size])) for i in range(nr)]
            # probe_function = lambda k,t,p,ke,o,l,pr,s,r,nr: 1
            # constant_props = [init, leave]
            # constant_parameters = [[ki], [0,kt]]
            # probe_parameters = 0
            # ribosome_props = [elong_step]
            # ribosome_parameters = [[1]]
            # reaction_ids =[0,1,2]                
        else:
            
            mRNA_length = len(mRNA_model.kelong)                  # get the length of the mRNA
            mRNA_model.generate_3frame_tags()                     # call to generate all open reading frames

            model = TranslationModel.TranslationModel(mRNA_model, 'default') # model object
            
            
            # Make the kelong mat (manually adding an extra location that is equal to zero, so particles dont run over the simulation)
            kelong_mat = np.zeros([3, mRNA_length+1])
            kelong_mat[0, :-1] = propf.get_k(mRNA.nt_seq, .033, 10, 10)[1:-1]
            kelong_mat[1, :-2] = propf.get_k(mRNA.nt_seq[1:-2], .1, 10, 10)[1:-1]
            kelong_mat[2, :-2] = propf.get_k(mRNA.nt_seq[2:-1], .1, 10, 10)[1:-1]
            kelong_mat[0, -2] = 0
            
            model._kelong_mat = kelong_mat #override the current kelongation mat
            
            footprint = 9
            parameters = [0.03, 10, 591]
            
            # DEFAULT STEPPING INITIATION
            init = lambda k,t,p,ke,o,l,pr,s,r,nr: (1-np.any(l[0:0+footprint]))*k[0]
            model.add_lattice_reaction(init, parameters, rxn_name = 'initiation', exclusion=1, frame=0, loc=0, dexist=1,)
            
            
            # DEFAULT TERMINATION 
            leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[k[2]-1]*k[1] #(lattice location 590 = 1) * parameter
            model.add_lattice_reaction(leave, parameters, rxn_name='termination', exclusion=0, frame=0, loc=590, dexist=-1,)
            
            
            # DEFAULT STEPPING OF ELONGATION USING THE ELONGATION MATRIX
            #default stepping
            elong_step = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+footprint])) for i in range(nr)]
            model.add_ribosome_reaction(elong_step, parameters, rxn_name='elongation', exclusion=1, dloc=1) 
            
            
            # finally specify which reactions are ribosome specific
            model._ribosome_reactions = [2]
            model._constant_reactions = [0,1,]
            model._lattice_arr0 = np.zeros([model._length+1], dtype=int)

            return model
        
            # #overwrite to a default model if just a protein obj was passed
            # particle_size = self.default_footprint_size 
            # n_colors = int(len(mRNA_model.tag_epitopes)) # number of color tags
            # n_states = 0
            # n_resources = 0

            # max_rib = int(mRNA_model.total_length/particle_size)+5 #total number of ribosomes possible
            

            # L = len(mRNA_model.kelong)
            # #build the default reaction matrix
            # rxn_mat = np.array([[  2.,   1.,   0.,   0.,   1.,   0.,   0.,   0.,   0.],  #initiation
            #                     [  2.,   0.,   0., L-1,  -1.,   0.,   0.,   0.,   0.],     #termination
            #                     [  0.,   1.,   0.,   0.,   0.,   0.,   1.,   0.,   0.]]).astype(np.int32) #stepping

            # n_rxns = int(rxn_mat.shape[0])
            # n_ribosome_reactions = np.sum(rxn_mat[:,0] == 0)
            # n_constant_reactions = len(rxn_mat) - n_ribosome_reactions
            # kelong_mat = mRNA_model.kelong_mat
            # kelong_mat[0,-1] = 0
            # probe_mat = mRNA_model.probe_mat
            
            
            # init = lambda k,t,p,ke,o,l,pr,s,r,nr: (1-np.any(l[0:0+int(particle_size)]))*k[0]
            # leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[L-1]*k[1]
            # elong_step = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+particle_size])) for i in range(nr)]
            # probe_function = lambda k,t,p,ke,o,l,pr,s,r,nr: 1
            # constant_props = [init, leave]
            # constant_parameters = [[mRNA_model.ki], [0,mRNA_model.kt]]
            # probe_parameters = 0
            # ribosome_props = [elong_step]
            # ribosome_parameters = [[mRNA_model.ke_mu]]
            # reaction_ids =[0,1,2]

    def __setup_python_sim(self, mRNA_model, probe_list, ki, kt):
        # initalize constants, flags, and initial state of the simulation
        footprint, max_rib, rxn_mat, n_colors, n_states, n_rxns, n_resources, n_constant_reactions, n_ribosome_reactions, L, kelong_mat, probe_mat, constant_props, constant_parameters, ribosome_props, ribosome_parameters, probe_function, probe_parameters, reaction_ids = self.__constants(mRNA_model, probe_list, ki, kt)
        #rib_arr, lattice_arr, state_arr, resource_arr, occupied, NR = self.__initalize_trajs(mRNA_model)    
            
        probe_fun = inspect.getsourcelines(probe_function)[0][0].split(':')[-1].replace('\n','').replace(' ','')
        if probe_fun == '1':
            use_probe_fun = False
        else:
            use_probe_fun = True
        
        
        return mRNA_model, footprint, max_rib, rxn_mat, n_colors, n_states, n_rxns, n_resources, n_constant_reactions, n_ribosome_reactions, L, kelong_mat, probe_mat, constant_props, constant_parameters,ribosome_props, ribosome_parameters, probe_function, use_probe_fun, probe_parameters, reaction_ids
        
    def __run(self, mRNA_model, footprint, max_rib, rxn_mat, n_colors, n_states, n_rxns,
            n_resources, n_constant_reactions, n_ribosome_reactions, L,
            kelong_mat, probe_mat, constant_props, constant_parameters,
            ribosome_props, ribosome_parameters, probe_function, use_probe_fun, 
            probe_parameters, reaction_ids, t, burnin=0, seed=None):
        
        rib_arr, lattice_arr, state_arr, resource_arr, occupied, NR = self.__initalize_trajs(mRNA_model, n_colors)    
        
        reaction_taken, dexist, rib_id, ribosome_moved, tindex = [0,]*5
        #rint(state_arr)
        if seed != None:
            np.random.seed(seed)

        def get_constant_props(*props):
            return [props[i](constant_parameters[i],tc, rib_arr, kelong_mat, occupied, lattice_arr, probe_mat, state_arr, resource_arr, NR) for i in range(len(props))]
        
        def get_ribosome_props(*props):
            rates = [props[i](ribosome_parameters[i],tc, rib_arr, kelong_mat, occupied, lattice_arr, probe_mat, state_arr, resource_arr, NR) for i in range(len(props))]
            return [item for sublist in rates for item in sublist]
        
        # Initalize arrays to store the trajectory
        #intensity_array = np.zeros([len(t), n_colors], dtype=int)
        ribosome_array = np.zeros([len(t), max_rib, 4+n_colors+n_constant_reactions+n_ribosome_reactions], dtype=int)
        state_array = np.zeros([len(t), n_states], dtype=int)
        resource_array = np.zeros([len(t), n_resources], dtype=int)
        
        
        tc = 0-burnin # current time
        tf = t[-1] # final time point
        max_lattice_node = len(lattice_arr)-1
        while tc < tf:

            # get propensities and where
            rates = get_constant_props(*constant_props) + get_ribosome_props(*ribosome_props)
            if sum(n < 0 for n in rates):
                msg = 'Negative Rate detected, double check model design.'\
                    ''
                raise custom_err.NegativeRateError(msg)

            # select propensity

            rate_sum = np.cumsum(rates)

            tc = (tc-np.log(np.random.rand())/rate_sum[-1]) # Update the time point randomly
            ro = rate_sum[-1]*np.random.rand()  #draw random number for reaction
            
            # record
            while tc >= t[tindex]:

                ribosome_array[tindex] = rib_arr
                state_array[tindex] = state_arr
                resource_array[tindex] = resource_arr
                #intensity_array[tindex] = np.sum(rib_arr[:,4:4+n_colors],axis=0)
                tindex += 1
                if tindex == len(t):
                    return ribosome_array, state_array, resource_array
                    break


            for i in range(len(rates)): #pick the next reaction rate
                if rate_sum[i] >= ro:
                    event = i #which raction happened
                    rid = i
                    break
            if event>=n_constant_reactions: #if its a ribosome reaction, on which ribosome did it occur

                rib_ind = (event-n_constant_reactions)%NR  #edit event to match reaction matrix
                event = int((event-n_constant_reactions)/NR) + n_constant_reactions

            event = reaction_ids.index(event) # map the event to the reactions (reaction matrix and propensities may not match)
            # do the reaction
            ribosome_moved = 0
            if rxn_mat[event][0] == 2: # lattice reaction
                # find the matching ribosome
                fr = rxn_mat[event][2]
                loc = rxn_mat[event][3]


                if rxn_mat[event][4] == 1: #ribosome arriving
                    rib_ind = NR
                    rib_id += 1
                    rib_arr[NR,0] = rib_id
                    rib_arr[NR,1] = 1
                    rib_arr[NR,2] += rxn_mat[event][2]
                    rib_arr[NR,3] += rxn_mat[event][3]
                    rib_arr[NR,4+n_colors+event] += 1
                    NR += 1
                    ribosome_moved = 1



                elif rxn_mat[event][4] == -1: #ribosome leaving
                    rib_ind = np.where(rib_arr[:,3] == loc)[0][0]
                    rib_arr[rib_ind,[1,2,3]] = rib_arr[rib_ind, [1,2,3]] + rxn_mat[event,4:7]
                    NR -= 1
                    rib_arr[rib_ind] = 0
                    rib_arr[rib_ind:-1] = rib_arr[rib_ind+1:]
                    rib_arr[-1] = 0
                    lattice_arr[:] = 0
                    lattice_arr[occupied[:NR]] = 1

                else: # ribosome moving
                    rib_ind = np.where(rib_arr[:,3] == loc)[0][0]
                    #rib_arr[rib_ind,[1,2,3]] = rib_arr[rib_ind, [1,2,3]] + rxn_mat[event,4:7]
                    rib_arr[rib_ind, 1:(n_colors+4+1)] = rib_arr[rib_ind, 1:(n_colors+4+1)] + rxn_mat[event,4:n_colors+8]
                    rib_arr[NR,4+n_colors+event] += 1
                #change states
                if n_states > 0:
                    state_arr = state_arr +  rxn_mat[event][n_colors+7:n_colors+7+n_states]
                    if sum(n < 0 for n in state_arr):
                        msg = 'Negative state detected, check model design.'\
                            ''
                        raise custom_err.StatesError(msg)
                    #print('***')
                    #print(rxn_mat[event][n_colors+7:n_colors+7+n_states]  )

                #change resources
                if n_resources > 0:
                    resource_arr = resource_arr +  rxn_mat[event][n_colors+7+n_states:]
                    if sum(n < 0 for n in resource_arr):
                        msg = 'Negative resource detected, check model design.'\
                            ''
                        raise custom_err.NegativeResourcesError(msg)


                if rxn_mat[event][6] !=0:
                    ribosome_moved = 1

            if rxn_mat[event][0] == 0: #ribosome reaction
                # change movement and colors
                # dexist, dframe, dloc, dprobe1... dprobeN
                rib_arr[rib_ind, 1:(n_colors+4)] = rib_arr[rib_ind, 1:(n_colors+4)] + rxn_mat[event,4:n_colors+7]

                #change states
                if n_states > 0:
                    state_arr = state_arr +  rxn_mat[event][n_colors+7:n_colors+7+n_states]
                    if sum(n < 0 for n in state_arr):
                        msg = 'Negative state detected, check model design.'\
                            ''
                        raise custom_err.StatesError(msg)

                #change resources
                if n_resources > 0:
                    resource_arr = resource_arr +  rxn_mat[event][n_colors+7+n_states:]
                    if sum(n < 0 for n in resource_arr):
                        msg = 'Negative resource detected, check model design.'\
                            ''
                        raise custom_err.NegativeResourcesError(msg)

                ribosome_moved = 0
                if rxn_mat[event][6] !=0:
                    ribosome_moved = 1

                if rxn_mat[event][4] == -1: #ribosome left
                    NR -= 1
                    rib_arr[rib_ind] = 0
                    rib_arr[rib_ind] = 0
                    rib_arr[rib_ind:-1] = rib_arr[rib_ind+1:]
                    rib_arr[-1] = 0
                    lattice_arr[:] = 0
                    lattice_arr[occupied[:NR]] = 1

                rib_arr[rib_ind,4+n_colors+event] += 1


            if rxn_mat[event][0] == 1: # state reaction
                state_arr = state_arr +  rxn_mat[event][n_colors+7:n_colors+7+n_states]
                if sum(n < 0 for n in state_arr):
                    msg = 'Negative state detected, check model design.'\
                        ''
                    raise custom_err.StatesError(msg)

            if rxn_mat[event][0] == 3: # resource reaction
                resource_arr = resource_arr +  rxn_mat[event][n_colors+7+n_states:]
                if sum(n < 0 for n in resource_arr):
                    msg = 'Negative resource detected, check model design.'\
                        ''
                    raise custom_err.NegativeResourcesError(msg)

            if rxn_mat[event][0] == 4: # probe reaction not implemented yet
                state_mat[rxn_mat[event][1]] += rxn_mat[event][2]
                state_mat[rxn_mat[event][3]] += rxn_mat[event][4]
                state_mat[rxn_mat[event][5]] += rxn_mat[event][6]


            # Check probes
            if ribosome_moved:
                # update occupied vector
                occupied = rib_arr[:,3]
                if max(occupied) > max_lattice_node:
                        msg = 'A particle stepped out of the maximum lattice node, check model design.'\
                            ''
                        raise custom_err.SteppedOutsideLattice(msg)
                if min(occupied) < 0:
                        msg = 'Negative lattice location stepped to, check model design.'\
                            ''
                        raise custom_err.SteppedOutsideLattice(msg)                    
                    
                # update the lattice vector
                lattice_arr[:] = 0
                lattice_arr[occupied[:NR]] = 1

                # check for probes
                pr = probe_mat[rib_arr[rib_ind,2], rib_arr[rib_ind,3]]
                if pr != 0:
                    if use_probe_fun:
                        if probe_function(ribosome_parameters[i],tc, rib_arr, kelong_mat, occupied, lattice_arr, probe_mat, state_arr, resource_arr, NR):
                            rib_arr[rib_ind, 3+pr] +=1 #add a probe if it passed a location and passes the probe function
                    else:
                        rib_arr[rib_ind, 3+pr] +=1
                ribosome_moved = 0
                

    
class CustomSSASoln:
    def __init__(self, mRNA_model, n_colors, rib_array, state_array, resource_array, t, burnin, n_traj, solve_time):
        self.ribosome_array = rib_array
        self.state_array = state_array
        self.resource_array = resource_array
        self.t = t
        self.burnin = burnin
        
        
        try:
            particle_size = mRNA_model.particle_size
            is_model_obj = True
        except:
            is_model_obj = False
            
        if is_model_obj:
            self.kelong_mat = mRNA_model._kelong_mat
            self.L = mRNA_model._kelong_mat.shape[1]
        else:
            if isinstance(mRNA_model, list):
                self.kelong_mat = np.zeros([3,len(mRNA_model)])
                self.kelong_mat[0,:] = mRNA_model
                self.L = len(mRNA_model)
            else:
                self.kelong_mat = mRNA_model.kelong_mat
                self.L = mRNA_model.kelong_mat.shape[1]

        #self.model_id = mRNA_model.model_id
        self.n_traj = n_traj
        if is_model_obj:
            self.n_colors = mRNA_model._n_colors
        else:
            if isinstance(mRNA_model, list):
                self.n_colors = n_colors
            else:
               self.n_colors = len(mRNA_model.tag_epitopes)
        self.solve_time = solve_time
        
    @property
    def lattice_arr(self):
        lattice_arr = np.zeros([self.n_traj, len(self.t), self.L ] ,dtype=int)
        for i in range(self.n_traj):
            for t in range(len(self.t)):
                rtraj = self.ribosome_array[i,t,:,3][self.ribosome_array[i,t,:,1]==1]
                if len(rtraj) > 0:
                    lattice_arr[i,t,rtraj] = 1
        return lattice_arr
    
    @property
    def intensity_arr(self):
        arr = np.zeros([self.n_traj, len(self.t), self.n_colors ] ,dtype=int)
        for i in range(self.n_traj):
            for t in range(len(self.t)):
                traj = np.sum(self.ribosome_array[i,t,:,4:4+self.n_colors],axis=0)
                arr[i,t,:] = traj
        return arr
    
    @property
    def I(self):
        return self.intensity_arr

    def condense_ribosome_array(self):
        '''
        Condenses the solution memory size by stripping zeros in the ribosome
        array.

        Returns
        -------
        None.

        '''
        max_ind = 0
        while self.ribosome_array[:,:,max_ind,:].sum() != 0 and max_ind < self.ribosome_array.shape[2]-1:
            max_ind += 1
        if max_ind < self.ribosome_array.shape[2]:
            self.ribosome_array = self.ribosome_array[:,:,:max_ind,:]
            
    def convert_to_aas(self, sequence):
        
        number_of_ribosomes_per_traj = []
        for i in range(len(self.ribosome_array)):
            number_of_ribosomes_per_traj.append(np.max(self.ribosome_array[i,:,:,0]))
        
        aa_over_traj = []
        for i in range(len(self.ribosome_array)):
            aa_over_time = []
            for j in range(1,len(number_of_ribosomes_per_traj)+1):
                
                codon = self.ribosome_array[i,:,:,3][np.where(self.ribosome_array[i,:,:,0] == j)]
                frame = self.ribosome_array[i,:,:,2][np.where(self.ribosome_array[i,:,:,0] == j)]
                time = np.where(self.ribosome_array[i,:,:,0] == j)
                
                aa_over_time.append(np.vstack([time, frame, codon]))
            aa_over_traj.append(aa_over_time)
            
        aas = []
        for i in range(len(self.ribosome_array)): #over trajectories
            aas2 = []
            for j in range(len(self.ribosome_array[0])): #over time
                aas2.append()
                
        return
        
    def save(self, filename, fmt ='.npy'):
        self.condense_ribosome_array() # compress the ribosome array to save memory
        
        constants = 5
        barcode = [constants, *self.ribosome_array.shape, *self.state_array.shape,
                   *self.resource_array.shape, *self.kelong_mat.shape, *self.t.shape]
        
        if fmt == '.npy':
            # compress everything into a singular linear array for saving with the barcode first
            linear_save = np.hstack([len(barcode)] + barcode + 
                                    [self.solve_time, self.n_traj, self.L, self.n_colors, self.burnin] + 
                                    [x.flatten() for x in [self.ribosome_array, self.state_array,
                                     self.resource_array, self.kelong_mat, self.t ]])
            np.save(filename, linear_save)
        
        if fmt == '.npz':
            # compress everything into npz format
            
            np.savez_compressed(filename, 
                                constants=[self.solve_time, self.n_traj, self.L, self.n_colors, self.burnin],
                                ribosome_array=self.ribosome_array,
                                state_array=self.state_array,
                                resource_array=self.resource_array,
                                kelong_mat = self.kelong_mat,
                                t=self.t, allow_pickle=False)
        
        return
    
    

class TranslationSolvers():
    '''
    Container class for the solvers
    '''
    def __init__(self, time=None,xi=None):
        self.k = None
        self.k_bind = None
        self.k_term = None
        self.multiframe = False
        self.additional_rxns = {}
        self.probe_locations = None
        self.colors = 1

        try:
            ssa_translation_lowmem.np
            self.cython_available = True
        except:
            self.cython_available = False


        self.default_conditions = {'low_mem':True,
                                   'perturb':[0,0,0],
                                   'leaky_probes':False,
                                   'bins':None,
                                   'kprobe':1,
                                   'footprint':9,
                                   'burnin':0,
                                   'record_stats':False,
                                   'kon':1.,
                                   'koff':1.,
                                   'n_traj':30,
                                   'bursting':False,
                                   'tRNA':False,
                                   'ribosome_preallocation':100,
                                   'dynamic_rib_prealloc':False
                                   }

        self.t = np.linspace(0,1000,1001)
        self.x0 = []
        self._poi = None

    @property
    def protein(self):
        return self._poi

    @protein.setter         #Auto detect amount of colors if protein is set
    def protein(self,newpoi):
        self._poi = newpoi
        try:
            self.colors = self._poi.probe_loc.shape[0]
        except:
            pass


    def solve_ballistic_model(self, ki,ke, nt_seq=None,
                              length=None, poi=None, tag= None):
        '''

        Given a iniation rate and set of elongation rates,
        predict analytically the intensity statistics
        and decorrelation time tau.

        **warning this ignores the effect of collisiosn**

        Parameters
        ----------
        ki : float
            initation rate of ribosomes.
        ke : float
            mean elongation rate along the transcript.
        poi : protein object, optional
            a particular protein object to pull geometry from.
            The default is None.
        tag : list, optional
            list of tag locations if geometry not given.
            The default is None.

        Returns
        -------
        tau_analyticals : list of floats
            analytically sovled decorrelation times tau.
        mean_analyticals : list of floats
            the visible ribosome means on the transcript.
        var_analyticals : list of floats
            variances of ribosomes on the transcript.

        '''

        if nt_seq == None:
            nt_seq = poi.nt_seq
            
        if tag == None:
            if poi ==None:
                tags = []
            else:
                colors = np.where((poi.probe_loc)== 1)[0]
                locs = np.where((poi.probe_loc)== 1)[1]
                tags = []
                for i in range(0,max(colors)+1):
    
                    tags.append(locs[np.where(colors == i)].tolist())
        
        #parse out stop codons
        if nt_seq[-3:].lower() in ['taa', 'tag', 'tga' ,'uga' ,'uaa', 'uag']:
            nt_seq = nt_seq[:-3]
        L = len(nt_seq)

        tau_analyticals = []
        mean_analyticals = []
        var_analyticals = []
        
        
        ke_analytical = L*ke / np.sum(self.__get_ui(nt_seq[:-3]))
        tau_analytical = (L )/ke_analytical  #analytical tau ie autocovariance time
        for tag in tags:


            L = poi.total_length #get the total length of the gene
            Lm = np.mean(tag)  #the mean location of the tag epitopes
            L_after_tag = L - tag[-1]
            L_tag = int((tag[-1] - tag[0]) / 2)

            mean_analytical = ki*tau_analytical * (1.-Lm/float(L)) # mean intensity
            var_analytical = ki*tau_analytical * (1.-Lm/float(L))**2  #var intensity
          
            mean_analyticals.append(mean_analytical)
            var_analyticals.append(var_analytical)

        return tau_analytical,mean_analyticals,var_analyticals


    def invert_ballistic(self,tau_measured, mu_I, poi= None):
        '''

        Given a measured decorrelation time tau from FCS,
        predict the elongation rates and initation rates

        Parameters
        ----------
        tau_measured : float, int, list, ndarray
            measured analytical decorrelation time(s).
        mu_I : float, int, list, ndarray
            average intensitie(s).
        poi : Poi object (optional)
            Protein of interest to pull geometry from.

        Returns
        -------
        kes : list
            analytical average k_elongation(s).
        kis : list
            analytical average k_initations(s).

        '''
        if poi == None:
            poi = self._poi
        #if tag == None:
        colors = np.where((poi.probe_loc)== 1)[0]
        locs = np.where((poi.probe_loc)== 1)[1]

        tags = []

        if isinstance(tau_measured,np.ndarray):
            tau_measured = tau_measured.tolist()

        if not isinstance(tau_measured,list):
            tau_measured = [tau_measured]


        if isinstance(mu_I,np.ndarray):
            mu_I = mu_I.tolist()

        if not isinstance(mu_I,list):
            mu_I = [mu_I]


        for i in range(0,max(colors)+1):

            tags.append(locs[np.where(colors == i)].tolist())

        kes = []
        kis = []

        for i in range(len(tags)):
            tag= tags[i]
            L = poi.total_length #get the total length of the gene
            Lm = np.mean(tag)  #the mean location of the tag epitopes
            L_after_tag = L - tag[-1]
            L_tag = int((tag[-1] - tag[0]) / 2)

            ke_analytical = (L)/tau_measured[i]

            ke = ke_analytical/(L)*np.sum(self.__get_ui(poi.nt_seq[:-3]))

            kes.append(ke)
            ki = (mu_I[i]/len(tag)) / ( (1.-Lm/float(L))*tau_measured[i])
            kis.append(ki)

        return kes, kis



    def __get_ui(self, nt_seq, cdict = None):
        '''
        return the ratio of average gene copy number / sequence codon copy number
        '''

        codon_dict = CodonDictionaries.CodonDictionaries()
        if cdict == None:
            cdict = codon_dict.human_codon_frequency_bias_nakamura
        mean_u =codon_dict.mean_genecopynumber
        ui = []
        for i in range(0, len(nt_seq), 3):
            ui.append(mean_u/ cdict[nt_seq[i:i+3]])
        return ui


    '''
    def solve_ssa_set_conditions(self):

        ssa_conditions = self.default_conditions
        kprobe = self.default_conditions['kprobe']
        kon = self.default_conditions['kon']
        koff = self.default_conditions['koff']

        if kprobe != 1:
            ssa_conditions['leaky_probes'] = True
            leaky_probes = True

        perturb = ssa_conditions['perturb']
        leaky_probes = ssa_conditions['leaky_probes']
        low_memory = ssa_conditions['low_mem']
        record_stats = ssa_conditions['record_stats']
        bins = ssa_conditions['bins']
        n_traj = ssa_conditions['n_traj']


        probe_vec = self.protein.probe_vec.astype(np.int32)
        t = self.t

        k = [self.protein.ki,] + self.protein.kelong + [self.protein.kt,]
        self.__check_rates(k)
        x0 = self.x0


        if self.cython_available:
            if low_memory:
                ssa_obj = self.__solve_ssa_lowmem_combined(
                    k,t,x0,n_traj,ssa_conditions=ssa_conditions,
                    kon=kon, koff=koff, kprobe=kprobe)
            else:

                if record_stats:

                    if leaky_probes == False:
                        # if low_memory:
                        #     ssa_obj = self.__solve_ssa_lowmem(k,t,x0,n_traj, ssa_conditions = ssa_conditions)
                        # else:
                        ssa_obj = self.__solve_ssa(k,t,x0,n_traj,ssa_conditions = ssa_conditions)
                    else:
                        # if low_memory:
                        #     ssa_obj = self.__solve_ssa_lowmem_leaky(k,t,x0,k_probe,n_traj, ssa_conditions = ssa_conditions)
                        # else:
                        ssa_obj = self.__solve_ssa_leaky(k,t,x0,kprobe,n_traj,ssa_conditions = ssa_conditions)

                else:
                    if leaky_probes == False:
                        # if low_memory:
                        #     ssa_obj = self.__solve_ssa_lowmem_nostats(k,t,x0,n_traj, ssa_conditions = ssa_conditions)
                        # else:
                        ssa_obj = self.__solve_ssa_nostats(k,t,x0,n_traj,ssa_conditions = ssa_conditions)
                    else:
                        # if low_memory:
                        #     ssa_obj = self.__solve_ssa_lowmem_leaky_nostats(k,t,x0,k_probe,n_traj, ssa_conditions = ssa_conditions)
                        # else:
                        ssa_obj = self.__solve_ssa_leaky_nostats(k,t,x0,kprobe,n_traj,ssa_conditions = ssa_conditions)

        else:
            ssa_obj = self.__solve_ssa_python(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)




        return ssa_obj
    '''

    def solve_ssa_trna(self,k_index, k_diffusion, k_bind, kelong,
                       k_compl, t,x0=[], k_trna = None,
                       perturb=[0,0,0],leaky_probes=False,
                       kprobe=np.ones(1),probe_vec = None, probe_loc=None,
                       kon=1,koff=1,bursting=False,n_traj=10,
                       connection_mat=None):
        
        self.__check_rates_trna(k_index)

        ssa_conditions = self.default_conditions
        if np.sum(kprobe != 1) !=0:
            ssa_conditions['leaky_probes'] = True
            leaky_probes = True

        ssa_conditions['perturb'] = perturb
        ssa_conditions['leaky_probes'] = leaky_probes


        ssa_conditions['bursting'] = bursting

        provided_probe = False
        try:
            probe_vec[0]
            provided_probe = True
        except:
            pass

        provided_protein = False
        try:
            self.protein.kelong[0]
            provided_protein = True
        except:
            pass

        try:
            k_trna[0]
        except:

            strGeneCopy = CodonDictionaries().human_codon_frequency_bias_nakamura
            strGeneCopy.pop('TAG')
            strGeneCopy.pop('TAA')
            strGeneCopy.pop('TGA')

            k_trna = np.array(list(CodonDictionaries().human_codon_frequency_bias_nakamura.values()))


        if not provided_probe:
            if provided_protein:
                probe_vec = self.protein.probe_vec.astype(np.int32)
                probe_loc = self.protein.probe_loc.astype(np.int32)
            else:
                probe_loc = np.zeros([1,len(kelong)], dtype=np.int32)
                probe_vec = np.zeros([1,len(kelong)], dtype=np.int32)
                warnings.warn('no provided probe vector, using a blank probe'\
                              '. This will result in zero sum intensity.')
        else:
            probe_vec = probe_vec


        ssa_conditions['probe_vec'] = probe_vec
        ssa_conditions['probe_loc'] = probe_loc


        try:
            connection_mat[0]
            ssa_obj = self.__solve_ssa_trna_connection_mat(k_index,k_trna, k_diffusion,k_bind,kelong, k_compl,t,x0,n_traj,ssa_conditions = ssa_conditions)

        except:
            ssa_obj = self.__solve_ssa_trna(k_index,k_trna, k_diffusion,k_bind,kelong, k_compl,t,x0,n_traj,ssa_conditions = ssa_conditions)

        return ssa_obj



    def solve_custom_model(self, model, parameters, kelong, t_array,
                           stoich_lattice, stoich_states,
                           xi_lattice, xi_states, n_traj=10,
                           probe_loc=None, poi=None):



        n_total_reactions = int(stoich_states.shape[0] + stoich_lattice.shape[0])
        length = len(kelong)

        try:
            probe_loc[0]
        except:
            probe_loc = poi.probe_loc


        #wash inputs:
        stoich_lattice = stoich_lattice.astype(np.int32)
        stoich_states = stoich_states.astype(np.int32)
        xi_lattice = xi_lattice.astype(np.int32)
        xi_states = xi_states.astype(np.int32)
        parameters = parameters.astype(np.float)
        t_array = t_array.astype(np.float)
        probe_loc = probe_loc.astype(np.int32)

        Ncolors = probe_loc.shape[0]
        max_rib = int(len(kelong)/9 + 5)
        Nt = len(t_array)

        seeds = np.random.randint(0,0x7fffff, n_traj)
        st = time.time()
        all_results = np.zeros([n_traj,  Nt, max_rib,],dtype=np.int32)
        all_intensities = np.zeros([n_traj,  Nt, Ncolors,],dtype=np.int32)
        all_states = np.zeros([n_traj,  Nt, max(xi_states.shape),],dtype=np.int32)

        for i in range(0,n_traj):
            result = np.zeros([Nt, max_rib,],dtype=np.int32)
            intensity = np.zeros([Nt, Ncolors,],dtype=np.int32)
            states = np.zeros([ Nt, max(xi_states.shape),],dtype=np.int32)
            a = model.run_ssa_cpp(result, intensity, states, stoich_states,
                                  stoich_lattice, parameters,
                                  np.array(kelong,dtype=np.float),  t_array,
                                  xi_lattice, xi_states,
                                  probe_loc,  length,
                                  seeds[i], n_total_reactions)
            all_results[i] = result.reshape(max_rib,Nt).T
            all_intensities[i] = intensity.reshape(Ncolors,Nt).T
            all_states[i] = states.reshape(max(xi_states.shape),Nt).T

        eval_time = time.time()-st



        rib_per_t = np.zeros((n_traj,Nt))


        validind = 0
        riblocs = []
        for i in range(len(all_results)):
            try:


                if np.where(np.sum(all_results[i].T,axis=1)!=0)[0][-1] > validind:

                    validind = np.where(np.sum(all_results[i].T,axis=1)!=0)[0][-1]

            #there was a blank
            except:
                pass
            
        if validind == 0:
            validind = 1
        all_results = all_results[:,:,:validind]

        for i in range(all_results.shape[0]):
            rib_loc = all_results[i,:,:]
            rt = np.count_nonzero(rib_loc.T,axis=0)
            rib_per_t[i,:] = rt

        rib_density_loc, rib_density_count = np.unique(all_results,return_counts=True)
        rib_density_loc = rib_density_loc[1:]
        rib_density_count = rib_density_count[1:]
        rib_density = rib_density_count / np.sum(rib_density_count)

        ssa_obj = SSA_Soln()


        ssa_obj.n_traj = n_traj
        ssa_obj.k = kelong
        # ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        ssa_obj.rib_per_t = rib_per_t
        ssa_obj.rib_mean = np.mean(rib_per_t.T)

        ssa_obj.rib_density = rib_density
        #ssa_obj.rib_means = ribosome_means
        ssa_obj.intensity_vec = all_intensities

        ssa_obj.time = t_array
        ssa_obj.time_rec = t_array
        ssa_obj.ribosome_locations = all_results
        ssa_obj.states = all_states
        ssa_obj.eval_time = eval_time

        return ssa_obj





    def solve_ssa_protease(self, ke, t, ki=.033, kt = 10, 
                  x0=[], n_traj=100, bins=None, protease_locations = [], protease_rates = [],
                  perturb=[0,0,0,0], leaky_probes=False, kprobe=np.ones(1),
                  record_stats=False, probe_vec=None, probe_loc=None,
                  kon=1, koff=1, bursting=False, rib_prealloc=200,
                  dynamic_prealloc=False, photobleaching=False, photobleaching_rate=0.0):
        '''

        Solve a stochastic simulation of ribosomes on an mRNA transcript
        given codon dependent rates (TASEP) with protease trimming

        Parameters
        ----------
        ke : list
            list of kelongation rates, can be accessed via a protein object by POI.kelong.
        t : list, arr
            a list/array of times to run the simulation for.
        ki : float, optional
            Initiation Rate. The default is .033 1/time.
        kt : float, optional
            Termination Rate. The default is 10 1/time.
        x0 : numpy array, optional
            1xRibsomes max intial condition. The default is []
            (no intial condition).
        n_traj : int, optional
            Number of trajectories to simulate. The default is 100.
        bins : numpy array, optional
            binning key array, for example [2,2,2,2,2] would bin a 10 length
            ke into 5 bins. The default is None.
        protease_locations : 1D array like, optional
            list of locations of where protease can trim nascent chains
        protease_rates : 1D array like, optional
            list of rates for each protease to trim available sites.
        low_memory : Bool, optional
            Use lower memory when simulating, meaning dont keep track of
            ribosome positions over time, just the Intensity. The default is True.
        perturb : list [bool, bool, float], optional
            Apply a pertubation [frap, inhibitor, application time], for
            example [1,0,50] would apply frap at time 50. The default is [0,0,0].
        leaky_probes : bool, optional
            Use leaky probes with a probability kprobe for binding when
            passing each epitope. The default is False.
        kprobe : list of floats, optional
            probabilities for each probe to bind, for example a two color
            system should have [.1,.9] for 10% chance of color 1 binding,
            and 90% of color two binding on each epitope.
            The default is np.ones(1).
        record_stats : bool, optional
            Record stats such as collisions or ribosome dwell times.
            The default is False.
        probe_vec : numpy array, optional
            probe cumsum to convert to Intensities Ncolor x L_transcript.
            If left blank will use POI. The default is None.
        probe_loc : numpy array, optional
            probe locations Ncolor x L_transcript with 1 for epitope
            locations, 0 otherwise. If left blank will use POI.
            The default is None.
        kon : float, optional
            on rate for bursting. The default is 1.
        koff : float, optional
            off rate for bursting. The default is 1.
        bursting : bool, optional
            use a on/off bursting dynamic for initiation. The default is False.
        rib_prealloc : int, optional
            how many ribosomes to consider binding at a time, for a
            transcript the maximum that can bind at a time is L/9.
            The default is 200.
        dynamic_prealloc : bool, optional
            calculate the expeted number of maximum ribosomes and use that
            to preallocate memory. The default is False.

        Returns
        -------
        ssa_obj : SSA_Soln_obj
            returns a container object containing all the stats from
            the simulation, such as intensity.

        '''


        self.__check_rates(ke)

        ssa_conditions = self.default_conditions
        if np.sum(kprobe != 1) !=0:
            ssa_conditions['leaky_probes'] = True
            leaky_probes = True

        if isinstance(protease_locations, list):
            protease_locations = np.array(protease_locations,dtype = np.int32)
            
        if isinstance(protease_rates, list):
            protease_rates = np.array(protease_rates,dtype = np.float64)

        ssa_conditions['perturb'] = perturb
        ssa_conditions['leaky_probes'] = leaky_probes
        ssa_conditions['low_mem'] = False
        ssa_conditions['record_stats'] = record_stats
        ssa_conditions['bins'] = bins
        ssa_conditions['bursting'] = bursting
        ssa_conditions['ribosome_preallocation'] = rib_prealloc
        ssa_conditions['dynamic_rib_prealloc'] = dynamic_prealloc
        ssa_conditions['photobleaching'] = photobleaching
        ssa_conditions['photobleaching_rate'] = photobleaching_rate
        ssa_conditions['protease_rates'] = protease_rates
        ssa_conditions['protease_locations'] = protease_locations
        provided_probe = False

        try:
            probe_vec[0]
            provided_probe = True
        except:
            pass

        provided_protein = False
        try:
            self.protein.kelong[0]
            provided_protein = True

        except:
            pass

        if provided_protein: #parse out stop codons if they are included
            includes_stop = False
            if self.protein.aa_seq[-1] == '*' or self.protein.nt_seq[-1].upper() in ['UAA','TAA','TAG','UAG', 'UGA','TGA']:
                if len(self.protein.aa_seq) == len(ke):
                    includes_stop = True



        if not provided_probe:
            if provided_protein:
                if includes_stop:
                    probe_vec = np.require(self.protein.probe_vec.astype(np.int32)[:,:-1],requirements=['C'])
                    probe_loc = np.require(self.protein.probe_loc.astype(np.int32)[:,:-1],requirements=['C'])
                    ke = ke[:-1]
                else:
                    probe_vec = self.protein.probe_vec.astype(np.int32)
                    probe_loc = self.protein.probe_loc.astype(np.int32)
            else:
                probe_loc = np.zeros([1,len(ke)], dtype=np.int32)
                probe_vec = np.zeros([1,len(ke)], dtype=np.int32)
                warnings.warn('no provided probe vector, using a blank probe'\
                              '. This will result in zero sum intensity.')
                
        else:
            probe_vec = probe_vec


        ssa_conditions['probe_vec'] = probe_vec
        ssa_conditions['probe_loc'] = probe_loc


        if dynamic_prealloc:

            ssa_conditions['ribosome_preallocation'] = int(len(ke)/ssa_conditions['footprint'])+5

        k = [ki,] + ke + [kt,]
        if self.cython_available:

            ssa_obj = self.__solve_ssa_replicon(k,t,x0,n_traj,protease_locations, protease_rates, ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)


        else:
            warnings.warn('Warning: Protease model is only included in the c++'\
              'This will return None until you install the translation c++.')


       #Generate metadata
        colors = ssa_obj.intensity_vec.shape[0]
        st = ssa_obj.start_time
        ft = ssa_obj.time_rec[-1]
        lt = len(ssa_obj.time_rec)
        if bins == None:
            bstr = 0
        else:
            bstr=1
        sid = 'L' + str(len(k)) + 'N' + str(n_traj) + 'T' + str(st) + '_' + str(ft) + '_' + str(lt)
        fstr = '-F' + str(int(self.cython_available)) + str(int(sum(perturb))) + str(int(leaky_probes)) + str(int(False))  + str(int(record_stats))  +  str(int(bstr)) + str(int(bursting))
        cstr = 'C' + str(int(colors))
        nprobe = np.sum(probe_loc,axis=1)
        for value in nprobe:
            cstr = cstr+ 'P' + str(int(value))
        sid = sid + cstr + fstr
        ssa_obj._SSA_Soln__meta['id'] = sid

        return ssa_obj


    def solve_ssa(self, ke, t, ki=.033, kt = 10,
                  x0=[], n_traj=100, bins=None, low_memory=True,
                  perturb=[0,0,0,0], leaky_probes=False, kprobe=np.ones(1),
                  record_stats=False, probe_vec=None, probe_loc=None,
                  kon=1, koff=1, bursting=False, rib_prealloc=200,
                  dynamic_prealloc=False, photobleaching=False, photobleaching_rate=0.0):
        '''

        Solve a stochastic simulation of ribosomes on an mRNA transcript
        given codon dependent rates (TASEP)

        Parameters
        ----------
        ke : list
            list of kelongation rates, can be accessed via a protein object by POI.kelong.
        t : list, arr
            a list/array of times to run the simulation for.
        ki : float, optional
            Initiation Rate. The default is .033 1/time.
        kt : float, optional
            Termination Rate. The default is 10 1/time.
        x0 : numpy array, optional
            1xRibsomes max intial condition. The default is []
            (no intial condition).
        n_traj : int, optional
            Number of trajectories to simulate. The default is 100.
        bins : numpy array, optional
            binning key array, for example [2,2,2,2,2] would bin a 10 length
            ke into 5 bins. The default is None.
        low_memory : Bool, optional
            Use lower memory when simulating, meaning dont keep track of
            ribosome positions over time, just the Intensity. The default is True.
        perturb : list [bool, bool, float], optional
            Apply a pertubation [frap, inhibitor, application time], for
            example [1,0,50] would apply frap at time 50. The default is [0,0,0].
        leaky_probes : bool, optional
            Use leaky probes with a probability kprobe for binding when
            passing each epitope. The default is False.
        kprobe : list of floats, optional
            probabilities for each probe to bind, for example a two color
            system should have [.1,.9] for 10% chance of color 1 binding,
            and 90% of color two binding on each epitope.
            The default is np.ones(1).
        record_stats : bool, optional
            Record stats such as collisions or ribosome dwell times.
            The default is False.
        probe_vec : numpy array, optional
            probe cumsum to convert to Intensities Ncolor x L_transcript.
            If left blank will use POI. The default is None.
        probe_loc : numpy array, optional
            probe locations Ncolor x L_transcript with 1 for epitope
            locations, 0 otherwise. If left blank will use POI.
            The default is None.
        kon : float, optional
            on rate for bursting. The default is 1.
        koff : float, optional
            off rate for bursting. The default is 1.
        bursting : bool, optional
            use a on/off bursting dynamic for initiation. The default is False.
        rib_prealloc : int, optional
            how many ribosomes to consider binding at a time, for a
            transcript the maximum that can bind at a time is L/9.
            The default is 200.
        dynamic_prealloc : bool, optional
            calculate the expeted number of maximum ribosomes and use that
            to preallocate memory. The default is False.

        Returns
        -------
        ssa_obj : SSA_Soln_obj
            returns a container object containing all the stats from
            the simulation, such as intensity.

        '''


        self.__check_rates(ke)

        ssa_conditions = self.default_conditions
        if np.sum(kprobe != 1) !=0:
            ssa_conditions['leaky_probes'] = True
            leaky_probes = True

        ssa_conditions['perturb'] = perturb
        ssa_conditions['leaky_probes'] = leaky_probes
        ssa_conditions['low_mem'] = low_memory
        ssa_conditions['record_stats'] = record_stats
        ssa_conditions['bins'] = bins
        ssa_conditions['bursting'] = bursting
        ssa_conditions['ribosome_preallocation'] = rib_prealloc
        ssa_conditions['dynamic_rib_prealloc'] = dynamic_prealloc
        ssa_conditions['photobleaching'] = photobleaching
        ssa_conditions['photobleaching_rate'] = photobleaching_rate
        provided_probe = False

        try:
            probe_vec[0]
            provided_probe = True
        except:
            pass

        provided_protein = False
        try:
            self.protein.kelong[0]
            provided_protein = True

        except:
            pass

        if provided_protein: #parse out stop codons if they are included
            includes_stop = False
            if self.protein.aa_seq[-1] == '*' or self.protein.nt_seq[-1].upper() in ['UAA','TAA','TAG','UAG', 'UGA','TGA']:
                if len(self.protein.aa_seq) == len(ke):
                    includes_stop = True



        if not provided_probe:
            if provided_protein:
                if includes_stop:
                    probe_vec = np.require(self.protein.probe_vec.astype(np.int32)[:,:-1],requirements=['C'])
                    probe_loc = np.require(self.protein.probe_loc.astype(np.int32)[:,:-1],requirements=['C'])
                    ke = ke[:-1]
                else:
                    probe_vec = self.protein.probe_vec.astype(np.int32)
                    probe_loc = self.protein.probe_loc.astype(np.int32)
            else:
                probe_loc = np.zeros([1,len(ke)], dtype=np.int32)
                probe_vec = np.zeros([1,len(ke)], dtype=np.int32)
                warnings.warn('no provided probe vector, using a blank probe'\
                              '. This will result in zero sum intensity.')
                
        else:
            probe_vec = probe_vec


        ssa_conditions['probe_vec'] = probe_vec
        ssa_conditions['probe_loc'] = probe_loc


        if dynamic_prealloc:

            ssa_conditions['ribosome_preallocation'] = int(len(ke)/ssa_conditions['footprint'])+5

        k = [ki,] + ke + [kt,]
        if self.cython_available:
            if low_memory:
                ssa_obj = self.__solve_ssa_lowmem_combined(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)
            else:
                ssa_obj = self.__solve_ssa_full_combined(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)


        else:
            ssa_obj = self.__solve_ssa_python(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)


       #Generate metadata
        colors = ssa_obj.intensity_vec.shape[0]
        st = ssa_obj.start_time
        ft = ssa_obj.time_rec[-1]
        lt = len(ssa_obj.time_rec)
        if bins == None:
            bstr = 0
        else:
            bstr=1
        sid = 'L' + str(len(k)) + 'N' + str(n_traj) + 'T' + str(st) + '_' + str(ft) + '_' + str(lt)
        fstr = '-F' + str(int(self.cython_available)) + str(int(sum(perturb))) + str(int(leaky_probes)) + str(int(low_memory))  + str(int(record_stats))  +  str(int(bstr)) + str(int(bursting))
        cstr = 'C' + str(int(colors))
        nprobe = np.sum(probe_loc,axis=1)
        for value in nprobe:
            cstr = cstr+ 'P' + str(int(value))
        sid = sid + cstr + fstr
        ssa_obj._SSA_Soln__meta['id'] = sid

        return ssa_obj



    def solve_ode(self,k,t,x0,ki,pl,bins=None, corr = False):
        '''
        Solve the system of odes describing the mRNA transcript by their forward rates.

        **WARNING** This cannot take into account collisions, if you want to
        simulate collisions you must use stochastic simulations.
        If you dont know if your system will encounter heavy collisions
        use the following rule of thumbs:
            ki ~ min(kelong)  or ki >=  L_transcript / ( 1/np.sum(1/kelong))

        these correspond to a single slow codon causing a jam or v_in >= v_out
        of the total system.

        Parameters
        ----------
        k : list
            elongation rates [ke...., kt].
        t : numpy array / list
            time vector.
        x0 : numpy array
            x0, 1xL_transcript with 1 for ribosomes occupying that location.
        ki : float
            initation rate.
        pl : numpy array
            probe location array, 0s and 1s.
        bins : numpy array, optional
            binning array. The default is None.
        corr : Bool, optional
            calculate the analytical correlation of intensity,
            **THIS WILL TAKE A LONG TIME**. The default is False.

        Returns
        -------
        ode_soln : ODEsoln OBJ
            ODE solution object.

        '''

        st = time.time()
        m_ode = models.TranslateODE()
        m_ode.N = len(k)
        m_ode.tf = t[-1]
        m_ode.ptimes = len(t)
        m_ode.ke = k
        m_ode.kb = ki
        m_ode.fi = 1
        m_ode.ti = t[0]
        m_ode.binary = pl

        #mu_intensity = m_ode.solve()

        m = models.TranslateCorrs()
        m.N = len(k)
        m.tf = t[-1]
        m.ptimes = len(t)
        m.ke = k
        m.kb = ki
        m.fi = 1
        m.ti = t[0]
        m.xi = x0
        m.binary = pl


        ode_soln = ODE_Soln()

        m.get_autonomous_matrix()
        #m.faster_solve()
        m.get_mean_SS()

        mean_ss = m.get_mean_SS()

        mean_I = m.map_to_fluorescence3(mean_ss)
        m.mu_ss = mean_ss

        var_ss = m.get_var_SS()
        var_I = m.map_to_fluorescence(var_ss)


        if corr:

            m.csolve()

            norm_acc = np.ravel((m.intensity)/var_I)
            intensity_acc = m.intensity
            ode_soln.intensity_acc = m.intensity
            ode_soln.intensity_acc_norm = norm_acc



        ode_soln.mu_state_ss = mean_ss
        ode_soln.var_state_ss = var_ss

        ode_soln.mu_It = m.solve()
        ode_soln.mu_I_ss = mean_I
        ode_soln.var_I_ss = var_I
        ode_soln.time = t
        ode_soln.bins = bins
        ode_soln.N = len(k)
        ode_soln.k = k
        ode_soln.x0 = x0
        ode_soln.fi = 1
        ode_soln.ki = ki
        ode_soln.prob_loc = pl
        ode_soln.solve_time = time.time()-st


        ode_soln.probvec = m.get_fluor_vec()
        if len(k) > 100:
            ode_soln.solve_method = 'expv'
        else:
            ode_soln.solve_method = 'odeint'


        return ode_soln



    '''
    def __solve_ssa(self,k,t,x0,n_traj,ssa_conditions=None):

        seeds = np.random.randint(0, 0x7FFFFFF, n_traj, dtype = np.int32)

        if ssa_conditions == None:
            ssa_conditions = self.default_conditions

        x0 = self.__check_x0(x0)


        rib_vec = []
        solutions = []
        solutionssave = []
        N_rib = 200
        colors = self.colors
        all_results,all_nribs,all_collisions,all_frapresults,all_ribtimes,all_col_points = self.__generate_mats(n_traj,k[0],t,N_rib,colors)
        footprint = ssa_conditions['footprint']
        evf = ssa_conditions['perturb'][0]
        evi = ssa_conditions['perturb'][1]
        intime = ssa_conditions['perturb'][2]
        non_consider_time = ssa_conditions['burnin']


        st = time.time()

        for i in range(n_traj):

            result,ribtimes,frapresult,coltimes,colpointsx,colpointst = self.__generate_vecs(k,t,N_rib,colors)
            nribs = np.array([0],dtype=np.int32)
            kelong = np.array(k[1:-1]).flatten()




            ssa_translation.run_SSA(result, ribtimes, coltimes, colpointsx,colpointst, kelong,frapresult, t, k[0], k[-1], evf, evi, intime, seeds[i],nribs,x0,footprint,200)
            #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)

            all_results[i, :] = result.T
            all_frapresults[i,:] = frapresult
            all_ribtimes[i,:] = ribtimes
            all_collisions[i,:] = coltimes
            all_nribs[i,:] = nribs

            endcolrec = np.where(colpointsx == 0)[0][0]

            colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
            all_col_points.append(colpoints.T)



        maxso = 0
        for i in range(n_traj):
            soln = all_results[i, :].reshape((N_rib, len(t)))


            validind = np.where(np.sum(soln,axis=1)!=0)[0]


            if np.max(validind) != N_rib-1:
                validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)

            so = soln[(validind,)]
            if so.shape[0] > maxso:
                maxso = so.shape[0]



        for i in range(n_traj):
            soln = all_results[i, :].reshape((N_rib, len(t)))


            validind = tuple([x for x in range(0,maxso)])


            so = soln[(validind,)]

            solutionssave.append(so)
            solutions.append(soln)


        collisions = np.array([[]])
        watched_ribs = []
        for i in range(n_traj):
            totalrib = all_nribs[i]

            if totalrib > all_collisions.shape[1]:
                collisions = np.append(collisions, all_collisions[i][:])
                watched_ribs.append(int(all_collisions.shape[1]))

            else:

                collisions = np.append(collisions, all_collisions[i][:int(totalrib[0])])
                watched_ribs.append(int(totalrib[0]))

        sttime = time.time() - st


        no_ribosomes = np.zeros((n_traj, (len(k)+1)))

        startindex = np.where(t >= non_consider_time)[0][0]

        #all_results = all_results[:,startindex*N_rib:]
        pv = self.protein.probe_vec
        I = np.zeros((colors,len(t), n_traj))

        for n in range(colors):
            for i in range(n_traj):
                traj = all_results[i,:].reshape((N_rib,len(t))).T
                for j in range(len(t)):
                    temp_output = traj[j,:]

                    I[n,j,i] = np.sum(pv[n][temp_output[temp_output>0]-1]  )

        for i in range(len(solutions)):
            for j in range(len(solutions[0][0][startindex:])):
                rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]

                no_ribosomes[i, rib_pos.astype(int)] += 1
        no_ribosomes = no_ribosomes[:, 1:]

        ribosome_means = np.mean(no_ribosomes, axis=0)
        ribosome_density = ribosome_means/len(k)

        no_ribosomes_per_mrna = np.mean(no_ribosomes)

        ssa_obj = SSA_Soln()
        ssa_obj.no_ribosomes = no_ribosomes
        ssa_obj.n_traj = n_traj
        ssa_obj.k = k
        ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        ssa_obj.rib_density = ribosome_density
        ssa_obj.rib_means = ribosome_means
        ssa_obj.rib_vec = rib_vec
        #ssa_obj.intensity_vec = intensity_vec
        ssa_obj.time_vec_fixed = t
        ssa_obj.time = t
        ssa_obj.time_rec = t[startindex:]
        ssa_obj.start_time = non_consider_time
        ssa_obj.watched_ribs = watched_ribs
        ssa_obj.intensity_vec = I
        ssa_obj.solutions = np.array(solutionssave)
        try:
            ssa_obj.col_points = all_col_points
        except:
            pass



        return ssa_obj

    '''


    def __map_to_intensity(self):

        return 1


    def __solve_ssa_trna(self, kindex, ktrna, kdiffusion, kbind, kelong,
                         kcompl, t, x0, n_traj, ssa_conditions=None, probe_vec=None ):

        seeds = np.random.randint(0, 0x7FFFFFF, n_traj, dtype=np.int32)

        if ssa_conditions == None:
            ssa_conditions = self.default_conditions

        x0 = self.__check_x0(x0)


        rib_vec = []
        solutions = []
        solutionssave = []
        N_rib = 200
        
        if probe_vec == None:
            pv = self.protein.probe_vec
        else:
            pv = probe_vec
            
        colors = int(pv.shape[0])

        n_trajectories = n_traj


        all_results,all_nribs,all_collisions,all_frapresults,all_ribtimes,all_col_points = self.__generate_mats(n_traj,kbind,t,N_rib,colors)
        all_trna_results = np.zeros((n_trajectories,61*len(t)), dtype=np.int32)

        footprint = ssa_conditions['footprint']
        evf = ssa_conditions['perturb'][0]
        evi = ssa_conditions['perturb'][1]
        intime = ssa_conditions['perturb'][2]
        non_consider_time = ssa_conditions['burnin']

        st = time.time()

        N_rib = 200


        result = np.zeros((len(t)*N_rib), dtype=np.int32  )

        #lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t<20)[0]))
        #all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)

        all_ribs = np.zeros((n_trajectories, 1))
        all_col_points = []

        #seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)
        x0 = np.zeros((N_rib),dtype=np.int32)


        for i in range(n_trajectories):



            trna_result = np.zeros((len(t)*61),dtype=np.int32)
            result,ribtimes,frapresult,coltimes,colpointsx,colpointst = self.__generate_vecs_trna(kbind,kindex,t,N_rib,colors)

            nribs = np.array([0],dtype=np.int32)

            if i == 0: #detect any int64
                inputs = [result,trna_result,ribtimes,coltimes,colpointsx,colpointst, kindex,ktrna,kdiffusion,frapresult,t,kbind,kcompl, seeds[i],nribs,x0,kelong]
                wash_inputs = self.__check_input_memview(inputs)

            if wash_inputs:
                #check memview so all given variables are in int32 if integer for C
                inputs = [result,trna_result,ribtimes,coltimes,colpointsx,colpointst, kindex,ktrna,kdiffusion,frapresult,t,kbind,kcompl, seeds[i],nribs,x0,kelong]

                result,trna_result,ribtimes,coltimes,colpointsx,colpointst, kindex,ktrna,kdiffusion,frapresult,t,kbind,kcompl, seed,nribs,x0,kelong = self.__check_memview(inputs)
            else:

                seed = seeds[i]

            ssa_translation_lowmem.run_SSA_trna_full(result,trna_result,ribtimes,coltimes,colpointsx,colpointst, kindex,ktrna,kdiffusion,frapresult,t,kbind,kcompl, 0,0,0, seed,nribs,x0,kelong)

            all_results[i, :] = result.T
            all_trna_results[i, :] = trna_result
           # all_frapresults[i, :] = frapresult
            all_ribtimes[i, :] = ribtimes
            all_collisions[i, :] = coltimes
            all_nribs[i, :] = nribs[0]

            endcolrec = np.where(colpointsx == 0)[0][0]

            colpoints = np.vstack((colpointsx[:endcolrec],
                                   colpointst[:endcolrec]))
            all_col_points.append(colpoints.T)

        maxso = 0
        for i in range(n_traj):
            soln = all_results[i, :].reshape((N_rib, len(t)))


            validind = np.where(np.sum(soln,axis=1)!=0)[0]


            if np.max(validind) != N_rib-1:
                validind = np.append(
                    np.where(np.sum(soln, axis=1)!=0)[0], np.max(validind)+1)

            so = soln[(validind,)]
            if so.shape[0] > maxso:
                maxso = so.shape[0]



        for i in range(n_traj):
            soln = all_results[i, :].reshape((N_rib, len(t)))


            validind = tuple([x for x in range(0,maxso)])
            so = soln[(validind,)]
            solutionssave.append(so)
            solutions.append(soln)

        collisions = np.array([[]])
        watched_ribs = []
        for i in range(n_traj):
            totalrib = all_nribs[i]

            if totalrib > all_collisions.shape[1]:
                collisions = np.append(collisions, all_collisions[i][:])
                watched_ribs.append(int(all_collisions.shape[1]))

            else:

                collisions = np.append(
                    collisions, all_collisions[i][:int(totalrib[0])])
                watched_ribs.append(int(totalrib[0]))

        sttime = time.time() - st


        no_ribosomes = np.zeros((n_traj, (len(kindex)+1)))

        startindex = np.where(t >= non_consider_time)[0][0]

        #all_results = all_results[:,startindex*N_rib:]
        
        I = np.zeros((colors, len(t), n_traj))

        for n in range(colors):
            for i in range(n_traj):
                traj = all_results[i, :].reshape((N_rib, len(t))).T
                for j in range(len(t)):
                    temp_output = traj[j, :]

                    I[n, j, i] = np.sum(pv[n][temp_output[temp_output>0]-1]  )

        for i in range(len(solutions)):
            for j in range(len(solutions[0][0][startindex:])):
                rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]

                no_ribosomes[i, rib_pos.astype(int)] += 1
        no_ribosomes = no_ribosomes[:, 1:]

        ribosome_means = np.mean(no_ribosomes, axis=0)
        ribosome_density = ribosome_means/len(kindex)

        no_ribosomes_per_mrna = np.mean(no_ribosomes)

        ssa_obj = SSA_Soln()
        ssa_obj.no_ribosomes = no_ribosomes
        ssa_obj.n_traj = n_traj
        ssa_obj.k = kindex
        ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        ssa_obj.rib_density = ribosome_density
        ssa_obj.rib_means = ribosome_means
        ssa_obj.rib_vec = rib_vec
        ssa_obj.I = I
        ssa_obj.eval_time = sttime
        ssa_obj.time_vec_fixed = t
        ssa_obj.time = t
        ssa_obj.time_rec = t[startindex:]
        ssa_obj.start_time = non_consider_time
        ssa_obj.watched_ribs = watched_ribs
        ssa_obj.intensity_vec = I
        ssa_obj.solutions = np.array(solutionssave)
        ssa_obj.all_trna_results = all_trna_results
        ssa_obj.ribtimes = all_ribtimes
        try:
            ssa_obj.col_points = all_col_points
        except:
            pass



        return ssa_obj

    def __solve_ssa_python(self, k, t, x0, n_traj, ssa_conditions=None,
                           kprobe=None, kon=None, koff=None, flags=[0,0,0]):

        seeds = np.random.randint(0, 0x7FFFFFF, n_traj, dtype=np.int32)

        if ssa_conditions == None:
            ssa_conditions = self.default_conditions

        x0 = self.__check_x0(x0)

        flags = [int(ssa_conditions['bursting']),
                 int(ssa_conditions['leaky_probes']),
                 int(ssa_conditions['record_stats'])]

        rib_vec = []
        solutions = []
        solutionssave = []
        N_rib = 200
        colors = self.colors

        all_results,all_nribs,all_collisions,all_frapresults,all_ribtimes,all_col_points = self.__generate_mats(n_traj,k[0],t,N_rib,colors)
        footprint = ssa_conditions['footprint']
        evf = ssa_conditions['perturb'][0]
        evi = ssa_conditions['perturb'][1]
        intime = ssa_conditions['perturb'][2]
        non_consider_time = ssa_conditions['burnin']


        st = time.time()


        warnings.warn('C++ extention missing or not compiled, using Python' \
                      ' Implementation. This will be much slower!' \
                      ' For C++ instillation instructions read the README or '\
                          'visit the github page: https://github.com/MunskyGroup/rSNAPsim')
        rib_vec = []

        solutions = []
        solutionssave = []
        N_rib = ssa_conditions['ribosome_preallocation']
        collisions = np.array([[]])
        all_results = np.zeros((n_traj, N_rib*len(t)), dtype=np.int32)
        I_internal = np.zeros((colors, len(t), n_traj))
        all_col_points = []
        watched_ribs = []
        for i in range(n_traj):

            soln,all_ribtimes,Ncol,col_points,intensity = self.__ssa_python(k, t, inhibit_time=intime+non_consider_time, FRAP=evf, Inhibitor=evi, flags=flags, kon=kon, kprobe=kprobe, koff=koff, ssa_conditions=ssa_conditions)
            #soln = soln.reshape((1, (len(time_vec_fixed)*N_rib)))

            collisions = np.append(collisions,Ncol)
            watched_ribs.append(int(len(collisions)))
            validind = np.where(np.sum(soln,axis=1) != 0)[0]
            all_col_points.append(np.array(col_points))
            if np.max(validind) != N_rib-1:
                validind = np.append(np.where(np.sum(soln, axis=1) != 0)[0],np.max(validind)+1)

            so = soln[(validind,)]

            #solutionssave.append(so)

            solutions.append(soln)

            result = soln.reshape((1, (len(t)*N_rib)))
            all_results[i, :] = result
            I_internal[:, : ,i] = intensity


        maxso = 0
        for i in range(n_traj):
            soln = all_results[i, :].reshape((N_rib, len(t)))

            validind = np.where(np.sum(soln,axis=1)!=0)[0]
            if np.max(validind) != N_rib-1:
                validind = np.append(np.where(np.sum(soln, axis=1)!=0)[0], np.max(validind)+1)

            so = soln[(validind,)]
            if so.shape[0] > maxso:
                maxso = so.shape[0]

        for i in range(n_traj):
            soln = all_results[i, :].reshape((N_rib, len(t)))

            validind = tuple([x for x in range(0, maxso)])
            so = soln[(validind,)]

            solutionssave.append(so)


        # for i in range(n_traj):
        #     # result,ribtimes,frapresult,coltimes,colpointsx,colpointst = self.__generate_vecs(k,t,N_rib,colors)
        #     # nribs = np.array([0],dtype=np.int32)
        #     # kelong = np.array(k[1:-1]).flatten()

        #     soln,all_ribtimes,Ncol,col_points  = self.__ssa_python(k,t,inhibit_time=intime+non_consider_time,FRAP=evf,Inhibitor=evi,flags=flags)

        #     all_results[i, :] = result.T
        #     all_frapresults[i,:] = frapresult
        #     all_ribtimes[i,:] = ribtimes
        #     all_collisions[i,:] = coltimes
        #     all_nribs[i,:] = nribs

        #     endcolrec = np.where(colpointsx == 0)[0][0]

        #     colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
        #     all_col_points.append(colpoints.T)


        # for i in range(n_traj):
        #     soln = all_results[i, :].reshape((N_rib, len(t)))


        #     validind = np.where(np.sum(soln,axis=1)!=0)[0]

        #     if np.max(validind) != N_rib-1:
        #         validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)

        #     so = soln[(validind,)]

        #     solutionssave.append(so)
        #     solutions.append(soln)

        collisions = np.array([[]])
        watched_ribs = []
        for i in range(n_traj):
            totalrib = all_nribs[i]

            if totalrib > all_collisions.shape[1]:
                collisions = np.append(collisions, all_collisions[i][:])
                watched_ribs.append(int(all_collisions.shape[1]))

            else:

                collisions = np.append(collisions, all_collisions[i][:int(totalrib[0])])
                watched_ribs.append(int(totalrib[0]))

        sttime = time.time() - st


        no_ribosomes = np.zeros((n_traj, (len(k)+1)))

        startindex = np.where(t >= non_consider_time)[0][0]

        #all_results = all_results[:,startindex*N_rib:]
        pv = self.protein.probe_vec
        # I = np.zeros((colors,len(t), n_traj))

        # for n in range(colors):
        #     for i in range(n_traj):
        #         traj = all_results[i,:].reshape((N_rib,len(t))).T
        #         for j in range(len(t)):
        #             temp_output = traj[j,:]

        #             I[n,j,i] = np.sum(pv[n][temp_output[temp_output>0]-1]  )

        for i in range(len(solutions)):
            for j in range(len(solutions[0][0][startindex:])):

                rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]
                if len(rib_pos.astype(int)) > 0:
                    no_ribosomes[i, rib_pos.astype(int)] += 1


        no_ribosomes = no_ribosomes[:, 1:]

        ribosome_means = np.mean(no_ribosomes, axis=0)
        ribosome_density = ribosome_means/len(k)

        no_ribosomes_per_mrna = np.mean(no_ribosomes)

        ssa_obj = SSA_Soln()
        ssa_obj.no_ribosomes = no_ribosomes
        ssa_obj.n_traj = n_traj
        ssa_obj.k = k
        ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        ssa_obj.rib_density = ribosome_density
        ssa_obj.rib_means = ribosome_means
        ssa_obj.rib_vec = rib_vec
        #ssa_obj.intensity_vec = intensity_vec
        ssa_obj.time_vec_fixed = t
        ssa_obj.time = t
        ssa_obj.time_rec = t[startindex:]
        ssa_obj.start_time = non_consider_time
        ssa_obj.watched_ribs = watched_ribs
        ssa_obj.intensity_vec = I_internal
        ssa_obj.I = I_internal
        ssa_obj.solutions = np.array(solutionssave)
        try:
            ssa_obj.col_points = all_col_points
        except:
            pass



        return ssa_obj



    def __solve_ssa_lowmem_combined(self, k, t, x0, n_traj,
                                    ssa_conditions=None, kon=1,
                                    koff=1, kprobe=[]):
        seeds = np.random.randint(0, 0x7FFFFFF, n_traj, dtype=np.int32)

        if isinstance(k, list):
            k = np.array(k).astype(np.float64)

        k = k.flatten()

        if kprobe == []:
            kprobe = np.ones(self.color)

        if isinstance(kprobe, list):
            kprobe = np.array(kprobe, dtype=np.float64)
        if isinstance(kprobe, int):
            kprobe = np.array([kprobe], dtype=np.float64)
        if isinstance(kprobe, float):
            kprobe = np.array([kprobe], dtype=np.float64)


        if ssa_conditions == None:
            ssa_conditions = self.default_conditions

        x0 = self.__check_x0(x0)

        probe_vec = ssa_conditions['probe_vec']

        colors = int(probe_vec.shape[0])
 

        if ssa_conditions['bursting'] == False:
            kon = 1
            koff = 1
        else:
            kon = kon
            koff = koff

        flags = np.array([int(ssa_conditions['bursting']),
                          int(ssa_conditions['leaky_probes']),
                          int(ssa_conditions['record_stats'])], dtype=np.int32)
        probe_loc = ssa_conditions['probe_loc']


        rib_vec = []
        solutions = []
        solutionssave = []


        N_rib = ssa_conditions['ribosome_preallocation']
        all_results,all_nribs,all_collisions,_,all_ribtimes,all_col_points = self.__generate_mats_lowmem(n_traj,k[0],t,N_rib,colors)
        footprint = ssa_conditions['footprint']
        evf = ssa_conditions['perturb'][0]
        evi = ssa_conditions['perturb'][1]
        intime = ssa_conditions['perturb'][2]
        intime_stop = ssa_conditions['perturb'][3]
        non_consider_time = ssa_conditions['burnin']
        photobleaching_rate = ssa_conditions['photobleaching_rate']
        photobleaching = ssa_conditions['photobleaching']

        st = time.time()

        for i in range(n_traj):

            result,ribtimes,_,coltimes,colpointsx,colpointst = self.__generate_vecs_lowmem(k,t,N_rib,colors)
            nribs = np.array([0], dtype=np.int32)


            if i == 0: #detect any int64
                inputs = [result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1], t, k[0], float(k[-1]), int(evf), int(evi), float(intime), float(intime_stop), seeds[i],nribs,x0,footprint, probe_vec ,int(colors), kon, koff, kprobe, probe_loc, flags,N_rib, photobleaching_rate, photobleaching]
                wash_inputs = self.__check_input_memview(inputs)

            if wash_inputs:
                #check memview so all given variables are in int32 if integer for C
                inputs = [result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1], t, k[0], float(k[-1]), int(evf), int(evi), float(intime),float(intime_stop), seeds[i],nribs,x0,footprint, probe_vec ,int(colors), kon, koff, kprobe, probe_loc, flags,N_rib, photobleaching_rate, photobleaching]
                result, ribtimes, coltimes, colpointsx,colpointst, kelong, t, ki, kt, evf, evi, intime,intime_stop, seed,nribs,x0,footprint, probe_vec ,colors, kon, koff, kprobe, probe_loc, flags, N_rib, photobleaching_rate, photobleaching = self.__check_memview(inputs)
            else:
                kelong = k[1:-1]
                ki = k[0]
                kt = float(k[-1])
                evf = int(evf)
                evi = int(evi)
                intime = float(intime)
                intime_stop = float(intime_stop)
                seed = seeds[i]
                colors = int(colors)
           
            ssa_translation_lowmem.run_SSA(result, ribtimes, coltimes, colpointsx,colpointst, kelong, t, ki, kt, evf, evi, intime,intime_stop, seed,nribs,x0,footprint, probe_vec ,colors, kon, koff, kprobe, probe_loc, flags, N_rib, photobleaching_rate, photobleaching)
            #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)

            all_results[i, :, :] = result.T
            #all_frapresults[i, :] = frapresult
            all_ribtimes[i, :] = ribtimes
            all_collisions[i, :] = coltimes
            all_nribs[i, :] = nribs

            endcolrec = np.where(colpointsx == 0)[0][0]

            colpoints = np.vstack((colpointsx[:endcolrec], colpointst[:endcolrec]))
            all_col_points.append(colpoints.T)


            # for i in range(n_traj):
            #     soln = all_results[i, :].reshape((N_rib, len(t),colors))

            #     so = soln
            #     solutionssave.append(so)
            #     solutions.append(soln)

            collisions = np.array([[]])
            watched_ribs = []
            for i in range(n_traj):
                totalrib = all_nribs[i]

                if totalrib > all_collisions.shape[1]:
                    collisions = np.append(collisions, all_collisions[i][:])
                    watched_ribs.append(int(all_collisions.shape[1]))

                else:

                    collisions = np.append(collisions, all_collisions[i][:int(totalrib[0])])
                    watched_ribs.append(int(totalrib[0]))

            sttime = time.time() - st


        no_ribosomes = np.zeros((n_traj, (len(k)+1)))

        startindex = np.where(t >= non_consider_time)[0][0]

        #all_results = all_results[:,startindex*N_rib:]

#        for i in range(len(solutions)):
#            for j in range(len(solutions[0][0][startindex:])):
#                rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]
#                print(rib_pos)
#
#                no_ribosomes[i, rib_pos.astype(int)] += 1
#        no_ribosomes = no_ribosomes[:, 1:]
#
#        ribosome_means = np.mean(no_ribosomes, axis=0)
#        ribosome_density = ribosome_means/len(k)
#
#        no_ribosomes_per_mrna = np.mean(no_ribosomes)

        ssa_obj = SSA_Soln()
        if ssa_conditions['record_stats']:
            ssa_obj.no_ribosomes = no_ribosomes
            ssa_obj.watched_ribs = watched_ribs
            ssa_obj.collisions = collisions

        ssa_obj.n_traj = n_traj
        ssa_obj.k = k
        #ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        #ssa_obj.rib_density = ribosome_density
        #ssa_obj.rib_means = ribosome_means

        ssa_obj.intensity_vec = all_results.T[:,startindex:,:]
        #ssa_obj.I = all_results.T[:,startindex:,:]

        ssa_obj.time = t
        ssa_obj.time_rec = t[startindex:]
        ssa_obj.start_time = non_consider_time
        ssa_obj.start_index = int(startindex)
        ssa_obj.photobleaching = photobleaching
        ssa_obj.photobleaching_rate = photobleaching_rate



        try:
            ssa_obj.col_points = all_col_points
        except:
            pass


        ssa_obj.eval_time = float(sttime)

        return ssa_obj


    def __solve_ssa_replicon(self,k, t, x0, n_traj, protease_locations, protease_rates, ssa_conditions=None,
                                  kon=1, koff=1, kprobe=[] ):
        seeds = np.random.randint(0, 0x7FFFFFF, n_traj, dtype=np.int32)

        if isinstance(k, list):
            k = np.array(k).astype(np.float64)

        k = k.flatten()

        if kprobe == []:
            kprobe = np.ones(self.color)

        if isinstance(kprobe, list):
            kprobe = np.array(kprobe, dtype=np.float64)
        if isinstance(kprobe, int):
            kprobe = np.array([kprobe], dtype=np.float64)
        if isinstance(kprobe, float):
            kprobe = np.array([kprobe], dtype=np.float64)


        if ssa_conditions == None:
            ssa_conditions = self.default_conditions

        x0 = self.__check_x0(x0)

        probe_vec = ssa_conditions['probe_vec']

        colors = self.colors


        if ssa_conditions['bursting'] == False:
            kon = 1
            koff = 1
        else:
            kon = kon
            koff = koff

        flags = np.array([int(ssa_conditions['bursting']),
                          int(ssa_conditions['leaky_probes']),
                          int(ssa_conditions['record_stats'])], dtype=np.int32)
        probe_loc = ssa_conditions['probe_loc']

        colors = int(probe_loc.shape[0])
        rib_vec = []
        solutions = []
        solutionssave = []


        N_rib = ssa_conditions['ribosome_preallocation']
        all_results,all_nribs,all_collisions,_,all_ribtimes,all_col_points = self.__generate_mats_lowmem(n_traj,k[0],t,N_rib,colors)
        all_rib_loc = np.zeros((n_traj, len(t), N_rib), dtype=np.int32)


        footprint = ssa_conditions['footprint']
        evf = ssa_conditions['perturb'][0]
        evi = ssa_conditions['perturb'][1]
        intime = ssa_conditions['perturb'][2]
        intime_stop = ssa_conditions['perturb'][3]
        non_consider_time = ssa_conditions['burnin']
        photobleaching_rate = ssa_conditions['photobleaching_rate']
        photobleaching = ssa_conditions['photobleaching']
        
         
        
        st = time.time()

        for i in range(n_traj):

            nribs = np.array([0], dtype=np.int32)
            result,ribtimes,_,coltimes,colpointsx,colpointst = self.__generate_vecs_lowmem(k,t,N_rib,colors)
            ribloc = np.zeros((N_rib,len(t)), dtype=np.int32)
            if i == 0: #detect any int64
                inputs = [ribloc,result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1], t, k[0], float(k[-1]), int(evf), int(evi), float(intime), float(intime_stop), seeds[i],nribs,x0,footprint, probe_vec ,int(colors), kon, koff, kprobe, probe_loc, flags,N_rib, photobleaching_rate, photobleaching, protease_locations, protease_locations]
                wash_inputs = self.__check_input_memview(inputs)

            if wash_inputs:
                #check memview so all given variables are in int32 if integer for C
                inputs = [ribloc,result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1], t, k[0], float(k[-1]), int(evf), int(evi), float(intime), float(intime_stop), seeds[i],nribs,x0,footprint, probe_vec ,int(colors), kon, koff, kprobe, probe_loc, flags,N_rib, photobleaching_rate, photobleaching, protease_locations, protease_locations]
                ribloc,result, ribtimes, coltimes, colpointsx,colpointst, kelong, t, ki, kt, evf, evi, intime,intime_stop, seed,nribs,x0,footprint, probe_vec ,colors, kon, koff, kprobe, probe_loc, flags, N_rib, photobleaching_rate, photobleaching, protease_locations, protease_locations = self.__check_memview(inputs)
            else:
                kelong = k[1:-1]
                ki = k[0]
                kt = float(k[-1])
                evf = int(evf)
                evi = int(evi)
                intime = float(intime)
                intime_stop = float(intime_stop)
                seed = seeds[i]
                colors = int(colors)

            
            ssa_translation_lowmem.run_SSA_replicon(ribloc, result, ribtimes, coltimes, colpointsx,colpointst, kelong, t, ki, kt, evf, evi, intime, intime_stop, seed,nribs,x0,footprint, probe_vec ,colors, kon, koff, kprobe, probe_loc, flags, photobleaching_rate, photobleaching,protease_locations, protease_rates)
            #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)

            all_results[i, :, :] = result.T
            #all_frapresults[i, :] = frapresult
            all_ribtimes[i, :] = ribtimes
            all_collisions[i, :] = coltimes
            all_nribs[i, :] = nribs
            all_rib_loc[i, :, :] = ribloc.T

            endcolrec = np.where(colpointsx == 0)[0][0]

            colpoints = np.vstack((colpointsx[:endcolrec],
                                   colpointst[:endcolrec]))
            all_col_points.append(colpoints.T)



            # for i in range(n_traj):
            #     soln = all_results[i,:, :]

            #     so = soln
            #     solutionssave.append(so)
            #     solutions.append(soln)

            collisions = np.array([[]])
            watched_ribs = []
            for i in range(n_traj):
                totalrib = all_nribs[i]

                if totalrib > all_collisions.shape[1]:
                    collisions = np.append(collisions, all_collisions[i][:])
                    watched_ribs.append(int(all_collisions.shape[1]))

                else:

                    collisions = np.append(collisions, all_collisions[i][:int(totalrib[0])])
                    watched_ribs.append(int(totalrib[0]))

            sttime = time.time() - st


        no_ribosomes = np.zeros((n_traj, (len(k)+1)))

        startindex = np.where(t >= non_consider_time)[0][0]

        all_col_points2 = []
        for i in range(len(all_col_points)):
            all_col_points2.append(all_col_points[i][ all_col_points[i][:, 1] > non_consider_time])


        rib_per_t = np.zeros((n_traj, len(t)))


        validind = 0
        riblocs = []
        for i in range(len(all_rib_loc)):
            try:


                if np.where(np.sum(all_rib_loc[i].T, axis=1)!=0)[0][-1] > validind:

                    validind = np.where(np.sum(all_rib_loc[i].T, axis=1)!=0)[0][-1]

                        #there was a blank
            except:
                pass


        all_rib_loc = all_rib_loc[:, :, :validind]

        for i in range(all_rib_loc.shape[0]):
            rib_loc = all_rib_loc[i, :, :]
            rt = np.count_nonzero(rib_loc.T,axis=0)
            rib_per_t[i, :] = rt

        rib_density_loc, rib_density_count = np.unique(all_rib_loc, return_counts=True)
        rib_density_loc = rib_density_loc[1:]
        rib_density_count = rib_density_count[1:]
        rib_density = rib_density_count / np.sum(rib_density_count)


        ssa_obj = SSA_Soln()
        ssa_obj.no_ribosomes = no_ribosomes
        ssa_obj.n_traj = n_traj
        ssa_obj.k = k
        # ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        ssa_obj.rib_per_t = rib_per_t.T[startindex:, :]
        ssa_obj.rib_mean = np.mean(rib_per_t.T[startindex:, :])

        ssa_obj.rib_density = rib_density
        # ssa_obj.rib_means = ribosome_means
        ssa_obj.intensity_vec = all_results.T[:, startindex:, :]

        ssa_obj.time = t
        ssa_obj.time_rec = t[startindex:]
        ssa_obj.start_time = non_consider_time
        ssa_obj.start_index = int(startindex)
        ssa_obj.watched_ribs = watched_ribs
        ssa_obj.collisions = collisions
        ssa_obj.ribosome_locations = all_rib_loc
        ssa_obj.ribtimes = all_ribtimes
        ssa_obj.photobleaching = photobleaching
        ssa_obj.photobleaching_rate = photobleaching_rate
        ssa_obj.protease_locations = protease_locations
        ssa_obj.protease_rates = protease_rates

        try:
            ssa_obj.col_points = all_col_points2
        except:
            pass


        ssa_obj.eval_time = sttime

        return ssa_obj


    def __solve_ssa_full_combined(self,k, t, x0, n_traj, ssa_conditions=None,
                                  kon=1, koff=1, kprobe=[] ):
        seeds = np.random.randint(0, 0x7FFFFFF, n_traj, dtype=np.int32)

        if isinstance(k, list):
            k = np.array(k).astype(np.float64)

        k = k.flatten()

        if kprobe == []:
            kprobe = np.ones(self.color)

        if isinstance(kprobe, list):
            kprobe = np.array(kprobe, dtype=np.float64)
        if isinstance(kprobe, int):
            kprobe = np.array([kprobe], dtype=np.float64)
        if isinstance(kprobe, float):
            kprobe = np.array([kprobe], dtype=np.float64)


        if ssa_conditions == None:
            ssa_conditions = self.default_conditions

        x0 = self.__check_x0(x0)

        probe_vec = ssa_conditions['probe_vec']

        colors = self.colors


        if ssa_conditions['bursting'] == False:
            kon = 1
            koff = 1
        else:
            kon = kon
            koff = koff

        flags = np.array([int(ssa_conditions['bursting']),
                          int(ssa_conditions['leaky_probes']),
                          int(ssa_conditions['record_stats'])], dtype=np.int32)
        probe_loc = ssa_conditions['probe_loc']

        colors = int(probe_loc.shape[0])
        rib_vec = []
        solutions = []
        solutionssave = []


        N_rib = ssa_conditions['ribosome_preallocation']
        all_results,all_nribs,all_collisions,_,all_ribtimes,all_col_points = self.__generate_mats_lowmem(n_traj,k[0],t,N_rib,colors)
        all_rib_loc = np.zeros((n_traj, len(t), N_rib), dtype=np.int32)


        footprint = ssa_conditions['footprint']
        evf = ssa_conditions['perturb'][0]
        evi = ssa_conditions['perturb'][1]
        intime = ssa_conditions['perturb'][2]
        intime_stop = ssa_conditions['perturb'][3]
        non_consider_time = ssa_conditions['burnin']
        photobleaching_rate = ssa_conditions['photobleaching_rate']
        photobleaching = ssa_conditions['photobleaching']
        
        st = time.time()

        for i in range(n_traj):

            nribs = np.array([0], dtype=np.int32)
            result,ribtimes,_,coltimes,colpointsx,colpointst = self.__generate_vecs_lowmem(k,t,N_rib,colors)
            ribloc = np.zeros((N_rib,len(t)), dtype=np.int32)
            if i == 0: #detect any int64
                inputs = [ribloc,result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1], t, k[0], float(k[-1]), int(evf), int(evi), float(intime), float(intime_stop), seeds[i],nribs,x0,footprint, probe_vec ,int(colors), kon, koff, kprobe, probe_loc, flags,N_rib, photobleaching_rate, photobleaching]
                wash_inputs = self.__check_input_memview(inputs)

            if wash_inputs:
                #check memview so all given variables are in int32 if integer for C
                inputs = [ribloc,result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1], t, k[0], float(k[-1]), int(evf), int(evi), float(intime), float(intime_stop), seeds[i],nribs,x0,footprint, probe_vec ,int(colors), kon, koff, kprobe, probe_loc, flags,N_rib, photobleaching_rate, photobleaching]
                ribloc,result, ribtimes, coltimes, colpointsx,colpointst, kelong, t, ki, kt, evf, evi, intime,intime_stop, seed,nribs,x0,footprint, probe_vec ,colors, kon, koff, kprobe, probe_loc, flags, N_rib, photobleaching_rate, photobleaching = self.__check_memview(inputs)
            else:
                kelong = k[1:-1]
                ki = k[0]
                kt = float(k[-1])
                evf = int(evf)
                evi = int(evi)
                intime = float(intime)
                intime_stop = float(intime_stop)
                seed = seeds[i]
                colors = int(colors)

            
            ssa_translation_lowmem.run_SSA_full(ribloc, result, ribtimes, coltimes, colpointsx,colpointst, kelong, t, ki, kt, evf, evi, intime, intime_stop, seed,nribs,x0,footprint, probe_vec ,colors, kon, koff, kprobe, probe_loc, flags, photobleaching_rate, photobleaching)
            #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)

            all_results[i, :, :] = result.T
            #all_frapresults[i, :] = frapresult
            all_ribtimes[i, :] = ribtimes
            all_collisions[i, :] = coltimes
            all_nribs[i, :] = nribs
            all_rib_loc[i, :, :] = ribloc.T

            endcolrec = np.where(colpointsx == 0)[0][0]

            colpoints = np.vstack((colpointsx[:endcolrec],
                                   colpointst[:endcolrec]))
            all_col_points.append(colpoints.T)



            # for i in range(n_traj):
            #     soln = all_results[i,:, :]

            #     so = soln
            #     solutionssave.append(so)
            #     solutions.append(soln)

            collisions = np.array([[]])
            watched_ribs = []
            for i in range(n_traj):
                totalrib = all_nribs[i]

                if totalrib > all_collisions.shape[1]:
                    collisions = np.append(collisions, all_collisions[i][:])
                    watched_ribs.append(int(all_collisions.shape[1]))

                else:

                    collisions = np.append(collisions, all_collisions[i][:int(totalrib[0])])
                    watched_ribs.append(int(totalrib[0]))

            sttime = time.time() - st


        no_ribosomes = np.zeros((n_traj, (len(k)+1)))

        startindex = np.where(t >= non_consider_time)[0][0]

        all_col_points2 = []
        for i in range(len(all_col_points)):
            all_col_points2.append(all_col_points[i][ all_col_points[i][:, 1] > non_consider_time])


        rib_per_t = np.zeros((n_traj, len(t)))


        validind = 0
        riblocs = []
        for i in range(len(all_rib_loc)):
            try:


                if np.where(np.sum(all_rib_loc[i].T, axis=1)!=0)[0][-1] > validind:

                    validind = np.where(np.sum(all_rib_loc[i].T, axis=1)!=0)[0][-1]

                        #there was a blank
            except:
                pass


        all_rib_loc = all_rib_loc[:, :, :validind]

        for i in range(all_rib_loc.shape[0]):
            rib_loc = all_rib_loc[i, :, :]
            rt = np.count_nonzero(rib_loc.T,axis=0)
            rib_per_t[i, :] = rt

        rib_density_loc, rib_density_count = np.unique(all_rib_loc, return_counts=True)
        rib_density_loc = rib_density_loc[1:]
        rib_density_count = rib_density_count[1:]
        rib_density = rib_density_count / np.sum(rib_density_count)


        ssa_obj = SSA_Soln()
        ssa_obj.no_ribosomes = no_ribosomes
        ssa_obj.n_traj = n_traj
        ssa_obj.k = k
        # ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        ssa_obj.rib_per_t = rib_per_t.T[startindex:, :]
        ssa_obj.rib_mean = np.mean(rib_per_t.T[startindex:, :])

        ssa_obj.rib_density = rib_density
        # ssa_obj.rib_means = ribosome_means
        ssa_obj.intensity_vec = all_results.T[:, startindex:, :]

        ssa_obj.time = t
        ssa_obj.time_rec = t[startindex:]
        ssa_obj.start_time = non_consider_time
        ssa_obj.start_index = int(startindex)
        ssa_obj.watched_ribs = watched_ribs
        ssa_obj.collisions = collisions
        ssa_obj.ribosome_locations = all_rib_loc
        ssa_obj.ribtimes = all_ribtimes
        ssa_obj.photobleaching = photobleaching
        ssa_obj.photobleaching_rate = photobleaching_rate


        try:
            ssa_obj.col_points = all_col_points2
        except:
            pass


        ssa_obj.eval_time = sttime

        return ssa_obj


    @classmethod
    def __get_ribosome_statistics(self, ssa_obj, result):


        return ssa_obj


    @classmethod
    def __generate_vecs(cls, k, t, N_rib, ncolor):
        tf = t[-1]
        ki = k[0]

        guessed_no_ribosomes = int(1.3*ki*tf)
        result = np.zeros((len(t)*N_rib), dtype=np.int32)
        ribtimes = np.zeros((guessed_no_ribosomes), dtype=np.float64)
        frapresult = np.zeros((len(t)*N_rib), dtype=np.int32)
        coltimes = np.zeros((guessed_no_ribosomes), dtype=np.int32)
        colpointsx = np.zeros(len(k[1:-1])*(guessed_no_ribosomes), dtype=np.int32)
        colpointst = np.zeros(len(k[1:-1])*(guessed_no_ribosomes), dtype=np.float64)
        return result, ribtimes, frapresult, coltimes, colpointsx, colpointst

    @classmethod
    def __generate_vecs_trna(cls, kbind, kind, t, N_rib, ncolor):
        tf = t[-1]
        ki = kbind

        guessed_no_ribosomes = int(1.3*ki*tf)
        result = np.zeros((len(t)*N_rib), dtype=np.int32)
        ribtimes = np.zeros((guessed_no_ribosomes), dtype=np.float64)
        frapresult = np.zeros((len(t)*N_rib), dtype=np.int32)
        coltimes = np.zeros((guessed_no_ribosomes), dtype=np.int32)
        colpointsx = np.zeros(len(kind)*(guessed_no_ribosomes), dtype=np.int32)
        colpointst = np.zeros(len(kind)*(guessed_no_ribosomes), dtype=np.float64)
        return result, ribtimes, frapresult, coltimes, colpointsx, colpointst

    @classmethod
    def __generate_vecs_lowmem(cls, k, t, N_rib, ncolor):
        tf = t[-1]
        ki = k[0]

        guessed_no_ribosomes = int(1.3*ki*tf)
        result = np.zeros((ncolor,len(t)), dtype=np.int32)
        ribtimes = np.zeros((guessed_no_ribosomes), dtype=np.float64)
        frapresult = np.zeros((len(t)*N_rib), dtype=np.int32)
        coltimes = np.zeros((guessed_no_ribosomes), dtype=np.int32)
        colpointsx = np.zeros(len(k[1:-1])*(guessed_no_ribosomes), dtype=np.int32)
        colpointst = np.zeros(len(k[1:-1])*(guessed_no_ribosomes), dtype=np.float64)
        return result, ribtimes, frapresult, coltimes, colpointsx, colpointst

    @classmethod
    def __generate_mats_lowmem(cls, ntraj, ki, t, N_rib, ncolor):
        tf = t[-1]
        guessed_no_ribosomes = int(1.3*ki*tf)
        all_results = np.zeros((ntraj, len(t),ncolor), dtype=np.int32)
        all_ribtimes = np.zeros((ntraj, guessed_no_ribosomes), dtype=np.float64)
        all_frapresults = np.zeros((ntraj, N_rib*len(t)),dtype=np.int32)
        all_collisions = np.zeros((ntraj, guessed_no_ribosomes), dtype=np.int32)
        all_nribs = np.zeros((ntraj, 1))
        all_col_points = []

        return all_results, all_nribs, all_collisions, all_frapresults, all_ribtimes, all_col_points


    @classmethod
    def __generate_mats(cls, ntraj, ki, t, N_rib, ncolor):
        tf = t[-1]
        guessed_no_ribosomes = int(1.3*ki*tf)
        all_results = np.zeros((ntraj, N_rib*len(t)), dtype=np.int32)
        all_ribtimes = np.zeros((ntraj, guessed_no_ribosomes), dtype=np.float64)
        all_frapresults = np.zeros((ntraj, N_rib*len(t)), dtype=np.int32)
        all_collisions = np.zeros((ntraj, guessed_no_ribosomes), dtype=np.int32)
        all_nribs = np.zeros((ntraj, 1))
        all_col_points = []




        return all_results, all_nribs, all_collisions, all_frapresults, all_ribtimes, all_col_points

    @classmethod
    def __generate_mats_lowmem_nostats(cls, ntraj, ki, t, N_rib, ncolor):

        all_results = np.zeros((ntraj, len(t),ncolor), dtype=np.int32)
        all_frapresults = np.zeros((ntraj, N_rib*len(t)), dtype=np.int32)
        return all_results, all_frapresults

    @classmethod
    def __generate_vecs_lowmem_nostats(cls, k, t, N_rib, ncolor):


        result = np.zeros((ncolor,len(t)), dtype=np.int32)

        frapresult = np.zeros((len(t)*N_rib), dtype=np.int32)

        return result, frapresult

    def __check_x0(self, x0):
        if len(x0) == 0:
            x0_clean = np.zeros((200), dtype=np.int32)
        else:
            if len(x0) >200:
                raise ValueError('Unrecognized initial condition, make sure the length of x0 is <= 200')

            x0_clean = np.zeros((200), dtype=np.int32)
            x0_clean[:len(x0)] = x0
        return x0_clean

    def __check_rates_trna(self, rates):
        if isinstance(rates, np.ndarray):
            pass
        else:
            rates = np.array(rates, dtype=np.int32)
        if rates.dtype !=int:
            raise ValueError('trna ID values are not integers, to run the trna ssa version, k_index needs to be an ID of 0-60 integer values')
        if len(np.where(rates < 0)[0]) or len(np.where(rates > 60)[0]):
            raise ValueError('trna ID values are out of bounds, the IDs are 0-60 for tRNA species')

    def __check_input_memview(self, inputs):
        check = False
        for item in inputs:
            try:
                item.dtype
                if item.dtype not in [np.int32]:
                    if 'int' in str(item.dtype):

                        check = True
            except:
                pass
        return check

    def __check_memview(self, inputs):
        #make sure everything is int 32, no int64
        returned_inputs = ()
        for item in inputs:
            try:
                item.dtype
                if item.dtype not in [np.int32]:
                    if 'int' in str(item.dtype):
                        item = item.astype(np.int32)
                returned_inputs = returned_inputs + (item,)

            except:
                returned_inputs = returned_inputs + (item,)
        return returned_inputs


    def __check_rates(self, rates):
        if isinstance(rates,np.ndarray):
            pass
        else:
            rates = np.array(rates)
        if len(np.where(rates < 0)[0]) > 0:
            raise ValueError('one or more model rates are negative, double check the provided rates')


    def __ssa_python(self, k, t_array,
                     inhibit_time=0, FRAP=False, Inhibitor=False,
                     flags=None, kon=1, koff=1, kprobe=1, ssa_conditions=None):
        '''
        mRNA Translation simulation python implementation

        given a propensity vector k, time array to record, 
        and inhibitory conditions, run a single trajectory
        of translation simulation

        The simulation is described here: [PUT LINK HERE TO PAPER]

        *args*

            **k**, propensity vector of size gene length + 2,
            [initiation rate,  Codon dependent rates,  completion rate / unbinding rate]
            for reference the codon dependent rates are refering to the time rate of a ribosome to move on to the next codon

            **t_array**, time points to record the ribosome posistions at

        *keyword args*

            **inhibit_time**, the time to start inhibition assays if FRAP or Inhibitor (harringtonine) == True

            **FRAP**, True or false to apply Fluorescence Recovery After Photobleaching (FRAP) https://en.wikipedia.org/wiki/Fluorescence_recovery_after_photobleaching

            **Inhibitor**, True or false to apply harringtonine at inhibit_time. Harringtonine acts as a protien translation initiation inhibitor

        '''

        #SSA params and propensities
        R = self.default_conditions['footprint'] #exclusion volume (ribosome footprint), ribosomes cant be less than 10 codons apart because of their physical size
        kelong = np.array([k[1:-1]]).T  #rates for ribosomes moving to the next codon, based on tRNA concentrations

        N = len(kelong)  #Number of codons in the mRNA
        kbind = k[0]   #rate for a ribosome to bind and start translation
        kcompl = k[-1]     #rate for a ribosome at the end of the mRNA to unbind
        X = np.array([0, 0], dtype=int)   #the updating ribosome posistion vector that is changed in the simulation

        bursting = flags[0]
        leaky = flags[1]


        if bursting:
            burst = np.random.rand() < (kon/(kon+koff) )
        else:
            burst = 1

        Ncol = np.zeros((1,0))

        N_rib = 200  #Maximum number of ribosomes on a single mRNA (hard limit for the simulation not a physical constant)


        colors = ssa_conditions['probe_loc'].shape[0]
        if leaky:
            leaky_probe_matrix = np.zeros((colors, N, N_rib))
        intensity = np.zeros((colors, len(t_array)))
        probe = np.array(np.where(ssa_conditions['probe_loc'] ==1))
        probevec = ssa_conditions['probe_vec']


        #example X arrays and how its formatted:
        # X = [423 30 10 0 ]  read from left to right theres a ribosome in position 423 30 and 10, with a 0 kept as a buffer for simulation

        t = t_array[0]  #time point
        Nt = len(t_array)  #number of time points to record over
        tf = t_array[-1]  #final time point
        col = np.zeros((1, N_rib))
        X_array = np.zeros((N_rib, Nt))  #recording array that records the ribosome posistions over time array points
        NR = 0  #number of ribosomes bound
        it = 1  #number of iterations
        Sn_p = np.eye(max(NR+1, 2), dtype=int) #stoichiometry for the SSA
        if bursting:
            wn_p = np.zeros((X.shape[0]+1, 1)) # propensities for the SSA
        else:
            wn_p = np.zeros((X.shape[0], 1)) # propensities for the SSA

        T = np.array([0, 0], dtype=float)
        ribtimes = np.array([[0,0]], dtype=float)
        col_points = []
        #wn_p = np.zeros((1,X.shape[0])).flatten()
        wshape = len(wn_p)
        Inhibit_condition = 1  #set up inhibitor flags
        while t < tf:


            if Inhibitor == True:
                if t >= inhibit_time:

                    Inhibit_condition = 0
                else:

                    Inhibit_condition = 1
            else:
                Inhibit_condition = 1
            if FRAP == True :   #if the Photobleaching is happening, "remove" ribosomes
                if t >= inhibit_time and t < inhibit_time + 20:
                    #X = np.array([0, 0])
                    a=1
                    #T = np.array([0,0])

            oldNR = NR
            NR = len(np.flatnonzero(X)) #each iteration get the number of ribosomes on the mRNA

            if X.shape[0] < NR+1:  #if the last reaction added a ribosome put a 0 on the end of X vec

                X = np.append(X, [0])
                T = np.append(T, [0])
                T[-2] = t



            X[-1] = 0
            T[-1] = 0

            X = X[0:max(NR, 1)+1]  #clear any additional 0's on the end
            T = T[0:max(NR, 1)+1]

            if oldNR != NR:     #if the number of ribosomes has changed reallocate the sizes of stoich and propensities

                if not bursting:
                    Sn_p = np.eye(max(NR+1, 2), dtype=int)
                    wn_p = np.zeros((X.shape[0], 1))
                else:
                    Sn_p = np.eye(max(NR+1, 2), dtype=int)
                    wn_p = np.zeros((X.shape[0]+1, 1))

                if leaky:
                    leaky_probe_matrix[probe[0], probe[1], NR-1] = (np.random.rand(len(probe[0])) < kprobe).astype(int)


                wshape = len(wn_p)


            Sn = Sn_p
            wn = wn_p


            #get indices of where X vecs are > 0 ie where the ribosome values are
            inds = X > 0

            if bursting:

                    wn[:-1][inds] = kelong[X[inds]-1]  #update propensities



            else:
                wn[inds] = kelong[X[inds]-1]


            if X[0] == N:  #if the ribosome in the 0 position is at the end of the mRNA set propensities to the reaction for completion

                Sn[:, 0] = (np.append(X[1:], np.array([0]))-X[0:])


                wn[0] = kcompl


            #if there are no ribosomes or when there is enough room for a new ribosome to bind add the propensity for binding
            if NR == 0:

                wn[NR] = kbind*Inhibit_condition*burst



            if NR > 0 and X[NR-1] > R:
                wn[NR] = kbind*Inhibit_condition*burst

            REST = np.less(X[1:]+10, X[0:-1])  #apply the footprint condition ie set any propensities where it violates the > 10 codons apart rule to 0

            if bursting:
                if NR > 1:
                    wn[1:-1] = (wn[1:-1].T*REST).T
                else:
                    wn[1:] = (wn[1:].T*REST).T
            else:
                wn[1:] = (wn[1:].T*REST).T  #apply that logical^ to propensities

            if bursting:
                if burst:
                    wn[-1] = koff
                else:
                    wn[-1] = kon

            w0 = sum(wn.flat)  #get the sum of propensities
            randnum = np.random.random_sample(2) #update time to point of next reaction (exponential waiting time distb)
            t = (t-np.log(randnum[0])/w0)

            while it < Nt and t > t_array[it]:  #record state if past timepoint
                X_array[0:len(X), it] = X

                if leaky:
                    int_tmp = np.zeros(colors)
                    validx = X[X>0]
                    c = leaky_probe_matrix[:, :, 0]
                    for ribind in validx:
                        int_tmp += np.sum(c[:, :(ribind-1)],axis=1)

                    intensity[:, it] = int_tmp

                else:
                    validx = X[X>0]
                    int_tmp = np.zeros(colors)
                    for ribind in validx:
                        int_tmp += probevec[:, (ribind-1)]
                    intensity[:, it] = int_tmp

                it += 1



            if t < tf:  #if still running simulation pick which reaction happened via random number and propensity sum
                r2 = w0*randnum[1]
                tmp = 0

                for i in range(wshape):

                    tmp = tmp + wn[i]
                    if tmp >= r2:
                        event = i
                        break


            if bursting:
                if event >= NR+1:
                    burst+=1
                    burst = burst%2
                else:
                    X = (X + Sn[:, event].T)

                    if np.sum(Sn[:, event]) < 0 :

                        ribtimes = np.vstack((ribtimes, [T[0],t]))
                        T[:-1] = T[1:]
                        Ncol = np.append(Ncol, col[0][0] )
                        col = np.atleast_2d(np.append(col[:,1:],[0]))

                    else:
                        if X[event-1] == X[event] + R:
                            col[0][event] +=1
                            col_points.append((X[event], t))

            else:
                X = (X + Sn[:, event].T)  #update X vector for new ribosome state


                if np.sum(Sn[:,event]) < 0 :

                    ribtimes = np.vstack((ribtimes, [T[0], t]))
                    T[:-1] = T[1:]
                    Ncol = np.append(Ncol, col[0][0])
                    col = np.atleast_2d(np.append(col[:, 1:],[0]))

                    if leaky: #remove the leaky probes
                        leaky_probe_matrix[probe[0], probe[1], 0] = leaky_probe_matrix[probe[0], probe[1], 1]


                else:
                    if X[event-1] == X[event] + R:
                        col[0][event] +=1
                        col_points.append((X[event], t))



        return X_array,ribtimes[1:, :], Ncol, col_points, intensity  #return the completed simulation

