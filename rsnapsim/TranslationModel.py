'''
import os
cwd = os.getcwd()
os.chdir('../../')
print(os.getcwd())
#import rsnapsim as rss
#from rsnapsim import seqmanip
#from rsnapsim import SSA_Soln
#from rsnapsim import GenericMetaData

import inspect
import numpy as np
import time
os.chdir(cwd)
import unittest
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

'''

import os
import numpy as np
from . import custom_errors as custom_err
from . import ModelMaker as mm
from itertools import product
import inspect
import importlib

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
import matplotlib.patches as mpatches
from matplotlib.patches import PathPatch
from matplotlib.path import Path

import pathlib

import json


class TranslationModel:
    def __init__(self, name, mRNA=None, particle_size=9):
        self.name = name
        self.blank(mRNA, particle_size)
        
        
    def blank(self, mRNA, particle_size):
        
        if mRNA is None:
            self._length=0
        else:
            self._length = mRNA.total_length
        
        self.particle_size = particle_size
        self.footprint = particle_size
        self._max_particles = int(self._length/particle_size+5)
        
        #self._particle_array = np.zeros([self._max_particles, 4 + self._n_colors], dtype=np.int32)
        self._state_mat = np.zeros([], dtype=np.int32)
        #self._kelong_mat = np.zeros([1, self._length+1])
        #self._kelong_mat[0,:-1] = mRNA.kelong
        self._probe_mat = np.zeros([3, self._length])
        #mRNA.generate_3frame_tags()
        k = 1
        used_tags = []
        if mRNA is None:
            pass
        else:
            for i in range(3):
                for key in mRNA.multiframe_epitopes[i]:
                    if key in used_tags:
                        self._probe_mat[i,mRNA.multiframe_epitopes[i][key]] = used_tags.index(key)+1
                    else:
                        self._probe_mat[i,mRNA.multiframe_epitopes[i][key]] = len(used_tags)+1    
                        used_tags.append(key)
                        
        self._n_colors= int(np.max(self._probe_mat))
        self._resource_mat = np.zeros([], dtype=np.int32)

        self._propensities = []
        self._parameters = []
        
        self._prop_ids = []
        
        self._constant_rxn_bool = [] #whether a reaction is per ribosome or not
        self._constant_reactions = []
        self._ribosome_reactions = []
        
        self._n_rxns = 0
        self._n_pars = 0
        self._n_kelong = 1
        self._n_frames = 1
        self._n_states = 0
        self._n_resources = 0
                                # rxn type, updates_states_resources?, frame, loc, dexist, dframe, dloc, 
                                # dprobe1, dprobe2... dprobeN, dstate1, dstate2... dstateN,
                                # dresource1, dresource2 ... dresourceN, 
        self._rxn_mat = np.zeros([0, 7 + self._n_colors+ self._n_states + self._n_resources])
        
        self.__rxn_size = self._rxn_mat.shape[1]
        
        self.__default_colors = ['#07BEB8', '#8F3985', '#5EDA6A', '#FA5E5E']
        self._state_names = []
        self._resource_names = []
        
        self._rib_arr0 = np.zeros([self._max_particles,  4+self._n_colors+self._n_rxns], dtype=int)
        self._lattice_arr0 = np.zeros([self._length], dtype=int)
        self._state_arr0 = np.zeros([self._n_states], dtype=int)
        self._resource_arr0 = np.zeros([self._n_resources], dtype=int)
        self._rxn_names = []
        
        
        self._codon_mat = np.zeros([3, self._length], dtype=int)
        
        if mRNA is None:
            pass
        else:
        
            codons  = [''.join(i) for i in product(['A','U','G','C'], repeat = 3)]
            codon_ids = [[ codons.index(line.upper().replace('T','U')[i:i+3]) for i in range(0, len(line), 3)] for line in mRNA.multiframe_nt_seq]
    
            
            
            self._codon_mat[0,:] = codon_ids[0]
            self._codon_mat[1,:-1] = codon_ids[1]
            self._codon_mat[2,:-1] = codon_ids[2]
            
        self._probe_function = lambda k,t,p,ke,o,l,pr,s,r,nr: 1
        self._probe_parameters = 0

    def add_states(self, n_states, state0=[], names=[]):
        self._rxn_mat = np.hstack([self._rxn_mat[:, :7+self._n_colors],
                                   np.zeros([self._rxn_mat.shape[0], n_states]),
                                   self._rxn_mat[:, 7+self._n_colors+self._n_states:]])

        
        if len(state0) == n_states:
            self._state_mat = np.array(state0).flatten().astype(int)
        else:
            print('warning: initializing states with a blank vector (provided initial state does not have correct number of states)')
            self._state_mat = np.zeros([n_states], dtype=int)
        self._state_names = names
        self.update()
        
    
    def add_resources(self, n_resources, resources0=[], names=[]):
        self._rxn_mat = np.hstack([self._rxn_mat[:, :7+self._n_colors+self._n_states],
                           np.zeros([self._rxn_mat.shape[0], n_resources]),
                           self._rxn_mat[:, 7+self._n_colors+self._n_states+self._n_resources:]])
        
        if len(resources0) == n_resources:
            self._resource_mat = np.array(resources0).flatten().astype(int)
        else:
            print('warning: initializing resources with a blank vector (provided initial state does not have correct number of states)')
            self._resource_mat = np.zeros([n_resources], dtype=int)
        self._resource_names = names
        self.update()
        
    

        
        
    ## Everything that can happen in the model
    def add_ribosome_reaction(self, propensity, rxn_name='', dexist = 0, dframe = 0, dloc = 0,
                              dprobes=[], dprobe_inds=[],
                              dstates=[], dstate_inds=[],
                              dresources=[], dresource_inds=[]):
        
        # update the flag for resources / state changes
        updates_states_resources = 0
        if np.any(dstates):
            updates_states_resources = 1
        if np.any(dresources):
            updates_states_resources = 1

        self.__make_rib_or_lattice_rxn(0, updates_states_resources, 0,0, dexist, dframe, dloc, dprobes, dprobe_inds, dstates, dstate_inds, dresources, dresource_inds)
        self.__add_propensity(propensity, rxn_name=rxn_name)
        
        
    def add_lattice_reaction(self, propensity, rxn_name='', frame=0, loc=0, dexist = 0, dframe = 0, dloc = 0,
                              dprobes=[], dprobe_inds=[],
                              dstates=[], dstate_inds=[],
                              dresources=[], dresource_inds=[]):
        
        # update the flag for resources / state changes
        updates_states_resources = 0
        if np.any(dstates):
            updates_states_resources = 1
        if np.any(dresources):
            updates_states_resources = 1
            
        self.__make_rib_or_lattice_rxn(2, updates_states_resources, frame,loc, dexist, dframe, dloc, dprobes, dprobe_inds, dstates, dstate_inds, dresources, dresource_inds)
        self.__add_propensity(propensity, rxn_name=rxn_name)
        

    def add_state_reaction(self, propensity, rxn_name='', inds=[], dstates=[]):
                # if the user did not provide indexes, pad with zeros for right shape
        if len(inds) < len(dstates):
            change = dstates + [0,]*(len(dstates)-len(self._n_states))
            self.__add_rxn(1, [0, 0, 0, 0, 0, 0,] +   [0]*self._n_colors +  change  )
            self.__add_propensity(propensity, rxn_name=rxn_name)
            return
            
        # if the user did  provide indexes, use them to build the reaction row
        change = [0]*self._n_states
        for i in range(self._n_states):
            if i in inds:
                change[i] = dstates[inds.index(i)]
        self.__add_rxn(1, [0, 0, 0, 0, 0, 0,] +  [0]*self._n_colors +  change   )
        self.__add_propensity(propensity, rxn_name=rxn_name)


    def _x0(self):
        return self._rib_arr0, self._lattice_arr0, self._state_arr0, self._resource_arr0

    def __make_rib_or_lattice_rxn(self, rtype, updates_states_resources, frame, loc, dexist, dframe, dloc, dprobes, dprobe_inds, dstates, dstate_inds, dresources, dresource_inds):
        change1 = [updates_states_resources, frame, loc, dexist, dframe, dloc,]
                ### calculate and add probe related changes
        if len(dprobe_inds) < len(dprobes):
            change = dprobes + [0,]*(len(dprobes)-len(self._n_colors))
        # if the user did  provide indexes, use them to build the reaction row
        else:
            change = [0]*self._n_colors
            for i in range(self._n_colors):
                if i in dprobe_inds:
                    change[i] = dprobes[dprobe_inds.index(i)]
        change1 = change1 + change
        
        ### calculate and add state related changes
        if len(dstate_inds) < len(dstates):
            change = dstates + [0,]*(len(dstates)-len(self._n_states))
        # if the user did  provide indexes, use them to build the reaction row
        else:
            change = [0]*self._n_states
            for i in range(self._n_states):
                if i in dstate_inds:
                    change[i] = dstates[dstate_inds.index(i)]
        change1 = change1 + change
        
        ### calculate and add state related changes

        if len(dresource_inds) < len(dresources):
            change = dstates + [0,]*(len(dresources)-len(self._n_resources))
        # if the user did  provide indexes, use them to build the reaction row
        else:
            change = [0]*self._n_resources
            for i in range(self._n_resources):
                if i in dresource_inds:
                    change[i] = dresources[dresource_inds.index(i)]

        change1 = change1 + change
                
        # finally add the reaction to the reaction matrix
        self.__add_rxn(rtype, change1)


    def add_resource_reaction(self, propensity, rxn_name='', dresources=[], inds=[] ):
        
        # if the user did not provide indexes, pad with zeros for right shape
        if len(inds) < len(dresources):
            change = dresources + [0,]*(len(dresources)-len(self._n_resources))
            self.__add_rxn(3, [0, 0, 0, 0, 0, 0,] + [0]*self._n_colors+ [0]*self._n_states +  change )
            self.__add_propensity(propensity, rxn_name=rxn_name)
            return
            
        # if the user did  provide indexes, use them to build the reaction row
        change = [0]*self._n_resources
        for i in range(self._n_resources):
            if i in inds:
                change[i] = dresources[inds.index(i)]
        
        self.__add_rxn(3, [0, 0, 0, 0, 0, 0,] + [0]*self._n_colors+ [0]*self._n_states +  change   )
        self.__add_propensity(propensity, rxn_name=rxn_name)

    def add_probe_reaction(self, propensity, rxn_name='', ind=0, dprobe=0):
        self.__add_rxn(4, [ind, dprobe],)
        self.__add_propensity(propensity, rxn_name=rxn_name)
    
    def add_probe_reaction(self, propensity, rxn_name='', inds=[], dprobes=[]):
                # if the user did not provide indexes, pad with zeros for right shape
        if len(inds) < len(dstates):
            change = dstates + [0,]*(len(dstates)-len(self._n_states))
            self.__add_rxn(3, [0, 0, 0, 0, 0, 0,] +  change + [0]*self._n_states  + + [0]*self._n_resources )
            self.__add_propensity(propensity, parameters, rxn_name=rxn_name)
            return
            
        # if the user did  provide indexes, use them to build the reaction row
        change = [0]*self._n_states
        for i in range(self._n_states):
            if i in inds:
                change[i] = dstates[inds.index(i)]
        self.__add_rxn(1, [0, 0, 0, 0, 0, 0,] +  [0]*self._n_colors +  change   )
        self.__add_propensity(propensity, parameters, rxn_name=rxn_name)
    
    def add_global_reaction(self, propensity, rxn_name=''):
        self.__add_rxn(5, [],)
        self.__add_propensity(propensity, rxn_name=rxn_name)
        
    def delete_reaction(self, row):
        self.__delete_rxn(row)
        
    def __add_propensity(self, propensity_fun, rxn_id=-1, rxn_name=''):
        if rxn_id == -1:
            try:
                rxn_id = max(self._prop_ids)
            except:
                rxn_id = 0
        self._propensities = self._propensities + [propensity_fun, ]
        #self._parameters = self._parameters + [parameters, ]
        self._prop_ids = self._prop_ids + [rxn_id, ]
        self._rxn_names = self._rxn_names + [rxn_name]
        
    
    def set_parameters(self,parameters):
        self._parameters = parameters 
        
    @property
    def _parameters_full(self):
        return [self._parameters for x in range(len(self._rxn_names))]
            
        
    def add_constant_propensities(self, constant_propensities, parameters, rxn_ids=None):
        self._constant_props = constant_propensities
        self._constant_parameters = parameters
        self._constant_ids = rxn_ids


    def add_ribosome_propensities(self, ribosome_propensities, parameters, rxn_ids=None):
        self._ribosome_props = ribosome_propensities
        self._ribosome_parameters = parameters
        self._ribosome_prop_ids = rxn_ids

    
    def add_reaction_matrix(self, matrix):
        return 1
    

    def __check_nframes_used(self):
        nframes = 1
        for i in range(len(self._rxn_mat)):
            if self._rxn_mat[i,0] == 2:
                nframes = max(nframes, self._rxn_mat[i,2])
        return nframes
    
    
    
    def __add_rxn(self, rxn_type, rxn_list,):
        row_list = [rxn_type] + rxn_list
        if (self.__rxn_size-len(row_list))%(self.__rxn_size) > 0:
            row_list = row_list + [0,]*((self.__rxn_size-len(row_list))%(self.__rxn_size))
        self._rxn_mat = np.vstack((self._rxn_mat, row_list))
        self._constant_rxn_bool = self._constant_rxn_bool + [rxn_type!=0]
        self.update()
        
    
    def generate_mats(self):
        return 1
        
    def check_validity(self):
        return 1

    def update(self):

        self._n_rxns = int(self._rxn_mat.shape[0])
        self.__rxn_size = int(self._rxn_mat.shape[1])
        self._n_pars = len(self._parameters)
        self._n_kelong = self._kelong_mat.shape[0]
        self._n_frames = int(self._kelong_mat.shape[0])
        try:           
            self._n_states = int(self._state_mat.shape[0])
        except:
            self._n_states = 0
        try:
            self._n_resources = int(self._resource_mat.shape[0])
        except:
            self._n_resources = 0
        self._n_colors= int(np.max(self._probe_mat))
        
        self._rib_arr0 = np.zeros([self._max_particles,  4+self._n_colors+ self._n_rxns], dtype=int)
        self._lattice_arr0 = np.zeros([self._length], dtype=int)
        self._state_arr0 = np.zeros([self._n_states], dtype=int)
        self._resource_arr0 = np.zeros([self._n_resources], dtype=int)
        
        
        self._constant_reactions =[]
        self._ribosome_reactions =[]
        for i in range(self._n_rxns):
            if self._constant_rxn_bool[i]:
                self._constant_reactions = self._constant_reactions + [i,]
            else:
                self._ribosome_reactions = self._ribosome_reactions + [i,]
    
    @property
    def model_id(self):
        hash_list = []
        
        id_list = list(sorted(self.__dict__.keys()))
        id_str = 'M'
        id_str += '_' + str(self._n_colors)
        id_str += '_' + str(self._n_states)
        id_str += '_' + str(self._n_resources)
        id_str += '_' + str(self._length)
        id_str += '_' + str(self._n_rxns)
        id_str += '_' + str(len(self._ribosome_reactions))
        id_str += '_' + str(self.particle_size)
        a = ''.join([str(x) for x in self._codon_mat.flatten()])  + ''.join([str(x) for x in self._rxn_mat.flatten()]) + ''.join([str(x) for x in self._kelong_mat.flatten()])
        id_str += '_' + hex(hash(a))
        
        id_str += '_' + hex(hash(''.join([str(inspect.getsourcelines(x)[0]) for x in self._propensities]) ))
            

        return id_str

    def __load_model_c(self,):
        self.cmodel = importlib.import_module('rsnapsim.models.%s.%s'%(self.name,self.name))
        #if self.model_id != cmodel.__model_id:
           # print('error')
        

    def compile_model_c(self, overwrite=True):
        '''
        propensities_str_list = [y[0].replace('\n','') for y in [inspect.getsourcelines(x)[0] for x in self._propensities]]
        propensity_names_list = [y[0].replace('\n','').split('=')[0].replace(' ','') for y in [inspect.getsourcelines(x)[0] for x in self._propensities]]
        propensity_function_list = [y[0].replace('\n','').split('=')[1].replace('lambda k,t,p,ke,o,l,pr,s,r,nr: ','') for y in [inspect.getsourcelines(x)[0] for x in self._propensities]]
        
        constant_propensities_strs = [propensities_str_list[x] for x in self._constant_reactions]
        ribosome_propensities_strs = [propensities_str_list[x] for x in self._ribosome_reactions]
        
        constant_propensities_c = RuleConverterLambda().make_c_propensities([self._propensities[x] for x in self._constant_reactions])
        ribosome_propensities_c = RuleConverterLambda().make_c_propensities([self._propensities[x] for x in self._ribosome_reactions], ribosome=1)
        original_rules = '/n'.join(propensities_str_list)
        
        #ModelFactory.compile_model(mRNA)
        
        '''
        
        return_code = mm.ModelFactory().compile_model(self.name, self.model_id, [self._propensities[x] for x in self._constant_reactions],
                                     [self._propensities[x] for x in self._ribosome_reactions], overwrite=overwrite )        
        
        if return_code == 0:
            self.cmodel = importlib.import_module('rsnapsim.models.%s.%s'%(self.name,self.name))
            
            # save a json copy of this model object along with the compiled model.

            self.save(os.path.join(os.path.dirname(__file__),'models',self.name, self.name + '.json'))

    
    def save(self, fname):
        
        save_dict = {}
        for i, entry in enumerate(self.__dict__.items()):
            key, value = entry
            
            if isinstance(self.__dict__[key], np.ndarray):
                save_dict[key] = self.__dict__[key].tolist()
                
            elif key == '_propensities':
                propensities_str_list = [y[0].replace('\n','') for y in [inspect.getsourcelines(x)[0] for x in self.__dict__[key]]]
                propensity_names_list = [y[0].replace('\n','').split('=')[0].replace(' ','') for y in [inspect.getsourcelines(x)[0] for x in self.__dict__[key]]]
                
                
                
                propensity_function_list = ['='.join(y[0].replace('\n','').split('=')[1:]).replace('lambda k,t,p,ke,o,l,pr,s,r,nr: ','') for y in [inspect.getsourcelines(x)[0] for x in self.__dict__[key]]]
                
                save_dict['_propensities_str_list'] = propensities_str_list
                save_dict['_propensity_names_list'] = propensity_names_list
                save_dict['_propensity_function_list'] = propensity_function_list
                
                pass
            
            elif key == '_probe_function':
                pf_str = [y[0].replace('\n','') for y in [inspect.getsourcelines(x)[0] for x in [self.__dict__[key] ] ]]
                pf_name = [y[0].replace('\n','').split('=')[0].replace(' ','') for y in [inspect.getsourcelines(x)[0] for x in [self.__dict__[key]]]]
                pf_function = ['='.join(y[0].replace('\n','').split('=')[1:]).replace('lambda k,t,p,ke,o,l,pr,s,r,nr: ','') for y in [inspect.getsourcelines(x)[0] for x in [self.__dict__[key]]]]
                
                save_dict['_probe_function_str'] = pf_str[0]
                save_dict['_probe_function_name'] = pf_name[0]
                save_dict['_probe_function_function'] = pf_function[0]
                
            elif key == 'cmodel':
                pass
                
            else:
                save_dict[key] = self.__dict__[key]
                
        for key in save_dict.keys():
            print(type(save_dict[key]))
        with open(fname, 'w') as fp:
            json.dump(save_dict, fp)
        
        return
    
    def load(self, fname,):
        
        with open(fname, 'r') as fp:
            load_dict = json.load(fp)
        
        propensities = []
        for i in range(len(load_dict['_propensity_function_list'])):
            propensities.append(load_dict['_propensity_function_list'][i])
        propensities = [eval('lambda k,t,p,ke,o,l,pr,s,r,nr:' + x) for x in propensities]
        
        load_dict['_probe_function'] = eval("lambda k,t,p,ke,o,l,pr,s,r,nr:" + load_dict['_probe_function_function'])
        
        load_dict['_propensites'] = propensities
        
        del load_dict['_probe_function_name']
        del load_dict['_probe_function_str']
        del load_dict['_probe_function_function']
        
        del load_dict['_propensities_str_list']
        del load_dict['_propensity_names_list']
        del load_dict['_propensity_function_list']
        
        convert_to_array_list_int32 = ['_state_mat', '_resource_mat', '_probe_mat',
                                       '_rxn_mat', '_rib_arr0', '_lattice_arr0', '_state_arr0',
                                       '_resource_arr0', '_codon_mat']
        for key in convert_to_array_list_int32:
            load_dict[key] = np.array(load_dict[key], dtype=np.int32)
        
        convert_to_array_list_float = [ '_kelong_mat']

        for key in convert_to_array_list_float:
            load_dict[key] = np.array(load_dict[key], dtype=float)
            
            
        
        self.__dict__ = load_dict
        
        
        return 
    

    def load_saved_cmodel(self, c_model_name):
        
        self.load(os.path.join(os.path.dirname(__file__),'models',c_model_name, c_model_name + '.json'))
        self.__load_model_c()
    
    def __get_mats(self):
        return self.kelong_mat, self.resource_mat, self.rxn_mat, self.state_mat, self.probe_mat


class ModelVisualizer():
    
    def __init__(self):
        pass
    
    
    def visualize(self, model, ax=None, show_speed=True, **kwargs):
        
        def movmean(a, window=3) :
            csum = np.cumsum(a, dtype=float)
            csum[window:] = csum[window:] - csum[:-window]
            return csum[window - 1:] / window

        if ax == None:
            fig, ax = plt.subplots(1,1,dpi=300, figsize=(5,6))
            
    
        rects = [mpatches.Rectangle([1.4,-.6], 18,.6, fc=None, ec='k', lw=2, fill=None),
                mpatches.Rectangle([1.4,1.2], 18,.6, fc=None, ec='k', lw=2, fill=None),
                mpatches.Rectangle([1.4,2.8], 18,.6,fc=None, ec='k', lw=2, fill=None)]
        ax.add_patch(rects[0])
        ax.add_patch(rects[1])
        ax.add_patch(rects[2])
        ax.set_xlim([0,20]); ax.set_ylim([-7,9])
        ax.text(.1,2.8+.3, '   0f', fontsize=8)
        ax.text(.1,1.2+.3, '+1f', fontsize=8)
        ax.text(.1,-.6+.3, ' -1f', fontsize=8)
        ax.text(.1, 4.2, 'Lattice')
        ax.text(.1, 8, 'States')
        ax.text(10, 8, 'Resources')
        #ax.text(.1, -2, 'Reactions')
                                
        plot_length = 18
        
        def convert_len(x, plot_length):
            return 1.4+x/self._length*plot_length
        
        # plot probe locations on the lattice
        
        # get a default or provided color list based on number of probes
        if self._n_colors> 4:
            try:
                colors = cm.get_cmap(kwargs['probe_cmap'])(np.linspace(0,1,self._n_colors))
            except:
                colors = cm.viridis(np.linspace(0,1,self._n_colors))
        else:
            try:
                colors = kwargs['colors']
            except:
                colors = self.__default_colors
            
        # for each color plot a line for each probe
        for i in range(1,self._n_colors+1):
            probe_locations = np.where(self._probe_mat == i)
            for j in range(len(probe_locations[0])):
                x = [convert_len(probe_locations[1][j],plot_length)]*2
                y = [[2.2, 2.6], [.6, 1.0], [-.8, -1.2], ][probe_locations[0][j]]
                ax.plot(x, [n+.7 for n in y], lw=2, color=colors[i-1] )
            
        # plot the elongation speeds per lattice
        if show_speed:
            try:
                ecmap = kwargs['elongation_cmap']
            except:
                ecmap = cm.coolwarm
            try:
                ebins = kwargs['elongation_bins']
            except:
                ebins = 20
            
            for i in range(self._n_frames):
               carray = movmean(self._kelong_mat[i,:], window=ebins)
               if i == 0:
                  # inset = fig.add_axes([.177, .486, .703, .1],  zorder=4)
                  ax.imshow(np.atleast_2d(carray), cmap=ecmap, aspect=1,
                             extent=(1.4,19.4,2.8,3.4), zorder=0, alpha=.5)
               if i == 1:
                  # inset = fig.add_axes([.177, .4, .703, .1], zorder=4) 
                  ax.imshow(np.atleast_2d(carray), cmap=ecmap, aspect=1,
                             extent=(1.4,19.4,1.2,1.8), zorder=0, alpha=.5)
               if i == 2:
                   ax.imshow(np.atleast_2d(carray), cmap=ecmap, aspect=1,
                             extent=(1.4,19.4,-.6,0), zorder=0, alpha=.5)
                #   inset = fig.add_axes([.177, .4-.098, .703, .1], zorder=4)
               #inset.axis('off')
              # inset.imshow(np.atleast_2d(carray), cmap=ecmap, aspect=10, zorder=4)
    

        # plot the enters, leaves
        for i in range(self._n_rxns):
            if self._rxn_mat[i,0] == 2 and self._rxn_mat[i,4] in [-1,1]:
                fr = [2.8, 1.2, -.6][int(self._rxn_mat[i,2])]+.6
                loc = convert_len(self._rxn_mat[i,3],plot_length)
                if self._rxn_mat[i,4] == 1:
                    arrow = mpatches.FancyArrowPatch((loc, fr+.2), (loc, fr),
                                  mutation_scale=15,shrinkA=10, shrinkB=10, fc='#00FF00')
                else:
                    arrow = mpatches.FancyArrowPatch((loc, fr-.6), (loc, fr+.7),
                                  mutation_scale=15, fc='#FF0000')
    
                ax.add_patch(arrow)
                
        ax.imshow(np.atleast_2d(np.linspace(0,1,50)), cmap=ecmap, aspect=1,
           extent=(.3,5,-4,-4.4), zorder=0, alpha=.5)
        ax.text(.2,-3.8,'Speed', fontsize=7)
        arrow = mpatches.FancyArrowPatch((1, -3.1+.2), (1, -3.1),
                      mutation_scale=10,shrinkA=10, shrinkB=10, fc='#00FF00')
        ax.add_patch(arrow)
        arrow = mpatches.FancyArrowPatch((3.4, -3.1+.2), (3.4, -3.1),
                                  mutation_scale=10, fc='#FF0000')
        ax.add_patch(arrow)
        ax.text(.3,-2.5,'Enters',fontsize=7)
        ax.text(3,-2.5,'Leaves',fontsize=7)
        ax.text (.3, -5.3, 'Tags', fontsize=7)
        ax.add_patch(mpatches.Rectangle([.1,-6.9],6.1,4.9,fc='None', lw=1, ec='k', fill=None ))
        
        jumpcolor = '#3EB879'
        def arrowup(ratio, start, m1, stop, cap, width, capsize):
            #w, h = fig.get_figwidth(), fig.get_figheight()
            lx = start[0] - stop[0]
            xwidth= width/ratio
            ywidth = width/ratio
            path = [[start[0]-xwidth, start[1]],
                    [start[0]+xwidth, start[1]],
                    [start[0]+xwidth, m1[1]-ywidth],
                    [stop[0]+xwidth,  m1[1]-ywidth],
                    [stop[0]+xwidth, stop[1]],
                    [stop[0]+xwidth+capsize, stop[1]],
                    [stop[0], cap[1]],
                    [stop[0]-xwidth-capsize, stop[1]],
                    [stop[0]-xwidth, stop[1]],
                    [stop[0]-xwidth, m1[1]+ywidth],
                    [start[0]-xwidth, m1[1]+ywidth],
                    [start[0]-xwidth, start[1]]]
            
            path = Path(path)
            patch = PathPatch(path, ec='k', lw=1, zorder=1, fill=True, facecolor=jumpcolor)
            return patch

        def arrowdown(ratio, start, m1, stop, cap, width, capsize):
            
            lx = start[0] - stop[0]
            xwidth= width/ratio
            ywidth = width/ratio
            path = [[start[0]-xwidth, start[1]-.07],
                    [start[0]+xwidth, start[1]-.07],
                    [start[0]+xwidth, m1[1]+ywidth],
                    [stop[0]+xwidth,  m1[1]+ywidth],
                    [stop[0]+xwidth, stop[1]],
                    [stop[0]+xwidth+capsize, stop[1]],
                    [stop[0], cap[1]],
                    [stop[0]-xwidth-capsize, stop[1]],
                    [stop[0]-xwidth, stop[1]],
                    [stop[0]-xwidth, m1[1]-ywidth],
                    [start[0]-xwidth, m1[1]-ywidth],
                    [start[0]-xwidth, start[1]-.05]]
            
            path = Path(path)
            patch = PathPatch(path, facecolor=jumpcolor, ec='k', lw=1, zorder=1,  fill=True,)
            return patch

        def arrowleft(ratio, start, m1, stop, cap, width, capsize):
            #w, h = fig.get_figwidth(), fig.get_figheight()
            lx = start[0] - stop[0]
            xwidth= width/ratio
            ywidth = width/ratio
            path = [[start[0]-xwidth, start[1]],
                    [start[0]+xwidth, start[1]],
                    [start[0]+xwidth, m1[1]-ywidth],
                    [stop[0]-xwidth,  m1[1]-ywidth],
                    [stop[0]-xwidth, stop[1]],
                    [stop[0]-xwidth-capsize, stop[1]],
                    [stop[0], cap[1]],
                    [stop[0]+xwidth+capsize, stop[1]],
                    [stop[0]+xwidth, stop[1]],
                    [stop[0]+xwidth, m1[1]+ywidth],
                    [start[0]-xwidth, m1[1]+ywidth],
                    [start[0]-xwidth, start[1]]]
            
            path = Path(path)
            patch = PathPatch(path, ec='k', lw=1, zorder=1, fill=True, facecolor=jumpcolor)
            return patch
        

        def arrowright(ratio, start, m1, stop, cap, width, capsize):
            #w, h = fig.get_figwidth(), fig.get_figheight()
            lx = start[0] - stop[0]
            xwidth= width/ratio
            ywidth = width/ratio
            path = [[start[0]+xwidth, start[1]],
                    [start[0]-xwidth, start[1]],
                    [start[0]-xwidth, m1[1]-ywidth],
                    [stop[0]+xwidth,  m1[1]-ywidth],
                    [stop[0]+xwidth, stop[1]],
                    [stop[0]+xwidth+capsize, stop[1]],
                    [stop[0], cap[1]],
                    [stop[0]-xwidth-capsize, stop[1]],
                    [stop[0]-xwidth, stop[1]],
                    [stop[0]-xwidth, m1[1]+ywidth],
                    [start[0]+xwidth, m1[1]+ywidth],
                    [start[0]+xwidth, start[1]]]
            
            path = Path(path)
            patch = PathPatch(path, ec='k', lw=1, zorder=1, fill=True, facecolor=jumpcolor)
            return patch


        def arrowlefttop(ratio, start, m1, stop, cap, width, capsize):
            #w, h = fig.get_figwidth(), fig.get_figheight()
            lx = start[0] - stop[0]
            xwidth= width/ratio
            ywidth = width/ratio
            path = [[start[0]-xwidth, start[1]],
                    [start[0]+xwidth, start[1]],
                    [start[0]+xwidth, m1[1]+ywidth],
                    [stop[0]-xwidth,  m1[1]+ywidth],
                    [stop[0]-xwidth, stop[1]],
                    [stop[0]-xwidth-capsize, stop[1]],
                    [stop[0], cap[1]],
                    [stop[0]+xwidth+capsize, stop[1]],
                    [stop[0]+xwidth, stop[1]],
                    [stop[0]+xwidth, m1[1]-ywidth],
                    [start[0]-xwidth, m1[1]-ywidth],
                    [start[0]-xwidth, start[1]]]
            
            path = Path(path)
            patch = PathPatch(path, ec='k', lw=1, zorder=1, fill=True, facecolor=jumpcolor)
            return patch

        def arrowrighttop(ratio, start, m1, stop, cap, width, capsize):
            #w, h = fig.get_figwidth(), fig.get_figheight()
            lx = start[0] - stop[0]
            xwidth= width/ratio
            ywidth = width/ratio
            path = [[start[0]+xwidth, start[1]],
                    [start[0]-xwidth, start[1]],
                    [start[0]-xwidth, m1[1]+ywidth],
                    [stop[0]+xwidth,  m1[1]+ywidth],
                    [stop[0]+xwidth, stop[1]],
                    [stop[0]+xwidth+capsize, stop[1]],
                    [stop[0], cap[1]],
                    [stop[0]-xwidth-capsize, stop[1]],
                    [stop[0]-xwidth, stop[1]],
                    [stop[0]-xwidth, m1[1]-ywidth],
                    [start[0]+xwidth, m1[1]-ywidth],
                    [start[0]+xwidth, start[1]]]
            
            path = Path(path)
            patch = PathPatch(path, ec='k', lw=1, zorder=1, fill=True, facecolor=jumpcolor)
            return patch
        
        def check_mid(midpoints, midpt, loc1, loc2, updown):
            if len(midpoints) == 0:
                return midpt
            mdpt_array = np.array(midpoints)
            pl1 = 1
            pl2 = 1
            redo = True
            midpt_original = midpt
            if midpt in mdpt_array[:,2]:
                while redo:

                    matches = np.where(mdpt_array[:,2] == midpt)[0]
                    redo=False

                    for match in matches:
                        if loc1 > mdpt_array[match,0] and loc1 > mdpt_array[match,1]:
                            pl1 = 2
                        if loc1 <= mdpt_array[match,0] and loc1 <= mdpt_array[match,1]:
                            pl1 = 0
                        if loc2 > mdpt_array[match,0] and loc2 > mdpt_array[match,1]:
                            pl2 = 2
                        if loc2 <= mdpt_array[match,0] and loc2 <= mdpt_array[match,1]:
                            pl2 = 0
 
                        if pl1 == 1 or pl2 == 1:
                            midpt = midpt+.15*updown
                            redo=True
                        if [pl1,pl2] in [[0,2], [2,0]]:
                            midpt = midpt+.15*updown
                            redo=True

    
                        if [pl1,pl2] in [[2,2], [0,0]]:
                            pass
                            print('break')

                        

            else:
                pass
            
            return midpt


        # plot the lattice jumps 
        # 0-0, 1-1, 2-2, 0-1, 1-2, 2-0, 1-0, 2-1, 0-2
        # fr1, loc1, fr2, loc2, midpoints,
        connection_mat = []
        used_mids = []
        ratio = 20/14
        for i in range(self._n_rxns):
            
            if self._rxn_mat[i,0] == 2 and abs(self._rxn_mat[i,5])+abs(self._rxn_mat[i,6]) != 0:
                print('lattice jump')
                # starting point 
                fr = [2.8, 1.2, -.6][int(self._rxn_mat[i,2])]+.3
                loc = convert_len(self._rxn_mat[i,3],plot_length)
                
                # stop point
                fr2 = [2.8, 1.2, -.6][int(self._rxn_mat[i,2] + self._rxn_mat[i,5])]+.3
                loc2 = convert_len(self._rxn_mat[i,3] + self._rxn_mat[i,6],plot_length)
                print(fr,fr2)
                if fr == fr2:
                    if int(self._rxn_mat[i,2]) != 0:
                        if loc > loc2:
                            fr2 = fr2-.6
                            fr = fr-.3
                            m1 = fr-.6
                            cap = (loc2, fr2+.3)                  
                            m1 = check_mid(used_mids, m1, loc, loc2, -1)
                            ax.add_patch(arrowleft(ratio, (loc,fr), (m1,m1), (loc2,fr2), cap, .15, .1), )
                            used_mids.append([loc,loc2,m1,])
                        else:
                            fr2 = fr2-.6
                            fr = fr-.3
                            m1 = fr-.6
                            cap = (loc2, fr2+.3)                        
                            m1 = check_mid(used_mids, m1, loc, loc2, -1)
                            ax.add_patch(arrowright(ratio, (loc,fr), (m1,m1), (loc2,fr2), cap, .15, .1), )
                            used_mids.append([loc,loc2,m1,])
                    else:
                        
                        if loc > loc2:
                            fr2 = fr2+.6
                            fr = fr+.3
                            m1 = fr+.6
                            cap = (loc2, fr2-.3)                        
                            m1 = check_mid(used_mids, m1, loc, loc2, 1)
                            ax.add_patch(arrowlefttop(ratio, (loc,fr), (m1,m1), (loc2,fr2), cap, .15, .1), )
                            used_mids.append([loc,loc2,m1,])
                        else:
                            fr2 = fr2+.6
                            fr = fr+.3
                            m1 = fr+.6
                            cap = (loc2, fr2-.3)                        
                            m1 = check_mid(used_mids, m1, loc, loc2, 1)
                            ax.add_patch(arrowrighttop(ratio, (loc,fr), (m1,m1), (loc2,fr2), cap, .15, .1), )
                            used_mids.append([loc,loc2,m1,])
                else:
                    if fr > fr2:
                        
                        fr2 = fr2 + .6
                        fr = fr -.3
                        m1 = fr-.45
                        cap = (loc2, fr2-.2)
                        m1 = check_mid(used_mids, m1, loc, loc2, -1)
                        ax.add_patch(arrowdown(ratio, (loc,fr), (m1,m1), (loc2,fr2), cap, .15, .1), )
                        used_mids.append([loc,loc2,m1,])
                    else:
                        fr2 = fr2-.6
                        fr = fr + .3
                        m1 = fr+.5
                        cap = (loc2, fr2+.3)                        
                        m1 = check_mid(used_mids, m1, loc, loc2, 1)
                        ax.add_patch(arrowup(ratio, (loc,fr), (m1,m1), (loc2,fr2), cap, .15, .1), )
                        used_mids.append([loc,loc2,m1,])
                        
                
                
        return fig, ax
    
        
    
    


