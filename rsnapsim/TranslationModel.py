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

class GenericMetaData():
    '''
    Class that generates some generic metadata and returns it in dictionary

    meta data currently returned:
        * user
        * id
        * datetime this function was called
        * python version
        * platform this is being run on

    '''

    def __init__(self):
        self.id = ''
        self.created_at = time.strftime('%Y-%m-%d %H:%M:%S',
                                        time.localtime(time.time()))
        self.user = os.path.expanduser("~")
        self.platform = platform.platform()
        self.python_version = sys.version

    def get(self):
        '''
        generate and return a metadata dictionary for a solver object

        Returns
        -------
        dict
            a dictionary of metadata such as solution id, time ran, user,
            platform and rss version.

        '''
        return self.__dict__


class CustomSSASoln:
    def __init__(self, mRNA_model, rib_array, state_array, resource_array, t, burnin, n_traj, solve_time):
        self.ribosome_array = rib_array
        self.state_array = state_array
        self.resource_array = resource_array
        self.t = t
        self.burnin = burnin
       # self.__meta = GenericMetaData.GenericMetaData().get()
        self.L = mRNA_model._length
        self.kelong_mat = mRNA_model._kelong_mat
        self.model_id = mRNA_model.ID
        self.n_traj = n_traj
        self.n_colors = mRNA_model._n_colors
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


class TranslationModel:
    def __init__(self, mRNA, name, particle_size=9):
        self.name = name
        self.blank(mRNA, particle_size)
        
        
    def blank(self, mRNA, particle_size):
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
                                # rxn type, exclusion?, frame, loc, dexist, dframe, dloc, 
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
        

        
        codons  = [''.join(i) for i in product(['A','U','G','C'], repeat = 3)]
        codon_ids = [[ codons.index(line.upper().replace('T','U')[i:i+3]) for i in range(0, len(line), 3)] for line in mRNA.multiframe_nt_seq]
        self._codon_mat = np.zeros([3, self._length], dtype=int)
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
    def add_ribosome_reaction(self, propensity, parameters, rxn_name='', exclusion=0, dexist = 0, dframe = 0, dloc = 0,
                              dprobes=[], dprobe_inds=[],
                              dstates=[], dstate_inds=[],
                              dresources=[], dresource_inds=[]):
        self.__make_rib_or_lattice_rxn(0, exclusion, 0,0, dexist, dframe, dloc, dprobes, dprobe_inds, dstates, dstate_inds, dresources, dresource_inds)
        self.__add_propensity(propensity, parameters, rxn_name=rxn_name)
        
        
    def add_lattice_reaction(self, propensity, parameters, rxn_name='', frame=0, loc=0, exclusion=0, dexist = 0, dframe = 0, dloc = 0,
                              dprobes=[], dprobe_inds=[],
                              dstates=[], dstate_inds=[],
                              dresources=[], dresource_inds=[]):
        self.__make_rib_or_lattice_rxn(2, exclusion,frame,loc, dexist, dframe, dloc, dprobes, dprobe_inds, dstates, dstate_inds, dresources, dresource_inds)
        self.__add_propensity(propensity, parameters, rxn_name=rxn_name)
        

    def add_state_reaction(self, propensity, parameters, rxn_name='', inds=[], dstates=[]):
                # if the user did not provide indexes, pad with zeros for right shape
        if len(inds) < len(dstates):
            change = dstates + [0,]*(len(dstates)-len(self._n_states))
            self.__add_rxn(1, [0, 0, 0, 0, 0, 0,] +   [0]*self._n_colors +  change  )
            self.__add_propensity(propensity, parameters, rxn_name=rxn_name)
            return
            
        # if the user did  provide indexes, use them to build the reaction row
        change = [0]*self._n_states
        for i in range(self._n_states):
            if i in inds:
                change[i] = dstates[inds.index(i)]
        self.__add_rxn(1, [0, 0, 0, 0, 0, 0,] +  [0]*self._n_colors +  change   )
        self.__add_propensity(propensity, parameters, rxn_name=rxn_name)


    def _x0(self):
        return self._rib_arr0, self._lattice_arr0, self._state_arr0, self._resource_arr0

    def __make_rib_or_lattice_rxn(self, rtype, exclusion, frame, loc, dexist, dframe, dloc, dprobes, dprobe_inds, dstates, dstate_inds, dresources, dresource_inds):
        change1 = [exclusion, frame, loc, dexist, dframe, dloc,]
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


    def add_resource_reaction(self, propensity, parameters, rxn_name='', dresources=[], inds=[] ):
        
        # if the user did not provide indexes, pad with zeros for right shape
        if len(inds) < len(dresources):
            change = dresources + [0,]*(len(dresources)-len(self._n_resources))
            self.__add_rxn(3, [0, 0, 0, 0, 0, 0,] + [0]*self._n_colors+ [0]*self._n_states +  change )
            self.__add_propensity(propensity, parameters, rxn_name=rxn_name)
            return
            
        # if the user did  provide indexes, use them to build the reaction row
        change = [0]*self._n_resources
        for i in range(self._n_resources):
            if i in inds:
                change[i] = dresources[inds.index(i)]
        
        self.__add_rxn(3, [0, 0, 0, 0, 0, 0,] + [0]*self._n_colors+ [0]*self._n_states +  change   )
        self.__add_propensity(propensity, parameters, rxn_name=rxn_name)

    def add_probe_reaction(self, propensity, parameters, rxn_name='', ind=0, dprobe=0):
        self.__add_rxn(4, [ind, dprobe],)
        self.__add_propensity(propensity, parameters, rxn_name=rxn_name)
    
    def add_probe_reaction(self, propensity, parameters, rxn_name='', inds=[], dprobes=[]):
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
    
    def add_global_reaction(self, propensity, parameters, rxn_name=''):
        self.__add_rxn(5, [],)
        self.__add_propensity(propensity, parameters, rxn_name=rxn_name)
        
    def delete_reaction(self, row):
        self.__delete_rxn(row)
        
    def __add_propensity(self, propensity_fun, parameters, rxn_id=-1, rxn_name=''):
        if rxn_id == -1:
            try:
                rxn_id = max(self._prop_ids)
            except:
                rxn_id = 0
        self._propensities = self._propensities + [propensity_fun, ]
        self._parameters = self._parameters + [parameters, ]
        self._prop_ids = self._prop_ids + [rxn_id, ]
        self._rxn_names = self._rxn_names + [rxn_name]
        
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

    def load_model_c(self, model_name):
        cmodel = importlib.import_module('rsnapsim.models.%s.%s'%(self.name,self.name))
        if self.model_id != cmodel.__model_id:
            print('error')
        

    def compile_model_c(self):
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
        
        mm.ModelFactory().compile_model(self.name, self.model_id, [self._propensities[x] for x in self._constant_reactions],
                                     [self._propensities[x] for x in self._ribosome_reactions], overwrite=True )        
        
        self.cmodel = importlib.import_module('rsnapsim.models.%s.%s'%(self.name,self.name))


    def visualize(self, ax=None, show_speed=True, **kwargs):
        
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
    
    
    def __get_mats(self):
        return kelong_mat, resource_mat, rxn_mat, state_mat, probe_mat


def arrow_3( start, stop, m1, m2, point):
    r1 = mpatches.Rectangle(start,1,m1[1]-start[1],ec='k',lw=2)
    r2 = mpatches.Rectangle(m1,m2[0]-m1[0],1,ec='k',lw=1)
    r3 = mpatches.Rectangle(m2,1,stop[1]-m2[1],ec='k',lw=2)
    
    
    return [r1,r2,r3]

from matplotlib.patches import PathPatch
from matplotlib.path import Path
def arrow3(fig, start, m1, stop, cap, width, capsize):
    w, h = fig.get_figwidth(), fig.get_figheight()
    lx = start[0] - stop[0]
    xwidth= width*h/w
    ywidth = width*1
    path = [[start[0]-xwidth, start[1]-ywidth],
            [start[0]+xwidth, start[1]-ywidth],
            [start[0]+xwidth, m1[1]-ywidth],
            [stop[0]+xwidth,  m1[1]-ywidth],
            [stop[0]+xwidth, stop[1]],
            [stop[0]+xwidth+capsize, stop[1]],
            [stop[0], cap[1]],
            [stop[0]-xwidth-capsize, stop[1]],
            [stop[0]-xwidth, stop[1]],
            [stop[0]-xwidth, m1[1]+ywidth],
            [start[0]-xwidth, m1[1]+ywidth],
            [start[0]-xwidth, start[1]-ywidth]]
            
    path = Path(path)
    patch = PathPatch(path, facecolor='none', ec='k', lw=2)
    return patch









class SSASolver:
    def __init__(self, NUMBER_OF_CORES=None):
        if NUMBER_OF_CORES == None:
            self.NUMBER_OF_CORES = int(multiprocessing.cpu_count()/2) #claim half the cores
        else:
            self.NUMBER_OF_CORES = 1
        

    def solve_ssa(self, mRNA_model, t, n_traj=1, burnin=0, seed=None, parallel=False):
        st = time.time()
        if seed == None:
            seeds = np.random.randint(0x7FFFFF, size=n_traj)
        
        if parallel:
            solns = Parallel(n_jobs = self.NUMBER_OF_CORES, prefer="threads")(delayed(self.__run)(*args, **kwargs) for args, kwargs in ([[(model, t), {'burnin':burnin, 'seed':seeds[i]}] for i in range(0,n_traj)]))

        else:
            solns = []
            for i in range(n_traj):
                print(seeds[i])
                solns.append(self.__run(mRNA_model, t, burnin=burnin, seed=seeds[i]))
        
        rib_array = np.array([solns[i][0] for i in range(len(solns))])
        resource_array = np.array([solns[i][2] for i in range(len(solns))])
        state_array = np.array([solns[i][1] for i in range(len(solns))])
        solve_time = time.time() - st
        soln = CustomSSASoln(mRNA_model, rib_array, state_array, resource_array, t, burnin, n_traj, solve_time) 
        
        return soln
        
        
    def __initalize_trajs(self, mRNA_model):
        rib_arr, lattice_arr, state_arr, resource_arr = mRNA_model._x0()
        NR = np.sum(lattice_arr)
        occupied = rib_arr[:,3]
        return np.copy(rib_arr), np.copy(lattice_arr), np.copy(state_arr), np.copy(resource_arr), np.copy(occupied), NR
        
    def __constants(self, mRNA_model):
        particle_size = mRNA_model.particle_size
        rxn_mat = mRNA_model._rxn_mat.astype(int)
        n_colors = int(np.max(mRNA_model._probe_mat))
        n_states = mRNA_model._n_states
        n_resources = mRNA_model._n_resources
        n_rxns = int(rxn_mat.shape[0])
        n_ribosome_reactions = np.sum(rxn_mat[:,0] == 0)
        max_rib = mRNA_model._max_particles
        n_constant_reactions = len(rxn_mat) - n_ribosome_reactions
        L = kelong_mat.shape[1]
        
        return particle_size, max_rib, rxn_mat, n_colors, n_states, n_rxns, n_resources, n_constant_reactions, n_ribosome_reactions, L
        


        
    def __run(self, mRNA_model, t, burnin=0, seed=None):
        
        # initalize constants, flags, and initial state of the simulation
        footprint, max_rib, rxn_mat, n_colors, n_states, n_rxns, n_resources, n_constant_reactions, n_ribosome_reactions, L = self.__constants(mRNA_model)
        reaction_taken, dexist, rib_id, ribosome_moved, tindex = [0,]*5
        rib_arr, lattice_arr, state_arr, resource_arr, occupied, NR = self.__initalize_trajs(mRNA_model)
        kelong_mat = mRNA_model._kelong_mat
        probe_mat = mRNA_model._probe_mat.astype(int)
        #rint(state_arr)
        if seed != None:
            np.random.seed(seed)
            
        # initalize propensities
        constant_props = [mRNA_model._propensities[i] for i in mRNA_model._constant_reactions]
        constant_parameters = [mRNA_model._parameters[i] for i in mRNA_model._constant_reactions]

        def get_constant_props(*props):
            return [props[i](constant_parameters[i],tc, rib_arr, kelong_mat, occupied, lattice_arr, probe_mat, state_arr, resource_arr, NR) for i in range(len(props))]
            
        ribosome_props = [mRNA_model._propensities[i] for i in mRNA_model._ribosome_reactions]
        ribosome_parameters = [mRNA_model._parameters[i] for i in mRNA_model._ribosome_reactions]
        
        def get_ribosome_props(*props):
            rates = [props[i](ribosome_parameters[i],tc, rib_arr, kelong_mat, occupied, lattice_arr, probe_mat, state_arr, resource_arr, NR) for i in range(len(props))]
            return [item for sublist in rates for item in sublist]
        
        probe_function = mRNA_model._probe_function
        probe_fun = inspect.getsourcelines(probe_function)[0][0].split(':')[-1].replace('\n','').replace(' ','')
        if probe_fun == '1':
            use_probe_fun = False
        else:
            use_probe_fun = True
        
        probe_parameters = mRNA_model._probe_parameters
        
        reaction_ids = mRNA_model._constant_reactions + mRNA_model._ribosome_reactions 
        
        # Initalize arrays to store the trajectory
        #intensity_array = np.zeros([len(t), n_colors], dtype=int)
        ribosome_array = np.zeros([len(t), max_rib, 4+n_colors+n_constant_reactions+n_ribosome_reactions], dtype=int)
        state_array = np.zeros([len(t), n_states], dtype=int)
        resource_array = np.zeros([len(t), n_resources], dtype=int)
        
        
        tc = 0-burnin # current time
        tf = t[-1] # final time point
        
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

            if rxn_mat[event][0] == 4: # probe reaction
                state_mat[rxn_mat[event][1]] += rxn_mat[event][2]
                state_mat[rxn_mat[event][3]] += rxn_mat[event][4]
                state_mat[rxn_mat[event][5]] += rxn_mat[event][6]


            # Check probes
            if ribosome_moved:
                # update occupied vector
                occupied = rib_arr[:,3]

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


