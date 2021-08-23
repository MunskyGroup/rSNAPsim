# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 11:21:05 2021

@author: willi
"""

import time
import numpy as np
import os
import importlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
import matplotlib.patches as mpatches

from . import rsnapsim_model_maker
ModelFactory = rsnapsim_model_maker.ModelFactory
RuleConverter = rsnapsim_model_maker.RuleConverter
from . import PropensityFactory, ProbeVectorFactory


class ModelBuilder():
    def __init__(self):
        self.poi = None
        self.n_total_reactions = 0
        self.n_state_reactions = 0
        self.max_ribs = 100
        self.parameters = None
        self.rules = None
        self.kin = 0.033
        self.kout = 10
        
    
    def get_max_ribs(self):
        self.max_rib = int(len(self.poi.kelong)/9 + 5)
        
    def build_lattice(self, N_reactions):
        forward_rates = self.poi.kelong
        self.stoich_lattice = np.zeros([N_reactions,
                                        len(forward_rates)], dtype=np.int32)
        self.xi_lattice = np.zeros([1,len(forward_rates)], dtype=np.int32)
        self.n_total_reactions += N_reactions
        
    def build_states(self, state_stoichiometry, state_inital_condition):
        self.stoich_state = state_stoichiometry.astype(np.int32)
        self.xi_state = state_inital_condition.astype(np.int32)
        
        self.n_total_reactions += int(state_stoichiometry.shape[0])
        self.n_state_reactions = int(state_stoichiometry.shape[0])
        
    
    def compile_model(self, name, overwrite=True, rules = None):
        mf = ModelFactory()
        mf.compile_model(name, overwrite=True, rules = rules)
        model = self.get_model(name)
        
    def get_model(self, name):
        model = importlib.import_module('rsnapsim.models.%s.%s'%(name,name))
        return model
    
    

            
                
                
        
               
        
        