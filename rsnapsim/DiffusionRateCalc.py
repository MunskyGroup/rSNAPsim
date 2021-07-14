# -*- coding: utf-8 -*-
"""
Created on Thu May 20 11:34:48 2021

@author: William Raymond
"""

import numpy as np
from . import SequenceManipMethods
smm = SequenceManipMethods.SequenceManipMethods

'''
This is a dictionary object that handles everything dealing with tag sequences and AA tables
'''


class DiffusionRateCalc():
    def __init__(self):
        
        # g/mol or Daltons
        self.mw_table = {'a': 267.2413,          #adenosine
                    'AMP': 267.2413,
                    'g': 363.22,   #guanosine
                    'GMP': 363.22,
                    'u': 244.2014,    #uridine
                    'UMP': 244.2014,
                    'c': 243.217,   #cytidine
                    'CMP': 243.217,
                    '7mG': 363.22 + (12.017 + 3.021)*2 + 3*(94.9714) ,  #7 methyl G cap
                    'mean 80S': 2e6,
                    'human 80S': 4.3e6,  #https://pubmed.ncbi.nlm.nih.gov/25901680/
                    'tRNA_base': 25000, #approximation of 76 - 90 nt        
                    'GFP': 26870,
                    'mCherry': 28000,
                    'Fab_blank': 48000,#approximation
                    'MS2CP': 13700,
                    'PCP': 16025,
            
            }
        
        # amino acid individual MWs
        self.aa_table_mw = {'A': 75,
                     'R':175,
                     'N': 132,
                     'D':133,
                     'B': 133,
                     'C': 121,
                     'Q':146,
                     'E': 147,
                     'Z':147,
                     'G':75,
                     'H':155,
                     'I':131,
                     'L':131,
                     'K':146,
                     'M':149,
                     'F':165,
                     'P':115,
                     'T':119,
                     'W':204,
                     'Y':181,
                     'V':117,
                     'S':105,
                     '*':0}
        
        
        
    def calculate_rna_strand_base_mw(self,seq, n_loops = 0, fluorophore = 'GFP', coat_protein = 'MS2CP'):
        '''
        Return the base molecular weight (no ribosomes) of a strand of RNA with Fluorophore FAB loops 

        Parameters
        ----------
        seq : str
            sequence string of the rna construct.
        n_loops : int, optional
            number of fluorophore tagging loops in the 3'UTR. The default is 0.
        fluorophore : str, optional
            type of fluorophore GFP or mCHerry. The default is 'GFP'.
        coat_protein : TYPE, optional
            type of coat protein for the MS2 or PP7 loops. The default is 'MS2CP'.
        Returns
        -------
        rna_weight : float
            Molecular weight in daltons of the construct.

        '''
    
        seq = seq.lower().replace('t','u')
        rna_weight = np.sum( [self.mw_table[x.lower()] for x in seq])
        rna_weight += self.mw_table['7mG']  #add weight of the cap
        rna_weight += self.mw_table[fluorophore]*n_loops  #add weight of MS2 or PP7 loops
        rna_weight += self.mw_table[coat_protein]*n_loops  #add weight MS2CP or PCP
        
        return rna_weight
        
        
    
    def calculate_base_diffusion_constant(self,molecular_weight, viscosity= 4.4e-2 , tempC = 37, r_fun = 'Default'):
        
        #reasonable range of diffusion: .1- 0.4 um^2/s
        
        #https://pubs.acs.org/doi/10.1021/nl2008218
        #https://pubs.acs.org/doi/10.1021/acs.jpclett.0c01748
        if r_fun == 'Default':
            radius = (molecular_weight**.8)/1000*1e-9   ## approximation for now  #meters
        else:
            radius = r_fun(molecular_weight)
        
        kbT = (tempC + 273)*1.380649e-23 # joules
        
    
        Dt = kbT/ (6*np.pi * viscosity*radius)   #meters^2/s
        Dr = kbT/ (8*np.pi * viscosity * radius**3)
        
        Dt = Dt*1e12 #convert to um^2/s
    
        return Dt
    
    def aa_sum_vector(self,aa_seq, probe_loc, ribosome = 'human 80S', fluorophore = 'GFP'):
        

        base_rib = self.mw_table[ribosome]
        
        probe_mw = probe_loc * np.array([self.mw_table[fluorophore]+ self.mw_table['Fab_blank']   ])
        print(probe_mw)
        aa_mw = np.array([self.aa_table_mw[x] for x in aa_seq])
        
        
        weight_vec = base_rib + np.cumsum( probe_mw + aa_mw )
        
        
        
        return weight_vec
    
    def convert_rib_pos_tensor(self,X, nt_seq, probe_loc, base_mw, ribosome = 'human 80S', fluorophore = 'GFP', hydraulic_radii_fun = 'Default'):
        
        aa_seq = smm.SequenceManipMethods().aa_seq(nt_seq)
        weight_vec = self.aa_sum_vector(aa_seq, probe_loc)
        
        mw_per_pos = np.zeros( X.shape[:-1] )
        
        for i in range(X.shape[0]):  # for each trajectory
        
        
        
        
        
        
        
        
        
        
        
        
        
   
    
    