# -*- coding: utf-8 -*-
"""
Created on Thu May 20 11:34:48 2021

@author: William Raymond
"""

import warnings
import numpy as np
from . import SequenceManipMethods
smm = SequenceManipMethods.SequenceManipMethods


class DiffusionRateCalc():

    '''
    This is a dictionary object that handles everything
     dealing with tag sequences and AA tables
    '''

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
                        '7mG': 363.22 + (12.017 + 3.021)*2 + 3*(94.9714),  #7 methyl G cap
                        'mean 80S': 2e6,
                        'human 80S': 4.3e6,  #https://pubmed.ncbi.nlm.nih.gov/25901680/
                        'tRNA_base': 25000, #approximation of 76 - 90 nt
                        'GFP': 26870,
                        'cy3':627.7,
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



    def calculate_rna_strand_base_mw(self, seq, n_loops=0,
                                     fluorophore='GFP', coat_protein='MS2CP'):
        '''
        Return the base molecular weight (no ribosomes) of a strand of
        RNA with Fluorophore FAB loops

        Parameters
        ----------
        seq : str
            sequence string of the rna construct.
        n_loops : int, optional
            number of fluorophore tagging loops in the 3'UTR. The default is 0.
        fluorophore : str, optional
            type of fluorophore GFP or mCHerry. The default is 'GFP'.
        coat_protein : TYPE, optional
            type of coat protein for the MS2 or PP7 loops.
            The default is 'MS2CP'.
        Returns
        -------
        rna_weight : float
            Molecular weight in daltons of the construct.

        '''
        if isinstance(fluorophore, str):
            fluorphore_weight = self.mw_table[fluorophore]
        else:
            fluorphore_weight = fluorophore

        if isinstance(coat_protein, str):
            coat_weight = self.mw_table[coat_protein]
        else:
            coat_weight = coat_protein

        seq = seq.lower().replace('t','u')
        rna_weight = np.sum([self.mw_table[x.lower()] for x in seq])
        rna_weight += self.mw_table['7mG']  #add weight of the cap
        rna_weight += fluorphore_weight*n_loops #add weight of MS2 or PP7 loops
        rna_weight += coat_weight*n_loops  #add weight MS2CP or PCP

        return rna_weight



    def calculate_diffusion_constant(self, molecular_weight,
                                          viscosity=4.4e-2,
                                          tempC=37,
                                          r_fun='Default'):
        '''


        Parameters
        ----------
        molecular_weight : float or ndarray of floats
            Molecular weight to convert to diffusion constant (over time or
            singular value).
        viscosity : float, optional
            viscosity in pascal seconds. The default is 4.4e-2 Pascal seconds.
        tempC : float, optional
            temperature in Celcius. The default is 37.
        r_fun : lambda, optional
            function to convert MW to hydraulic radius.
            The default is Rh = MW^(1/3).

        Returns
        -------
        diffusion_coeff : float or ndarray of floats
            Diffusion coeff of the mrna strand over time.

        '''


        #reasonable range of diffusion: .6- 0.04 um^2/s

        #https://pubs.acs.org/doi/10.1021/nl2008218
        #https://pubs.acs.org/doi/10.1021/acs.jpclett.0c01748
        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1831683/
        if r_fun == 'Default':
             ## approximation for now  #meters
            radius = (molecular_weight**.33)*1e-9
        else:
            radius = r_fun(molecular_weight)

        kbT = (tempC + 273.15)*1.380649e-23 # joules

        diffusion_coeff = kbT/ (6*np.pi * viscosity*radius)   #meters^2/s
        #Dr = kbT/ (8*np.pi * viscosity * radius**3)
        diffusion_coeff = diffusion_coeff*1e12 #convert to um^2/s
        return diffusion_coeff

    def calculate_single_rib_mw(self, aa_seq, probe_loc,
                                ribosome='human 80S', fluorophore='GFP',
                                fab='Fab_blank'):
        '''
        calculate the change in molecular weight for a single ribosome as it
        moves along codon position.

        Parameters
        ----------
        aa_seq : str
            Amino Acid sequence of the protein being translated.
        probe_loc : 2d numpy array
            probe locations, binary 0,1 for locations of probes.
        ribosome : str or float, optional
            key for the .mw_table or a molecular weight value in daltons.
            The default is 'human 80S'.
        fluorophore : str or float, optional
            key for the .mw_table or a molecular weight value in daltons.
            The default is 'GFP'.
        fab : TYPE, optional
            key for the .mw_table or a molecular weight value in daltons.
            The default is 'Fab_blank'.

        Returns
        -------
        weight_vec : 1d Numpy array
            molecular weight over position for a ribosome.

        '''

        if isinstance(ribosome, str):
            base_rib = self.mw_table[ribosome]
        else:
            base_rib = ribosome

        if isinstance(fab, str):
            fab_weight = self.mw_table['Fab_blank']
        else:
            fab_weight = fab

        if isinstance(fluorophore, str):
            fluorophore_weight = self.mw_table[fluorophore]
            fluorophore_weight = [fluorophore_weight]*probe_loc.shape[0]
        else:
            fluorophore_weight = fluorophore
            if not isinstance(fluorophore_weight, list):
                fluorophore_weight = [fluorophore_weight]

            if len(fluorophore_weight) >= probe_loc.shape[0]:
                fluorophore_weight = fluorophore_weight[:probe_loc.shape[0]]
            else:
                msg = 'warning: Provided probe location vector implies there'\
                      ' is more than one color, however only one fluorophore'\
                      ' was provided. Defaulting to repeating the one '\
                      'fluorophore MW multiple times. '
                warnings.warn(msg)
                fluorophore_weight = [fluorophore_weight[0]]*probe_loc.shape[0]

        probe_mw = probe_loc * np.array([np.sum(fluorophore_weight)
                                         + fab_weight])


        aa_mw = np.array([self.aa_table_mw[x] for x in aa_seq])
        weight_vec = base_rib + np.cumsum( probe_mw + aa_mw )
        return weight_vec

    def mw_over_time(self, ribosome_position_tensor, nt_seq, probe_loc,
                     base_mw,
                     ribosome='human 80S',
                     fluorophore='GFP',
                     fab='Fab_blank'):
        '''
        Calculate molecular weight over time for a ribosome position trajectory

        Parameters
        ----------
        X : 3d numpy array
            ribosome position tensor of size:
                [ntraj, n_time_points, n_ribs (at a time)].
        nt_seq : TYPE
            nucleotide sequence of the mRNA species.
        probe_loc : 2d Ndarray of probe locations
            binary vector describing where probe locations are.
        base_mw : float
            base molecular weight of the tagged RNA with no ribosomes actively
            translating.
        ribosome : str or float, optional
            key for the .mw_table or a molecular weight value in daltons.
            The default is 'human 80S'.
        fluorophore : str or float, optional
            key for the .mw_table or a molecular weight value in daltons.
            The default is 'GFP'.
        fab : TYPE, optional
            key for the .mw_table or a molecular weight value in daltons.
            The default is 'Fab_blank'.


        Returns
        -------
        1D numpy array
            molecular weight over time for a translating mRNA.

        '''

        if isinstance(fluorophore, str):
            fluorophore = self.mw_table[fluorophore]

        if isinstance(ribosome, str):
            ribosome = self.mw_table[ribosome]

        if isinstance(fab, str):
            fab = self.mw_table['Fab_blank']

        aa_seq = smm().nt2aa(nt_seq)
        weight_vec = self.calculate_single_rib_mw(aa_seq, probe_loc,
                                                  ribosome=ribosome,
                                                  fluorophore=fluorophore,
                                                  fab=fab)

        mw_per_pos = np.zeros(ribosome_position_tensor.shape[:-1])

        for i in range(ribosome_position_tensor.shape[0]):
            mw_per_pos[i,:] = np.sum(weight_vec[
                ribosome_position_tensor[i,:,:]], axis=1)

        return mw_per_pos + base_mw
