# -*- coding: utf-8 -*--
"""
Created on Tue Oct 23 09:42:24 2018

@author: William
"""

import re #import regex

import FileParser
import warnings
import sys
import os
sys.path.append("C:\\Users\\willi\\Documents\\Github\\rSNAPsim\\generalized_cpp")
sys.path.append("C:\\Users\\willi\\Documents\\Github\\rSNAPsim\\ssa_cpp")
sys.path.append("C:\\Users\\willi\\Documents\\Github\\rSNAPsim\\trna_ssa")

path_to_cpp = ''
path_to_gen = ''
path_to_trna =''
#OS walk to find the cpp compilation
for root, dirs, files in os.walk(".", topdown=False):
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
        
        import ssa_translation
        import ssa_translation_lowmem
        os.chdir(cwd)
    except:
        os.chdir(cwd)
    

if path_to_gen != '':
    try:
        cwd = os.getcwd()
        
        os.chdir(path_to_gen)
        print('importing C++ models')
        import ssa_translation_generic
        import ssa_translation_generic_lowmem
        print('c++ models loaded successfully')
        os.chdir(cwd)
    except:
        os.chdir(cwd)

if path_to_trna !='':
    try:
        cwd = os.getcwd()
        
        os.chdir(path_to_gen)
        print('importing C++ tRNA models')
        import ssa_trna
        print('c++ models loaded successfully')
        os.chdir(cwd)
    except:
        os.chdir(cwd)   
    
    
try:
    
    from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord
except:
    pass

import time
import json, codecs

from scipy import sparse
from scipy.stats import pearsonr


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import matplotlib.patches as mpatches
import matplotlib.animation as animation
from matplotlib.collections import PatchCollection
from matplotlib import cm
from matplotlib import gridspec
from matplotlib.patches import Ellipse
#import scipy.stats.trim_mean as tmean

from scipy.stats import kde
import scipy as sci

import pandas as pd
import copy
try:
    from Bio import SeqIO
    from Bio import Entrez
except:
    
    print('BioPython is not installed, polling genbank will not be possible')
    pass


import translation_models as models

import platform
from scipy.stats import multivariate_normal,chi2
import pickle

class rSNAPsim():

    """
    The Single Molecule Simulator (SMS) provides a python class for running
    single molecule mRNA translation simulations

    When presented with a valid protein sequence the SMS can find open reading frames
    and simulate intensity trajectories from translation of the protein with given fluorescent tags.

    *model description*

        link to paper here / image


    *main functions*

        -open_seq_file(filepath), opens a txt or .gb file and gets the sequence

        -get_orfs(nt_sequence, min_codons), returns open reading frames of a given
        sequence and a minimum codon length per protein

        -get_temporal_proteins(), gets the proteins after get_orfs

        -analyze_poi(aa_seq,nt_seq), analyzes the proteins of intrest for
        codon sensitivity and elongation rates

        -__.poi(), class to contain proteins of intrest after analyzed

        -run_default(), runs get_orfs, get_temporal proteins, and analyze_poi
        with the first protien found in the sequence



    *attributes*

        **gene_sequence_str** = string of the nucleotide sequence
        **tag_dict** = dictionary with various types of fluorescent tag epitopes
        
        **tag_full** = dictionary of full tag sequences
        
        **aa_keys** = amino acid single letter keys
        
        **codon_types** = flag dictionary of which amino acids are set to Wild-type, fast, or slow
        
        **aa_table** = dictionary of amino acids
        
        **aa_table_r** = reverse dictionary (amino acid letters are the keys)
        
        **strGeneCopy** = dictionary of wild-type tRNA copy numbers
        
        **strGeneCopy_fast** = dictionary of fast tRNA copy numbers
        
        **strGeneCopy_slow** = dictionary of slow tRNA copy numbers
        
        **slow_codons_value** = list of slowest codon tRNA copy numbers
        
        **fast_codons_value** = list of fastest codon tRNA copy numbers
        
        **sensitivity_fast_slow** = list of sensitivity for amino acids
        
        **poi** = Class container for proteins of intrest
        
        **orfs** = dictionary of open reading frames with keys 1,2,3
        
        **seq_str** = sequence string
        
        **proteins** = dictionary of proteins detected in the sequence by ORF
        
        **tagged_proteins** = dictionary of proteins that were detected and tagged


    *POI*

        Protein of intrest has the following attributes:

        **aa_seq** = amino acid sequence
        
        **nt_seq** = nucleotide sequence
        
        **gene_length** = length of the gene
        
        **tag_length** = length of the tags
        
        **total_length** = total length of the full amino acid sequence
        
        **name** = name of the gene
        
        **tag_types** = what types of tags does the protien have
        
        **tag_epitopes** = type of tags and epitope lists per tag
        
        **codon_sensitivity** = how sensitive is the protein per amino acid sequence?
        
        **CAI** = codon activation index
        
        **CAI_codons** = means of the codon activation

    *ssa*

        The ssa container class has the following attributes:


        **no_ribosomes** = number of ribosomes
        
        **n_traj** = number of trajectories
        
        **k** = all kelongation rates (calculated from codon sequence)
        
        **no_rib_per_mrna** = number of ribosomes per mRNA strand on average
        
        **rib_density** = ribosome density
        
        **rib_means** = ribosome means
        
        **rib_vec** = raw ribosome location matrix for each trajectory
        
        **intensity_vec** = fluorescence intensities
        
        **time_vec_fixed** = the time vector
        
        **start_time** = the time the simulation was started

        **evaluating_inhibitor** = was there an inhibitor present?
        
        **evaluating_frap** = was the simulation subjected to a FRAP test
        
        **time_inhibit** = the time of the perturbation

        **autocorr_vec** = autocorrelation vector of intensities
        
        **mean_autocorr** = the average autocorrelations, averaged over trajectories
        
        **error_autocorr** = the standard deviation of the autocorrelation
        
        **dwell_time** = how long do the ribosomes stay on the mRNA strand calculated by the simulation
        
        **ke_sim** = the calculated average elongation rate from the simulations


    """

    def __init__(self):
        #self.gene_sequence_str = lambda: self.sequence_str()
        self.__version__ = '1.0a-dev'
        self.codon_dicts = CodonDictionaries()
    
    @property
    def gene_sequence_str(self):
        return self.sequence_str
    @gene_sequence_str.setter
    def gene_sequence_str(self,value):
        self.sequence_str = value


    def __update_sensitivity(self):

        """
        updates sensitivities for the GUI implementation call
        """

        self.fast_codons_value = []
        for key in self.aa_keys:
            values = []
            codons = self.aa_table_r[key]
            for codon in codons:
                values.append(self.strGeneCopy[codon])

            self.fast_codons_value.append(max(values))

            for codon in codons:
                self.strGeneCopy_fast[codon] = max(values)



        self.slow_codons_value = []
        for key in self.aa_keys:
            values = []
            codons = self.aa_table_r[key]
            for codon in codons:
                values.append(self.strGeneCopy_slow[codon])

            self.slow_codons_value.append(min(values))

            for codon in codons:
                self.strGeneCopy_slow[codon] = min(values)

        codonkeys = ['GCT', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGT', 'CAT', 'ATT',
                     'TTA', 'AAA', 'ATG', 'TTT', 'CCT', 'TCT', 'ACT', 'TGG', 'TAT', 'GTT', 'TAA']

        self.sensitivity_fast_slow = []
        for i in range(len(codonkeys)):
            self.sensitivity_fast_slow.append(self.strGeneCopy_fast[codonkeys[i]] / self.strGeneCopy_slow[codonkeys[i]])





    def get_k_construct(self, nt_seq, k_init, k_elong_mean, codon_types=None):
        '''
        Returns the k_elongation rates of a given nucleotide sequence under constructed conditions
        given some sort of key describing which amino acids are slow, fast or natural

        *args*

            **nt_seq**, nucleotide sequence to get the propensities of

            **k_init**, initiation rate of starting translation

            **k_elong_mean**, average rate of elongation for the protein translation


        *keyword args*

            **codon_types**, a dictonary or identifier determining which amino acids are slow, fast or natural

                self.codon_types is an example dictionary for the user to change / utilize, if codon_types is left blank
                get_k_construct uses this internal dictonary

                ex: codon_types = 'slow' or 'rare'  all amino acids set to slow
                    codon_types = 'fast' or 'common'  all amino acids set to fast
                    codon_types = 'natural' all amino acids set to fast

                    codon_types = {'A':[0], 'T':[2]}  A set to slow, T set to fast
                    codon_types = {'rare':['A','R'],'common':['L']}  A and R set to slow, L set to fast


        '''



        if codon_types == None:
            codon_types = self.codon_types
        else:
            all_natural = dict(zip(self.aa_keys, np.ones((1, 20)).flatten().astype(int).tolist()))

            if isinstance(codon_types, str):
                if codon_types == 'rare' or codon_types == 'slow':
                    all_natural = dict(zip(self.aa_keys, np.zeros((1, 20)).flatten().astype(int).tolist()))
                if codon_types == 'common' or codon_types == 'fast':
                    all_natural = dict(zip(self.aa_keys, (2*np.ones((1, 20))).flatten().astype(int).tolist()))
            if isinstance(codon_types, dict):
                for key in codon_types.keys():
                    if isinstance(key, str):
                        if key.lower() not in ['rare', 'common', 'natural']:
                            if key.upper() in self.aa_keys:
                                if codon_types[key] in [0, 1, 2]:
                                    all_natural[key] = key
                                if codon_types[key] in ['rare', 'common', 'natural']:
                                    if codon_types[key] == 'rare':
                                        all_natural[key] = 0
                                    if codon_types[key] == 'common':
                                        all_natural[key] = 2
                                    if codon_types[key] == 'natural':
                                        all_natural[key] = 1
                        else:
                            newkeys = codon_types[key]
                            for newkey in newkeys:
                                if newkey.upper() in self.aa_keys:
                                    if key.lower() == 'rare':
                                        all_natural[newkey.upper()] = 0
                                    if key.lower() == 'common':
                                        all_natural[newkey.upper()] = 2
                                    if key.lower() == 'natural':
                                        all_natural[newkey.upper()] = 1


                    if isinstance(key, int):
                        newkeys = codon_types[key]
                        for newkey in newkeys:
                            all_natural[newkey] = key




            codon_types = all_natural


        aa_seq = self.nt2aa(nt_seq)



        tRNA_design = np.zeros((1, len(aa_seq)))
        tRNA_norm = np.zeros((1, len(aa_seq)))

        seperated_codons = [nt_seq[i:i+3] for i in range(0, len(nt_seq), 3)] #split codons by 3


        for i in range(len(seperated_codons)):
            tRNA_norm[0, i] = self.strGeneCopy[seperated_codons[i]]





        for i in range(len(self.aa_keys)-1):

            fs = codon_types[self.aa_keys[i]]
            indexes = [m.start() for m in re.finditer(self.aa_keys[i], aa_seq)]
            for index in indexes:

                if fs == 0:
                    tRNA_design[0, index] = self.slow_codons_value[i]
                if fs == 2:
                    tRNA_design[0, index] = self.fast_codons_value[i]
                if fs == 1:
                    tRNA_design[0, index] = tRNA_norm[0, index]


        tRNA_design[0, -1] = tRNA_norm[0, -1]




        mean_tRNA_copynumber = np.mean(list(self.strGeneCopy_single.values()))



        k_elongation_design = (tRNA_design / mean_tRNA_copynumber) * k_elong_mean

        all_k_design = [k_init] + k_elongation_design.flatten().tolist() + [k_elong_mean]

        return all_k_design
    
        
        


    def get_temporal_proteins(self):
        '''
        gets all the temporal proteins after getting the ORFs

        __.tagged_proteins = dictionary with keys of tag types and a list of proteins
        __.pois = list of proteins of intrest
        __.pois_seq = list of nucleotide sequences of proteins of sequences
        __.proteins = dictonary with keys of 1 2 or 3 orfs

        '''


        self.proteins = {'1':[], '2':[], '3':[]}
        self.tagged_proteins = {a:[] for a in self.tag_dict.keys()}
        self.tagged_protein_seq = {a:[] for a in self.tag_dict.keys()}

        for i in range(len(self.orfs)):
            for j in range(len(self.orfs[str(i+1)])):
                pro = self.nt2aa(self.sequence_str[self.orfs[str(i+1)][j][0]:self.orfs[str(i+1)][j][1]+3])
                nt_seq = self.sequence_str[self.orfs[str(i+1)][j][0]:self.orfs[str(i+1)][j][1]+3]
                self.proteins[str(i+1)].append(pro)
                for tag in self.tag_dict.keys():
                    if self.tag_dict[tag] in pro:
                        self.tagged_protein_seq[tag].append(nt_seq)
                        self.tagged_proteins[tag].append(pro)

        tags = 0
        for key in self.tagged_proteins.keys():
            tags += len(self.tagged_proteins[key])


        self.pois = []
        self.pois_seq = []
        for tag in self.tag_dict.keys():
            for i in range(len(self.tagged_proteins[tag])):
                if self.tagged_proteins[tag][i] not in self.pois:
                    self.pois.append(self.tagged_proteins[tag][i])
                    self.pois_seq.append(self.tagged_protein_seq[tag][i])

        if len(self.pois) == 0:

            POIs = []
            pois_s = []
            pois_nt = []
            for i in range(len(self.gb_obj.features)):

                try:

                    self.gb_obj.features[i].qualifiers['translation']

                    if tags == 0:

                        POIs.append(self.gb_obj.features[i])
                        pois_s.append(self.nt2aa(self.tag_full['T_Flag']) + self.gb_obj.features[i].qualifiers['translation'][0])
                        pois_nt.append(self.tag_full['T_Flag'] + str(self.gb_obj.seq)[int(self.gb_obj.features[i].location.start):int(self.gb_obj.features[i].location.end)])
                    else:

                        POIs.append(self.gb_obj.features[i])
                        pois_s.append(self.gb_obj.features[i].qualifiers['translation'][0])
                        pois_nt.append(str(self.gb_obj.seq)[int(self.gb_obj.features[i].location.start):int(self.gb_obj.features[i].location.end)])

                except:
                    pass


            self.pois = pois_s
            self.pois_seq = pois_nt


    def seq_to_protein_obj(self, sequence_str, min_codons=80):
        
        smm = SequenceManipMethods(sequence_str)
        
        orfs = smm.get_orfs(sequence_str,min_codons=min_codons)
        protein_strs,proteins, tagged_proteins  = smm.get_proteins(orfs,sequence_str)
        
        return proteins
                

    def open_seq_file(self, seqfile,min_codons=80):
        
        '''
        Reads a sequence file, either a .txt file or a .gb genbank file

        *args*

            **seqfile**, sequence file either in txt, gb, gbk format
        '''

        fp = FileParser.FileParser()
        self.sequence_name = fp.get_name(seqfile)
        self.sequence_description = fp.get_description(seqfile)
        self.sequence_str= fp.get_sequence(seqfile).upper()
        
        smm = SequenceManipMethods(self.sequence_str)
        
        self.orfs = smm.get_orfs(self.sequence_str,min_codons=min_codons)
        protein_strs,self.proteins, self.tagged_proteins  = smm.get_proteins(self.orfs,self.sequence_str)
        
        


    def multitau_acc(self, ivec, n, sampling_rate, sample_rate_seconds):
        '''
        Multi-tau acc
        '''
        sigmas = 3
        acc = np.array([[]])
        for i in range(0, n):
            tempdata = ivec[i, :].flatten()
            tempdata[np.where(tempdata > tmean(tempdata, 10)) + sigmas*np.std(tempdata)] = 0
            tempdata[np.where(tempdata < tmean(tempdata, 10)) - sigmas*np.std(tempdata)] = 0

            if np.isnan(tempdata[0]):
                tempdata = tempdata[1:]
            if np.isnan(tempdata[-1]):
                tempdata = tempdata[:-1]

            outliers = np.where(tempdata == 0)[0]
            if outliers[-1] == len(tempdata)-1:
                outliers = outliers[:-1]
            if outliers[0] == 0:
                outliers = outliers[1:]

            tempdata[outliers] = 1/2*(tempdata[outliers-1] + tempdata[outliers+1])
            tempdata = tempdata-np.mean(tempdata)

            preacc = self.get_acc2(tempdata)
            if i == 0:
                acc = preacc
            else:
                acc = np.hstack((acc, preacc))
        for i in range(0, n):
            data = acc[i]
            data[0:sample_rate_seconds] = []

            binnedData_1 = data



    def analyze_seq_file(self, filename):
        '''
        General catch all to run all functions necessary before a SSA and store the first POI found from any given sequence

        *args*

            **filename** a txt or gb file to be read and analyzed

        '''


        self.open_seq_file(filename)
        self.get_orfs(self.sequence_str, min_codons=80)
        self.get_temporal_proteins()
        self.analyze_poi(self.pois[0], self.pois_seq[0])
        self.POI.k = self.get_k(self.POI.nt_seq, .03, 10)
        probe_vec,probe_loc = self.get_probvec()
        self.POI.probe_vec = probe_vec
        self.POI.probe_loc = probe_loc

    def run_default(self):
        '''
        this is now depreceated but I am keeping it to not break codes
        '''
        return
#        self.open_seq_file()
#        
#        self.get_orfs(self.sequence_str, min_codons=80)
#        self.get_temporal_proteins()
#        self.analyze_poi(self.pois[0], self.pois_seq[0])



    def get_gb_file(self, accession_number, savetofile=False):
        '''
        A function to poll genbank given an accession number and pull the relevant gb file

        *args*

            **accession_number**, the accession number of the sequence to find.
            http://www.nslc.wustl.edu/elgin/genomics/bio4342/1archives/2006/AccReference.pdf

        *keyword args*

            **savetofile**, true or false to save the gb file in the same directory as sms for future use



        '''



        Entrez.email = "wsraymon@rams.colostate.edu"
        Entrez.tool = 'SingleMoleculeSimulator'

        er = False

        try:
            handle =  Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=accession_number)
            gb_record = SeqIO.read(handle, "genbank") #using "gb" as an alias for "genbank"
            handle.close()
        except:
            er = True


        time.sleep(2)

        if er == True:
            print('HTTP Error: Could not find specified ascession ID')

            return


        self.gb_rec = gb_record
        self.gb_obj = gb_record

        self.sequence_str = str(gb_record.seq)
        self.sequence_name = gb_record.name

        if savetofile:
            filename = self.sequence_name
            f = open(filename, 'w')


            f.write(self.gb_rec.format('gb'))

            f.close()



    def tau_plot(self,ssa_obj,t,tau,plot_type='contour', plot_all = False):

        
        
        stime = ssa_obj.time_rec-ssa_obj.start_time
        idx_t = (np.abs(stime - t)).argmin()
        idx_tau = (np.abs(stime - tau)).argmin()
        
        diff = idx_tau - idx_t
        difftime = t-tau
        
        if plot_type == 'Average':

            fig,ax= plt.subplots()
            for i in range(len(stime)-idx_tau,0,-4):
                idx_tau = (np.abs(stime- (stime[i]+difftime ))).argmin() 
                Itau = ssa_obj.intensity_vec[:,idx_tau]
                x,y = np.mean(ssa_obj.intensity_vec[:,idx_tau]/np.sum(ssa_obj.probe)),np.mean(ssa_obj.intensity_vec[:,idx_tau+diff]/np.sum(ssa_obj.probe))

        if plot_type == 'window':
            
            minx = 10000000
            maxx = 0
            
            miny = 10000000
            maxy = 0
            
            fig,ax= plt.subplots()
            for i in range(len(stime)-idx_tau,0,-10):
                idx_tau = (np.abs(stime - (idx_t+i))).argmin()  
                Itau = ssa_obj.intensity_vec[:,idx_tau]
                x,y = np.mean(ssa_obj.intensity_vec[:,idx_tau]/np.sum(ssa_obj.probe)),np.mean(ssa_obj.intensity_vec[:,idx_tau+diff]/np.sum(ssa_obj.probe))
                minx = min(np.min(x),minx)
                miny = min(np.min(y),miny)
                maxx = max(np.max(x),maxx)
                maxy = max(np.max(y),maxy)
                
                ax.scatter(x, y,zorder=3,color= cm.viridis_r(1.*i/len(stime)))


            c_map_ax = fig.add_axes([.95, 0.1, 0.1, 0.8])
          

            c_map_ax.axes.get_xaxis().set_visible(False)

            cbar = mpl.colorbar.ColorbarBase(c_map_ax, cmap=cm.viridis_r, orientation = 'vertical')
            
            
            cbar.ax.set_yticklabels(np.linspace(idx_t,stime[-1],6).astype(int) )
            cbar.ax.set_title('t')
            
            ax.plot([min(minx,miny),max(maxx,maxy)],[min(minx,miny),max(maxx,maxy)], color='red',ls='--')
            
            ax.set_ylabel(('<I(t=' + 't + tau'+')>'))
            ax.set_xlabel(('<I(t=' +'t'+')>'))
            ax.set_title(( 'Average I(t) vs Average I(t+tau) for tau = ' + str(diff) ) )
            
            
        if plot_type == 'density':
            fig,ax= plt.subplots()
            nbins = int(np.max(ssa_obj.intensity_vec/np.sum(ssa_obj.probe)))+2
            x, y = ssa_obj.intensity_vec[:,idx_t]/np.sum(ssa_obj.probe),ssa_obj.intensity_vec[:,idx_tau]/np.sum(ssa_obj.probe)
            k = kde.gaussian_kde([x,y])
            xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(), yi.flatten()])) 
            
            R = pearsonr(x,y)[0]
            ax.set_title(('Density Plot' + ' R = ' + str(np.round(R,3))))
            ax.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.viridis)
            ax.contour(xi, yi, zi.reshape(xi.shape) )   
            ax.set_ylabel(('I(t=' + str(tau)+')'))
            ax.set_xlabel(('I(t=' + str(t)+')'))
            fig.show()    
            
            
        if plot_type == 'set_tau':
            fig,ax= plt.subplots()
            for i in range(len(stime)-diff-idx_t):
                    idx_tau = (np.abs(stime - (idx_t+i))).argmin()                
                    plt.scatter(ssa_obj.intensity_vec[:,i]/np.sum(ssa_obj.probe), ssa_obj.intensity_vec[:,i+diff]/np.sum(ssa_obj.probe),c= cm.viridis(1.*i/len(stime)),alpha=.5  )
            plt.ylabel('I(t + s)')
            plt.xlabel(('I(t)'))
            plt.title(('Set tau, all times s = ' + str(diff) ))            

            c_map_ax = fig.add_axes([.95, 0.1, 0.1, 0.8])
          

            c_map_ax.axes.get_xaxis().set_visible(False)

            cbar = mpl.colorbar.ColorbarBase(c_map_ax, cmap=cm.viridis, orientation = 'vertical')
            
            cbar.ax.set_yticklabels(np.linspace(idx_t,stime[-1],6).astype(int) )
        
        if plot_type == 'scatter':
            if not plot_all:
                
                plt.scatter(ssa_obj.intensity_vec[:,idx_t]/np.sum(ssa_obj.probe), ssa_obj.intensity_vec[:,idx_tau]/np.sum(ssa_obj.probe) )
                plt.ylabel(('I(t=' + str(tau)+')'))
                
                
            else:
               
                for i in range(idx_t,len(stime)):
                    idx_tau = (np.abs(stime - (idx_t+i))).argmin()                
                    plt.scatter(ssa_obj.intensity_vec[:,idx_t]/np.sum(ssa_obj.probe), ssa_obj.intensity_vec[:,idx_tau]/np.sum(ssa_obj.probe),c= cm.viridis(1.*i/len(stime)),alpha=.1  )
                    plt.ylabel('I(tau)')
            plt.xlabel(('I(t=' + str(t)+')'))
            
        if plot_type == 'contour':
            fig,ax= plt.subplots()
            if not plot_all:
                It = ssa_obj.intensity_vec[:,idx_t]
                Itau = ssa_obj.intensity_vec[:,idx_tau]
                
                cov = np.cov(It,Itau)
                
                eigs, v = np.linalg.eig(cov)
                eigs = np.sqrt(eigs)
                plt.ylabel(('I(t=' + str(tau)+')'))
                colors = [cm.viridis(1.0),cm.viridis(.5),cm.viridis(0.0),cm.viridis(0.0)]
      
                
                for j in xrange(3, 0,-1):
                   
                    ell_artist = Ellipse(xy=(np.mean(It), np.mean(Itau)),
                                  width=eigs[0]*j*2, height=eigs[1]*j*2,
                                  angle=np.rad2deg(np.arccos(v[0, 0])))
                    
                    ell_artist.set_linewidth(2)
                    ell_artist.set_edgecolor(colors[j-1])
                    ell_artist.set_color(colors[j-1])
                    ax.add_patch(ell_artist)
                    
                ax.autoscale()      
                ax.set_xlim(0)
                ax.set_ylim(0)
                ax.scatter(It, Itau,zorder=3,alpha=0.3,color='red',marker='.')
                fig.show()
            else:
                plt.ylabel('I(tau)')
                It = ssa_obj.intensity_vec[:,idx_t]
                for i in range(len(stime)-idx_t,0,-10):
                    idx_tau = (np.abs(stime - (idx_t+i))).argmin()  
                    Itau = ssa_obj.intensity_vec[:,idx_tau]
                   
                    cov = np.cov(It,Itau)
                    
                    eigs, v = np.linalg.eig(cov)
                    eigs = np.sqrt(eigs)
                    
                    
                    j = 3
                    ell_artist = Ellipse(xy=(np.mean(It), np.mean(Itau)),
                                  width=eigs[0]*j*2, height=eigs[1]*j*2,
                                  angle=np.rad2deg(np.arccos(v[0, 0])))
                    
                    ell_artist.set_linewidth(2)
                    ell_artist.set_edgecolor( cm.viridis_r(1.*i/len(stime)))
                    ell_artist.set_color( cm.viridis_r(1.*i/len(stime)))
                    ax.autoscale()    
                    ax.add_patch(ell_artist)
                    ax.figure.canvas.draw()
            
                
                    
                plt.xlabel(('I(t=' + str(t)+')'))
                ax.set_xlim(0)
                ax.set_ylim(0)
                c_map_ax = fig.add_axes([.95, 0.1, 0.1, 0.8])
              

                c_map_ax.axes.get_xaxis().set_visible(False)

                cbar = mpl.colorbar.ColorbarBase(c_map_ax, cmap=cm.viridis_r, orientation = 'vertical')
                
                cbar.ax.set_yticklabels(np.linspace(idx_t,stime[-1],6).astype(int) )
                
                fig.show()       
            
    

    def kymograph(self,ssa_obj,n_traj,bg_intense=True,show_intense = True,tag = 0, show_col=True,col_size = 1.5, custom_fig = None, facecolor='black', *args,**kwargs):
        '''
        Constructs a kymograph of ribosome locations
        '''
        startfrags = 0
        for i in range(n_traj):
            startfrags += ssa_obj.frag_per_traj[i]
            
        endfrags = startfrags + ssa_obj.frag_per_traj[n_traj]
        fragments = ssa_obj.fragments[startfrags:endfrags]

        

        time = ssa_obj.time#[0:len(ssa_obj.time_rec)-1]
        

        
        
        
        
        if len(ssa_obj.intensity_vec.shape) ==3:
            ivec = ssa_obj.intensity_vec[tag][n_traj]
        else:        
            ivec = ssa_obj.intensity_vec[n_traj]
        ftimes = ssa_obj.fragtimes[startfrags:startfrags+endfrags]


        nfrag = fragments.shape[0]
        maxlen= fragments.shape[1]

        
        #plt.figure(figsize=(5,10))
        if show_intense == True:
            gs = gridspec.GridSpec(1, 2, custom_fig, width_ratios=[3, 1]) 
        else:
            gs = gridspec.GridSpec(1, 1)
        
        plt.subplot(gs[0])
        lenplot = np.max(fragments)
        maxin = np.max(ivec)
        ax = plt.gca()
        ax.set_facecolor(facecolor)

        if bg_intense == True:
            for i in range(len(time)):
                plt.plot([0,lenplot],[time[i],time[i]],color = cm.summer(1.*ivec[i]/maxin),lw=1)
            
        for i in range(nfrag):
            
            
            
            if maxlen <= np.where(fragments[i] > 0 )[0][-1]:       
                timeseg = time[ftimes[i]:ftimes[i]+maxlen]
                
                plt.plot(fragments[i][0:len(timeseg)] ,timeseg[::-1] )
                
            else:
                timeseg = time[ftimes[i]:]
                
                stop = np.where(fragments[i] > 0 )[0][-1]
                timelen = len(fragments[i][0:stop]) 

                plt.plot(fragments[i][0:stop]   ,timeseg[0:timelen],**kwargs )

        plt.xlabel('Ribosome position')
        plt.ylabel('Time (sec)')
        segtime = ssa_obj.time[0:len(ssa_obj.time_rec)]
        plt.ylim(ssa_obj.time_rec[-1], ssa_obj.time_rec[0])
                
        if show_col == True:
            try:
                col = ssa_obj.col_points[n_traj]

                plt.plot(col[:,0],col[:,1],color='#00ff00',markersize=col_size,linestyle='none',marker='o')
            except:
                pass

        if show_intense == True:
            plt.subplot(gs[1])
            ax = plt.gca()
            ax.set_facecolor(facecolor)
 
            plt.plot(ivec.T/ np.sum(ssa_obj.probe),segtime,**kwargs)
            plt.xlabel('Intensity (ump)')
            plt.xlim(0,30)
            plt.ylim(segtime[-1], segtime[0])
            
        
            
            plt.tight_layout()
                

            



   
    def generate_additional_ks(self,k_enters,k_pauses,k_jumps,k_stops,L):
    
        max_enter = 0
        max_pause = 0
        max_stop = 0
        max_jump = 0
        
        if k_enters != []:
            k_enters[:,0] = k_enters[:,0]+L*k_enters[:,1]
            k_enters[:,1] = k_enters[:,2]    
            k_enters = k_enters[:,0:2]
            max_enter = np.max( k_enters[:,0])
    
        if k_pauses != []:
            k_pauses[:,0] = k_pauses[:,0]+ L*k_pauses[:,1]
            k_pauses[:,1] = k_pauses[:,2]
            k_pauses = k_pauses[:,0:2]
            max_pause = np.max( k_pauses[:,0])
    
        if k_stops != []:
            k_stops[:,0] = k_stops[:,0]+L*k_stops[:,1]
            k_stops[:,1] = k_stops[:,2]    
            k_stops = k_stops[:,0:2]
            max_stop = np.max( k_stops[:,0])
        
        if k_jumps != []:
            k_jumps[:,0] = k_jumps[:,0]+ L*k_jumps[:,1]
            
            k_jumps[:,1] = k_jumps[:,2]+ L*k_jumps[:,3]
            k_jumps[:,2] = k_jumps[:,4]
            k_jumps = k_jumps[:,0:3]
            
            max_jump = max([np.max( k_jumps[:,0]),np.max( k_jumps[:,1])])
            
        max_loc = max(max_jump,max_stop,max_pause,max_enter)
        
        if max_loc <=L: 
            frames_used = 0
        if max_loc > L:
            frames_used = 1
        if max_loc > 2*L :
            frames_used = 1
        
        return k_enters, k_pauses, k_stops, k_jumps, frames_used
        
    
        
        
    
    def get_all_autocovariances(self,intensity_vec,time_vec,geneLength,shotnoise=True):
        '''
        Get all autocovariances for all 4 routines of normalization / means
        '''
        
        not_equal = False
        
        
    
        firstlen = len(intensity_vec[0])
        for traj in intensity_vec:
            if len(traj) != firstlen:
                not_equal = True
                    
            
        if not_equal == True:
            if shotnoise == True:
                new_ivec = []
                for traj in intensity_vec:
                    new_ivec.append(traj)
                
                ivec = new_ivec
            else:
                ivec = intensity_vec
            
            ivecflat = []
            
            maxtraj = 0
            
            for traj in ivec:
                ivecflat = ivecflat + traj.tolist()
                if len(traj) > maxtraj:
                    maxtraj = len(traj)
                
                
                
            ivecflat = np.array(ivecflat)
            
            ug = np.mean(ivecflat)
     
            varg = np.var(ivecflat)      
            ivecflat = 0
            
            autocorr_ui = np.zeros( (len(ivec), maxtraj))
            for i in range(len(ivec)):
              
                autocorr_ui[i,:len(ivec[i])] = self.get_acc2( (ivec[i]-np.mean(ivec[i]))/np.std(ivec[i]) ) 
                
            autocorr_ug = np.zeros((autocorr_ui.shape))
            for i in range(len(ivec)):
                autocorr_ug[i,:len(ivec[i])] = self.get_acc2(ivec[i]-ug)


            

            mean_autocorr_ug = np.mean(autocorr_ug.T, axis=1)
            mean_autocorr_ui = np.mean(autocorr_ui.T, axis=1)
            
            

            
            mean_autocorr_ug_norm = np.mean(autocorr_ug.T/varg, axis=1)
            
            
            autocorr_ui_norm = np.zeros((autocorr_ui.shape))
            for i in range(len(ivec)):
                autocorr_ui_norm[i,:len(ivec[i])] = self.get_acc2((ivec[i]-np.mean(ivec[i]))/np.std(ivec[i])) 
                
 
            if shotnoise == True:
                X = [1,2,3,4]
                V = mean_autocorr_ui[X]
                G0 = np.interp(0,X,V)
                
                
            else:
                for i in range(len(ivec)):
                    autocorr_ui_norm[i,:len(ivec[i])] = autocorr_ui_norm[i,:len(ivec[i])]/np.var(ivec[i])
                    
            
            if shotnoise == True: 
             
                mean_autocorr_ui_norm = np.mean(autocorr_ui_norm.T, axis=1)/G0  
            
            
            
            ntraj = len(ivec)                        
            
        else:
            if shotnoise == True:
                ivec = intensity_vec[:,1:]
            else:
                ivec = intensity_vec
            ug = np.mean(ivec.flatten())
     
            varg = np.var(ivec.flatten())
            
            ntraj = ivec.shape[0]
            
            autocorr_ui = np.zeros((ivec.shape))
            for i in range(ivec.shape[0]):
                autocorr_ui[i,:len(ivec[i])] = self.get_acc2( (ivec[i]-np.mean(ivec[i]))/np.std(ivec[i]) ) 
                
            autocorr_ug = np.zeros((ivec.shape))
            for i in range(ivec.shape[0]):
                autocorr_ug[i,:] = self.get_acc2(ivec[i]-ug)
                
                
            
            mean_autocorr_ug = np.mean(autocorr_ug.T, axis=1)
            mean_autocorr_ui = np.mean(autocorr_ui.T, axis=1)
            
            mean_autocorr_ug_norm = np.mean(autocorr_ug.T/varg, axis=1)
            
            
            autocorr_ui_norm = np.zeros((autocorr_ui.shape))
            for i in range(len(ivec)):
                autocorr_ui_norm[i,:len(ivec[i])] = self.get_acc2((ivec[i]-np.mean(ivec[i]))/np.std(ivec[i])) 
                
 
            if shotnoise == True:
                X = [1,2,3,4]
                V = mean_autocorr_ui[X]
                G0 = np.interp(0,X,V)
                
                
            else:
                for i in range(len(ivec)):
                    autocorr_ui_norm[i,:len(ivec[i])] = autocorr_ui_norm[i,:len(ivec[i])]/np.var(ivec[i])
                    
            
            if shotnoise == True: 
                
                mean_autocorr_ui_norm = np.mean(autocorr_ui_norm.T, axis=1)/G0  
            
            
                
            mean_autocorr_ui_norm = np.mean(autocorr_ui_norm.T, axis=1)
        
        sem_autocorr_ui_norm = 1.0/np.sqrt(ntraj)*np.std(autocorr_ui_norm.T,ddof=1,axis=1)
        if shotnoise==True:
            sem_autocorr_ui_norm = 1.0/np.sqrt(ntraj)*np.std(autocorr_ui_norm.T,ddof=1,axis=1)/G0**2
        
        sem_autocorr_ug_norm = 1.0/np.sqrt(ntraj)*np.std(autocorr_ug.T/varg,ddof=1,axis=1)
        
        sem_autocorr_ui = 1.0/np.sqrt(ntraj)*np.std(autocorr_ui.T,ddof=1,axis=1)
        sem_autocorr_ug = 1.0/np.sqrt(ntraj)*np.std(autocorr_ug.T,ddof=1,axis=1)
        
        
        dwelltime_ui = None
        try:
            dwelltime_ui = time_vec[np.where(mean_autocorr_ui < .01)[0][0]]
        except:
            try:
                dwelltime_ui = time_vec[np.where(mean_autocorr_ui < .05)[0][0]]
            except:
                dwelltime_ui = 1
        ke_exp_ui = np.round(geneLength/dwelltime_ui ,1)

        dwelltime_ug = None
        try:
            dwelltime_ug = time_vec[np.where(mean_autocorr_ug < .01)[0][0]]
        except:
            try:
                dwelltime_ug = time_vec[np.where(mean_autocorr_ug < .05)[0][0]]
            except:
                dwelltime_ug = 1
        ke_exp_ug = np.round(geneLength/dwelltime_ug ,1)
        
        dwelltime_ui_norm = None
        try:
            dwelltime_ui_norm = time_vec[np.where(mean_autocorr_ui_norm < .01)[0][0]]
        except:
            try:
                dwelltime_ui_norm = time_vec[np.where(mean_autocorr_ui_norm < .05)[0][0]]
            except:
                dwelltime_ui_norm = 1
        ke_exp_ui_norm = np.round(geneLength/dwelltime_ui_norm ,1)

        dwelltime_ug_norm = None
        try:
            dwelltime_ug_norm = time_vec[np.where(mean_autocorr_ug_norm < .01)[0][0]]
        except:
            try:
                dwelltime_ug_norm = time_vec[np.where(mean_autocorr_ug_norm < .05)[0][0]]
            except:
                dwelltime_ug_norm = 1
        ke_exp_ug_norm = np.round(geneLength/dwelltime_ug_norm ,1)        
        
        nacov = {'global': {'sem': sem_autocorr_ug_norm, 'mean':mean_autocorr_ug_norm,'traj':autocorr_ug/varg,'ke': ke_exp_ug_norm, 'dwelltime':dwelltime_ug_norm},
                 'indiv': {'sem': sem_autocorr_ui_norm, 'mean':mean_autocorr_ui_norm,'traj':autocorr_ui_norm,'ke': ke_exp_ui_norm, 'dwelltime':dwelltime_ui_norm}}
        acov =  {'global': {'sem': sem_autocorr_ug, 'mean':mean_autocorr_ug,'traj':autocorr_ug,'ke': ke_exp_ug, 'dwelltime':dwelltime_ug },
                 'indiv': {'sem': sem_autocorr_ui, 'mean':mean_autocorr_ui,'traj':autocorr_ui,'ke': ke_exp_ui, 'dwelltime':dwelltime_ui}}
        
        return nacov,acov



    

        
    
    

           


