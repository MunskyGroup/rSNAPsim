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





class SequenceManipMethods():
    '''
    class that handles anything dealing with sequences
    '''
    def __init__(self,sequence):
        self.sequence = sequence
        self.codon_dicts = CodonDictionaries()
        pass    
    
    
    def nt2aa(self,nt_seq):
        '''
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.

        Returns
        -------
        aa : str
            amino acid sequence.

        '''
        aa = ''
        nt_seq = nt_seq.upper()
        for i in range(0, len(nt_seq), 3):
            aa += self.codon_dicts.aa_table[nt_seq[i:i+3]]
        return aa    

    def get_orfs(self, nt_seq='', min_codons=80):

        '''
        Returns open reading frames of the nucleotide sequence given

        orfs = {'1':[proteins],
                '2':[proteins],
                '3':[proteins]}

        *keyword args*

            **nt_seq**, nucleotide sequence as a string. If left blank uses
            the self.sequence_str

            **min_codons**, minimum amount of codons to be considered
            a protein in the open reading frame
        '''

        if nt_seq == '':
            nt_seq = self.sequence.upper()
        nt_seq = nt_seq.upper()
        allstarts = np.array([m.start() for m in re.finditer('(?=A[TU]G((?:.{3})+?)[TU](?:AG|AA|GA))', nt_seq)])
        
 
        #allsegments = re.findall('(?=A[TU]G((?:.{3})+?)[TU](?:AG|AA|GA))',self.sequence_str)
        allstops = np.array([m.start() for m in re.finditer('(?=[TU](?:AG|AA|GA))', nt_seq)])
        start_frames = allstarts%3
        stop_frames = allstops%3
        min_len = min_codons*3
        orf1_starts = allstarts[np.where(start_frames == 0)]
        orf2_starts = allstarts[np.where(start_frames == 1)]
        orf3_starts = allstarts[np.where(start_frames == 2)]

        orf1_stops = allstops[np.where(stop_frames == 0)]
        orf2_stops = allstops[np.where(stop_frames == 1)]
        orf3_stops = allstops[np.where(stop_frames == 2)]

        #self.starts = [orf1_starts, orf2_starts, orf3_starts]
        #self.stops = [orf1_stops, orf2_stops, orf3_stops]
        orfs = {'1':[], '2':[], '3':[]}


        orfs = {'1':[], '2':[], '3':[]}
        laststop = 0
        for start in orf1_starts:
            nextstop = orf1_stops[np.where(orf1_stops > start)[0][0]]
            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    orfs['1'].append((start, nextstop))

                    laststop = nextstop

        laststop = 0
        for start in orf2_starts:
            nextstop = orf2_stops[np.where(orf2_stops > start)[0][0]]
            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    orfs['2'].append((start, nextstop))
                    laststop = nextstop

        laststop = 0
        for start in orf3_starts:
            nextstop = orf3_stops[np.where(orf3_stops > start)[0][0]]

            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    orfs['3'].append((start, nextstop))
                    laststop = nextstop
                    
        
                    
        return orfs
    
    def codon_usage(self, nt_seq):
        '''
        Analyzes codon useage from the nucleotide sequence

        *args*

            **nt_seq**,  nucleotide sequence as a string

        *returns*

            **codon_sensitivity**, a list of codon sensitivity for the nucleotide sequence

            **cai**, cai value

        '''


        codon_usage = np.zeros((1, 21))
        gene_len = len(nt_seq)/3
        aa_seq = self.nt2aa(nt_seq)

        for i in range(len(self.codon_dicts.aa_keys)-1):

            codon_usage[0, i] = len(re.findall(self.codon_dicts.aa_keys[i], aa_seq))
        codon_usage[0, 20] = len(re.findall('\*', aa_seq))
        codon_norm = codon_usage/gene_len
        codon_sensitivity = np.round(codon_norm*self.codon_dicts.sensitivity_fast_slow, 2)

        cai_codons = []
        for i in range(0, len(nt_seq), 3):
            cai_codons.append(self.codon_dicts.strGeneCopy[nt_seq[i:i+3]] / self.codon_dicts.strGeneCopy_fast[nt_seq[i:i+3]])

        cai = self.geomean(cai_codons)

        return codon_sensitivity, cai, cai_codons
    
    
    
    def get_proteins(self,orfs,seq):
        '''
        Parameters
        ----------
        orfs : dict
            dictionary of open reading frames.
            {'1': [[starts],[stops] ],'2': [[starts],[stops] ],'3': [[starts],[stops] ] }
        seq : str
            nucleotide sequence.

        Returns
        -------
        proteins_strs : dict
            aa strings of all proteins found in the given orfs.
        protein_objs : dict
            container objects for proteins found in the given orf.
        proteins_w_tags : dict
            conatiner objects for any proteins with detected tags.

        '''
        cd = self.codon_dicts
        proteins_strs = {'1':[], '2':[], '3':[]}
        protein_objs = {'1':[], '2':[], '3':[]}
        proteins_w_tags = {'1':[], '2':[], '3':[]}
        
        
        tagged_proteins = {a:[] for a in cd.tag_dict.keys()}
        tagged_protein_seq = {a:[] for a in cd.tag_dict.keys()}

        for i in range(len(orfs)):
            for j in range(len(orfs[str(i+1)])):
                
                protein = poi()
                
                pro = self.nt2aa(seq[orfs[str(i+1)][j][0]:orfs[str(i+1)][j][1]+3])
                nt_seq = seq[orfs[str(i+1)][j][0]:orfs[str(i+1)][j][1]+3]
                if pro[-1] == '*':
                    pro = pro[:-1]
                    nt_seq = nt_seq[:-3]
                
                
                
                
                
                protein.aa_seq = pro
                protein.nt_seq = nt_seq
                
                proteins_strs[str(i+1)].append(pro)
                
                protein.gene_length =  len(pro) #length of the gene
                protein.tag_length = 0   #length of the tags
                protein.total_length = len(pro)  #total length of the full amino acid sequence
                protein.source_seq = seq
                protein.orf = i
                protein.loc = (orfs[str(i+1)][j][0],orfs[str(i+1)][j][1]+3)      
                protein.tags = []
                
                protein_objs[str(i+1)].append(protein)
                
                
                
                
        for i in range(len(orfs)):  
            for pr in protein_objs[str(i+1)]:
                
                tag_detected = False
                
                for tag in cd.tag_dict.keys():
                   
                    if cd.tag_dict[tag] in pr.aa_seq:
                        tag_detected = True
                        
                        
                if tag_detected:
                    self.analyze_protein_w_tags(pr)
                    pr.tag_added = False
                    proteins_w_tags[str(i+1)].append(protein)
                else:
                    self.add_tag_to_protein(pr)
                    pr.tag_added = True
                

                        
        return proteins_strs, protein_objs, proteins_w_tags
        
    def add_tag_to_protein(self,POI,tag_type='T_Flag'):
        '''
        Parameters
        ----------
        POI : poi object
            protein of interest object.
        tag_type : str, optional
            What kind of tag to append onto the protein object. The default is 'T_Flag'.

        Returns
        -------
        None.

        '''
        cd = self.codon_dicts
        
        POI.nt_seq = cd.tag_full[tag_type] + POI.nt_seq
        POI.aa_seq = self.nt2aa(POI.nt_seq)
        
        self.analyze_protein_w_tags(POI)
        
        
    
    
    def analyze_protein_w_tags(self,POI,epitope_loc='front'):
        cd = self.codon_dicts
        nt_seq = POI.nt_seq
        aa_seq = POI.aa_seq
        #self.POI.name = self.sequence_name
        total_length = len(POI.aa_seq)

        '''
        for key in self.tagged_proteins:
            if protein in self.tagged_proteins[key]:
                self.POI.tag_types.append(key)
        '''
        POI.tag_types = []
        for tag in cd.tag_dict.keys():
            if cd.tag_dict[tag] in aa_seq:
                POI.tag_types.append(tag)

                #''.join(sms.poi[0].split('DYKDDDDK')

        POI.tag_epitopes = {a:[] for a in POI.tag_types}
        gs = POI.aa_seq


        for i in range(len(POI.tag_types)):
            
            try:
                nt_tag = cd.tag_full[POI.tag_types[i]]
                aa_tag = self.nt2aa(nt_tag)
            except:
                epi = cd.tag_dict[POI.tag_types[i]]
                firstep = POI.aa_seq.find(epi) 
                lastep = len(POI.aa_seq) - POI.aa_seq[::-1].find(epi[::-1])                
                aa_tag = POI.aa_seq[firstep:lastep]
                nt_tag = POI.nt_seq[3*firstep:3*lastep]
                
            if epitope_loc == 'front':
                offset = 0
            if epitope_loc == 'middle':
                offset = int(len(cd.tag_dict[POI.tag_types[i]])/2)
            if epitope_loc == 'back':
                offset = len(cd.tag_dict[POI.tag_types[i]])
                
            POI.tag_epitopes[POI.tag_types[i]] = [m.start()+1+offset for m in re.finditer(cd.tag_dict[POI.tag_types[i]], POI.aa_seq)]

            gs = gs.replace(aa_tag, '')

        POI.gene_seq = gs
        POI.gene_length = len(gs)
        POI.total_length = total_length
        POI.tag_seq = aa_tag 
        POI.tag_length = len(aa_tag)
        
        codons = []
        for i in range(0, len(nt_seq), 3):
            codons.append(nt_seq[i:i+3])
        #POI.codons = codons

        #POI.codon_sensitivity, POI.CAI, POI.CAI_codons = self.codon_usage(POI.nt_seq)   
        POI.ki = .03
        POI.ke = 10
        POI.kt = 10
        
    def get_tag_loc(self,aa_seq,tag,epitope_loc='front'):
        cd = self.codon_dicts
        
        if epitope_loc == 'front':
            offset = 0
        if epitope_loc == 'middle':
            offset = int(len(tag)/2)
        if epitope_loc == 'back':
            offset = len(tag)
                
        return [m.start()+1+offset for m in re.finditer(tag, aa_seq)]

        
    @staticmethod
    def geomean(iterable):
        '''geometric mean used for codon sensitivity calculations
        '''
        a = np.array(iterable)
        return a.prod()**(1.0/len(a))    


class PropensityFactory():
    '''
    factory class for the k's
    '''
    def __init__(self):
        self.codon_dicts = CodonDictionaries()
        pass    

    @staticmethod
    def bin_k(k,inds):
        try: 
            k = k.tolist()
        except:
            pass
        k_binned = np.zeros(len(inds)-1)
        binned_ks = []
        for i in range(0,len(inds)-1):
            binned_ks = binned_ks +   [k[inds[i]:inds[i+1]],]  
        
         
        for i in range(0,len(inds)-1):            
            k_binned[i] = 1/ np.sum( 1/np.array(binned_ks[i]))
        return k_binned


    def get_trna_ids(self,nt_seq):
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] 
        return [self.codon_dicts.trna_dict[x] for x in seperated_codons]
        
        

    def get_k(self, nt_seq, k_init, k_elong_mean,k_end):
        '''
        returns all propensities for a given nucleotide sequence

        *args*

            **nt_seq**,  nucleotide sequence as a string

            **k_initiation**, initiation rate of ribosome binding

            **k_elong_mean**, average rate of elgonation experimentally found


        '''
        codons = nt_seq.upper()
        
        genelength = int(len(codons)/3)
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
        k_elongation = np.zeros((1, genelength))
        tRNA_copynumber = np.zeros((1, genelength))

     
        for i in range(len(seperated_codons)):
            tRNA_copynumber[0, i] = self.codon_dicts.strGeneCopy[seperated_codons[i]]

        mean_tRNA_copynumber = np.mean(list(self.codon_dicts.strGeneCopy_single.values()))

        k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean
        all_k = [k_init] + k_elongation.flatten().tolist() + [k_end]
        
        return all_k
    
    
    def get_k_3_frame(self,nt_seq,k_elong_mean):
        '''
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.
        k_elong_mean : float
            mean elongation rate aa/s.

        Returns
        -------
        kelongs : list
            elongation rates for all 3 frames of a nucleotide sequence stacked in:
                [0+frame L, 1+frame L-1, 2+frame L-1] .

        '''
        kelongs = []
        
        for n in range(3):
            if n !=0:
                codons = nt_seq[n:-(3-n)]
            else:
                codons = nt_seq
            genelength = int(len(codons)/3)
            seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
            k_elongation = np.zeros((1, genelength))
            tRNA_copynumber = np.zeros((1, genelength))

     
            for i in range(len(seperated_codons)):
                tRNA_copynumber[0, i] = self.codon_dicts.strGeneCopy[seperated_codons[i]]
    
            mean_tRNA_copynumber = np.mean(list(self.codon_dicts.strGeneCopy_single.values()))
    
            k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean
        
           # k_elongation.flatten().tolist()[:-1]
        
            kelongs = kelongs + k_elongation.flatten().tolist()
        
        return kelongs
        
    
    @staticmethod
    def get_binned_k(k,bins):
        '''
        evenly bins elongation rates as best it can.
        
        '''
        binsize = int(np.floor(len(k)/bins))
        binned_ks = []
        
        k_binned = np.zeros(bins)
        k_lens = np.ones(bins)*binsize
        
        to_redistribute = len(k)%bins

        k_lens[-to_redistribute:] = binsize+1
        
        inds = np.hstack(([0.], np.cumsum(k_lens))).astype(int)     

        
        for i in range(0,bins):
            binned_ks = binned_ks +   [k[inds[i]:inds[i+1]].tolist(),]  
         
        for i in range(0,bins):            
            k_binned[i] = np.mean(binned_ks[i])/len(binned_ks[i])

        return k_binned,k_lens
    
    
    @staticmethod
    def intellegent_bin(pl,nbins ,min_bin = 1):
        '''
        Function to do intellegent binning, focuses resolution on the areas
        defined in the probe location vector
        
        Note if you pass it a minium bin that when min_bin*nbins > length of your sequence
        this function will fail
        '''
        
        if min_bin*nbins > pl.shape[1]:
            warnings.warn('Desired minimum binsize and desired number of bins is not possible with the length of the probe vector, returning best guess')

        pl_inds = np.where(pl == 1)[1]
       
        if 0 not in pl_inds:
            pl_inds = np.hstack((np.array([0]), pl_inds))
        if pl.shape[1]+1 not in pl_inds:
            pl_inds = np.hstack(( pl_inds,np.array(pl.shape[1]+1) ))
     
        used_bins = len(pl_inds)
        k = len(pl_inds)-1
        j = 0
        to_add = []
        while used_bins < nbins+1:
            if j == k:
                
                j = 0

                prev_pl_inds = pl_inds
                pl_inds = np.hstack(( pl_inds,np.array(to_add)) )
                pl_inds = np.sort(pl_inds)
             
                if np.array_equal(prev_pl_inds,pl_inds):
                    break
                k =  len(pl_inds)-1
                to_add = []

            newbin = int( pl_inds[j]   +  (pl_inds[j+1]-pl_inds[j]  )/2)

            if newbin not in pl_inds:
                if not (np.abs(pl_inds - newbin) <= min_bin).any():
                    to_add.append(newbin)
                    used_bins+=1
            j+=1    
            
        pl_inds = np.hstack(( pl_inds,np.array(to_add)) )
        pl_inds = np.sort(pl_inds) 
        
        return pl_inds

    @staticmethod
    def even_bin(length,nbins):
        '''
        evenly bins a length over a given amount of bins as best it can
        
        '''
        binsize = int(np.floor(length/nbins))
        
        k_lens = np.ones(nbins)*binsize
        
        to_redistribute = length%nbins

        k_lens[-to_redistribute:] = binsize+1
        
        inds = np.hstack(([0.], np.cumsum(k_lens))).astype(int)     
        
        return inds
    
class ProbeVectorFactory():
    def __init__(self):
        pass        
    
    def get_probe_vec(self, tag_epitope_dict, length):
        
        pv = np.zeros( (len(list(tag_epitope_dict)), length))
        for i in range(len(list(tag_epitope_dict))):
            pv[i,[tag_epitope_dict[list(tag_epitope_dict.keys())[i]]]] = 1
        pv = np.cumsum(pv,axis=1)        
        return pv

    def get_probe_loc(self, tag_epitope_dict, length):
        
        pv = np.zeros( (len(list(tag_epitope_dict)), length))
        for i in range(len(list(tag_epitope_dict))):
            pv[i,[tag_epitope_dict[list(tag_epitope_dict.keys())[i]]]] = 1
        return pv

    @staticmethod
    def bin_probe_vecs(probe_loc,inds):

        probeloc_binned = np.zeros((probe_loc.shape[0],   len(inds)-1 ) )
        for i in range(0,len(inds)-1):
            probeloc_binned[:,i] = np.sum(probe_loc[:,inds[i]:inds[i+1]],axis=1)
            
            
        #probeloc_binned[:,-1] = np.sum(probe_loc[:,inds[-1]:],axis=1)
        probevec_binned = np.cumsum(probeloc_binned,axis=1)
        
        return probeloc_binned.astype(int),probevec_binned.astype(int)
        
    

    @staticmethod
    def intellegent_bin(pl,nbins ,min_bin = 1):
        '''
        Function to do intellegent binning, focuses resolution on the areas
        defined in the probe location vector
        
        Note if you pass it a minium bin that when min_bin*nbins > length of your sequence
        this function will fail
        '''
        
        if min_bin*nbins > pl.shape[1]:
            warnings.warn('Desired minimum binsize and desired number of bins is not possible with the length of the probe vector, returning best guess')

        pl_inds = np.where(pl == 1)[1]
        if 0 not in pl_inds:
            pl_inds = np.hstack((np.array([0]), pl_inds))
        if pl.shape[1] not in pl_inds:
            pl_inds = np.hstack(( pl_inds,np.array(pl.shape[1]) ))
     
        used_bins = len(pl_inds)
        k = len(pl_inds)-1
        j = 0
        to_add = []
        while used_bins < nbins+1:
            if j == k:
                
                j = 0

                prev_pl_inds = pl_inds
                pl_inds = np.hstack(( pl_inds,np.array(to_add)) )
                pl_inds = np.sort(pl_inds)
             
                if np.array_equal(prev_pl_inds,pl_inds):
                    break
                k =  pl.shape[1]-1
                to_add = []
                
            newbin = int( pl_inds[j]   +  (pl_inds[j+1]-pl_inds[j]  )/2)

            if newbin not in pl_inds:
                if not (np.abs(pl_inds - newbin) <= min_bin).any():
                    to_add.append(newbin)
                    used_bins+=1
            j+=1    
            
        pl_inds = np.hstack(( pl_inds,np.array(to_add)) )
        pl_inds = np.sort(pl_inds) 
        
        return pl_inds
        
    @staticmethod
    def even_bin(length,nbins):
        '''
        Parameters
        ----------
        length : int
            Length of the vector to bin.
        nbins : int
            Number of desired bins.

        Returns
        -------
        inds : TYPE
            n bin locations over the vector length given.

        '''
        
        binsize = int(np.floor(length/nbins))
        
        k_lens = np.ones(nbins)*binsize
        
        to_redistribute = length%nbins

        k_lens[-to_redistribute:] = binsize+1
        
        inds = np.hstack(([0.], np.cumsum(k_lens))).astype(int)     
        
        return inds            


class CodonOptimizer():
    def __init__(self):
        self.codon_dict = CodonDictionaries()
    
    def optimize_ntseq(self,nt_seq,dictionary = None):
        if dictionary == None:
            dictionary = self.codon_dict.strGeneCopy
            
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3    
        aa = [ self.codon_dict.aa_table[x] for x in seperated_codons ]     
        opt_seq = ''
        for i in range(0,len(aa)):
            ind = np.argmax([self.codon_dict.strGeneCopy[x] for x in self.codon_dict.aa_table_r[aa[i]]])
            opt_codon = self.codon_dict.aa_table_r[aa[i]][ind]
            opt_seq = opt_seq + opt_codon
        return opt_seq
            
        
    def deoptimize_ntseq(self,nt_seq,dictionary = None):
        if dictionary == None:
            dictionary = self.codon_dict.strGeneCopy
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3    
        aa = [ self.codon_dict.aa_table[x] for x in seperated_codons ]     
        opt_seq = ''
        for i in range(0,len(aa)):
            ind = np.argmin([self.codon_dict.strGeneCopy[x] for x in self.codon_dict.aa_table_r[aa[i]]])
            opt_codon = self.codon_dict.aa_table_r[aa[i]][ind]
            opt_seq = opt_seq + opt_codon
        return opt_seq        



class CodonDictionaries():
    
    def __init__(self):
        self.tag_dict = {'T_SunTag':'EELLSKNYHLENEVARLKK',
                         'T_Flag':'DYKDDDDK',
                         'T_Hemagglutinin':'YPYDVPDYA'}
        
        self.tag_colors = {'T_SunTag':'green',
                         'T_Flag':'blue',
                         'T_Hemagglutinin':'blue'}
        
        self.tag_full = {'T_Flag':('ATGGACTACAAGGACGACGACGACAAAGGTGAC'
                                   'TACAAAGATGATGACGATAAAGGCGACTATA'
                                   'AGGACGATGACGACAAGGGCGGAAACTCACTGA'
                                   'TCAAGGAAAACATGCGGATGAAGGTGGTGAT'
                                   'GGAGGGCTCCGTGAATGGTCACCAGTTCAAGTG'
                                   'CACCGGAGAGGGAGAGGGAAACCCGTACATG'
                                   'GGAACTCAGACCATGCGCATTAAGGTCATCGAA'
                                   'GGAGGTCCGCTGCCGTTCGCTTTCGATATCC'
                                   'TGGCCACTTCGTTCGGAGGAGGGTCGCGCACGTTC'
                                   'ATCAAGTACCCGAAGGGAATCCCGGACTT'
                                   'CTTTAAGCAGTCATTCCCGGAAGGATTCACTTGGG'
                                   'AACGGGTGACCCGGTATGAAGATGGAGGT'
                                   'GTGGTGACTGTCATGCAAGATACTTCGCTGGAGGATGGG'
                                   'TGCCTCGTGTACCACGTCCAAGTCC'
                                   'GCGGAGTGAATTTCCCGTCCAACGGACCAGTGATGCAG'
                                   'AAAAAGACGAAGGGTTGGGAACCTAA'
                                   'TACTGAAATGATGTACCCCGCAGACGGAGGGCTGAGGG'
                                   'GCTACACCCACATGGCGCTGAAGGTC'
                                   'GACGGAGGAGATTACAAGGATGACGACGATAAGCAACAA'
                                   'GATTACAAAGACGATGATGACAAGG'
                                   'GCCAGCAGGGCGACTACAAGGACGACGACGACAAGCAG'
                                   'CAGGACTACAAAGATGACGATGATAA'
                                   'AGGAGGAGGACATCTGTCCTGTTCGTTCGTGACCACCT'
                                   'ACAGATCAAAGAAAACCGTGGGAAAC'
                                   'ATCAAGATGCCGGGCATTCATGCCGTCGACCACCGCCT'
                                   'GGAGCGGCTCGAAGAATCAGACAATG'
                                   'AGATGTTCGTCGTGCAAAGAGAACATGCCGTGGCCAAGTT'
                                   'CGCGGGACTGGGAGGCGGTGGAGG'
                                   'CGATTACAAAGACGATGATGACAAGGGTGACTATAAAGA'
                                   'CGACGATGACAAAGGGGATTACAAG'
                                   'GATGATGATGATAAGGGAGGCGGTGGATCAGGTGGAG'
                                   'GAGGTTCACTGCAG')}

        self.aa_keys = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
                       'P', 'S', 'T', 'W', 'Y', 'V', '*']

        self.codon_types = dict(zip(self.aa_keys, np.ones((1, 21)).flatten().astype(int).tolist()))

        self.aa_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
            
            'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
            'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
            'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
            'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
            'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
            'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
            'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
            'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
            'UAC':'Y', 'UAU':'Y', 'UAA':'*', 'UAG':'*',
            'UGC':'C', 'UGU':'C', 'UGA':'*', 'UGG':'W',}

        self.aa_table_r = {'A':['GCA', 'GCC', 'GCG', 'GCT','GCU'],
                          'R':['CGA', 'CGC', 'CGG', 'CGT','AGG','AGA','CGU'],
                          'N':['AAC', 'AAT','AAU'],
                          'D':['GAC', 'GAT','GAU'],
                          'C':['TGC', 'TGT','UGC','UGU'],
                          'Q':['CAA', 'CAG'],
                          'E':['GAA', 'GAG'],
                          'G':['GGT', 'GGC', 'GGA', 'GGC','GGU'],
                          'H':['CAC', 'CAT','CAU'],
                          'I':['ATT', 'ATC', 'ATA','AUU','AUC','AUA'],
                          'L':['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG','CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'],
                          'K':['AAA', 'AAG'],
                          'M':['ATG','AUG'],
                          'F':['TTC', 'TTT','UUC','UUU'],
                          'P':['CCT', 'CCC', 'CCG', 'CCA','CCU'],
                          'S':['TCA', 'TCC', 'TCG', 'TCT','AGC','AGT','UCA','UCC','UCG'],
                          'T':['ACA', 'ACC', 'ACG', 'ACT','ACU'],
                          'W':['TGG','UGG'],
                          'Y':['TAT', 'TAC','UAC','UAU'],
                          'V':['GTA', 'GTC', 'GTT','GTG','GUG','GUU','GUC','GUA'],
                          '*':['TGA', 'TAG', 'TAA','UGA','UAG','UAA']
                         }


        self.strGeneCopy = {'TTT': 17.6, 'TCT': 15.2, 'TAT': 12.2, 'TGT': 10.6, 'TTC': 20.3,
                            'TCC': 17.7, 'TAC': 15.3, 'TGC': 12.6, 'TTA': 7.7, 'TCA': 12.2,
                            'TAA': 1.0, 'TGA': 1.6, 'TTG': 12.9, 'TCG':  4.4, 'TAG': 0.8,
                            'TGG': 13.2, 'CTT': 13.2, 'CCT': 17.5, 'CAT': 10.9, 'CGT': 4.5,
                            'CTC': 19.6, 'CCC': 19.8, 'CAC': 15.1, 'CGC': 10.4, 'CTA':  7.2,
                            'CCA': 16.9, 'CAA': 12.3, 'CGA':  6.2, 'CTG': 39.6, 'CCG':  6.9,
                            'CAG': 34.2, 'CGG': 11.4, 'ATT': 16.0, 'ACT': 13.1, 'AAT': 17.0,
                            'AGT': 12.1, 'ATC': 20.8, 'ACC': 18.9, 'AAC': 19.1, 'AGC': 19.5,
                            'ATA':  7.5, 'ACA': 15.1, 'AAA': 24.4, 'AGA': 12.2, 'ATG': 22.0,
                            'ACG': 6.1, 'AAG': 31.9, 'AGG': 12.0, 'GTT': 11.0, 'GCT': 18.4,
                            'GAT': 21.8, 'GGT': 10.8, 'GTC': 14.5, 'GCC': 27.7, 'GAC': 25.1,
                            'GGC': 22.2, 'GTA':  7.1, 'GCA': 15.8, 'GAA': 29.0, 'GGA': 16.5,
                            'GTG': 28.1, 'GCG': 7.4, 'GAG': 39.6, 'GGG': 16.5}
        
        self.strGeneCopy_single = {'TTT': 17.6, 'TCT': 15.2, 'TAT': 12.2, 'TGT': 10.6, 'TTC': 20.3,
                            'TCC': 17.7, 'TAC': 15.3, 'TGC': 12.6, 'TTA': 7.7, 'TCA': 12.2,
                            'TAA': 1.0, 'TGA': 1.6, 'TTG': 12.9, 'TCG':  4.4, 'TAG': 0.8,
                            'TGG': 13.2, 'CTT': 13.2, 'CCT': 17.5, 'CAT': 10.9, 'CGT': 4.5,
                            'CTC': 19.6, 'CCC': 19.8, 'CAC': 15.1, 'CGC': 10.4, 'CTA':  7.2,
                            'CCA': 16.9, 'CAA': 12.3, 'CGA':  6.2, 'CTG': 39.6, 'CCG':  6.9,
                            'CAG': 34.2, 'CGG': 11.4, 'ATT': 16.0, 'ACT': 13.1, 'AAT': 17.0,
                            'AGT': 12.1, 'ATC': 20.8, 'ACC': 18.9, 'AAC': 19.1, 'AGC': 19.5,
                            'ATA':  7.5, 'ACA': 15.1, 'AAA': 24.4, 'AGA': 12.2, 'ATG': 22.0,
                            'ACG': 6.1, 'AAG': 31.9, 'AGG': 12.0, 'GTT': 11.0, 'GCT': 18.4,
                            'GAT': 21.8, 'GGT': 10.8, 'GTC': 14.5, 'GCC': 27.7, 'GAC': 25.1,
                            'GGC': 22.2, 'GTA':  7.1, 'GCA': 15.8, 'GAA': 29.0, 'GGA': 16.5,
                            'GTG': 28.1, 'GCG': 7.4, 'GAG': 39.6, 'GGG': 16.5}
        
        self.trna_ids = ['TTT','TCT','TAT','TGT','TTC','TCC','TAC','TGC','TTA','TCA','TAA','TGA','TTG',
                        'TCG','TAG','TGG', 'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CTA',
                        'CCA', 'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG', 'ATT', 'ACT', 'AAT', 'AGT', 'ATC',
                        'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ATG', 'ACG', 'AAG', 'AGG', 'GTT',
                        'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA', 'GCA', 'GAA', 'GGA', 'GTG',
                        'GCG', 'GAG', 'GGG']
        
        self.trna_ids.remove('TAG')
        self.trna_ids.remove('TAA')
        self.trna_ids.remove('TGA')
        
        trna_temp = [x.replace('T','U') for x in self.trna_ids]
        self.trna_ids = self.trna_ids + trna_temp
        
        
        
        
        self.trna_ids_vals = np.linspace(0,60,61).astype(int).tolist() +  np.linspace(0,60,61).astype(int).tolist()
        
        self.trna_dict = dict(zip(self.trna_ids, self.trna_ids_vals))
        self.id_to_trna = dict(map(reversed, self.trna_dict.items()))

        
        
        # add the U codons
        for key in list(self.strGeneCopy.keys()):
            if 'T' in key:
                val = self.strGeneCopy[key]
                newkey = key.replace('T','U')
                self.strGeneCopy[newkey] = val



        self.strGeneCopy_fast = {'GCT': 27.7, 'GCC': 27.7, 'GCA': 27.7, 'GCG': 27.7,  #A
                                 'CGT': 12.2, 'CGC': 12.2, 'CGA': 12.2, 'CGG': 12.2,
                                 'AGA': 12.2, 'AGG': 12.2,   # R
                                 'AAT': 19.1, 'AAC': 19.1,   #N
                                 'GAT': 25.1, 'GAC': 25.1,   # D
                                 'TGT': 12.6, 'TGC': 12.6,  # C
                                 'CAA': 34.2, 'CAG': 34.2,  # Q
                                 'GAA': 39.6, 'GAG': 39.6,  #E
                                 'GGT': 22.2, 'GGC': 22.2, 'GGA': 22.2, 'GGG': 22.2,  # G
                                 'CAT': 15.1, 'CAC': 15.1,  # H
                                 'ATT': 20.8, 'ATC': 20.8, 'ATA': 20.8,  # I
                                 'TTA': 39.6, 'TTG': 39.6, 'CTT': 39.6, 'CTC': 39.6,
                                 'CTA': 39.6, 'CTG': 39.6, # L
                                 'AAA': 31.9, 'AAG': 31.9,  # K
                                 'ATG': 22.0,   #M
                                 'TTT': 20.3, 'TTC': 20.3,    # F
                                 'CCT': 19.8, 'CCC': 19.8, 'CCA': 19.8, 'CCG': 19.8,  # P
                                 'TCT': 19.5, 'TCC': 19.5, 'TCA': 19.5, 'TCG': 19.5,
                                 'AGT': 19.5, 'AGC': 19.5,  # S
                                 'ACT': 18.9, 'ACC': 18.9, 'ACA': 18.9, 'ACG': 18.9, # T
                                 'TGG': 13.2,   #W
                                 'TAT': 15.3, 'TAC': 15.3,  # Y
                                 'GTT': 28.1, 'GTC': 28.1, 'GTA':28.1, 'GTG': 28.1,  # V
                                 'TAA': 1.6, 'TAG': 1.6, 'TGA':1.6 #STOP
                                }


        for key in list(self.strGeneCopy_fast.keys()):
            if 'T' in key:
                val = self.strGeneCopy_fast[key]
                newkey = key.replace('T','U')
                self.strGeneCopy_fast[newkey] = val

        self.strGeneCopy_slow = {'GCT': 7.4, 'GCC': 7.4, 'GCA': 7.4, 'GCG': 7.4,  #A
                                 'CGT': 4.5, 'CGC': 4.5, 'CGA': 4.5, 'CGG': 4.5,
                                 'AGA':4.5, 'AGG':4.5,   #R
                                 'AAT': 17.0, 'AAC':17.0,  #%N
                                 'GAT': 21.8, 'GAC': 21.8,  #D
                                 'TGT': 10.6, 'TGC':10.6,  #C
                                 'CAA': 12.3, 'CAG': 12.3,  #Q
                                 'GAA': 29.0, 'GAG': 29.0,  #E
                                 'GGT': 10.8, 'GGC': 10.8, 'GGA': 10.8, 'GGG': 10.8,  #G
                                 'CAT': 10.9, 'CAC':10.9,  #H
                                 'ATT': 7.5, 'ATC': 7.5, 'ATA': 7.5, #I
                                 'TTA': 7.2, 'TTG':7.2, 'CTT': 7.2, 'CTC': 7.2,
                                 'CTA': 7.2, 'CTG': 7.2, #L
                                 'AAA': 24.4, 'AAG': 24.4, #K
                                 'ATG': 22.0, #M
                                 'TTT': 17.6, 'TTC': 17.6, #F
                                 'CCT': 6.9, 'CCC': 6.9, 'CCA': 6.9, 'CCG': 6.9, #P
                                 'TCT': 4.4, 'TCC': 4.4, 'TCA': 4.4, 'TCG': 4.4,
                                 'AGT': 4.4, 'AGC': 4.4, #S
                                 'ACT': 6.1, 'ACC': 6.1, 'ACA': 6.1, 'ACG': 6.1,#T
                                 'TGG': 13.2, #W
                                 'TAT': 12.2, 'TAC': 12.2, #Y
                                 'GTT': 7.1, 'GTC':7.1, 'GTA': 7.1, 'GTG': 7.1, # V
                                 'TAA': 0.8, 'TAG': 0.8, 'TGA': 0.8 #STOP CODON}
                                }
        
        for key in list(self.strGeneCopy_slow.keys()):
            if 'T' in key:
                val = self.strGeneCopy_slow[key]
                newkey = key.replace('T','U')
                self.strGeneCopy_slow[newkey] = val


        self.fast_codons_value = [27.7, 12.2, 19.1, 25.1, 12.6, 34.2, 39.6, 22.2, 15.1,
                                  20.8, 39.6, 31.9, 22, 20.3, 19.8, 19.5,
                                  18.9, 13.2, 15.3, 28.1, 1.6]


        self.slow_codons_value = [7.4, 4.5, 17, 21.8, 10.6, 12.3, 29, 10.8, 10.9, 7.5, 7.2,
                                  24.4, 22, 17.6, 6.9, 4.4, 6.1, 13.2, 12.2, 7.1, .8]

        fullcodonkeys = ['GCT', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGT', 'CAT',
                     'ATT', 'TTA', 'AAA', 'ATG', 'TTT', 'CCT', 'TCT',
                     'ACT', 'TGG', 'TAT', 'GTT', 'TAA',
                     'GCU', 'CGU', 'AAU', 'GAU', 'UGU', 'CAA', 'GAA', 'GGU', 'CAU',
                     'AUU', 'UUA', 'AAA', 'AUG', 'UUU', 'CCU', 'TCU',
                     'ACU', 'UGG', 'UAU', 'GUU', 'UAA',                     ]

        codonkeys = ['GCT', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGT', 'CAT',
                     'ATT', 'TTA', 'AAA', 'ATG', 'TTT', 'CCT', 'TCT',
                     'ACT', 'TGG', 'TAT', 'GTT', 'TAA']
        self.sensitivity_fast_slow = []
        for i in range(len(codonkeys)):
            self.sensitivity_fast_slow.append(self.strGeneCopy_fast[codonkeys[i]] / self.strGeneCopy_slow[codonkeys[i]])

        self.load_tags()



    def load_tags(self):
        try:
            f= open("custom_tags.txt","r")
        except:
            return
         
        raw = f.readlines()
        previous_tags = []
    
        for line in raw:
            if line != '\n':
                previous_tags.append(line)
                
        for line in previous_tags:
          
            custom_tag = line.strip('\n').split('---')
   
            if custom_tag[0] not in self.tag_dict.keys():
                self.tag_dict[custom_tag[0]] = custom_tag[2]
              
                self.tag_full[custom_tag[0]] = custom_tag[1]
        f.close()
                

    def add_tag(self,nt_seq,name):
        '''
        add a custom tag sequence
        '''

        f= open("custom_tags.txt","r")
         
        raw = f.readlines()
        previous_tags = []
    
        for line in raw:
            if line != '\n':
                previous_tags.append(line)
            
        if not set(nt_seq.lower()).issubset(  set(['a','t','c','g','u'])):
            print('invalid NT sequence')
            f.close()
            return

        
        aa_seq = ''
        for i in range(0, len(nt_seq), 3):
            aa_seq += self.aa_table[nt_seq[i:i+3]]
        

        newtag =name+'---'+ nt_seq.lower() + '---'+ aa_seq.upper()+'\n'  
        
        if newtag not in previous_tags:
            previous_tags.append(newtag)
        f.close()
        
        f= open("custom_tags.txt","w+")
        
        for item in previous_tags:
            f.write('%s' % item)
            
        f.close()    
    
    
class FragmentSeperator():
    '''
    Class to manage ribosomal movement and kymograph seperation
    
    TODO: Comment and explain this code, its very convoluted and unclear how it works
    '''
    def __init__(self):
        pass
    
    def get_fragments(self,position_tensor, total_length = None):
        '''
        attempts to get the individual trajectories back from ribosomal posistion,
        this is by no means a perfect process
        '''
        if total_length == None:
            total_length = np.max(position_tensor)
            
        fragmented_trajectories = []
        fragtimes = []
        maxlen = 0
    
        fragmentspertraj= []
        solutions = [position_tensor]
        nsteps = position_tensor.shape[0]
        k = 0

        ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
        changes = ind[1:] - ind[:-1]
        addindexes = np.where(changes > 0)[0]
        subindexes = np.where(changes < 0)[0]
        
        sub = solutions[k][:,1:] - solutions[k][:,:-1]
        neutralindexes = np.unique(np.where(sub < 0)[1])
        neutralindexes = np.setxor1d(neutralindexes, subindexes)
        
        for index in neutralindexes:
            pre = solutions[k][:,index]
            post = solutions[k][:,index+1]
            changecount = 0
            while len(np.where(post - pre < 0)[0]) > 0:

                post = np.append([total_length],post)
                pre = np.append(pre,0)
                
                changecount+=1
            
            for i in range(changecount):
                addindexes = np.sort(np.append(addindexes,index))
                subindexes = np.sort(np.append(subindexes,index))
                
            changes[index] = -changecount
            ind[index] += changecount
         
            
        for index in np.where(np.abs(changes)>1)[0]:
            if changes[index] < 0:
                for i in range(np.abs(changes[index])-1):
                    subindexes = np.sort(np.append(subindexes,index))
            else:
                for i in range(np.abs(changes[index])-1):
                    addindexes = np.sort(np.append(addindexes,index))   
            
        truefrags = len(subindexes)
 
            
    
       
        if len(subindexes) < len(addindexes):
            subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(nsteps-1)).astype(int))
            
        
        fragmentspertraj.append(len(subindexes))
        
        for m in range(min(len(subindexes),len(addindexes))):
            traj = solutions[k][:, addindexes[m]:subindexes[m]+1]
            traj_ind = changes[addindexes[m]:subindexes[m]+1]
            
            startind = ind[addindexes[m]]
            minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
            fragment = np.array([])
        
                
            
            iterind = startind
            
            if subindexes[m]-addindexes[m] > 0:
                if len(minusloc) > 1:
                    if m <= truefrags:
                        for n in range(len(minusloc)-1):
                            iterind = iterind + min(0,traj_ind[minusloc[n]])
                            fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                            
                            
                            
              
            
                  
                        
                        fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                        
                    else:
                        for n in range(len(minusloc)-1):

                            iterind = iterind + min(0,traj_ind[minusloc[n]])
                            
                            fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
              
                            
                        fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
      
                    
                
                else:

                    fragment = solutions[k][startind][addindexes[m]:subindexes[m]+1].flatten()
               
            
                
                fragtimes.append(addindexes[m]+1)
                   
                
                fragmented_trajectories.append(fragment)
                #if m <= truefrags:
                    #kes.append(genelength/truetime[len(fragment)])
        
                if len(fragment) > maxlen:
                    maxlen = len(fragment)
                

        fragarray = np.zeros((len(fragmented_trajectories), maxlen))
        for i in range(len(fragmented_trajectories)):
            fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
        
        return fragtimes,fragarray
            
            
        
    def get_frags_ssa_obj(self,ssa_obj, total_length = None):
        
        fragmented_trajectories = []
        fragtimes = []
        maxlen = 0
    
        fragmentspertraj= []
        n_traj = ssa_obj.n_traj
        
        solutions = ssa_obj.solutions
        if total_length == None:
            total_length = 0
            for k in range(n_traj):
                if total_length < np.max(solutions[k]):
                    total_length = np.max(solutions[k])
                
        nsteps = solutions[0].shape[1]
        for k in range(n_traj):
            ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
            changes = ind[1:] - ind[:-1]
            addindexes = np.where(changes > 0)[0]
            subindexes = np.where(changes < 0)[0]
            
            sub = solutions[k][:,1:] - solutions[k][:,:-1]
            neutralindexes = np.unique(np.where(sub < 0)[1])
            neutralindexes = np.setxor1d(neutralindexes, subindexes)
            
            for index in neutralindexes:
                pre = solutions[k][:,index]
                post = solutions[k][:,index+1]
                changecount = 0
                while len(np.where(post - pre < 0)[0]) > 0:
    
                    post = np.append([total_length],post)
                    pre = np.append(pre,0)
                    
                    changecount+=1
                
                for i in range(changecount):
                    addindexes = np.sort(np.append(addindexes,index))
                    subindexes = np.sort(np.append(subindexes,index))
                    
                changes[index] = -changecount
                ind[index] += changecount
             
                
            for index in np.where(np.abs(changes)>1)[0]:
                if changes[index] < 0:
                    for i in range(np.abs(changes[index])-1):
                        subindexes = np.sort(np.append(subindexes,index))
                else:
                    for i in range(np.abs(changes[index])-1):
                        addindexes = np.sort(np.append(addindexes,index))   
                
            truefrags = len(subindexes)
     
                
        
           
            if len(subindexes) < len(addindexes):
                subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(nsteps-1)).astype(int))
                
            
            fragmentspertraj.append(len(subindexes))
            
            for m in range(min(len(subindexes),len(addindexes))):
                traj = solutions[k][:, addindexes[m]:subindexes[m]+1]
                traj_ind = changes[addindexes[m]:subindexes[m]+1]
                
                startind = ind[addindexes[m]]
                minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
                fragment = np.array([])
            
                    
                
                iterind = startind
                
                if subindexes[m]-addindexes[m] > 0:
                    if len(minusloc) > 1:
                        if m <= truefrags:
                            for n in range(len(minusloc)-1):
                                iterind = iterind + min(0,traj_ind[minusloc[n]])
                                fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                                
                                
                                
                  
                
                      
                            
                            fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                            
                        else:
                            for n in range(len(minusloc)-1):

                                iterind = iterind + min(0,traj_ind[minusloc[n]])
                                
                                fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                  
                                
                            fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
          
                        
                    
                    else:

                        fragment = solutions[k][startind][addindexes[m]:subindexes[m]+1].flatten()
                   
                
                    
                    fragtimes.append(addindexes[m]+1)
                       
                    
                    fragmented_trajectories.append(fragment)
                    #if m <= truefrags:
                        #kes.append(genelength/truetime[len(fragment)])
            
                    if len(fragment) > maxlen:
                        maxlen = len(fragment)
                    
    
            fragarray = np.zeros((len(fragmented_trajectories), maxlen))
            for i in range(len(fragmented_trajectories)):
                fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
            
       
        ssa_obj.fragments = fragarray
        ssa_obj.fragtimes = fragtimes
        ssa_obj.frag_per_traj = fragmentspertraj
        ssa_obj.full_frags = truefrags
        
    
    
    
class IntensityAnalysesRagged():
    def __init__(self):
        pass
    
    @staticmethod
    def get_acc2(data, trunc=False):
        '''
        Get autocorrelation function

        *NOT* multi-tau
        '''
        N = len(data)
        fvi = np.fft.fft(data, n=2*N)
        acf = fvi*np.conjugate(fvi)
        acf = np.fft.ifft(acf)
        acf = np.real(acf[:N])/float(N)
        if trunc:
            acf[acf < 0]=0
            for i in range(1, len(acf)):
                if acf[i] > acf[i-1]:
                    acf[i] = acf[i-1]
        return acf       
    
    def get_g0(self,covariance, mode = 'interp'):
        '''
        
        '''
        if mode.lower() in ['interp','inter','extrapolate','interpolate']:
            X = [1,2,3,4]
            G0 = np.zeros((covariance.shape[0], covariance.shape[2]))
            for i in range(covariance.shape[0]):
                for n in range(covariance.shape[0]):
                    V = covariance[i,X,n]
                    print(V.shape)
                    G0[i,n] = np.interp(0,X,V)
           
            
        if mode.lower() in ['g1','1']:
            G0 = covariance[:,1,:]
            
        if mode.lower() in ['g0','0']:
            G0 = covariance[:,0,:]

        if mode.lower() in ['max','maximum']:
            G0 = np.max(covariance,axis=1)
            
        return G0
    
    def get_autocorr(self,autocov,g0='G0'):
        '''
        normalize the autocovariance over g0
        '''
        autocorr = np.copy(autocov)
        n_traj = autocorr.shape[-1]
        autocov_error =  1.0/np.sqrt(n_traj)*np.std(autocov,ddof=1,axis=2)
        
        g0 = self.get_g0(autocov,g0)
        for n in range(autocov.shape[0]):
            autocorr[n] = autocorr[n]/g0[n]
        
        err_autocorr =  1.0/np.sqrt(n_traj)*np.std(autocorr,ddof=1,axis=2)
        return autocorr,err_autocorr
                

    def get_autocov(self,intensity_vec,max_lag,norm='raw'):
        
        colors, n_traj = self.__get_ivec_dims(intensity_vec)
        
        autocorr_vec = np.zeros((colors, max_lag, n_traj    ) )
        autocorr_err = np.zeros((colors, max_lag, n_traj  ) )
        

        
        
        for n in range(colors):
            if norm in [ 'Individual','I','individual','ind']:
                for i in range(n_traj):
                    ivec = intensity_vec[i][n]
                    print(ivec.shape)
                    trace_len = min(max_lag, len(ivec) )
                    autocorr_vec[n,:trace_len,i] = self.get_acc2( (ivec-np.mean(ivec))/np.var(ivec)  )[:trace_len]
                    
            elif norm in ['global','Global','g','G']:
                means = []
                var = []
                for i in range(n_traj):
                    ivec = intensity_vec[i][n]
                    means.append( np.mean(ivec) )
                    var.append(np.var(ivec))         
                    
                global_mean = np.mean(means)
                global_var = np.var(var)
                
                for i in range(n_traj):
                    ivec = intensity_vec[i][n]
                    trace_len = min(max_lag, len(ivec) )
                    autocorr_vec[n,:trace_len,i] = self.get_acc2((ivec-global_mean)/global_var )[:trace_len]
            elif norm in ['raw','Raw']:
                for i in range(n_traj):
                    ivec = intensity_vec[i][n]
                    trace_len = min(max_lag, len(ivec) )
                    autocorr_vec[n,:trace_len,i] = self.get_acc2(ivec)[:trace_len]  
                    
            # elif norm in ['max_ind','Max','maximum','Maximum','max']:
            #     for i in range(n_traj):
            #         ivec = intensity_vec[i][n]
            #         trace_len = min(max_lag, len(ivec) )
            #         autocorr_vec[n,:trace_len,i] = self.get_acc2(ivec)[:trace_len]
            #         autocorr_vec[n,:trace_len,i] = autocorr_vec[n,:trace_len,i]/np.max(autocorr_vec[n,:trace_len,i])
                    
            else:
                print('unrecognized normalization, please use individual, global, or none')
                return

        autocorr_err =  1.0/np.sqrt(n_traj)*np.std(autocorr_vec,ddof=1,axis=2)
                
        return autocorr_vec, autocorr_err
    
    
    def get_intensity_histogram(self,intensity_vec,n_bins,scale_factor,time_slice=1):
        colors, n_traj = self.__get_ivec_dims(intensity_vec)
        hist_list = []
        hist_edges_list = []
        for n in range(colors):
            for i in range(n_traj):
                if i == 0:
                    ivec = intensity_vec[i][n][::time_slice]
                    
                else:
                    tmp_vec = intensity_vec[i][n][::time_slice]
                    
                    ivec = np.hstack((ivec, tmp_vec))
                    
                    
            exp_int_dist_ump = np.divide(ivec,scale_factor)
            exp_int_dist_ump[exp_int_dist_ump<0] = 0
            exp_int_dist_ump[np.isnan(exp_int_dist_ump)] = 0
            
            temp_hist = np.histogram(exp_int_dist_ump, bins=n_bins)
            hist_list.append(temp_hist[1])
            hist_edges_list.append(temp_hist[0])
            
        hist_data = np.array(hist_list)
        hist_bins = np.array(hist_edges_list)   
        
        return hist_bins, hist_data
            
    
    def __get_ivec_dims(self,ivec):
        colors = ivec[0].shape[0]
        n_traj = len(ivec)           
        return colors,n_traj
     
        
    
class IntensityAnalyses():
    def __init__(self):
        pass
    
    @staticmethod
    def get_acc2(data, trunc=False):
        '''
        Get autocorrelation function

        *NOT* multi-tau
        '''
        N = len(data)
        fvi = np.fft.fft(data, n=2*N)
        acf = fvi*np.conjugate(fvi)
        acf = np.fft.ifft(acf)
        acf = np.real(acf[:N])/float(N)
        if trunc:
            acf[acf < 0]=0
            for i in range(1, len(acf)):
                if acf[i] > acf[i-1]:
                    acf[i] = acf[i-1]
        return acf    
        
    
    def get_g0(self,correlation, mode = 'interp'):
        '''
        
        '''
        if mode.lower() in ['interp','inter','extrapolate','interpolate']:
            X = [1,2,3,4]
            V = correlation[:,X,:]
            G0 = np.interp(0,X,V)      
            
        if mode.lower() in ['g1','1']:
            G0 = correlation[:,1,:]
            
        if mode.lower() in ['g0','0']:
            G0 = correlation[:,0,:]
            
        if mode.lower() in ['max','maximum']:
            G0 = np.max(correlation,axis=1)
        return G0
    
    def normalize_cc(self, correlation,mode='max'):
        '''
        normalize cc via either center or maximum. 
        '''
        
        if mode.lower() in ['max','maximum']:
            norm_cor = correlation/np.max(correlation,1)
            
        if mode.lower() in ['center','middle']:
            centerpoint = int((correlation.shape[1]+1)/2)-1
            norm_cor = correlation/(correlation[:,centerpoint])
            
        return norm_cor 
        
    
    def get_crosscorr(self,intensity_vecs):
        
        ncolors = intensity_vecs.shape[0]
        time_pts = intensity_vecs.shape[1]
        traj = intensity_vecs.shape[2]
        cross_corr = np.zeros( (ncolors**2,time_pts*2-1,traj ))
        
        i = 0
        k = 0
        inds = []
        for n in range(intensity_vecs.shape[0]):
            for m in range(intensity_vecs.shape[0]):
                iv1 = intensity_vecs[n].T
                iv2 = intensity_vecs[m].T
                inds.append((n,m))
                
                # slen = np.correlate(iv1[i]-np.mean(iv1[i]),iv2[i]-np.mean(iv2[i]),'full').shape[0]
                # crosscorr_vec = np.zeros((iv1.shape[0],slen))
                
                for i in range(traj):
                    cross_corr[k,:,i] = np.correlate(iv1[i,:]-np.mean(iv1[i,:]),iv2[i,:]-np.mean(iv2[i,:]),'full')/time_pts
                k +=1
                
        return cross_corr  ,inds     

    def get_crosscorr2(self, iv1,iv2):
        '''
        returns the autocorrelations
        '''

        i = 0
        slen = np.correlate(iv1[i]-np.mean(iv1[i]),iv2[i]-np.mean(iv2[i]),'full').shape[0]
        crosscorr_vec = np.zeros((iv1.shape[0],slen))
        
        for i in range(iv1.shape[0]):
            crosscorr_vec[i,:] = np.correlate(iv1[i]-np.mean(iv1[i]),iv2[i]-np.mean(iv2[i]),'full')/len(iv1)

        normalized_autocorr = crosscorr_vec.T/ crosscorr_vec[:,len(iv1[i])-1]
        mean_autocorr = np.mean(normalized_autocorr, axis=1)

        return crosscorr_vec, mean_autocorr
        

    def get_autocorr_norm(self, intensity_vec, time_vec, totalSimulationTime, geneLength,normalization= 'Individual'):
        '''
        returns the autocorrelations
        '''

        autocorr_vec = np.zeros((intensity_vec.shape))


        if normalization in [ 'Individual','I','individual','ind']:
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-np.mean(intensity_vec[i]))
        elif normalization in ['global','Global','g','G']:
            global_mean = np.mean(intensity_vec)
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-global_mean)
            
        else:   
            print('unrecognized normalization, using indivdual means')
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-np.mean(intensity_vec[i]))
                
                


        normalized_autocorr = autocorr_vec.T/ autocorr_vec[:,0]
        mean_autocorr = np.mean(normalized_autocorr, axis=1)
        
        error_autocorr = np.std(normalized_autocorr, axis=1)/np.sqrt(intensity_vec.shape[0])
 
        dwelltime = None

        try:
            dwelltime = time_vec[np.where(mean_autocorr < .01)[0][0]]
            
        except:
            try:
                dwelltime = time_vec[np.where(mean_autocorr < .05)[0][0]]
                
            except:
                dwelltime = 1

        
        try:
            zeroind = np.where(mean_autocorr<0)[0][0]
            length = int(.3*len(mean_autocorr))
            zeromean = np.mean(mean_autocorr[-length:])
            zeromean2 = np.mean(mean_autocorr[zeroind:])
    
            normalized_autocorr = normalized_autocorr-zeromean2
            mean_autocorr = np.mean(normalized_autocorr, axis=1)
            
            error_autocorr = np.std(normalized_autocorr, axis=1)/np.sqrt(intensity_vec.shape[0])
     
        except:
            pass
        


        ke_exp = np.round(geneLength/dwelltime ,1)

        return normalized_autocorr, mean_autocorr, error_autocorr, dwelltime, ke_exp
    
    
    
    def get_autocov(self,intensity_vec,norm='raw'):
        autocorr_vec = np.zeros((intensity_vec.shape))
        autocorr_err = np.zeros((intensity_vec.shape))
        colors = intensity_vec.shape[0]
        n_traj = intensity_vec.shape[2]
        
        for n in range(colors):
            if norm in [ 'Individual','I','individual','ind']:
                for i in range(intensity_vec.shape[2]):
                    ivec = intensity_vec[n,:,i]
                    autocorr_vec[n,:,i] = self.get_acc2( (ivec-np.mean(ivec))/np.var(ivec)  )
                    
            elif norm in ['global','Global','g','G']:
                global_mean = np.mean(intensity_vec[n])
                global_var = np.var(intensity_vec[n])
                for i in range(intensity_vec.shape[2]):
                    autocorr_vec[n,:,i] = self.get_acc2((intensity_vec[n,:,i]-global_mean)/global_var )
            elif norm in ['raw','Raw']:
                for i in range(intensity_vec.shape[2]):
                    autocorr_vec[n,:,i] = self.get_acc2(intensity_vec[n,:,i])     
            else:
                print('unrecognized normalization, please use individual, global, or none')
                return

        autocorr_err =  1.0/np.sqrt(n_traj)*np.std(autocorr_vec,ddof=1,axis=2)
                
        return autocorr_vec, autocorr_err
    
    
    def get_autocorr(self,autocov,g0='G0'):
        '''
        normalize the autocovariance over g0
        '''
        autocorr = np.copy(autocov)
        n_traj = autocorr.shape[-1]
        autocov_error =  1.0/np.sqrt(n_traj)*np.std(autocov,ddof=1,axis=2)
        
        g0 = self.get_g0(autocov,g0)
        for n in range(autocov.shape[0]):
            autocorr[n] = autocorr[n]/g0[n]
        
        err_autocorr =  1.0/np.sqrt(n_traj)*np.std(autocorr,ddof=1,axis=2)
        return autocorr,err_autocorr
                
                
                
    def get_autocorr2(self, intensity_vec, time_vec, totalSimulationTime, geneLength, normalization='Individual'):
        '''
        returns the autocorrelations
        '''

        autocorr_vec = np.zeros((intensity_vec.shape))
        
        
        if normalization in [ 'Individual','I','individual','ind']:
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-np.mean(intensity_vec[i]))
        elif normalization in ['global','Global','g','G']:
            global_mean = np.mean(intensity_vec)
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-global_mean)
            
        else:   
            print('unrecognized normalization, using indivdual means')
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-np.mean(intensity_vec[i]))
        autocorr = autocorr_vec.T
        mean_autocorr = np.mean(autocorr, axis=1)
        error_autocorr = np.std(autocorr, axis=1)/np.sqrt(intensity_vec.shape[0])
        
 
        dwelltime = None

        try:
            dwelltime = time_vec[np.where(mean_autocorr < .01)[0][0]]

        except:
            try:
                dwelltime = time_vec[np.where(mean_autocorr < .05)[0][0]]
            except:
                dwelltime = 1




        ke_exp = np.round(geneLength/dwelltime ,1)

        return autocorr, mean_autocorr, error_autocorr, dwelltime, ke_exp  



class ModelBuilder():
    '''
    Container class for the solvers
    '''
    def __init__(self, time=None,xi=None):
        self.k = None
        self.k_bind = None
        self.k_term = None
        self.multiframe = True
        self.additional_rxns = {}
        self.probe_locations = None
        self.colors = 1
        
        self.default_conditions = {'low_mem':True,
                                   'perturb':[0,0,0],
                                   'leaky_probes':False,
                                   'bins':None,
                                   'k_probe':0,
                                   'footprint':9,
                                   'burnin':500,
                                   'record_stats':False}
        
        self._poi = None
        self.k_enters = np.array([])
        self.k_stops = np.array([])
        self.k_jumps = np.array([])
        self.k_pauses = np.array([])
        self.frames_used = 1
        self.probe = None
        
        
    @property    
    def n_extra_rxns(self):
        return len(self.k_enters) + len(self.k_stops) + len(self.k_jumps) + len(self.k_pauses) 
        
    def add_k(self, poi, kelong_mean = 10):
        self._poi = poi
        ks = PropensityFactory().get_k_3_frame(poi.nt_seq, kelong_mean)        
        self.k = ks

        L = len(poi.aa_seq)
        self.k[L-1] = 0
        self.k[L*2-1  -1] = 0
        self.k[L*3-2 -1] = 0
        self.k_stops = np.array([[L,0,10],[L-1,1,10],[L-1,2,10]],dtype=np.float64)
        
    
    def add_jumps(self, locs, frames, dest_locs, dest_frames, rates):
        if isinstance(locs,int):
            locs = [locs]
        if isinstance(frames,int):
            frames = [frames]
        if isinstance(dest_locs,int):
            dest_locs = [dest_locs]
        if isinstance(dest_frames,int):
            dest_frames = [dest_frames]
        if isinstance(rates,float):
            rates = [rates]
        self.__check_frames(frames)
        self.k_jumps = np.array((locs,frames,dest_locs,dest_frames,rates )).T.astype(np.float64)
        
    
    def add_enters(self,locs, frames, rates):
        if isinstance(locs,int):
            locs = [locs]
        if isinstance(frames,int):
            frames = [frames]
        if isinstance(rates,float):
            rates = [rates]        
        self.__check_frames(frames)
        self.k_enters = np.array((locs,frames,rates )).T.astype(np.float64)
    
    def add_stops(self,locs, frames, rates):
        if isinstance(locs,int):
            locs = [locs]
        if isinstance(frames,int):
            frames = [frames]
        if not isinstance(rates,list):
            rates = [rates]        
        self.__check_frames(frames)

        self.k_stops = np.array((locs,frames,rates )).T.astype(np.float64)
        
    
    def clear_model(self):
        self._poi = None
        self.k_enters = np.array([])
        self.k_stops = np.array([])
        self.k_jumps = np.array([])
        self.k_pauses = np.array([])
        self.frames_used = 1        
        self.probe = None
    
    def add_pauses(self,locs, frames, rates):
        if isinstance(locs,int):
            locs = [locs]
        if isinstance(frames,int):
            frames = [frames]
        if isinstance(rates,float):
            rates = [rates]        
            
        self.__check_frames(frames)
        self.k_pauses = np.array((locs,frames,rates )).T.astype(np.float64)
    
    def __check_frames(self,frames):
        for frame in frames:
            if frame > self.frames_used:
                self.frames_used = frame
        
    def visualize_transcript(self):
        
        probe = self._poi.probe_loc
        #self._poi.generate_3frame_tags()
        
        keys = []
        for dics in self._poi.multiframe_epitopes:
            for key in dics.keys():
                if key not in keys:
                    keys.append(key)
        
        
        fig,ax = plt.subplots(1,1,dpi=300)
        N = len(self._poi.aa_seq)
        ncolors = len(keys)
        xs = [.1, 2.1, 4.1][::-1]
        fr = ["+0","+1","+2"]
        
        
        opt = {'head_width': .1, 'head_length': .1,
        'length_includes_head': True}
        
        
        cmap = cm.get_cmap('viridis')
        colors = cmap(np.linspace(.01,.95, ncolors))
        
        go_color = '#2A9D8F'
        stop_color = '#E76F51'
        jump_color = '#843b62'
        pause_color = '#e9c46a'
        
        
        
        if ncolors <4:
            colors = ['#64ff21','#08ffff','#fb00ff']
        
        for i in range(3):
            

            
            rectangle =  mpatches.Rectangle((0,xs[i]), N ,.7,linewidth=1,edgecolor='k',facecolor='darkgray')
    
            ax.add_patch(rectangle)
            
            color = []
            color_dict = self._poi.multiframe_epitopes[i]
            location = []
            for n in range(len(keys)):
                if keys[n] in color_dict.keys():
                
                    location.append(color_dict[keys[n]])
                    color.append( [n for x in color_dict[keys[n]]])
                
          
            
            colorlabels = ['Color %d'% j for j in range(ncolors)    ]
            
        
            color = [j for i in color for j in i]
            location = [j for i in location for j in i]
            
            for c in range(ncolors):
                ax.plot([-10,-10],[-4,-4],color = colors[c]  )  #fix the legend colors
            
            
            for c,loc in zip(color,location):
                ax.plot([loc,loc],[xs[i] ,xs[i]+.8],color = 'k',lw=2 )
                ax.plot([loc,loc],[xs[i] ,xs[i]+.8],color = colors[c]  )
                
            ax.set_ylim([-.1,8])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.set_xlabel('codon')
            
            ax.axes.get_yaxis().set_visible(False)
            if i == 0:
                leg1 = ax.axes.legend(colorlabels,loc=9)
                ax.text(-20,7.5,'Transcript Name: %s' % self._poi.name)
                ax.text(-20,6.7,'Total Length: %d aa' % self._poi.total_length)
                ax.text(-20,5.9,'Seq: %s ...' % self._poi.aa_seq[:10])
                

                
            ax.text(-N*.07 ,xs[i] + .4,fr[i])
                
        for j in range(self.k_enters.shape[0]):
            yloc = xs[int(self.k_enters[j][1])]
            loc = int(self.k_enters[j][0])
            ax.plot([loc,loc],[yloc, yloc+.8],color = 'k',lw=3  )
            ax.plot([loc,loc],[yloc, yloc+.8],color = go_color,lw=2  )

            ax.annotate("",
                    xy=(loc, yloc+.85), xycoords='data',
                    xytext=(loc-1, yloc+.85), textcoords='data',
                    arrowprops=dict(arrowstyle="-|>",lw=1.5,facecolor='k',edgecolor='k',
                                    connectionstyle="arc3"),
                    )            
            
            ax.annotate("",
                    xy=(loc, yloc+.85), xycoords='data',
                    xytext=(loc-1, yloc+.85), textcoords='data',
                    arrowprops=dict(arrowstyle="-|>",lw=1,facecolor=go_color,edgecolor=go_color,
                                    connectionstyle="arc3"),
                    )
            

        for j in range(self.k_stops.shape[0]):
            yloc = xs[int(self.k_stops[j][1])]
            loc = int(self.k_stops[j][0])
            ax.plot([loc,loc],[yloc, yloc+.8],color = 'k',lw=3  )
            ax.plot([loc,loc],[yloc, yloc+.8],color = stop_color,lw=2  )
            
        for j in range(self.k_jumps.shape[0]):
            
            yloc1 = xs[int(self.k_jumps[j][1])]
            loc1 = int(self.k_jumps[j][0])
            yloc2 = xs[int(self.k_jumps[j][3])]
            loc2 = int(self.k_jumps[j][2])   

            ax.plot([loc1,loc1],[yloc1, yloc1+.8],color = 'k',lw=3  )
            ax.plot([loc2,loc2],[yloc2, yloc2+.8],color = 'k',lw=3  )
            
            ax.plot([loc1,loc1],[yloc1, yloc1+.8],color = jump_color,lw=2  )
            ax.plot([loc2,loc2],[yloc2, yloc2+.8],color = jump_color,lw=2  )
            
            if yloc2 == yloc1:
                ax.plot([loc1,loc1],[yloc1, yloc1-.3],color='k')
                ax.plot([loc1,loc2],[yloc1-.3, yloc1-.3],color='k')
                
                ax.annotate("",
                        xy=(loc2, yloc2), xycoords='data',
                        xytext=(loc2, yloc1-.3), textcoords='data',
                        arrowprops=dict(arrowstyle="-|>",lw=2,
                                        connectionstyle="arc3"),
                        )
                
            else:
                dx = loc2 - loc1
            
                dy = yloc2- yloc1 
                
                if yloc1 > yloc2:
                    dx = loc2 - loc1
            
                    dy = yloc2+.8- yloc1 
                
                   # ax.arrow(loc1,yloc1,dx,dy,color = 'k',lw=2)#,  head_width=5, head_length=.2)
                    ax.annotate("",
                        xy=(loc2, yloc2+.7), xycoords='data',
                        xytext=(loc1, yloc1), textcoords='data',
                        arrowprops=dict(arrowstyle="-|>",lw=2,
                                        connectionstyle="arc3"),
                        )
                else:
                    dx = loc2 - loc1
            
                    dy = yloc2- yloc1+.8
                    #ax.arrow(loc1,yloc1,dx,dy,color = 'k',lw=1)
                    
                    ax.annotate("",
                        xy=(loc2, yloc2), xycoords='data',
                        xytext=(loc1, yloc1+.8), textcoords='data',
                        
                        
                        arrowprops=dict(arrowstyle="-|>",lw=2,
                                        connectionstyle="arc3"),
                        )

        for j in range(self.k_pauses.shape[0]):
            yloc = xs[int(self.k_pauses[j][1])]
            loc = int(self.k_pauses[j][0])
            ax.plot([loc,loc],[yloc, yloc+.8],color = 'k',lw=3  )
            ax.plot([loc,loc],[yloc, yloc+.8],color = pause_color,lw=2  )         
            
            
            
            
        custom_lines = [Line2D([0], [0], color=go_color, lw=2),
                Line2D([0], [0], color=stop_color, lw=2),
                Line2D([0], [0], color=pause_color, lw=2), 
                Line2D([0], [0], color=jump_color, lw=2)]
        
        ax.legend(custom_lines, ['Enter', 'Stop', 'Pause','Jump'])
        ax.add_artist(leg1)
        #ax.set_facecolor('darkgray')
        fig.show()    
        
        

    def __generate_additional_ks(self,enters,pauses,jumps,stops,L):
        
        def frame_check_1(L,arr):        
            return (L- arr[:,1]+1)*(arr[:,1]>0) + L*(arr[:,1]>1)
        
        def frame_check_3(L,arr):        
            return (L- arr[:,3]+1)*(arr[:,3]>0) + L*(arr[:,3]>1)            
                    
        def gen_ks_1_loc(L,arr):
            arr[:,0] = arr[:,0]+frame_check_1(L,arr)
            arr[:,1] = arr[:,2]    
            arr = arr[:,0:2]
            max_arr = np.max( arr[:,0])     
            return arr,max_arr
        
        def gen_ks_3_loc(L,arr):
            arr[:,0] = arr[:,0]+ frame_check_1(L,arr)     
            arr[:,1] = arr[:,2]+ frame_check_3(L,arr)
            arr[:,2] = arr[:,4]
            arr = arr[:,0:3]
            max_arr = max([np.max( arr[:,0]),np.max( arr[:,1])])
            return arr,max_arr
    
        max_enter = 0
        max_pause = 0
        max_stop = 0
        max_jump = 0
        k_jumps = np.copy(jumps)
        k_pauses = np.copy(pauses)
        k_stops = np.copy(stops)
        k_enters = np.copy(enters)
        if len(k_enters) != 0:
            k_enters,max_enter = gen_ks_1_loc(L,k_enters)
    
        if len(k_pauses) != 0:
            k_pauses,max_pause = gen_ks_1_loc(L,k_pauses)
    
        if len(k_stops) != 0:
            k_stops,max_stop = gen_ks_1_loc(L,k_stops)
        
        if len(k_jumps) != 0:
            k_jumps,max_jump = gen_ks_3_loc(L,k_jumps)
            
        max_loc = max(max_jump,max_stop,max_pause,max_enter)
        
        if max_loc <=L: 
            frames_used = 0
        if max_loc > L:
            frames_used = 1
        if max_loc > 2*L-1 :
            frames_used = 2
        
        return k_enters, k_pauses, k_stops, k_jumps, frames_used            

    def generate_probe(self):
        keys = []
        for frame in self._poi.multiframe_epitopes:
            for key in frame.keys():
                if key not in keys:
                    keys.append(key)
                    
        ncolors = len(keys)
        
        N = len(''.join(self._poi.multiframe_aa_seq))
        Ns = [0,len(self._poi.multiframe_aa_seq[0]),len(self._poi.multiframe_aa_seq[1])+len(self._poi.multiframe_aa_seq[0])   ]
        probe = np.zeros((ncolors,N) )
        
        for j in range(3):
            dic = self._poi.multiframe_epitopes[j]
            for i in range(len(keys)):
                if keys[i] in dic.keys():
                    epis = dic[keys[i]]
                    
                    probe[i,Ns[j] + np.array(epis)] = 1
        return probe

    def run_ssa(self,t, n_traj=10,low_mem = True ):
        
        if self.probe == None:
            self.probe = self.generate_probe()
        
        probe = self.probe.astype(int).copy(order='C')
        t_array = t
        t0 = t[0]
        ncolors = probe.shape[0]

        N_rib = 200
        #result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
        #kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
        n_trajectories = n_traj
        start = time.time()
        
        lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))
        
        all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
        
        if low_mem == False:
            all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
            all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
            all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
            nribs = np.array([0],dtype=np.int32)
            all_ribs = np.zeros((n_trajectories,1))
            
        else:
            all_results = np.zeros((n_trajectories,ncolors,len(t_array)),dtype=np.int32)
        
        seeds = np.random.randint(0,0x7FFFFFF,n_trajectories, dtype= np.int32)
        
        k_enters,k_pauses,k_stops,k_jumps,frames_used = self.__generate_additional_ks(self.k_enters,self.k_pauses,self.k_jumps,self.k_stops,self._poi.total_length)
        
        #k_add = np.hstack((self.k_enters.flatten(),self.k_pauses.flatten(),self.k_stops.flatten(),self.k_jumps.flatten() ))
        k_add = np.hstack((k_enters.flatten(),k_pauses.flatten(),k_stops.flatten(),k_jumps.flatten() ))
        
        

        max_ribs = 200
        
        

        n_enters = self.k_enters.shape[0]
        n_pauses = self.k_pauses.shape[0]
        n_stops = self.k_stops.shape[0]
        n_jumps = self.k_jumps.shape[0]
        

                
        if low_mem == False:
            for i in range(n_trajectories):
                result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
                frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
                
                ribtimes = np.zeros((400),dtype=np.float64)
                coltimes = np.zeros((400),dtype=np.int32)
                ssa_translation_generic.run_SSA_generic(result,ribtimes,coltimes, self.k,frapresult,t_array, np.array([0,0,0],dtype=np.float64), seeds[i],nribs, k_add.flatten()  ,n_enters,n_pauses,n_stops,n_jumps)
                all_results[i,:] = result
                all_frapresults[i,:] = frapresult
                all_coltimes[i,:] = coltimes
                all_ribtimes[i,:] = ribtimes
                all_ribs[i,:] = nribs[0]
            
        else:
            
            
        
            kelong = np.array(self.k).astype(np.float64)
            
            for i in range(n_trajectories):
                result = np.zeros((ncolors,len(t_array)),dtype=np.int32)    
                frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
                
                

                ssa_translation_generic_lowmem.run_SSA_generic(result,kelong ,frapresult,t_array, np.array([0,0,0],dtype=np.float64), seeds[i], k_add.flatten()  ,n_enters,n_pauses,n_stops,n_jumps, probe, ncolors ,max_ribs )
                all_results[i,:,:] = result
                #all_frapresults[i,:] = frapresult
        

      
        #traj = all_results[0,:].reshape((N_rib,len(t_array))).T
        return all_results
    
    
    

class TranslationOptimization():
    '''
    optimization class adjusted from Micheal May's GenericOpt
    '''        
    def __init__(self):
        
        self.solver_obj = None
        self.data_obj = None
        self.int_a = IntensityAnalyses()
        
        self.parnames = []
        self.params = np.array([])
        self.initial_params = np.array([])
        
        self.pff = PropensityFactory()
        self.pvf = ProbeVectorFactory()
        
        
        self._int_a = IntensityAnalyses()
        self.methods={'met_hast':self.methast,\
                      'methast':self.methast,\
                      'MH':self.methast,\
                      'lin_opt':self.lin_opt,\
                      'linopt':self.lin_opt,\
                      'genetic':self.genetic,\
                      'sim_anneal':self.sim_anneal,\
                      'simanneal':self.sim_anneal}
            
        self.objective_funs={'sse':self.sse,\
                   'SSE':self.sse,\
                   'chi':self.chisq,\
                   'chisq':self.chisq,\
                   'LL_acorr':self.get_loglikelihood_autocorr,\
                   'I_mu_sse':self.get_intensity_sse,\
                   'LL_I_distb':self.get_loglikelihood_intensity_distribution,\
                   'LL_acorr_ode':self.get_loglikelihood_autocorr_ode}
            
        self.opts={'bounds':([0,10],)*len(self.params)}  
        
        self.args = {'LL_acorr': (20,'raw','G0'), 'I_mu_sse':(), 'combined_objective':[],'LL_I_distb':(30,True)}
        
        self._obj_weights = self.args.copy()
        for key in self._obj_weights.keys():
            self._obj_weights[key] = 1
        
     
        
        self.last_run = None
        self.chain = None
        
    def print_args(objfun):
        x=1
        
    
    def add_bounds(self,lower_bounds, upper_bounds):
        '''
        Adds bounds to the optimizer object

        Parameters
        ----------
        lower_bounds : List or ndarray
            lower bound of each parameter.
        upper_bounds : List or ndarray
            upper bound of each parameter.

        Returns
        -------
        None.

        '''
        if isinstance(lower_bounds,np.ndarray):
            lower_bounds = lower_bounds.tolist()
        if isinstance(upper_bounds,np.ndarray):
            upper_bounds = upper_bounds.tolist() 
            
        self.opts['bounds'] = tuple(list(x) for x in zip(lower_bounds,upper_bounds))
        
        
    def intensity_fun(self,x):
        self.solver_obj._poi.ke_mu = x[1]
        self.solver_obj._poi.ki = x[0]
        ssa_soln = self.solver_obj.solve_ssa_set_conditions()
        return ssa_soln.intensity_vec
    
    def autocovariance_fun(self,intensity,norm='ind'):
        acov,err_acov = self._int_a.get_autocov(intensity,norm=norm)        
        return acov,err_acov

    def autocorrelation_fun(self,intensity,norm='ind',g0='G0'):
        acov,err_acov = self._int_a.get_autocov(intensity,norm=norm)
        acorr,err_acorr = self._int_a.get_autocorr(acov,g0=g0)
        return acorr,err_acorr    

    def intensity_distribution(self,intensity,bins=None,density=True,norm=1):
        
       
        int_dist = np.histogram(intensity/norm, bins=bins, density=density,) 
        int_dist_bins = int_dist[1]
        int_dist_heights = int_dist[0]
        return int_dist_heights,int_dist_bins
    
    def analytical_autocorrelation(self,x,bins=None,bin_method='intellegent'):
        
        self.solver_obj._poi.ke_mu = x[1]
        self.solver_obj._poi.ki = x[0]       
        
        if not isinstance(bins,None):
            
            
            inds = self.pff.intellegent_bin(np.atleast_2d(self.solver_obj._poi.probe_loc),100)
            bpl,bpv = self.pvf.bin_probe_vecs(self.solver_obj._poi.probe_loc,inds)
            k_bin = self.pff.bin_k(self.solver_obj._poi.kelong, inds)
            x0 = np.zeros((k_bin.shape[0],1))
            t = self.solver_obj.t
            ode_soln_bin = self.solver_obj.solve_ode(k_bin,   t, x0, self.solver_obj._poi.ki, bpl,corr=True)
        
        
    def genetic(self,objfun,**kwargs):
        
        return sci.optimize.differential_evolution(objfun,self.opts['bounds'],**kwargs)
    
    def sim_anneal(self,objfun,**kwargs):
        
        minimizer_kwargs = {"bounds":self.opts['bounds'] }
        return sci.optimize.basinhopping(objfun,self.initial_params,minimizer_kwargs =minimizer_kwargs,**kwargs)
    
    
    def methast(self, optfun, optfun_type, niter=1000, burnin=100, stepsize=None,mut_rate=.3, disp=False, logspace=True, proposal = None):
        '''

        Parameters
        ----------
        optfun : function
            objective function to optimize.
        optfun_type : str
            the name of the objective function for __.args to access its arguments
        niter : int, optional
            Number of iterations. The default is 1000.
        burnin : int, optional
            Number of iterations to burn (disregard from the final chain). The default is 100.
        stepsize : list, ndarray, optional
            The stepsize of the proposal distribution. The default is None.
        mut_rate : float 0-1, optional
            mutation probability. The default is .3.
        disp : boolean, optional
            return a display of the current progress. The default is False.
        logspace : boolean, optional
            use logspace to search. The default is True.
        proposal : function, optional
            function to draw proposal distributions from. The default is None.

        Returns
        -------
        OptimizeResult
            scipy result object.

        '''
        if logspace:
            bounds = tuple([np.log10(x).tolist() for x in self.opts['bounds']])
            initial_par = np.log10(self.initial_params)
        else:
            bounds = self.opts['bounds']
            initial_par = self.initial_params
        
        if stepsize is None:
            stepsize = (initial_par/10).tolist()
            
        if proposal is None:
            proposal = np.random.normal
            
        
        def evolvepars(p,stepsize):
            new=np.copy(p)
            for j in range(len(new)):
                    if np.random.rand() < mut_rate:
                        while 1:
                            new[j]=p[j]+proposal(0,stepsize[j])
                    #print p
                    #print new
                            if (bounds[j][0]<new[j]<bounds[j][1]):
                                break
            return new
                        
        objective_fun_args = self.args[optfun_type]
        

        if logspace:
            oldpars=initial_par
            f_old=f_best=optfun(10**oldpars, *objective_fun_args)             
            bestpars=initial_par
        else:
            oldpars=initial_par
            f_old=f_best=optfun(oldpars, *objective_fun_args) 
            bestpars=initial_par
            
        if np.isnan(f_old):
            f_old = np.inf
        
        
    
        if disp:
            print('Burning in....')
        
        for i in range(-burnin,niter):
            if disp:
                if i > 0 and i%(niter/10) == 0:
                    self.__mh_print_report(i,bestpars,f_best,niter)
            
            newpars=evolvepars(oldpars,stepsize)

            if logspace:
                tmp_pars = copy.deepcopy(newpars)
                f_new=optfun(10**tmp_pars,*objective_fun_args)
            else:
                f_new=optfun(newpars,*objective_fun_args)
            #print("newpars = "+str(newpars)+ "     oldpars = "+str(oldpars))
            #print("fnew="+str(f_new)+", fold="+str(f_old))
            
            if np.isnan(f_new):
                f_new = np.inf
     
            if f_new<f_old or np.random.rand()<np.exp(f_old-f_new):

                #raw_input("Press Enter to continue...")
                f_old=f_new
                oldpars=newpars
                self.__update_mh_chain(newpars,f_new)
                if f_new<f_best:
                    f_best=f_new
                    bestpars=newpars
                    
        result = sci.optimize.OptimizeResult()
        
        if logspace:
            result.x = 10**bestpars
        else:
            result.x = bestpars
        
        result.fun = f_best
        result.success = True
        result.nit = niter
        result.logspace = logspace
        return result
    
    def __mh_print_report(self,iteration,bestpars,fbest,niter):
        print('current iteration: %d out of %d | best_parameters: %s | best evaulation: %f' % (iteration,niter,''.join(str(bestpars.tolist())),fbest )  )
    
    def __update_mh_chain(self,pars,funeval):
       
        self.chain.accepted_parchain = np.vstack( (self.chain.accepted_parchain, pars) )
        self.chain.accepted_evalchain = np.vstack((self.chain.accepted_evalchain,np.sum(funeval)))
        #self.chain.accepted_objfunchain  = np.vstack((self.chain.accepted_objfunchain,np.atleast_1d(funeval)))
        self.chain.accepted = self.chain.accepted + 1
        

            
    def lin_opt(self):
        return 1
    
    
    def combined_objective(self,x,objfun_list,intensity_fun):
        obj_sum = 0
        intensity = intensity_fun(x)
        obj_fun_evals = np.zeros(len(objfun_list))
        k = 0
        for i in range(len(objfun_list)):
            objective = objfun_list[i]
            obj_args = self.args[objective]
            
            obj_fun_evals[i] = self._obj_weights[objective]*self.objective_funs[objective](intensity,*obj_args)
         
            #obj_sum += self._obj_weights[objective]*self.objective_funs[objective](intensity,*obj_args)  
        obj_sum = np.sum(obj_fun_evals)
        self.__update_chain(x, obj_fun_evals )
        
        return obj_sum
    
    def run_optimization(self, objective_fun_list, method ,model = None, data = None, intensity_fun = None,**kwargs):
        if isinstance(objective_fun_list,str):
            objective_fun_list = [objective_fun_list]
        
        if model == None:
            model = self.solver_obj
        if data == None:
            data = self.data_obj
        if intensity_fun == None:
            intensity_fun = self.intensity_fun
                        
        obj_fun = self.combined_objective
                
        method_fun = self.methods[method]
        
        self.chain = OptChain()
        
        self.chain.parchain = self.initial_params
        self.chain.initial_params = self.initial_params
        self.chain.iterations = 0
        self.chain.parnames = self.parnames
        
        self.chain.bestpar = None
        self.chain.besteval = None
        self.chain.opt_method = method
        

        self.chain.evalchain  = np.array([11110])
        self.chain.objfunchain = np.zeros((1,len(objective_fun_list)))
        
        starttime = time.time()
        
        args = (objective_fun_list,intensity_fun)
        
        if method in ['met_haste','methaste','MH']:
            self.args['combined_objective'] = [objective_fun_list,intensity_fun]
            self.chain.accepted = 0
            self.chain.accepted_evalchain = np.array([11110])
            self.chain.accepted_objfunchain =  np.zeros((1,len(objective_fun_list)))
            self.chain.accepted_parchain = self.initial_params
            result = method_fun(obj_fun,'combined_objective', **kwargs)
            self.chain.accepted_evalchain = self.chain.accepted_evalchain[1:]
            #self.chain.accepted_objfunchain =  self.chain.accepted_objfunchain[1:,:]
            self.chain.accepted_parchain = self.chain.accepted_parchain[1:,:]     
            
        else:
            kwargs['args'] = args
            result =  method_fun(obj_fun,**kwargs)
            
        self.chain.runtime = time.time()-starttime
        
        #result =  method_fun(obj_fun,**kwargs)
        
        self.chain.bestpar = result.x
        self.chain.besteval = result.fun
        
        self.chain.parchain = self.chain.parchain[1:,:]
        self.chain.objfunchain = self.chain.objfunchain[1:,:]
        self.chain.evalchain  = self.chain.evalchain[1:]
        

    def __update_chain(self,pars, funeval):   

  
        self.chain.parchain = np.vstack( (self.chain.parchain, pars) )

        self.chain.evalchain = np.vstack((self.chain.evalchain,np.sum(funeval)))
        self.chain.objfunchain  = np.vstack((self.chain.objfunchain,np.atleast_1d(funeval)))
        self.chain.iterations = self.chain.iterations + 1
        
        
    def __intensity_generator(self,pars,objective_fun_list):
        intensity = self.intensity_fun(pars)
        
        return intensity
    
    def get_loglikelihood_autocorr(self, intensity, n_points,norm,g0):
        
        model_acorr,model_acorr_err = self.autocorrelation_fun(intensity,norm=norm,g0=g0)

        total_n_spots = self.solver_obj.n_traj
        data_autocorrelation = self.data_obj.acorr
        data_acc_err = self.data_obj.acorr_err
        LL = self.loglikelihood_acc(model_acorr[:,:n_points,:],  data_autocorrelation[:,:n_points,:], data_acc_err[:,:n_points], total_n_spots)
    
        
        
        return LL

    def get_loglikelihood_autocorr_ode(self, intensity, n_points,norm,g0):
        
        model_acorr,model_acorr_err = self.autocorrelation_fun(intensity,norm=norm,g0=g0)

        total_n_spots = self.solver_obj.n_traj
        data_autocorrelation = self.data_obj.acorr
        data_acc_err = self.data_obj.acorr_err
        LL = self.loglikelihood_acc(model_acorr[:,:n_points,:],  data_autocorrelation[:,:n_points,:], data_acc_err[:,:n_points], total_n_spots)
    
       
        
        return LL

    def get_intensity_sse(self,intensity):
        
        return (np.mean(intensity) - np.mean(self.data_obj.I_mu))**2
    
    def get_loglikelihood_intensity_distribution(self,intensity,norm):
        
        dist_sim_data = self.intensity_distribution(intensity,bins=self.data_obj.histogram_bins,density=True,norm=norm)[0]
        dist_sim_data[dist_sim_data==0] = 1e-7
        
        LL = -np.dot(self.data_obj.histogram,np.log(dist_sim_data))
        
        if LL == -np.inf:
            LL = np.inf
        return LL
        
    @staticmethod
    def sse(model_data,data):        
        return np.sum((model_data-data)**2)
    
    @staticmethod
    def chisq(model_data,model_var, data):
        return np.sum(((model_data-data)**2)/model_var)
    
    
    @staticmethod
    def loglikelihood_acc( model, data,data_err, nspots):
        '''
        Parameters
        ----------
        model : ndarray
            autocorrelation model generated data array.
        data : ndarray
            autocorrelation data array.
        data_err : ndarray
            error array (SEM, STD) of the data.
        nspots : int
            number of spots.

        Returns
        -------
        Loglikelihood (float)
            returns the logliklihood from the formula:
                
                log L(G|M) = Const. - 1/Nt * sum( (G_D(t_i) - G_M(t_i))^2 / sigma_G_D(t_i) )

        '''
        
        d = np.mean(data,axis=-1)
        m = np.mean(model,axis=-1)
       

    
        return (nspots/2) * np.sum(( d[:,2:] - m[:,2:])**2/ data_err[:,2:])
    
    @staticmethod
    def loglikelihood_distb( model_intensity, data_intensity,nbins=30,norm=1):
        '''
    
        Parameters
        ----------
        model_intensity : ndarray
            array of model simulated intensity over time.
        data_intensity : ndarray
            array of data intensity over time.
        nbins : int, optional
            Number of bins. The default is 20.

        Returns
        -------
        float
            LogLikelihood of the intensity distributions.

        '''
        dist_sim_data = np.histogram(model_intensity/norm, bins=nbins, density=True)[1]
        hist_exp_data = np.histogram(data_intensity, bins=nbins)
        return -np.dot(hist_exp_data,np.log(dist_sim_data))
    
    
class IntensityData():
    def __init__(self):
      
        self.ragged = False
        self.head = None
        self.ssa_obj = None    
        
    
    def add_data(self, t, intensity_vec):
        '''

        Parameters
        ----------
        intensity_vec : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.intensity_vec = intensity_vec
        if isinstance(intensity_vec,np.ndarray):
            self.ragged=False
            
        else:
            self.ragged=True
        self.times = t
        
         
    def get_stats(self):
        if self.ragged == False:
            self.I_mu = np.mean(self.intensity_vec,axis=2)
            self.I_var = np.mean(self.intensity_vec,axis=2)
        else:
            spot_no = len(self.intensity_vec)
            def autocorr(x):
                result = np.correlate(x, x, mode='full')
                return result[result.size // 2:]        
            # Calculating the autocovariance
            n_lags = 300
            n_selected_Particles = np.max(spot_no)
            ac_array = np.zeros((n_selected_Particles , n_lags ))
            norm_ac_array = np.zeros((n_selected_Particles , n_lags ))
            counter = 0
            # Calculating Autocorrelation.
            for i in range(1,n_selected_Particles+1):
                intensity_green = self.intensity_vec[i][0]
                temp_ac = autocorr(intensity_green)
                size_ac_temp  = len(temp_ac)
                ac_array[counter, 0:size_ac_temp] = autocorr(intensity_green)
                norm_ac_array[counter, 0:size_ac_temp] = autocorr(intensity_green)/ float(ac_array[counter, :] .max())
                counter += 1
            # Plotting mean autocovariance
            lag_time = [*range(0, self.max_lag_output, 1)]
            lag_time = [element * self.sampling_rate for element in lag_time]
            mean_ac_data = norm_ac_array.mean(0)
            std_ac_data = norm_ac_array.std(0)
            # normalized autocovariance, removing shot noise
            
            if self.remove_shotnoise ==1:
                #mean_ac_data_norm = mean_ac_data[1:-1]/mean_ac_data[1]
                return lag_time[1:self.max_lag_output], mean_ac_data[1:self.max_lag_output], std_ac_data[1:self.max_lag_output]
            else:
                #mean_ac_data_norm = mean_ac_data[1:-1]/mean_ac_data[0]
                return lag_time[0:self.max_lag_output], mean_ac_data[0:self.max_lag_output], std_ac_data[0:self.max_lag_output]
            
        
    def load_data(self, file):
        extension = file.split('.')[-1]
        if extension in ['xls', 'xlsx']:
            df = pd.read_excel(file)
            
            self.head = df.head()
            n_cells = max(df['Cell_#'])
            max_spots =  max(df['Spot_#'])
            total_n_spots = 0
            data_struct = []
            for i in range(0,n_cells):
                i_g = []
                i_r = []
                intensitys = []
                times = []
                for j in range(0,max_spots):
                    temp_len_spot =  df[(df['Cell_#'] == i) & (df['Spot_#'] ==j) ]
                    if len(temp_len_spot) >0:
                        time = df[(df['Cell_#'] ==i) & (df['Spot_#'] ==j)]['Time_(sec)'].values
                        intensity_red = df[(df['Cell_#'] ==i) & (df['Spot_#'] ==j)]['R_GaussianInt'].values
                        intensity_green = df[(df['Cell_#'] ==i) & (df['Spot_#'] ==j)]['G_GaussianInt'].values
            
                        i_g.append(intensity_green)
                        i_r.append(intensity_red)
                        
                        intensity = np.vstack((intensity_red,intensity_green))
                        intensitys.append(intensity)
                        times.append(time)
                data_struct.append([i_g,i_r,intensitys,times ])
                
            time_arrays = []
            intensity_arrays = []
            cells = []
            for n in range(len(data_struct)):
                cells = cells + [n for x in range( len(data_struct[n][3])  )]
                time_arrays = time_arrays + data_struct[n][3]
            
            for n in range(len(data_struct)):
                intensity_arrays = intensity_arrays + data_struct[n][2]
            
            len_t = len(time_arrays[0])
            self.ragged = False
            for n in range(len(time_arrays)):
                if len(time_arrays[n]) != len_t:
                    self.ragged = True
                    
            self.intensity_vec = intensitys
            self.times = time_arrays
            self.cells = cells
            
            
    def __guess_headers(self,dataframe):
        columns = dataframe.columns.tolist()
        
        cell_header = None; spot_header = None; time_header = None;
        R_header=None; G_header = None; B_header=None;
    
        kwargs = {}
        for col in columns:
            if 'spot' in col.lower():
                spot_header = col
                kwargs['spot_header'] = spot_header
            if np.sum([x in col.lower() for x in ['green','g_']]) >0:
                green_header = col
                kwargs['green_header'] = green_header
            if np.sum([x in col.lower() for x in ['red','r_']]) >0:
                red_header = col
                kwargs['red_header'] = red_header
            if np.sum([x in col.lower() for x in ['blue','b_']]) >0:
                blue_header = col  
                kwargs['blue_header'] = blue_header
            if 'cell' in col.lower():
                cell_header = col
                kwargs['cell_header'] = cell_header
            if np.sum([x in col.lower() for x in ['time','sec']]) >0:
                time_header = col   
                kwargs['time_header'] = time_header
        
        return kwargs
        
    
    def load_dataframe(self,dataframe):
        kwargs = self.__guess_headers(dataframe)
        print(kwargs)
        self.__load_dataframe(dataframe, **kwargs)
        
            
    def __load_dataframe(self, dataframe, cell_header = 'Cell_No', spot_header = 'Spot_No',time_header = 'Time_sec', red_header='Red_Int', green_header = 'Green_int', blue_header='Blue_int' ):
            df = dataframe
            
            self.head = df.head()
            n_cells = max(df[cell_header])+1
            max_spots =  max(df[spot_header])+1
            total_n_spots = 0
            data_struct = []
            

            for i in range(0,n_cells):
                i_g = []
                i_r = []
                i_b = []
                intensitys = []
                times = []
                for j in range(0,max_spots):
                    temp_len_spot =  df[(df[cell_header] == i) & (df[spot_header] ==j) ]
    
                    if len(temp_len_spot) >0:
                        time = df[(df[cell_header] ==i) & (df[spot_header] ==j)][time_header].values
                        intensity_red = df[(df[cell_header] ==i) & (df[spot_header] ==j)][red_header].values
                        intensity_green = df[(df[cell_header] ==i) & (df[spot_header] ==j)][green_header].values
                        intensity_blue = df[(df[cell_header] ==i) & (df[spot_header] ==j)][blue_header].values
                        i_g.append(intensity_green)
                        i_r.append(intensity_red)
                        i_b.append(intensity_blue)
                        
                        intensity = np.vstack((intensity_red,intensity_green,intensity_blue))
                        intensitys.append(intensity)
                        times.append(time)
                data_struct.append([i_r,i_g,i_b,intensitys,times ])
         
            time_arrays = []
            intensity_arrays = []
            cells = []
            for n in range(len(data_struct)):
                cells = cells + [n for x in range( len(data_struct[n][4])  )]
                time_arrays = time_arrays + data_struct[n][4]

            len_t = len(time_arrays[0])
            self.ragged = False
            for n in range(len(time_arrays)):
                if len(time_arrays[n]) != len_t:
                    self.ragged = True
                    
            self.intensity_vec = intensitys
            self.data_struct = data_struct
            self.times = time_arrays
            self.cells = cells               


class OptChain():
    '''
    Container class for parameter optimization chains
    '''
    def __init__(self):
        self.parchain = np.array([])
        self.iterations = 0
        self.parnames = np.array([])
        self.evalchain  = np.array([])
        self.objfunchain = np.array([])
        self.bestpar = None
        self.besteval = None
        self.opt_method = None
        self.opt_args = None
        self.objective_args = None
        self.objective_fun_list  = None
        self.intensity_fun = None
        self.logspace = False
        
    def report(self):
        print('=====================')
        print('Optimizer: %s ran for %d iterations ' % (self.opt_method,self.iterations))
        print('Optimizer arguments: ' + str(self.opt_args))
        print('Objective function: ' + str(self.objective_fun_list))
        print('Objective arguments: ' + str(self.objective_args))
        print('_____________________')
        print('Best Parameter Set: %s, feval: %d'%  (''.join(str(self.bestpar.tolist())),self.besteval ) )
        print('=====================')
        
    def check_parameter_convergence(self):
        
        def get_acc2(data, trunc=False):
            '''
            Get autocorrelation function
    
            *NOT* multi-tau
            '''
            N = len(data)
            fvi = np.fft.fft(data, n=2*N)
            acf = fvi*np.conjugate(fvi)
            acf = np.fft.ifft(acf)
            acf = np.real(acf[:N])/float(N)
            if trunc:
                acf[acf < 0]=0
                for i in range(1, len(acf)):
                    if acf[i] > acf[i-1]:
                        acf[i] = acf[i-1]
            return acf    
        
        for i in range(0,self.parchain.shape[1] ):
            acc = get_acc2(self.parchain[:,i] )       
            if i == 0:
                self.par_acc = acc
            else:
                self.par_acc = np.vstack( (self.par_acc, acc))
        self.par_acc = self.par_acc.T
            
    def check_objfun_convergence(self):
        
        def get_acc2(data, trunc=False):
            '''
            Get autocorrelation function
    
            *NOT* multi-tau
            '''
            N = len(data)
            fvi = np.fft.fft(data, n=2*N)
            acf = fvi*np.conjugate(fvi)
            acf = np.fft.ifft(acf)
            acf = np.real(acf[:N])/float(N)
            if trunc:
                acf[acf < 0]=0
                for i in range(1, len(acf)):
                    if acf[i] > acf[i-1]:
                        acf[i] = acf[i-1]
            return acf    
        
        for i in range(0,self.objfunchain.shape[1] ):
            acc = get_acc2(self.objfunchain[:,i] )       
            if i == 0:
                self.objfun_acc = acc
            else:
                self.objfun_acc = np.vstack( (self.objfun_acc, acc))
                
        self.objfun_acc = self.objfun_acc.T
        
    
    def __clear_invalid_values(self):
        evalchain = self.evalchain
        parchain = self.parchain
        trimmed_evals = evalchain[~np.isnan(evalchain)]        
        final_evals = trimmed_evals[np.isfinite(trimmed_evals)]
        trimmed_parchains = parchain[~np.isnan(evalchain)]
        final_parchains = trimmed_parchains[np.isfinite(trimmed_evals)]
        
        
        return final_evals, final_parchains
        
        
    def parplot(self,ellipse=True,logspace=False):
        n_par = len(self.bestpar)
       
        fig, ax = plt.subplots( n_par,n_par,dpi=300)
        
        plotnum = np.arange(n_par**2).reshape((n_par,n_par)) + 1
        triangle_inds = np.tril(plotnum)[np.where(np.tril(plotnum) !=0 )]
        nplots = len(triangle_inds)
        
        
        if np.sum(np.isnan(self.evalchain))>0 or np.sum(~np.isfinite(self.evalchain)) >0 :
            print('Warning: NaN or Infinite values detected within function evaluation chain, these parameter sets are left out of the plot ')
        eval_chain,par_chain = self.__clear_invalid_values()
         
        
        viridis = cm.get_cmap('viridis', int(np.ceil(np.max(eval_chain))))
        colors = eval_chain
        covariances = np.cov(par_chain.T)
        used_pairs = []
        for i in range(n_par-1,-1,-1):
            for j in range(n_par-1,-1,-1):
                if set([i,j]) not in used_pairs:
                    if i != j:
                        a = ax[i][j].scatter(par_chain[:,i],par_chain[:,j], marker='.',c= colors)
                        if ellipse == True:
                            self.__get_ellipse([np.mean(par_chain[:,i]),np.mean(par_chain[:,j])], covariances,ax=ax[i][j] )
                                         
                    else:
                        b = ax[i][j].hist(par_chain[:,i],bins=40,density=True)
                        
           
                    used_pairs.append(set([i,j]))
                else: 
                    ax[i][j].axis('off')
                
                if len(self.parnames) >0:
                    if i == n_par-1:
                        ax[i][j].set_xlabel(self.parnames[i])
                    if j == 0:
                        ax[i][j].set_ylabel(self.parnames[j])
                        
                if logspace:
                    if i != j:
                        ax[i][j].set_yscale('log')
                    ax[i][j].set_xscale('log')
 
        fig.tight_layout()
        fig.colorbar(a, ax=ax)
        fig.show()
        
    def __get_ellipse(self,mu,cov,ax=None,crosshairs=False,*args,**kwargs):
        ''' Command to plot the fit instance 
        '''
        cmap= cm.viridis
        cmap = cm.coolwarm
        ci = .95
        vals,vecs = np.linalg.eig(cov)
        theta = float(abs((360/(2*np.pi))*np.arctan(vecs[1,0]/vecs[0,0])))
        # able to change CI now.
        scale = chi2.ppf(ci,2)
        w = np.sqrt(vals[0]*scale)*2
        h = np.sqrt(vals[1]*scale)*2
        

        e = Ellipse(xy = tuple(mu) ,width = w, height = h,angle=90-theta,**kwargs )

        ax.add_artist(e)
    
        e.set_clip_box(ax.bbox)
        # e.set_alpha(.75)
        e.set_edgecolor(('red'))
        e.set_linestyle(('-'))
        e.set_facecolor(('none'))
        return ax
    
    def save(self,filename):
        ext = filename.split('.')[-1]
        if ext == 'txt':
            x=1
        if ext in ['p','pickle']:
            pickle.dump(self,open('filename','wb'))
        if ext == 'csv':
            x=1
        if ext == 'npz':
            x=1
            
    def load(self,filename):
        ext = filename.split('.')[-1]
        if ext == 'txt':
            x=1
        if ext in ['p','pickle']:
            tmp_obj = pickle.load(open('filename','wb'))
            for item in tmp_obj.__dict__.keys():
                self.__dict__[item] = tmp_obj.__dict__[item]
            
        if ext == 'csv':
            x=1
        if ext == 'npz':
            x=1
                    

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
                                   'tRNA':False
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
        

    def solve_ballistic_model(self, ki,ke, poi=None, tag= None):
        '''
        

        Parameters
        ----------
        ki : float
            initation rate of ribosomes.
        ke : float
            mean elongation rate along the transcript.
        poi : protein object, optional
            a particular protein object to pull geometry from. The default is None.
        tag : list, optional
            list of tag locations if geometry not given. The default is None.

        Returns
        -------
        tau_analyticals : list of floats
            analytically sovled decorrelation times tau.
        mean_analyticals : list of floats
            the visible ribosome means on the transcript.
        var_analyticals : list of floats
            variances of ribosomes on the transcript.

        '''
   
        if poi == None:
            poi = self._poi
        if tag == None:
            colors = np.where((poi.probe_loc)== 1)[0]
            locs = np.where((poi.probe_loc)== 1)[1]
            tags = []
            for i in range(0,max(colors)+1):
                
                tags.append(locs[np.where(colors == i)].tolist())
            
            
        tau_analyticals = []
        mean_analyticals = []
        var_analyticals = []
        for tag in tags:

            
            L = poi.total_length #get the total length of the gene
            Lm = np.mean(tag)  #the mean location of the tag epitopes
            L_after_tag = L - tag[-1]
            L_tag = int((tag[-1] - tag[0]) / 2)

            ke_analytical = L*ke / np.sum(self.__get_ui(poi.nt_seq[:-3]))

            tau_analytical = (L )/ke_analytical  #analytical tau ie autocovariance time 
            mean_analytical = ki*tau_analytical * (1.-Lm/float(L)) # mean intensity
            var_analytical = ki*tau_analytical * (1.-Lm/float(L))**2  #var intensity
            
            tau_analyticals.append(tau_analytical)
            mean_analyticals.append(mean_analytical)
            var_analyticals.append(var_analytical)
        
        return tau_analyticals,mean_analyticals,var_analyticals
    
    
    def invert_ballistic(self,tau_measured, mu_I, poi= None):
        '''

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
            print(tag)
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


    
    def __get_ui(self, nt_seq):
        '''
        return the ratio of average gene copy number / sequence codon copy number
        '''
        
        codon_dict = CodonDictionaries()
        mean_u = np.mean(list(codon_dict.strGeneCopy_single.values()) )
        ui = []
        for i in range(0, len(nt_seq), 3):
            ui.append(mean_u/ codon_dict.strGeneCopy_single[nt_seq[i:i+3]])
        return ui
        

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
                ssa_obj = self.__solve_ssa_lowmem_combined(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)
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
            

    def solve_ssa_trna(self,k_index, k_diffusion, k_bind, kelong, k_compl, t,x0=[], k_trna = None, perturb=[0,0,0],leaky_probes=False,kprobe=np.ones(1),probe_vec = None, probe_loc=None, kon=1,koff=1,bursting=False,n_traj=10   ):
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
            
            strGeneCopy = CodonDictionaries().strGeneCopy_single
            strGeneCopy.pop('TAG')
            strGeneCopy.pop('TAA')
            strGeneCopy.pop('TGA')

            k_trna = np.array(list(CodonDictionaries().strGeneCopy_single.values()))
        
        
        if not provided_probe:
            if provided_protein:
                probe_vec = self.protein.probe_vec.astype(np.int32)
                probe_loc = self.protein.probe_loc.astype(np.int32)
            else:
                print("no provided probe vector, please set the solver.protein with a protein object or provide a probe vector")
                raise 
        else:
            probe_vec = probe_vec
            
            
        ssa_conditions['probe_vec'] = probe_vec
        ssa_conditions['probe_loc'] = probe_loc
        
        
     
        ssa_obj = self.__solve_ssa_trna(k_index,k_trna, k_diffusion,k_bind,kelong, k_compl,t,x0,n_traj,ssa_conditions = ssa_conditions)
                        
        return ssa_obj
    
    def solve_ssa(self,k,t,x0=[],n_traj=100,bins=None,low_memory=True,perturb=[0,0,0],leaky_probes=False,kprobe=np.ones(1),record_stats=False,probe_vec = None, probe_loc=None, kon=1,koff=1,bursting=False):
        
        self.__check_rates(k)
        
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
            
            
        if not provided_probe:
            if provided_protein:
                probe_vec = self.protein.probe_vec.astype(np.int32)
                probe_loc = self.protein.probe_loc.astype(np.int32)
            else:
                print("no provided probe vector, please set the solver.protein with a protein object or provide a probe vector")
                raise 
        else:
            probe_vec = probe_vec
            
            
        ssa_conditions['probe_vec'] = probe_vec
        ssa_conditions['probe_loc'] = probe_loc
        
        
        if self.cython_available:
            if low_memory:
                ssa_obj = self.__solve_ssa_lowmem_combined(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)
            else:
                ssa_obj = self.__solve_ssa_lowmem_combined(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)
 
                
                # if record_stats:
                
                #     if leaky_probes == False:
                #         # if low_memory:
                #         #     ssa_obj = self.__solve_ssa_lowmem(k,t,x0,n_traj, ssa_conditions = ssa_conditions)
                #         # else:
                #         ssa_obj = self.__solve_ssa_lowmem_combined(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)

                #         #ssa_obj = self.__solve_ssa(k,t,x0,n_traj,ssa_conditions = ssa_conditions)
                #     else:
                #         # if low_memory:
                #         #     ssa_obj = self.__solve_ssa_lowmem_leaky(k,t,x0,k_probe,n_traj, ssa_conditions = ssa_conditions)
                #         # else:
                #         ssa_obj = self.__solve_ssa_lowmem_combined(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)


                #         #ssa_obj = self.__solve_ssa_leaky(k,t,x0,kprobe,n_traj,ssa_conditions = ssa_conditions)         
                   
                # else:
                #     if leaky_probes == False:
                #         # if low_memory:
                #         #     ssa_obj = self.__solve_ssa_lowmem_nostats(k,t,x0,n_traj, ssa_conditions = ssa_conditions)
                #         # else:
                            
                #         ssa_obj = self.__solve_ssa_lowmem_combined(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)

                #         #ssa_obj = self.__solve_ssa_nostats(k,t,x0,n_traj,ssa_conditions = ssa_conditions)
                #     else:
                #         # if low_memory:
                #         #     ssa_obj = self.__solve_ssa_lowmem_leaky_nostats(k,t,x0,k_probe,n_traj, ssa_conditions = ssa_conditions)
                #         # else:
                #         ssa_obj = self.__solve_ssa_lowmem_combined(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)


                #         #ssa_obj = self.__solve_ssa_leaky_nostats(k,t,x0,kprobe,n_traj,ssa_conditions = ssa_conditions)         
                           
        else:
            ssa_obj = self.__solve_ssa_python(k,t,x0,n_traj,ssa_conditions = ssa_conditions, kon=kon, koff=koff, kprobe=kprobe)

       
        
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
    
    

    def solve_ode(self,k,t,x0,kbind,pl,bins=None, corr = False):
        st = time.time()
        m_ode = models.TranslateODE()
        m_ode.N = len(k)
        m_ode.tf = t[-1]
        m_ode.ptimes = len(t)
        m_ode.ke = k
        m_ode.kb = kbind 
        m_ode.fi = 1
        m_ode.ti = t[0]
        m_ode.binary = pl
        
        #mu_intensity = m_ode.solve()
        
        m = models.TranslateCorrs()
        m.N = len(k)
        m.tf = t[-1]
        m.ptimes = len(t)
        m.ke = k
        m.kb = kbind 
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
        
        ode_soln.N = len(k)
        ode_soln.k = k
        ode_soln.x0 = x0
        ode_soln.fi = 1
        ode_soln.kb = kbind
        ode_soln.prob_loc = pl
        ode_soln.solve_time = time.time()-st
        
        if len(k) > 100:
            ode_soln.solve_method = 'expv'
        else:
            ode_soln.solve_method = 'odeint'
        
        
        return ode_soln
    
    
    
    
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
        
        print(t)
        
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
    

    
                
    def __map_to_intensity(self): 
        
        return 1           


    def __solve_ssa_trna(self,kindex,ktrna,kdiffusion,kbind,kelong,kcompl, t, x0, n_traj, ssa_conditions=None ):

        

        seeds = np.random.randint(0, 0x7FFFFFF, n_traj, dtype = np.int32)
    
        if ssa_conditions == None:
            ssa_conditions = self.default_conditions
        
        x0 = self.__check_x0(x0)
   
       
        rib_vec = []
        solutions = []            
        solutionssave = []       
        N_rib = 200
        colors = self.colors
        
        n_trajectories = n_traj
        

        all_results,all_nribs,all_collisions,all_frapresults,all_ribtimes,all_col_points = self.__generate_mats(n_traj,kbind,t,N_rib,colors)
        all_trna_results = np.zeros((n_trajectories,61*len(t)),dtype=np.int32)
        
        footprint = ssa_conditions['footprint']
        evf = ssa_conditions['perturb'][0]
        evi = ssa_conditions['perturb'][1]
        intime = ssa_conditions['perturb'][2]
        non_consider_time = ssa_conditions['burnin']
        
        st = time.time()
        
        N_rib = 200
        
       
        result = np.zeros((len(t)*N_rib),dtype=np.int32  )
        
        #lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t<20)[0]))
        #all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)

        all_ribs = np.zeros((n_trajectories,1))
        all_col_points = []
    
        #seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)
        x0 = np.zeros((N_rib),dtype=np.int32)
     
        
        for i in range(n_trajectories):
            
            trna_result = np.zeros((len(t)*61),dtype=np.int32)    
            result,ribtimes,frapresult,coltimes,colpointsx,colpointst = self.__generate_vecs_trna(kbind,kindex,t,N_rib,colors)

            nribs = np.array([0],dtype=np.int32)
            
            ssa_trna.run_SSA(result,trna_result,ribtimes,coltimes,colpointsx,colpointst, kindex,ktrna,kdiffusion,frapresult,t,kbind,kcompl, 0,0,0, seeds[i],nribs,x0,kelong)
                  
            all_results[i, :] = result.T
            all_trna_results[i,:] = trna_result
            all_frapresults[i,:] = frapresult
            all_ribtimes[i,:] = ribtimes
            all_collisions[i,:] = coltimes
            all_nribs[i,:] = nribs[0]
         
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


        no_ribosomes = np.zeros((n_traj, (len(kindex)+1)))
        
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
        try:
            ssa_obj.col_points = all_col_points
        except:
            pass
        
        
        
        return ssa_obj                

    def __solve_ssa_python(self,k,t,x0,n_traj,ssa_conditions=None, kprobe=None, kon=None, koff=None, flags=[0,0,0]):
        
        seeds = np.random.randint(0, 0x7FFFFFF, n_traj, dtype=np.int32)
        
        if ssa_conditions == None:
            ssa_conditions = self.default_conditions
        
        x0 = self.__check_x0(x0)
   
        flags = [int(ssa_conditions['bursting']), int(ssa_conditions['leaky_probes']), int(ssa_conditions['record_stats'])]

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
        
        
        print('C++ library failed, Using Python Implementation')
        rib_vec = []

        solutions = []            
        solutionssave = []
        N_rib = 200
        collisions = np.array([[]])
        all_results = np.zeros((n_traj, N_rib*len(t)), dtype=np.int32)
        I_internal = np.zeros((colors,len(t), n_traj)) 
        all_col_points = []
        watched_ribs = []
        for i in range(n_traj):
            
            soln,all_ribtimes,Ncol,col_points,intensity = self.__ssa_python(k, t, inhibit_time=intime+non_consider_time, FRAP=evf, Inhibitor=evi, flags=flags, kon=kon, kprobe=kprobe, koff=koff, ssa_conditions=ssa_conditions)
            #soln = soln.reshape((1, (len(time_vec_fixed)*N_rib)))
            
            collisions = np.append(collisions,Ncol)
            watched_ribs.append(int(len(collisions)))
            validind = np.where(np.sum(soln,axis=1)!=0)[0]
            all_col_points.append(np.array(col_points))
            if np.max(validind) != N_rib-1:
                validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)
            
            so = soln[(validind,)]
          
            #solutionssave.append(so)

            solutions.append(soln)
        
            result = soln.reshape((1, (len(t)*N_rib)))
            all_results[i, :] = result
            I_internal[:,:,i] = intensity
        
                        
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
        
    def __solve_ssa_lowmem_nostats(self,k,t,x0,n_traj,ssa_conditions=None):
        seeds = np.random.randint(0, 0x7FFFFFF, n_traj, dtype=np.int32)
        
        if isinstance(k,list):
            k = np.array(k).astype(np.float64)
        
        k = k.flatten()

        
        
        if ssa_conditions == None:
            ssa_conditions = self.default_conditions
        
        x0 = self.__check_x0(x0)
        
        
        pl = ssa_conditions['probe_vec']
        
        
        colors = self.colors
 
        
        rib_vec = []
        solutions = []            
        solutionssave = []       
        N_rib = 1
        all_results,all_frapresults = self.__generate_mats_lowmem_nostats(n_traj,k[0],t,N_rib,colors)
        footprint = ssa_conditions['footprint']
        evf = ssa_conditions['perturb'][0]
        evi = ssa_conditions['perturb'][1]
        intime = ssa_conditions['perturb'][2]
        non_consider_time = ssa_conditions['burnin']
        
        st = time.time()
        
        for i in range(n_traj):
    
            result,frapresult = self.__generate_vecs_lowmem_nostats(k,t,N_rib,colors)
            nribs = np.array([0],dtype=np.int32)


            ssa_translation_lowmem_nostats.run_SSA(result, k[1:-1],frapresult,t,k[0], float(k[-1]), evf, evi, float(intime), seeds[i],x0,footprint, pl,colors)

            #ssa_translation_lowmem_nostats.run_SSA(result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1],frapresult, t, k[0], float(k[-1]), evf, evi, float(intime), seeds[i],nribs,x0,footprint, pl,colors)
            #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
            all_results[i, :] = result.T
            all_frapresults[i,:] = frapresult

            for i in range(n_traj):
                soln = all_results[i, :].reshape((N_rib, len(t),colors))

                so = soln
                solutionssave.append(so)
                solutions.append(soln)
            

            
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
        ssa_obj.no_ribosomes = None
        ssa_obj.n_traj = n_traj
        ssa_obj.k = k
        #ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        #ssa_obj.rib_density = ribosome_density
        #ssa_obj.rib_means = ribosome_means
        ssa_obj.rib_vec = rib_vec
        ssa_obj.intensity_vec = all_results.T[:,startindex:,:]
        ssa_obj.I = all_results.T[:,startindex:,:]
        ssa_obj.time_vec_fixed = t
        ssa_obj.time = t[startindex:]
        ssa_obj.time_rec = t[startindex:]
        ssa_obj.start_time = non_consider_time
        ssa_obj.watched_ribs = None
        try:
            ssa_obj.col_points = None
        except:
            pass
        
        
        ssa_obj.eval_time = sttime
        
        return ssa_obj
    
 
    def __solve_ssa_lowmem_combined(self,k,t,x0,n_traj,ssa_conditions=None, kon=1, koff=1, kprobe=[] ):
        seeds = np.random.randint(0, 0x7FFFFFF, n_traj, dtype=np.int32)
        
        if isinstance(k,list):
            k = np.array(k).astype(np.float64)
        
        k = k.flatten()

        if kprobe == []:
            kprobe = np.ones(self.color)
            
        if isinstance(kprobe,list):
            kprobe = np.array(kprobe,dtype=np.float64)
        if isinstance(kprobe,int):
            kprobe = np.array([kprobe],dtype=np.float64)
        if isinstance(kprobe,float):
            kprobe = np.array([kprobe],dtype=np.float64)
     
        
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
            
        flags = np.array([int(ssa_conditions['bursting']), int(ssa_conditions['leaky_probes']), int(ssa_conditions['record_stats'])],dtype=np.int32)
        probe_loc = ssa_conditions['probe_loc']
        
        
        rib_vec = []
        solutions = []            
        solutionssave = []       
        N_rib = 1
        all_results,all_nribs,all_collisions,all_frapresults,all_ribtimes,all_col_points = self.__generate_mats_lowmem(n_traj,k[0],t,N_rib,colors)
        footprint = ssa_conditions['footprint']
        evf = ssa_conditions['perturb'][0]
        evi = ssa_conditions['perturb'][1]
        intime = ssa_conditions['perturb'][2]
        non_consider_time = ssa_conditions['burnin']
        
        st = time.time()
        
        for i in range(n_traj):
            
            result,ribtimes,frapresult,coltimes,colpointsx,colpointst = self.__generate_vecs_lowmem(k,t,N_rib,colors)
            nribs = np.array([0],dtype=np.int32)
            
            ssa_translation_lowmem.run_SSA(result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1],frapresult, t, k[0], float(k[-1]), int(evf), int(evi), float(intime), seeds[i],nribs,x0,footprint, probe_vec ,int(colors), kon, koff, kprobe, probe_loc, flags)
            #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
            all_results[i, :] = result.T
            all_frapresults[i,:] = frapresult
            all_ribtimes[i,:] = ribtimes
            all_collisions[i,:] = coltimes
            all_nribs[i,:] = nribs
            
            endcolrec = np.where(colpointsx == 0)[0][0]
            
            colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
            all_col_points.append(colpoints.T)
                
    
            for i in range(n_traj):
                soln = all_results[i, :].reshape((N_rib, len(t),colors))

                so = soln
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
        ssa_obj.no_ribosomes = no_ribosomes
        ssa_obj.n_traj = n_traj
        ssa_obj.k = k
        #ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        #ssa_obj.rib_density = ribosome_density
        #ssa_obj.rib_means = ribosome_means
        ssa_obj.rib_vec = rib_vec
        ssa_obj.intensity_vec = all_results.T
        ssa_obj.I = all_results.T
        ssa_obj.time_vec_fixed = t
        ssa_obj.time = t
        ssa_obj.time_rec = t[startindex:]
        ssa_obj.start_time = non_consider_time
        ssa_obj.watched_ribs = watched_ribs
        ssa_obj.collisions = collisions
        
        
        try:
            ssa_obj.col_points = all_col_points
        except:
            pass
        
        
        ssa_obj.eval_time = sttime
        
        return ssa_obj        
 
    
#     def __solve_ssa_lowmem(self,k,t,x0,n_traj,ssa_conditions=None):
#         seeds = np.random.randint(0, 0x7FFFFFF, n_traj)
        
#         if isinstance(k,list):
#             k = np.array(k).astype(np.float64)
        
#         k = k.flatten()

        
        
#         if ssa_conditions == None:
#             ssa_conditions = self.default_conditions
        
#         x0 = self.__check_x0(x0)
        
#         pl = ssa_conditions['probe_vec']
#         colors = self.colors
        
#         rib_vec = []
#         solutions = []            
#         solutionssave = []       
#         N_rib = 1
#         all_results,all_nribs,all_collisions,all_frapresults,all_ribtimes,all_col_points = self.__generate_mats_lowmem(n_traj,k[0],t,N_rib,colors)
#         footprint = ssa_conditions['footprint']
#         evf = ssa_conditions['perturb'][0]
#         evi = ssa_conditions['perturb'][1]
#         intime = ssa_conditions['perturb'][2]
#         non_consider_time = ssa_conditions['burnin']
        
#         st = time.time()
        
#         for i in range(n_traj):
            
#             result,ribtimes,frapresult,coltimes,colpointsx,colpointst = self.__generate_vecs_lowmem(k,t,N_rib,colors)
#             nribs = np.array([0],dtype=np.int32)

#             ssa_translation_lowmem.run_SSA(result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1],frapresult, t, k[0], float(k[-1]), evf, evi, float(intime), seeds[i],nribs,x0,footprint, pl,colors)
#             #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
#             all_results[i, :] = result.T
#             all_frapresults[i,:] = frapresult
#             all_ribtimes[i,:] = ribtimes
#             all_collisions[i,:] = coltimes
#             all_nribs[i,:] = nribs
            
#             endcolrec = np.where(colpointsx == 0)[0][0]
            
#             colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
#             all_col_points.append(colpoints.T)
                
    
#             for i in range(n_traj):
#                 soln = all_results[i, :].reshape((N_rib, len(t),colors))

#                 so = soln
#                 solutionssave.append(so)
#                 solutions.append(soln)
            
#             collisions = np.array([[]])
#             watched_ribs = []
#             for i in range(n_traj):
#                 totalrib = all_nribs[i]
            
#                 if totalrib > all_collisions.shape[1]:
#                     collisions = np.append(collisions, all_collisions[i][:])
#                     watched_ribs.append(int(all_collisions.shape[1]))
            
#                 else:
                   
#                     collisions = np.append(collisions, all_collisions[i][:int(totalrib[0])])
#                     watched_ribs.append(int(totalrib[0]))
            
#             sttime = time.time() - st


#         no_ribosomes = np.zeros((n_traj, (len(k)+1)))
        
#         startindex = np.where(t >= non_consider_time)[0][0]
        
#         #all_results = all_results[:,startindex*N_rib:]

# #        for i in range(len(solutions)):
# #            for j in range(len(solutions[0][0][startindex:])):
# #                rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]
# #                print(rib_pos)
# #               
# #                no_ribosomes[i, rib_pos.astype(int)] += 1
# #        no_ribosomes = no_ribosomes[:, 1:]
# #
# #        ribosome_means = np.mean(no_ribosomes, axis=0)
# #        ribosome_density = ribosome_means/len(k)
# #
# #        no_ribosomes_per_mrna = np.mean(no_ribosomes)
        
#         ssa_obj = SSA_Soln()
#         ssa_obj.no_ribosomes = no_ribosomes
#         ssa_obj.n_traj = n_traj
#         ssa_obj.k = k
#         #ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
#         #ssa_obj.rib_density = ribosome_density
#         #ssa_obj.rib_means = ribosome_means
#         ssa_obj.rib_vec = rib_vec
#         ssa_obj.intensity_vec = all_results.T
#         ssa_obj.I = all_results.T
#         ssa_obj.time_vec_fixed = t
#         ssa_obj.time = t
#         ssa_obj.time_rec = t[startindex:]
#         ssa_obj.start_time = non_consider_time
#         ssa_obj.watched_ribs = watched_ribs
#         ssa_obj.collisions = collisions
        
        
#         try:
#             ssa_obj.col_points = all_col_points
#         except:
#             pass
        
        
#         ssa_obj.eval_time = sttime
        
#         return ssa_obj
            

#     def __solve_ssa_lowmem_leaky_nostats(self,k,t,x0,k_probe,n_traj,ssa_conditions=None):
#         seeds = np.random.randint(0, 0x7FFFFFF, n_traj)
#         k = np.array(k).astype(np.float64)
#         if ssa_conditions == None:
#             ssa_conditions = self.default_conditions
        
#         x0 = self.__check_x0(x0)
        
#         pl = self.protein.probe_loc.astype(int)
        
        
#         if isinstance(k_probe,list):
#             k_probe = np.array(k_probe)
        
#         pv = ssa_conditions['probe_vec']
#         colors = self.colors
        
#         rib_vec = []
#         solutions = []            
#         solutionssave = []       
#         N_rib = 1
#         all_results,all_frapresults = self.__generate_mats_lowmem_nostats(n_traj,k[0],t,N_rib,colors)
#         footprint = ssa_conditions['footprint']
#         evf = ssa_conditions['perturb'][0]
#         evi = ssa_conditions['perturb'][1]
#         intime = ssa_conditions['perturb'][2]
#         non_consider_time = ssa_conditions['burnin']
        
#         st = time.time()
        
#         for i in range(n_traj):
            
#             result,frapresult = self.__generate_vecs_lowmem_nostats(k,t,N_rib,colors)
#             nribs = np.array([0],dtype=np.int32)

#             ssa_translation_lowmem_leaky.run_SSA(result, k[1:-1],frapresult, t, k[0], k[-1], evf, evi, intime, seeds[i],x0,footprint, pv,k_probe,pl, colors)
#             #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
#             all_results[i, :] = result.T
#             all_frapresults[i,:] = frapresult


#             for i in range(n_traj):
#                 soln = all_results[i, :].reshape((N_rib, len(t),colors))

#                 so = soln
#                 solutionssave.append(so)
#                 solutions.append(soln)
            
#             collisions = np.array([[]])
#             watched_ribs = []

            
#             sttime = time.time() - st


#         no_ribosomes = np.zeros((n_traj, (len(k)+1)))
        
#         startindex = np.where(t >= non_consider_time)[0][0]
        
#         #all_results = all_results[:,startindex*N_rib:]

# #        for i in range(len(solutions)):
# #            for j in range(len(solutions[0][0][startindex:])):
# #                rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]
# #                print(rib_pos)
# #               
# #                no_ribosomes[i, rib_pos.astype(int)] += 1
# #        no_ribosomes = no_ribosomes[:, 1:]
# #
# #        ribosome_means = np.mean(no_ribosomes, axis=0)
# #        ribosome_density = ribosome_means/len(k)
# #
# #        no_ribosomes_per_mrna = np.mean(no_ribosomes)
        
#         ssa_obj = SSA_Soln()
#         ssa_obj.no_ribosomes = no_ribosomes
#         ssa_obj.n_traj = n_traj
#         ssa_obj.k = k
#         #ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
#         #ssa_obj.rib_density = ribosome_density
#         #ssa_obj.rib_means = ribosome_means
#         ssa_obj.rib_vec = rib_vec
#         ssa_obj.intensity_vec = all_results.T
#         ssa_obj.I = all_results.T
#         ssa_obj.time_vec_fixed = t
#         ssa_obj.time = t
#         ssa_obj.time_rec = t[startindex:]
#         ssa_obj.start_time = non_consider_time
#         ssa_obj.watched_ribs = watched_ribs

#         ssa_obj.eval_time = sttime
#         return ssa_obj       
    


#     def __solve_ssa_lowmem_leaky(self,k,t,x0,k_probe,n_traj,ssa_conditions=None):
#         seeds = np.random.randint(0, 0x7FFFFFF, n_traj)
#         k = np.array(k).astype(np.float64)
#         if ssa_conditions == None:
#             ssa_conditions = self.default_conditions
        
#         x0 = self.__check_x0(x0)
        
#         pl = self.protein.probe_loc.astype(int)
        
        
#         if isinstance(k_probe,list):
#             k_probe = np.array(k_probe)
        
#         pv = ssa_conditions['probe_vec']
#         colors = self.colors
        
#         rib_vec = []
#         solutions = []            
#         solutionssave = []       
#         N_rib = 1
#         all_results,all_nribs,all_collisions,all_frapresults,all_ribtimes,all_col_points = self.__generate_mats_lowmem(n_traj,k[0],t,N_rib,colors)
#         footprint = ssa_conditions['footprint']
#         evf = ssa_conditions['perturb'][0]
#         evi = ssa_conditions['perturb'][1]
#         intime = ssa_conditions['perturb'][2]
#         non_consider_time = ssa_conditions['burnin']
        
#         st = time.time()
        
#         for i in range(n_traj):
            
#             result,ribtimes,frapresult,coltimes,colpointsx,colpointst = self.__generate_vecs_lowmem(k,t,N_rib,colors)
#             nribs = np.array([0],dtype=np.int32)

#             ssa_translation_lowmem_leaky.run_SSA(result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1],frapresult, t, k[0], k[-1], evf, evi, intime, seeds[i],nribs,x0,footprint, pv,k_probe,pl, colors)
#             #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
#             all_results[i, :] = result.T
#             all_frapresults[i,:] = frapresult
#             all_ribtimes[i,:] = ribtimes
#             all_collisions[i,:] = coltimes
#             all_nribs[i,:] = nribs
            
#             endcolrec = np.where(colpointsx == 0)[0][0]
            
#             colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
#             all_col_points.append(colpoints.T)
                
    
#             for i in range(n_traj):
#                 soln = all_results[i, :].reshape((N_rib, len(t),colors))

#                 so = soln
#                 solutionssave.append(so)
#                 solutions.append(soln)
            
#             collisions = np.array([[]])
#             watched_ribs = []
#             for i in range(n_traj):
#                 totalrib = all_nribs[i]
            
#                 if totalrib > all_collisions.shape[1]:
#                     collisions = np.append(collisions, all_collisions[i][:])
#                     watched_ribs.append(int(all_collisions.shape[1]))
            
#                 else:
                   
#                     collisions = np.append(collisions, all_collisions[i][:int(totalrib[0])])
#                     watched_ribs.append(int(totalrib[0]))
            
#             sttime = time.time() - st


#         no_ribosomes = np.zeros((n_traj, (len(k)+1)))
        
#         startindex = np.where(t >= non_consider_time)[0][0]
        
#         #all_results = all_results[:,startindex*N_rib:]

# #        for i in range(len(solutions)):
# #            for j in range(len(solutions[0][0][startindex:])):
# #                rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]
# #                print(rib_pos)
# #               
# #                no_ribosomes[i, rib_pos.astype(int)] += 1
# #        no_ribosomes = no_ribosomes[:, 1:]
# #
# #        ribosome_means = np.mean(no_ribosomes, axis=0)
# #        ribosome_density = ribosome_means/len(k)
# #
# #        no_ribosomes_per_mrna = np.mean(no_ribosomes)
        
#         ssa_obj = SSA_Soln()
#         ssa_obj.no_ribosomes = no_ribosomes
#         ssa_obj.n_traj = n_traj
#         ssa_obj.k = k
#         #ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
#         #ssa_obj.rib_density = ribosome_density
#         #ssa_obj.rib_means = ribosome_means
#         ssa_obj.rib_vec = rib_vec
#         ssa_obj.intensity_vec = all_results.T
#         ssa_obj.I = all_results.T
#         ssa_obj.time_vec_fixed = t
#         ssa_obj.time = t
#         ssa_obj.time_rec = t[startindex:]
#         ssa_obj.start_time = non_consider_time
#         ssa_obj.watched_ribs = watched_ribs
#         try:
#             ssa_obj.col_points = all_col_points
#         except:
#             pass
        
        
#         ssa_obj.eval_time = sttime
#         return ssa_obj       
    
    
    @classmethod
    def __get_ribosome_statistics(self,ssa_obj,result):
        
        
        return ssa_obj
    
        
    @classmethod
    def __generate_vecs(cls,k,t,N_rib,ncolor):
        tf = t[-1]
        ki = k[0]
        
        guessed_no_ribosomes = int(1.3*ki*tf)
        result = np.zeros((len(t)*N_rib), dtype=np.int32)
        ribtimes = np.zeros((guessed_no_ribosomes),dtype=np.float64)
        frapresult = np.zeros((len(t)*N_rib),dtype=np.int32)
        coltimes = np.zeros((guessed_no_ribosomes),dtype=np.int32)
        colpointsx = np.zeros(len(k[1:-1])*(guessed_no_ribosomes),dtype=np.int32)
        colpointst = np.zeros(len(k[1:-1])*(guessed_no_ribosomes),dtype=np.float64)
        return result,ribtimes,frapresult,coltimes,colpointsx,colpointst

    @classmethod
    def __generate_vecs_trna(cls,kbind,kind,t,N_rib,ncolor):
        tf = t[-1]
        ki = kbind
        
        guessed_no_ribosomes = int(1.3*ki*tf)
        result = np.zeros((len(t)*N_rib), dtype=np.int32)
        ribtimes = np.zeros((guessed_no_ribosomes),dtype=np.float64)
        frapresult = np.zeros((len(t)*N_rib),dtype=np.int32)
        coltimes = np.zeros((guessed_no_ribosomes),dtype=np.int32)
        colpointsx = np.zeros(len(kind)*(guessed_no_ribosomes),dtype=np.int32)
        colpointst = np.zeros(len(kind)*(guessed_no_ribosomes),dtype=np.float64)
        return result,ribtimes,frapresult,coltimes,colpointsx,colpointst
    
    @classmethod
    def __generate_vecs_lowmem(cls,k,t,N_rib,ncolor):
        tf = t[-1]
        ki = k[0]
        
        guessed_no_ribosomes = int(1.3*ki*tf)
        result = np.zeros((ncolor,len(t)), dtype=np.int32)
        ribtimes = np.zeros((guessed_no_ribosomes),dtype=np.float64)
        frapresult = np.zeros((len(t)*N_rib),dtype=np.int32)
        coltimes = np.zeros((guessed_no_ribosomes),dtype=np.int32)
        colpointsx = np.zeros(len(k[1:-1])*(guessed_no_ribosomes),dtype=np.int32)
        colpointst = np.zeros(len(k[1:-1])*(guessed_no_ribosomes),dtype=np.float64)
        return result,ribtimes,frapresult,coltimes,colpointsx,colpointst

    @classmethod
    def __generate_mats_lowmem(cls,ntraj, ki,t,N_rib,ncolor):
        tf = t[-1]
        guessed_no_ribosomes = int(1.3*ki*tf)           
        all_results = np.zeros((ntraj, len(t),ncolor), dtype=np.int32)
        all_ribtimes = np.zeros((ntraj,guessed_no_ribosomes),dtype=np.float64)
        all_frapresults = np.zeros((ntraj,N_rib*len(t)),dtype=np.int32)
        all_collisions = np.zeros((ntraj,guessed_no_ribosomes),dtype=np.int32)
        all_nribs = np.zeros((ntraj,1))
        all_col_points = []
    
        return all_results,all_nribs,all_collisions,all_frapresults,all_ribtimes,all_col_points
        
        
    @classmethod
    def __generate_mats(cls,ntraj, ki,t,N_rib,ncolor):
        tf = t[-1]
        guessed_no_ribosomes = int(1.3*ki*tf)           
        all_results = np.zeros((ntraj, N_rib*len(t)), dtype=np.int32)
        all_ribtimes = np.zeros((ntraj,guessed_no_ribosomes),dtype=np.float64)
        all_frapresults = np.zeros((ntraj,N_rib*len(t)),dtype=np.int32)
        all_collisions = np.zeros((ntraj,guessed_no_ribosomes),dtype=np.int32)
        all_nribs = np.zeros((ntraj,1))
        all_col_points = []
        
        

    
        return all_results,all_nribs,all_collisions,all_frapresults,all_ribtimes,all_col_points
 
    @classmethod
    def __generate_mats_lowmem_nostats(cls,ntraj, ki,t,N_rib,ncolor):

        all_results = np.zeros((ntraj, len(t),ncolor), dtype=np.int32)
        all_frapresults = np.zeros((ntraj,N_rib*len(t)),dtype=np.int32)
        return all_results,all_frapresults
        
    @classmethod
    def __generate_vecs_lowmem_nostats(cls,k,t,N_rib,ncolor):


        result = np.zeros((ncolor,len(t)), dtype=np.int32)

        frapresult = np.zeros((len(t)*N_rib),dtype=np.int32)

        return result,frapresult
    
    def __check_x0(self,x0):
        if len(x0) == 0:
            x0_clean = np.zeros((200),dtype=np.int32)  
        else:
            if len(x0) >200:
                raise ValueError('Unrecognized initial condition, make sure the length of x0 is <= 200')
            
            x0_clean = np.zeros((200),dtype=np.int32)  
            x0_clean[:len(x0)] = x0
        return x0_clean

    def __check_rates_trna(self,rates):
        if isinstance(rates,np.ndarray):
            pass
        else:
            rates = np.array(rates,dtype=np.int32)
        if rates.dtype !=int:
            raise ValueError('trna ID values are not integers, to run the trna ssa version, k_index needs to be an ID of 0-60 integer values')
        if len(np.where(rates < 0)[0]) or len(np.where(rates > 60)[0]):
            raise ValueError('trna ID values are out of bounds, the IDs are 0-60 for tRNA species')             
        
    def __check_rates(self,rates):
        if isinstance(rates,np.ndarray):
            pass
        else:
            rates = np.array(rates)
        if len(np.where(rates < 0)[0]) > 0:
            raise ValueError('one or more model rates are negative, double check the provided rates')
       

    def __ssa_python(self, k, t_array, inhibit_time=0, FRAP=False, Inhibitor=False, flags=None, kon=1, koff=1, kprobe=1, ssa_conditions=None):
        '''
        mRNA Translation simulation python implementation

        given a propensity vector k, time array to record, and inhibitory conditions, run a single trajectory of translation simulation

        The simulation is described here: [PUT LINK HERE TO PAPER]

        *args*

            **k**, propensity vector of size gene length + 2, [initiation rate,  Codon dependent rates,  completion rate / unbinding rate]
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
            leaky_probe_matrix = np.zeros((colors,N,N_rib))
        intensity = np.zeros((colors,len(t_array)))
        probe = np.array(np.where(ssa_conditions['probe_loc'] ==1))
        probevec = ssa_conditions['probe_vec']
        
            
        #example X arrays and how its formatted:
        # X = [423 30 10 0 ]  read from left to right theres a ribosome in position 423 30 and 10, with a 0 kept as a buffer for simulation

        t = t_array[0]  #time point
        Nt = len(t_array)  #number of time points to record over
        tf = t_array[-1]  #final time point
        col = np.zeros((1,N_rib))
        X_array = np.zeros((N_rib, Nt))  #recording array that records the ribosome posistions over time array points
        NR = 0  #number of ribosomes bound
        it = 1  #number of iterations
        Sn_p = np.eye(max(NR+1, 2), dtype=int) #stoichiometry for the SSA
        if bursting:
            wn_p = np.zeros((X.shape[0]+1, 1)) # propensities for the SSA
        else:
            wn_p = np.zeros((X.shape[0], 1)) # propensities for the SSA
        
        T = np.array([0, 0], dtype=float)
        ribtimes = np.array([[0,0]],dtype=float)
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
                    leaky_probe_matrix[ probe[0] , probe[1] ,NR-1 ] = (np.random.rand(len(probe[0])) < kprobe).astype(int)
                
                    
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
                    c = leaky_probe_matrix[:,:,0]
                    for ribind in validx:
                        int_tmp += np.sum(c[:,:(ribind-1)],axis=1)
                        
                    intensity[:,it] = int_tmp
                    
                else:
                    validx = X[X>0]
                    int_tmp = np.zeros(colors)
                    for ribind in validx:
                        int_tmp += probevec[:,(ribind-1)]
                    intensity[:,it] = int_tmp
                    
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
            
                    if np.sum(Sn[:,event]) < 0 :
                        
                        ribtimes = np.vstack((ribtimes,[T[0],t]))
                        T[:-1] = T[1:]
                        Ncol = np.append(Ncol,col[0][0] )
                        col = np.atleast_2d(np.append(col[:,1:],[0]))
                        
                    else:
                        if X[event-1] == X[event] + R:
                            col[0][event] +=1
                            col_points.append( (X[event],t) )                   

            else:
                X = (X + Sn[:, event].T)  #update X vector for new ribosome state
            
            
                if np.sum(Sn[:,event]) < 0 :
                    
                    ribtimes = np.vstack((ribtimes,[T[0],t]))
                    T[:-1] = T[1:]
                    Ncol = np.append(Ncol,col[0][0] )
                    col = np.atleast_2d(np.append(col[:,1:],[0]))
                    
                    if leaky: #remove the leaky probes
                        leaky_probe_matrix[ probe[0] , probe[1] ,0 ] = leaky_probe_matrix[ probe[0] , probe[1] ,1 ]
                
                    
                else:
                    if X[event-1] == X[event] + R:
                        col[0][event] +=1
                        col_points.append( (X[event],t) )
                    
                
            
        return X_array,ribtimes[1:,:],Ncol,col_points,intensity  #return the completed simulation


class suite():
    '''
    a coterie of constructs
    '''

    def __init__(self):
        self.pois = []
        self.discernable = False
        self.models = []
        self.combo_list = []
        
        


class poi():
    '''
    Protein of Intrest class

    Holds all the information for analyzed proteins
    '''
    def __init__(self):
        self.aa_seq = ''    #amino sequence
        self.nt_seq = ''    #nucleotide sequence
        self.gene_length = 0   #length of the gene
        self.tag_length = 0   #length of the tags
        self.total_length = 0  #total length of the full amino acid sequence
        self.name = ''         #name of the gene
        self.tag_types = []
        self.tag_epitopes = {}  #type of tags and epitope lists per tag
        self.ki = .03
        self.ke_mu = 10
        self.kt = 10

    @property
    def CAI(self):
        cs, CAI, cc = SequenceManipMethods().codon_usage(self.nt_seq)
        return CAI
    
    @property
    def codon_sensitivity(self):
        cs, CAI, cc = SequenceManipMethods().codon_usage(self.nt_seq)
        return cs
       
        

    @property
    def codons(self):
        codons = self.nt_seq.upper()
        return [codons[i:i+3] for i in range(0, len(codons), 3)]

    @property
    def ktrna_id(self):
        return PropensityFactory().get_trna_ids(self.nt_seq)
        
    @property    
    def kelong(self):
        return PropensityFactory().get_k(self.nt_seq,self.ki,self.ke_mu,self.kt)[1:-1]
    
    @property
    def probe_vec(self):
        pv = np.zeros( (len(list(self.tag_epitopes)), self.total_length))
        for i in range(len(list(self.tag_epitopes))):
            pv[i,[self.tag_epitopes[list(self.tag_epitopes.keys())[i]]]] = 1
        pv = np.cumsum(pv,axis=1)
        return pv
    
    @property
    def probe_loc(self):
        pv = np.zeros( (len(list(self.tag_epitopes)), self.total_length))
        for i in range(len(list(self.tag_epitopes))):
            pv[i,[self.tag_epitopes[list(self.tag_epitopes.keys())[i]]]] = 1      
        return pv        


    
    @property
    def all_k(self):
        return PropensityFactory().get_k(self.nt_seq,self.ki,self.ke_mu,self.kt)   
    
    
    def get_binned_vectors(self,nbins,strategy='intellegent',min_binsize=3):
        if strategy == 'intellegent':    
            inds = PropensityFactory().intellegent_bin(self.probe_loc,nbins,min_bin=min_binsize)
            
        if strategy == 'even':
            inds = PropensityFactory().even_bin(len(self.kelong),nbins)
            
        bin_k = PropensityFactory().bin_k(self.kelong,inds)  
        probeloc_binned,probevec_binned = ProbeVectorFactory().bin_probe_vecs(self.probe_loc,inds)
                
        return inds,  probeloc_binned, probevec_binned,bin_k
    
    
    def generate_3frame_tags(self):
        
        #multiframe_epitopes = {}
        codons_seq = ''
        self.multiframe_epitopes = []
        self.multiframe_nt_seq = []
        self.multiframe_aa_seq = []
        
        for n in range(3):
            if n == 0:
                codons_seq = codons_seq + self.nt_seq
                self.multiframe_nt_seq.append(self.nt_seq)
            else:
                codons_seq = codons_seq + self.nt_seq[n:-(3-n)]
                self.multiframe_nt_seq.append(self.nt_seq[n:-(3-n)])
        
        cd = CodonDictionaries()
        smm = SequenceManipMethods(codons_seq)
        codons_aa_seq = smm.nt2aa(codons_seq)
        
        for frame in self.multiframe_nt_seq:
            self.multiframe_aa_seq.append(smm.nt2aa(frame))
        
        for frame in self.multiframe_aa_seq:
            multiframe_epitopes = {}
            for tag in cd.tag_dict.keys():
          
                if cd.tag_dict[tag] in frame:
                  
                    tag_detected = True
                
                    epi = smm.get_tag_loc(frame,cd.tag_dict[tag] )
                    
                    multiframe_epitopes[tag] = epi
                
            self.multiframe_epitopes.append(multiframe_epitopes)

                
                

    
    def visualize_probe(self):
        probe = self.probe_loc
        fig,ax = plt.subplots(1)
        N = len(self.kelong)
        ncolors = probe.shape[0]
        
        cmap = cm.get_cmap('gist_rainbow')
        colors = cmap(np.linspace(.01,.95, ncolors))
        
        rectangle =  mpatches.Rectangle((0,.1), N ,.8,linewidth=1,edgecolor='k',facecolor='lightgray')

        ax.add_patch(rectangle)
        
        color = np.where(probe == 1)[0]
        location = np.where(probe == 1)[1]
        
        colorlabels = ['Color %d'% i for i in range(ncolors)    ]
        
        for c in range(ncolors):
            ax.plot([-10,-10],[-4,-4],color = colors[c]  )  #fix the legend colors
        
        
        for c,loc in zip(color,location):
            ax.plot([loc,loc],[.1,.9],color = colors[c]  )
            
        ax.set_ylim([-.1,10])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xlabel('codon')
        ax.axes.legend(colorlabels,loc=7)
        ax.axes.get_yaxis().set_visible(False)
        
        ax.text(0,5,'Transcript Name: %s' % self.name)
        ax.text(0,4,'Total Length: %d codons' % self.total_length)
        ax.text(0,3,'Seq: %s ...' % self.aa_seq[:10])
        
        fig.show()    
        
        
        
    


class ODE_Soln():
    def __init__(self):
        
        self.__arrays = [    'probe_loc',
                             'x0',
                             'k',
                             'N',
                             'time',
                             'mu_It',
                             'intensity_acc',
                             'intensity_acc_norm',]
        
        self.__arrays_2d = ['var_state_ss']
        
        self.__vals = ['fi','kb','N','solve_time','mu_I_ss','var_I_ss']
        self.__meta = GenericMetaData().get()
        
        
        pass
    
    def save(self,filename):
        ext = filename.split('.')[-1]
        
        if 'txt' == ext:
            self.__save_txt(filename)
        if 'json' == ext:
            self.__save_json(filename)
            

    def load(self,filename):
        ext = filename.split('.')[-1]
        
        if 'txt' == ext:
            self.__load_from_txt(filename)
        if 'json' == ext:
            self.__load_from_json(filename)
            
    def __save_txt(self,filename):

        if '.txt' not in filename:    
            filename = filename + '.txt'


      
        f = open(filename,'w')
        for key in self.__dict__.keys():
            
            if key in self.__arrays + self.__arrays_2d + self.__vals:

                f.write((key + '\r\n'))
                np.savetxt(f, np.atleast_2d(self.__dict__[key]), delimiter=',', fmt='%s')
                f.write(('\r\n'))
                
        f.close()

    def __load_from_txt(self, filename):
        if '.txt' in filename:
            ode_obj = np.loadtxt(filename, dtype=str,delimiter='\n')
          
            solutions = []
            for i in range(0,len(ode_obj)-1):
                label = ode_obj[i]                
                
                if label in self.__arrays:
                    
                    array = np.fromstring(ode_obj[i+1], dtype=float, sep=',')
                    exec(('self.'+label+ '=array'))  
                
                if label in self.__vals:
                    value = np.fromstring(ode_obj[i+1], dtype=float, sep=',')
                    exec(('self.'+label+'=value'))
                    
                if label in self.__arrays_2d:
                    j = 1
                    array2d = np.array([[]])
                    while ode_obj[i+j] not in self.__arrays + self.__vals:
                        if j == 1:
                            array2d = np.array([np.fromstring(ode_obj[i+j], dtype=float, sep=',')   ])
                        else:
                            str_2d = np.fromstring(ode_obj[i+j], dtype=float, sep=',')    
                            array2d = np.vstack((array2d,str_2d))
                        
                        j+=1
                        
                   
                    exec(('self.'+label+'=array2d'))
                        
                


    def __save_from_json(self, filename):

        if '.json' not in filename:
            filename =  filename + '.json'

        odedict = {}
        for key in self.__dict__.keys():
            if key in self.__arrays:
                
                odedict[key] = self.ssa_harr.__dict__[key].tolist()
            else:
                odedict[key] = self.ssa_harr.__dict__[key]

        json.dump(odedict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)


    def __load_json(self,filename):
        if '.json' in filename:

            obj_text = codecs.open(filename, 'r', encoding='utf-8').read()
            ode_dict = json.loads(obj_text)

            for key in ode_dict.keys():
                if key in self.__arrays:

                    self.__dict__[key] = np.array(ode_dict[key])
                else:
                    self.__dict__[key] = ode_dict[key]
           


class SSA_Soln():
    '''
    SSA container class

    holds intensity / ribosome data as well as the propensities used

    __.n_traj = number of trajectories
    __.k = propensities used for the simulation
    __.rib_density = ribosome density per mRNA strand
    __.ribosome_means


    '''
    def __init__(self):
        self.n_traj = 0   #number trajectories
        self.k = []       #propensities
        self.no_rib_per_mrna = 0    #number of ribosomes per mrna strand
        self.rib_density = 0      #ribosome density
        self.ribosome_means = 0  #mean ribosomes
        self.rib_vec = 0          #actual vector of ribosome locations
        self.intensity_vec = []   #intensity vectors per SSA trajectory
        self.time_vec_fixed = []   #time array
        self.__meta = GenericMetaData().get()



    def save(self,filename):
        ext = filename.split('.')[-1]
        
        if 'txt' == ext:
            self.__save_txt(filename)
        if 'json' == ext:
            self.__save_json(filename)
            

    def load(self,filename):
        ext = filename.split('.')[-1]
        
        if 'txt' == ext:
            self.__load_from_txt(filename)
        if 'json' == ext:
            self.__load_from_json(filename)

    def __save_txt(self,filename):

        if '.txt' in filename:
            f = open(filename, 'a')
            for key in self.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    f.write((key + '\r\n'))
                    np.savetxt(f, np.atleast_2d(self.__dict__[key]), delimiter=',', fmt='%s')
                    f.write(('\r\n'))

        else:
            filename = filename + '.txt'
            f = open(filename,'a')
            for key in self.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    f.write((key + '\r\n'))
                    np.savetxt(f, np.atleast_2d(self.__dict__[key]), delimiter=',', fmt='%s')
                    f.write(('\r\n'))
        f.close()
        

    def __load_from_txt(self, filename):
        if '.txt' in filename:
            ssa_obj = np.loadtxt(filename, dtype=str,delimiter='\n')
            solutions = []
            for i in range(0,len(ssa_obj)-1):
                label = ssa_obj[i]
                
                
                if label in ['rib_means',
                             'rib_vec',
                             'n_traj',
                             'start_time',
                             'k',
                             'time_vec_fixed',
                             'dwelltime',
                             'mean_autocorr',
                             'no_rib_per_mrna',
                             'ke_sim',
                             'autocorr_vec',
                             'ribosome_means',
                             'error_autocorr',
                             'rib_density',
                             'time',
                             'ke_true',
                             'evaluating_inhibitor',
                             'time_inhibit',
                             'evaluating_frap']:

                    if label in ['start_time', 'no_rib_per_mrna', 'ke_sim', 'dwelltime','ke_true','time_inhibit']:
                        
                        array = np.fromstring(ssa_obj[i+1], dtype=float, sep=',')[0]
                        exec(('self.'+label+ '=array'))
                    elif label in ['n_traj']:
                        array = int(np.fromstring(ssa_obj[i+1], dtype=float, sep=',')[0])
                        exec(('self.'+label+ '=array'))
                    else:
                        array = np.fromstring(ssa_obj[i+1], dtype=float, sep=',')
                        exec(('self.'+label+ '=array'))

                if label in ['evaluating_inhibitor','evaluating_frap']:
                    if 'False' in ssa_obj[i+1]:                         
                        exec(('self.'+label+ '=False'))
                    if 'True' in ssa_obj[i+1]:                         
                        exec(('self.'+label+ '=True'))


            for i in range(0,len(ssa_obj)-1):
                label = ssa_obj[i]                    
                    
                if label == 'intensity_vec':

                    tvec = self.time_vec_fixed[np.where(self.time_vec_fixed >= self.start_time)]
                    i_vec = np.zeros((self.n_traj, len(self.time)))

                    for j in range(self.n_traj):
                        array = np.fromstring(ssa_obj[i+j+1], dtype=float,sep=',')
                        i_vec[j] = array
                        
                    exec(('self.'+label+ '=i_vec'))
                     
                if label == 'solutions':    
                    for j in range(self.n_traj):
                        array = np.fromstring(ssa_obj[i+j+1], dtype=float,sep=',')
                        solutions.append(array)
                        
                    exec(('self.'+label+ '=solutions'))
                        
                  
                    
                    
    def make_dict(self):
        ssadict = {}
        for key in self.__dict__.keys():
            print(key)
            if key != 'rib_vec' and key != 'ribosome_means':
                try:
                    ssadict[key] = self.__dict__[key].tolist()
                except:
                    ssadict[key] = self.__dict__[key]
                
            if key == 'col_points':
                col_pt = [x.tolist() for x in self.__dict__[key] ] 
                ssadict[key] = col_pt                
                    
                    
        return ssadict          
                     
                     

    def __save_json(self, filename):

        if '.json' in filename:

            ssadict = {}
            for key in self.__dict__.keys():
               
                #if key != 'rib_vec' and key != 'ribosome_means':
                try:
                    ssadict[key] = self.__dict__[key].tolist()
                except:
                    ssadict[key] = self.__dict__[key]
                    
                if key == 'col_points':
                    col_pt = [x.tolist() for x in self.__dict__[key] ] 
                    ssadict[key] = col_pt
                    
            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

        else:
            filename =  filename + '.json'

            ssadict = {}
            for key in self.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa_harr.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa_harr.__dict__[key]

            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)


    def __load_from_json(self,filename):
        if '.json' in filename:

            obj_text = codecs.open(filename, 'r', encoding='utf-8').read()
            ssadict = json.loads(obj_text)

            for key in ssadict.keys():
                if key in ['all_trna_results','rib_means','time_vec_fixed','mean_autocorr','autocorr_vec','error_autocorr','rib_density','intensity_vec','I']:

                    self.__dict__[key] = np.array(ssadict[key])
                    
                elif key in ['colpoints']: 
                    
                    cpts = [np.array(x) for x in ssadict[key]]
                    self.__dict__[key] = cpts
                else:
                    self.__dict__[key] = ssadict[key]



class GenericMetaData():
    def __init__(self):
        self.id = ''
        self.created_at = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime( time.time() ))
        self.user = os.path.expanduser("~")
        self.rss_version = rSNAPsim().__version__
        self.platform = platform.platform()
        self.python_version = sys.version
        
    def get(self):
        return self.__dict__
        
        

class GenericSSA():

    '''
    Generic SSA solver - used for the simulated cell animations
    '''

    def __init__(self,type='linear'):

        self.time_variant = False
        self.xi=np.array([])
        self.ti= None
        self.tf=None
        self.S=np.array([])
        self.type=type
        self.ptimes=100
        self.params={}
        if type=='linear':
            #self.fast_rxn = 0.5
            self.W0=np.array([])
            self.W1=np.array([])
        if type == 'nonlinear':
            #self.fast_rxn = 0.5
            self.P=lambda x,t:None



    def gettvec(self):
        return np.linspace(self.ti,self.tf,self.ptimes)

    def updateparams(self):
        for k, v in self.params.iteritems():
            setattr(self, k, v)

    def _run_trajectory(self):#!!!!!!!!!!!!!!!!!renamed run to solve(big deal)
        x=self.xi
        t=self.ti
        __n=len(x)
        self.time=self.gettvec()
        data=np.zeros((len(self.xi),self.ptimes))
        ip=0

        if self.type=='linear':
            if self.time_variant == False:
                while t<self.tf:
                    rate=np.atleast_2d(np.dot(self.W1,x))+self.W0
                    rate=np.cumsum(rate)
                    with np.errstate(divide='ignore', invalid='ignore'):
                        t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()
                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            data[:,ip]=x.reshape(__n)
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i
                            break
                    x=x+np.atleast_2d(self.S[:,event]).T

            else:


                x = np.concatenate((np.array([0]),self.xi.flatten()))

                #x = self.xi
                __n=len(x)
                self.time=self.gettvec()
                data=np.zeros((len(self.xi),self.ptimes))
                a,b = self.S.shape
                S = np.vstack((np.zeros(b),self.S))
                S = np.hstack((np.zeros((a+1,1)),S))
                while t<self.tf:
                    __n=len(x)
                    self.time=self.gettvec()
                    data=np.zeros((len(self.xi),self.ptimes))
                    a,b = self.S.shape
                    S = np.vstack((np.zeros(b),self.S))
                    S = np.hstack((np.zeros((a+1,1)),S))
                    while t<self.tf:
                        trate=self.get_P(x[1:],t)

                        rate = np.concatenate((np.array([self.fast_rxn]),trate))
                        rate=np.cumsum(rate)


                        t=(t-np.log(np.random.rand(1))/rate[-1])
                        ro=rate[-1]*np.random.rand()

                        while t>self.time[ip]:
                            if t>self.tf:
                                b = len(self.time[ip:])
                                fill = np.repeat(x[1:],b)
                                data[:,ip:]=fill.reshape(__n-1,b)
                                return data
                            else:
                                #data[:,ip]=x.reshape(__n)
                                data[:,ip]=x[1:]
                                ip=ip+1
                        for i in range(len(rate)):
                            if rate[i]>=ro:
                                event=i

                                break

                        x=x+S[:,event].ravel()
                    '''
                    rate=np.atleast_2d(np.dot(self.W1(t),x))+self.W0(t)

                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    print(t)
                    ro=rate[-1]*np.random.rand()
                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            data[:,ip]=x.reshape(__n)
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i
                            break
                    x=x+np.atleast_2d(self.S[:,event]).T

                    '''



        elif self.type=='nonlinear':
            if self.time_variant == True:  #if time variant use fast reaction
                x = np.concatenate((np.array([0]),self.xi.flatten()))

                #x = self.xi
                __n=len(x)
                self.time = self.gettvec()
                data = np.zeros((len(self.xi), self.ptimes))
                a, b = self.S.shape
                S = np.vstack((np.zeros(b), self.S))
                S = np.hstack((np.zeros((a+1, 1)), S))
                while t < self.tf:
                    trate=self.get_P(x[1:],t)
                    rate = np.concatenate((np.array([self.fast_rxn]),trate))
                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()

                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x[1:],b)
                            data[:,ip:]=fill.reshape(__n-1,b)
                            return data
                        else:
                            #data[:,ip]=x.reshape(__n)
                            data[:,ip]=x[1:]
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i

                            break

                    x=x+S[:,event].ravel()

            else:   #if not time variant ignore fast reaction

                x = self.xi.flatten()

                #x = self.xi
                __n=len(x)
                self.time=self.gettvec()

                while t<self.tf:
                    rate=self.get_P(x,t)

                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()

                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            #data[:,ip]=x.reshape(__n)
                            data[:,ip]=x
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i

                            break

                    x=x+self.S[:,event].ravel()



        else:
            'Error'
        self.data=data
        return data

    def _solve(self,n):
        __data=np.zeros((len(self.xi),self.ptimes,n))

        for i in range(n):
            __d=self._run_trajectory()
            __data[:,:,i]=__d
        self.data = __data

        return __data

    def setpar(self,key,val):
        self.params[key]=val

    def get_dist(self,specID=0):
        '''
        build distribution (non-normalized and pdf)
        of rna for the model)
        '''
        n_specs, n_times, n_traj = self.data.shape
        n_specs = int(n_specs)
        n_times = int(n_times)
        n_traj = int(n_traj)
        specID = int(specID)
        max_rna = int(np.max(self.data[specID,:,:]))
        self.pdf = np.zeros((n_times,max_rna+1))
        self.fdist = np.zeros((n_times,max_rna+1))
        for i in range(n_times):
            ind = int(i)
            for j in range(n_traj):
                jnd = int(j)
                self.fdist[ind,int(self.data[int(specID),ind,jnd])] +=1
            self.pdf[ind,:] = self.fdist[ind,:] / np.sum(self.fdist[ind,:])

    def get_traj(self,specID=0,ntraj='all'):
        n_specs, n_times, n_traj = self.data.shape
        n_specs = int(n_specs)
        n_times = int(n_times)




        if isinstance(specID,int):
            if ntraj == 'all':
                return self.data[specID],ntraj
            else:
                try:
                    ntraj = ntraj.flatten().astype(int).tolist()
                except:
                    ntraj = int(ntraj)
                    pass

                return self.data[specID][:,ntraj],ntraj
        else:
            if specID == 'all':
                if ntraj == 'all':
                    return self.data,ntraj
                else:

                    try:
                        ntraj = ntraj.flatten().astype(int).tolist()
                    except:
                        pass

                    return self.data,ntraj

    def get_means(self,specID=0):
        '''
        get the first moment.
        '''
        n_specs, n_times, n_traj = self.data.shape
        max_rna = int(np.max(self.data[int(specID),:,:])+1)
        self.means = np.zeros(n_times)
        for i in range(n_times):
            self.means[i] = np.sum(np.arange(max_rna)*self.pdf[i,:])

    def get_variances(self,specID=0):
        '''
        get the second moment.
        '''
        self.get_dist()
        n_specs, n_times, n_traj = self.data.shape
        max_rna = int(np.max(self.data[specID,:,:])+1)
        self.variances = np.zeros(n_times)
        self.covariances = np.zeros((n_specs,n_specs,n_times))
        for i in range(n_times):
            self.variances[i] = np.sum((np.arange(max_rna)**2)*self.pdf[i,:])-(np.sum(np.arange(max_rna)*self.pdf[i,:])**2)
            self.covariances[:,:,i] = np.cov(self.data[:,i,:])

    def return_all_var(self):
        all_members = self.__dict__.keys()

        return [ (item, self.__dict__[item]) for item in all_members if not item.startswith("_")]

    def return_names(self):
        all_members = self.__dict__.keys()
        return [ item for item in all_members if not item.startswith("_")]





#     def ssa_binned(self,nt_seq=None, bins = 50,all_k=None, k_elong_mean=10, k_initiation=.03, probePosition=[], n_traj=100, tf=1000, start_time=0, tstep=1000, time_inhibit=0, evaluating_frap=False, evaluating_inhibitor=False,force_python = False):
        
#         if nt_seq == None:                  #get sequence if none was passed
#             nt_seq = self.POI.nt_seq
#         genelength = int(len(nt_seq)/3)
        
#         if len(probePosition) == 0:
#             pv,probePosition = self.get_probvec()
                    

#         if all_k == None:          # build the k vector if one was not provided

#             codons = nt_seq
#             genelength = int(len(codons)/3)
#             seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
#             k_elongation = np.zeros((1, genelength))
#             tRNA_copynumber = np.zeros((1, genelength))

#             for i in range(len(seperated_codons)):
#                 tRNA_copynumber[0, i] = self.strGeneCopy[seperated_codons[i]]

#             mean_tRNA_copynumber = np.mean(list(self.strGeneCopy.values()))

#             k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean
#             all_k = [k_initiation] + k_elongation.flatten().tolist()[:-1] + [10]        
        

#         kbin,klen = self.get_binned_k(k_elongation.flatten()[:-1],bins)
#         all_k =     [k_initiation] + kbin.flatten().tolist()  # 
            
            
#         pv,probePosition = self.get_binned_probe_vec(probePosition,bins)

        
#         footprint = 0
#         if isinstance(probePosition,list):
#             probePosition = np.array([probePosition]).astype(int)
            
#         ssa_obj = self.__solve_ssa(genelength, all_k,pv,probePosition,n_traj, tf, start_time, tstep, time_inhibit, evaluating_frap, evaluating_inhibitor,force_python,footprint)
#         return ssa_obj




#     def ssa_solver(self, nt_seq=None, all_k=None, k_elong_mean=10, k_initiation=.03, probePosition=[], n_traj=100, tf=1000, start_time=0, tstep=1000, time_inhibit=0, evaluating_frap=False, evaluating_inhibitor=False,force_python = False):
#         '''
#         Solve stochastic simulation algorithms (SSA) for the translation simulation.

#         *keyword args*

#             **nt_seq**, nucleotide sequence to simulate

#             **all_k**, the propensity rates for each codon location (obtained via get_k)

#             **k_elong_mean**, average elongation rate to normalize by

#             **k_initiation**, rate of mRNA translation initiation

#             **probePosition**, binary vector of probe positions, i.e. where the tag epitopes start by codon position

#             **n_traj**, number of trajectories

#             **tf**, final time point

#             **tstep**, number of time steps to record from 0 to tf

#             **time_inhibit**, inhibition time of translation either, harringtonine assay or FRAP

#             **evaluating_frap**, true or false for evaluating frap assay at time_inhibit

#             **evaluating_inhibitor**, true or false for evaluating harringtonine at time_inhibit

#         *returns*

#             **ssa_obj**, a ssa() class containing the raw ribosome posistions simulated and statistics such as intensity vectors from the SSA trajectory group

#         '''

#         if len(probePosition) == 0:
#             '''
#             try:
#                 probePosition = []
#                 for key in self.POI.tag_epitopes.keys():
#                     probePosition = probePosition + self.POI.tag_epitopes[key]
#                 probePosition = np.unique(probePosition).tolist()
#             except:
#                 print('No POI found')
#                 #nt_seq = self.tag_full['T_flag'] + nt_seq
#             '''
            
#             pv,probePosition = self.get_probvec()
            
        

#         if nt_seq == None:
#           nt_seq = self.POI.nt_seq
#         genelength = int(len(nt_seq)/3)
        

#         if all_k == None:


#             codons = nt_seq
#             genelength = int(len(codons)/3)
#             seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
#             k_elongation = np.zeros((1, genelength))
#             tRNA_copynumber = np.zeros((1, genelength))

#             for i in range(len(seperated_codons)):
#                 tRNA_copynumber[0, i] = self.strGeneCopy[seperated_codons[i]]

#             mean_tRNA_copynumber = np.mean(list(self.strGeneCopy.values()))

#             k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean
#             all_k = [k_initiation] + k_elongation.flatten().tolist()[:-1] + [10]

        
#         if isinstance(probePosition,list):
#             probePosition = np.array([probePosition]).astype(int)
            
#         footprint = 9
#         ssa_obj = self.__solve_ssa(genelength,  all_k,pv,probePosition,n_traj, tf, start_time, tstep, time_inhibit, evaluating_frap, evaluating_inhibitor,force_python, footprint)
#         return ssa_obj       




            
            
#     def __solve_ssa(self,genelength,all_k,pv,probePosition,n_traj, tf, start_time, tstep, time_inhibit, evaluating_frap, evaluating_inhibitor,force_python,footprint):
     

#         non_consider_time = start_time
      
#         '''
#         if probePosition.shape[0] <= 1:
#             pv = np.zeros((1, genelength+1)).astype(int).flatten()
            
#             for i in range(len(probePosition[0])):
#                 pv[probePosition[0][i]:] = i+1
#         else:
#             pv = np.zeros((probePosition.shape[0], genelength+1)).astype(int)
#             for j in range(probePosition.shape[0]):
#                 for i in range(len(probePosition)):
#                     pv[j][probePosition[j][i]:] = i+1      
#         '''

#         npoints = tstep #non_consider_time + tstep
        
#         time_vec_fixed = np.linspace(0, npoints-1, npoints, dtype=np.float64)
#         truetime = np.linspace(0, tf, tstep, dtype=np.float64)


#         rib_vec = []

#         solutions = []
        


        

#         evf = int(evaluating_frap)
#         evi = int(evaluating_inhibitor)
#         try:
#             intime = float(time_inhibit)
#         except:
#             intime = 0

# #        if evaluating_frap == True or evaluating_inhibitor == True:
# #            for i in range(nRepetitions):
# #
# #                soln = self.SSA(all_k,time_vec_fixed,inhibit_time=time_inhibit+non_consider_time,FRAP=evaluating_frap,Inhibitor=evaluating_inhibitor)
# #                solutions.append(soln)
# #        else:

#         solutionssave = []
        
#         st = time.time() 
        
#         try:
#             if force_python == True:
#                 st[0]
                
#             rib_vec = []
    
#             solutions = []            
#             solutionssave = []
#             N_rib = 200
#             all_results = np.zeros((n_traj, N_rib*len(time_vec_fixed)), dtype=np.int32)
#             all_ribtimes = np.zeros((n_traj,int(1.3*all_k[0]*truetime[-1])),dtype=np.float64)
#             result = np.zeros((len(time_vec_fixed)*N_rib), dtype=np.int32)
#             nribs = np.array([0],dtype=np.int32)
#             k = np.array(all_k)
#             seeds = np.random.randint(0, 0x7FFFFFF, n_traj)
#             all_frapresults = np.zeros((n_traj,N_rib*len(time_vec_fixed)),dtype=np.int32)
#             all_collisions = np.zeros((n_traj,int(1.3*all_k[0]*truetime[-1])),dtype=np.int32)
#             all_nribs = np.zeros((n_traj,1))
#             all_col_points = []
#             x0 = np.zeros((N_rib),dtype=np.int32)
#             for i in range(n_traj):
#                 result = np.zeros((len(time_vec_fixed)*N_rib), dtype=np.int32)
#                 ribtimes = np.zeros((int(1.3*k[0]*truetime[-1])),dtype=np.float64)
#                 frapresult = np.zeros((len(time_vec_fixed)*N_rib),dtype=np.int32)
#                 coltimes = np.zeros((int(1.3*k[0]*truetime[-1])),dtype=np.int32)
#                 colpointsx = np.zeros(len(k[1:-1])*(int(1.3*k[0]*truetime[-1])),dtype=np.int32)
#                 colpointst = np.zeros(len(k[1:-1])*(int(1.3*k[0]*truetime[-1])),dtype=np.float64)
                
#                 nribs = np.array([0],dtype=np.int32)
               
#                 ssa_translation.run_SSA(result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs,x0,footprint)
#                 #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
#                 all_results[i, :] = result
#                 all_frapresults[i,:] = frapresult
#                 all_ribtimes[i,:] = ribtimes
#                 all_collisions[i,:] = coltimes
#                 all_nribs[i,:] = nribs
                
#                 endcolrec = np.where(colpointsx == 0)[0][0]
                
#                 colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
#                 all_col_points.append(colpoints.T)
                
                
    
#             for i in range(n_traj):
#                 soln = all_results[i, :].reshape((N_rib, len(time_vec_fixed)))
#                 validind = np.where(np.sum(soln,axis=1)!=0)[0]
#                 if np.max(validind) != N_rib-1:
#                     validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)
            
#                 so = soln[(validind,)]
                
#                 solutionssave.append(so)
#                 solutions.append(soln)
            
#             collisions = np.array([[]])
#             watched_ribs = []
#             for i in range(n_traj):
#                 totalrib = all_nribs[i]
            
#                 if totalrib > all_collisions.shape[1]:
#                     collisions = np.append(collisions, all_collisions[i][:])
#                     watched_ribs.append(int(all_collisions.shape[1]))
            
#                 else:
                   
#                     collisions = np.append(collisions, all_collisions[i][:int(totalrib[0])])
#                     watched_ribs.append(int(totalrib[0]))
            
#             sttime = time.time() - st

#         except:
            
#             print('C++ library failed, Using Python Implementation')
#             rib_vec = []
    
#             solutions = []            
#             solutionssave = []
#             N_rib = 200
#             collisions = np.array([[]])
#             all_results = np.zeros((n_traj, N_rib*len(time_vec_fixed)), dtype=np.int32)
#             all_col_points = []
#             watched_ribs = []
#             for i in range(n_traj):
                
#                 soln,all_ribtimes,Ncol,col_points = self.SSA(all_k, truetime, inhibit_time=time_inhibit+non_consider_time, FRAP=evaluating_frap, Inhibitor=evaluating_inhibitor)
#                 #soln = soln.reshape((1, (len(time_vec_fixed)*N_rib)))
                
#                 collisions = np.append(collisions,Ncol)
#                 watched_ribs.append(int(len(collisions)))
#                 validind = np.where(np.sum(soln,axis=1)!=0)[0]
#                 all_col_points.append(np.array(col_points))
#                 if np.max(validind) != N_rib-1:
#                     validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)
                
#                 so = soln[(validind,)]
               
#                 solutionssave.append(so)

#                 solutions.append(soln)
            
#                 result = soln.reshape((1, (len(time_vec_fixed)*N_rib)))
#                 all_results[i, :] = result
            
#             sttime = time.time() - st


#                 #rb = sparse.lil_matrix((len(time_vec_fixed),genelength),dtype=int)
#                 #for j in range(soln.shape[1]):

#                     #if len(np.where(soln[:,j]!=0)[0]) !=0:
#                     #print(np.where(soln[:,j]!=0)[0])


#                     #rb[j,np.where(soln[:,j]!=0)[0]] = 1


#                         #for value in soln[:,j][np.where(soln[:,j]!=0)[0]].astype(int):

#                             #rb[j, value-1] = 1

#                 #rib_vec.append(rb)

        



#         no_ribosomes = np.zeros((n_traj, (genelength+1)))
        
#         startindex = np.where(truetime >= non_consider_time)[0][0]
        
#         #all_results = all_results[:,startindex*N_rib:]

#         for i in range(len(solutions)):
#             for j in range(len(solutions[0][0][startindex:])):
#                 rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]
               
#                 no_ribosomes[i, rib_pos.astype(int)] += 1
#         no_ribosomes = no_ribosomes[:, 1:]

#         ribosome_means = np.mean(no_ribosomes, axis=0)
#         ribosome_density = ribosome_means/npoints

#         no_ribosomes_per_mrna = np.mean(no_ribosomes)
        
 

#         if probePosition.shape[0] <=1:
#             I = np.zeros((n_traj, len(time_vec_fixed[startindex:])))
         
            
#         else:
#             I = np.zeros((int(probePosition.shape[0]),n_traj, len(time_vec_fixed[startindex:])))
         

#         #I = np.zeros((1,tstep+1))
        
#         if evaluating_frap == False:
#             if probePosition.shape[0] <=1:
#                 for i in range(n_traj):
        
#                     traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T

#                     I[i, :] = np.sum(np.multiply(pv.flatten()[traj], traj>0), axis=1)[startindex:].T
#             else:
#                 for j in range(probePosition.shape[0]):
#                     for i in range(n_traj):
            
#                         traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
            
#                         I[j,i, :] = np.sum(pv[j][traj], axis=1)[startindex:].T                
    
    
#             intensity_vec = I
        
#         else:
#             fraptime = time_inhibit
         
            
#             inds = np.where(truetime > fraptime)

#             inds2 = np.where(truetime  < fraptime+20)
#             inds = np.intersect1d(inds,inds2)
#             endfrap = inds[-1]-1
            

            
#             for i in range(n_traj):
    
#                 traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
                
#                 nribs = np.sum(solutionssave[i][:,endfrap]!=0)
             
#                 #ribloc = solutionssave[i][:,endfrap]
                
#                 #adj_pv = pv[solutionssave[i][:,inds[-1]][:nribs]]
#                 frap_app = 20

#                 revI = self.get_negative_intensity(traj,genelength,pv,truetime,fraptime+start_time,fraptime+start_time+frap_app)
                

#                 I[i, :] = np.sum(pv[traj], axis=1)[startindex:].T
                             
#                 I[i,inds[0]:inds[0]+20] = 0
#                 #I[i,endfrap-startindex:] = np.sum(pv[traj],axis=1)[endfrap-startindex:].T

#                 I[i,inds[0]+frap_app:len(revI)+inds[0]+frap_app] = I[i,inds[0]+frap_app:len(revI)+inds[0]+frap_app] + revI
                
                
      
                
                
    
    
#             intensity_vec = I




#         ssa_obj = SSA_soln()
#         ssa_obj.no_ribosomes = no_ribosomes
#         ssa_obj.n_traj = n_traj
#         ssa_obj.k = all_k
#         ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
#         ssa_obj.rib_density = ribosome_density
#         ssa_obj.rib_means = ribosome_means
#         ssa_obj.rib_vec = rib_vec
#         ssa_obj.intensity_vec = intensity_vec
#         ssa_obj.time_vec_fixed = time_vec_fixed
#         ssa_obj.time = truetime
#         ssa_obj.time_rec = truetime[startindex:]
#         ssa_obj.start_time = non_consider_time
#         ssa_obj.watched_ribs = watched_ribs
#         try:
#             ssa_obj.col_points = all_col_points
#         except:
#             pass


#         ssa_obj.evaluating_inhibitor = evaluating_inhibitor
#         ssa_obj.evaluating_frap = evaluating_frap
#         ssa_obj.time_inhibit = time_inhibit
#         ssa_obj.solutions = solutionssave
#         ssa_obj.solvetime = sttime
#         ssa_obj.collisions = collisions
        
        
#         try:
#             ssa_obj.ribtimes = all_ribtimes[np.where(all_ribtimes > 0)]
#         except:
#             pass


#         #solt = solutions.T

#         fragmented_trajectories = []
#         fragtimes = []
#         maxlen = 0
    
#         fragmentspertraj= []
#         for k in range(n_traj):
#             ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
#             changes = ind[1:] - ind[:-1]
#             addindexes = np.where(changes > 0)[0]
#             subindexes = np.where(changes < 0)[0]
            
#             sub = solutions[k][:,1:] - solutions[k][:,:-1]
#             neutralindexes = np.unique(np.where(sub < 0)[1])
#             neutralindexes = np.setxor1d(neutralindexes, subindexes)
            
#             for index in neutralindexes:
#                 pre = solutions[k][:,index]
#                 post = solutions[k][:,index+1]
#                 changecount = 0
#                 while len(np.where(post - pre < 0)[0]) > 0:
    
#                     post = np.append([genelength],post)
#                     pre = np.append(pre,0)
                    
#                     changecount+=1
                
#                 for i in range(changecount):
#                     addindexes = np.sort(np.append(addindexes,index))
#                     subindexes = np.sort(np.append(subindexes,index))
                    
#                 changes[index] = -changecount
#                 ind[index] += changecount
             
                
#             for index in np.where(np.abs(changes)>1)[0]:
#                 if changes[index] < 0:
#                     for i in range(np.abs(changes[index])-1):
#                         subindexes = np.sort(np.append(subindexes,index))
#                 else:
#                     for i in range(np.abs(changes[index])-1):
#                         addindexes = np.sort(np.append(addindexes,index))   
                
#             truefrags = len(subindexes)
     
                
        
           
#             if len(subindexes) < len(addindexes):
#                 subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(len(truetime)-1)).astype(int))
                
            
#             fragmentspertraj.append(len(subindexes))
            
#             for m in range(min(len(subindexes),len(addindexes))):
#                 traj = solutions[k][:, addindexes[m]:subindexes[m]+1]
#                 traj_ind = changes[addindexes[m]:subindexes[m]+1]
                
#                 startind = ind[addindexes[m]]
#                 minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
#                 fragment = np.array([])
            
                    
                
#                 iterind = startind
                
#                 if subindexes[m]-addindexes[m] > 0:
#                     if len(minusloc) > 1:
#                         if m <= truefrags:
#                             for n in range(len(minusloc)-1):
#                                 iterind = iterind + min(0,traj_ind[minusloc[n]])
#                                 fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                                
                            
#                             fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                            
#                         else:
#                             for n in range(len(minusloc)-1):

#                                 iterind = iterind + min(0,traj_ind[minusloc[n]])
                                
#                                 fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                  
                                
#                             fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
          
                        
                    
#                     else:

#                         fragment = solutions[k][startind][addindexes[m]:subindexes[m]+1].flatten()
                   
                
                    
#                     fragtimes.append(addindexes[m]+1)
                       
                    
#                     fragmented_trajectories.append(fragment)
#                     #if m <= truefrags:
#                         #kes.append(genelength/truetime[len(fragment)])
            
#                     if len(fragment) > maxlen:
#                         maxlen = len(fragment)
                    
    
#             fragarray = np.zeros((len(fragmented_trajectories), maxlen))
#             for i in range(len(fragmented_trajectories)):
#                 fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
            
#         ssa_obj.fragments = fragarray
#         ssa_obj.fragtimes = fragtimes
#         ssa_obj.frag_per_traj = fragmentspertraj
#         ssa_obj.full_frags = truefrags
#         ssa_obj.all_results = all_results
        
#         if probePosition.shape[0] > 1:
#             for i in range(probePosition.shape[0]):
#                 if i > 0:
#                     autocorr_vec2, mean_autocorr2, error_autocorr2, dwelltime2, ke_sim2  = self.get_autocorr(intensity_vec[i], truetime, 0, genelength)
#                     autocorr_vec = np.vstack((autocorr_vec,autocorr_vec2))
#                     mean_autocorr = np.vstack((mean_autocorr,mean_autocorr2))
#                     error_autocorr = np.vstack((error_autocorr,error_autocorr2))
#                     dwelltime.append(dwelltime2)
#                     ke_sim.append(ke_sim2)
#                 else:
#                     autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.get_autocorr(intensity_vec[i], truetime, 0, genelength)
#                     autocorr_vec_norm, mean_autocorr_norm, error_autocorr_norm, dwelltime, ke_sim = self.get_autocorr_norm(intensity_vec[i], truetime, 0, genelength)
#                     dwelltime = [dwelltime]
#                     ke_sim = [ke_sim]
            
#         else:
#             autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.get_autocorr(intensity_vec, truetime, 0, genelength)
#             autocorr_vec_norm, mean_autocorr_norm, error_autocorr_norm, dwelltime, ke_sim = self.get_autocorr_norm(intensity_vec, truetime, 0, genelength)
            
#             acov,nacov = self.get_all_autocovariances(intensity_vec,truetime,genelength )              
        
#         ssa_obj.autocorr_vec = autocorr_vec
#         ssa_obj.mean_autocorr = mean_autocorr
#         ssa_obj.error_autocorr = error_autocorr
#         ssa_obj.autocorr_vec_norm = autocorr_vec_norm
#         ssa_obj.mean_autocorr_norm = mean_autocorr_norm
#         ssa_obj.error_autocorr_norm = error_autocorr_norm
#         ssa_obj.dwelltime = dwelltime
#         ssa_obj.ke_sim = ke_sim
#         ssa_obj.ke_true = float(genelength)/np.mean(ssa_obj.ribtimes)
#         ssa_obj.probe = probePosition
        
        
#         try:
#             ssa_obj.autocovariance_dict  = acov
#             ssa_obj.autocovariance_norm_dict = nacov
#         except:
#             pass

#         return ssa_obj
        
    # def get_negative_intensity(self,solution,gene_length,pv,tvec,ti,stop_frap):
        
    #     startindex = np.where(tvec >= ti)[0][0]
    #     stop_frap = np.where(tvec >= stop_frap)[0][0]
      
            
    #     solution = solution.T
    #     fragmented_trajectories = []
    #     fragtimes = []
    #     endfragtimes = []
    #     maxlen = 0
        
    #     fragmentspertraj= []
  
    #     ind = np.array([next(j for j in range(0,solution.shape[0]) if int(solution[j, i]) == 0 or int(solution[j, i]) == -1) for i in range(0, solution.shape[1])])
    #     changes = ind[1:] - ind[:-1]
    #     addindexes = np.where(changes > 0)[0]
    #     subindexes = np.where(changes < 0)[0]
        
    #     sub = solution[:,1:] - solution[:,:-1]
    #     neutralindexes = np.unique(np.where(sub < 0)[1])
    #     neutralindexes = np.setxor1d(neutralindexes, subindexes)
        
    #     for index in neutralindexes:
    #         pre = solution[:,index]
    #         post = solution[:,index+1]
    #         changecount = 0
    #         while len(np.where(post - pre < 0)[0]) > 0:
    
    #             post = np.append([gene_length],post)
    #             pre = np.append(pre,0)
                
    #             changecount+=1
            
    #         for i in range(changecount):
    #             addindexes = np.sort(np.append(addindexes,index))
    #             subindexes = np.sort(np.append(subindexes,index))
                
    #         changes[index] = -changecount
    #         ind[index] += changecount
         
            
    #     for index in np.where(np.abs(changes)>1)[0]:
    #         if changes[index] < 0:
    #             for i in range(np.abs(changes[index])-1):
    #                 subindexes = np.sort(np.append(subindexes,index))
    #         else:
    #             for i in range(np.abs(changes[index])-1):
    #                 addindexes = np.sort(np.append(addindexes,index))   
            
    #     truefrags = len(subindexes)
     
            
    
       
    #     if len(subindexes) < len(addindexes):
    #         subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(len(tvec)-1)).astype(int))
            
        
    #     fragmentspertraj.append(len(subindexes))
        
    #     for m in range(min(len(subindexes),len(addindexes))):
    #         traj = solution[:, addindexes[m]:subindexes[m]+1]
    #         traj_ind = changes[addindexes[m]:subindexes[m]+1]
            
    #         startind = ind[addindexes[m]]
    #         minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
    #         fragment = np.array([])
        
                
            
    #         iterind = startind
            
    #         if subindexes[m]-addindexes[m] > 0:
    #             if len(minusloc) > 1:
    #                 if m <= truefrags:
    #                     for n in range(len(minusloc)-1):
    #                         iterind = iterind + min(0,traj_ind[minusloc[n]])
    #                         fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                            
                            
                            
              
            
                  
                        
    #                     fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                        
    #                 else:
    #                     for n in range(len(minusloc)-1):
    
    #                         iterind = iterind + min(0,traj_ind[minusloc[n]])
                            
    #                         fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
              
                            
    #                     fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
      
                    
                
    #             else:
    
    #                 fragment = solution[startind][addindexes[m]:subindexes[m]+1].flatten()
               
            
                
    #             fragtimes.append(addindexes[m]+1)
    #             if addindexes[m]+1  + len(fragment) > len(tvec):
    #                 endfragtimes.append(len(tvec))
    #             else:
    #                 endfragtimes.append(addindexes[m]+1  + len(fragment))
                   
                
    #             fragmented_trajectories.append(fragment)
    #             #if m <= truefrags:
    #                 #kes.append(genelength/truetime[len(fragment)])
        
    #             if len(fragment) > maxlen:
    #                 maxlen = len(fragment)
                
    
    #     fragarray = np.zeros((len(fragmented_trajectories), maxlen))
    #     for i in range(len(fragmented_trajectories)):
    #         fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
            
            
    #     affected_frags = []
    #     fragindexes = []
        
        
    #     for i in range(len(fragtimes)):
        
    #        if  np.sum([fragtimes[i]> np.array([startindex, stop_frap]), endfragtimes[i] > np.array([startindex, stop_frap])]) in [1,2,3]:
    #            affected_frags.append(i)
    #            fragindexes.append([fragtimes[i],endfragtimes[i]])
            
        
    #     #affected_frags = np.intersect1d(np.where(np.array(fragtimes) >=  startindex), np.where(np.array(fragtimes)<= stop_frap))
    #     if len(fragindexes)> 0:
    #         findexes = np.array(fragindexes)
  
    #         frange = findexes[:,1]-stop_frap
    #         afterfrapribs = findexes[np.where(frange > 0 )]
            
            
    #         relevantfrags = np.array(affected_frags)[np.where(frange > 0 )]
    #         if len(relevantfrags) > 0:
    #             cooked_ribs = 0#(len(affected_frags) - len(relevantfrags))*max(pv)
      
    #             stopfrapindex = stop_frap - afterfrapribs[:,0]
                
    #             rfrags = fragarray[relevantfrags]
    #             np.diag(rfrags[:,stopfrapindex])
    #             laglen = afterfrapribs[:,1] - stop_frap
    #             posistions_at_end_of_FRAP = np.diag(rfrags[:,stopfrapindex])
             
    #             offset = pv[posistions_at_end_of_FRAP.astype(int)]
             
    #             trailing_intensity = np.zeros((max(laglen)))
                
    #             for i in range(len(laglen)):
    #                 trailing_intensity[:laglen[i]] -= offset[i] 
                    
    #             trailing_intensity= trailing_intensity-cooked_ribs
    #         else:
    #             trailing_intensity = np.array([0])
    #     else:
    #         trailing_intensity = np.array([0])
        
    #     return trailing_intensity        



#     def ssa_solver_append(self, ssa_obj, n=100):

#         nRepetitions = ssa_obj.n_traj
#         all_k = ssa_obj.k
#         no_ribosomes_per_mrna = ssa_obj.no_rib_per_mrna
#         ribosome_density = ssa_obj.rib_density
#         ribosome_means = ssa_obj.rib_means
#         rib_vec = ssa_obj.rib_vec
#         intensity_vec = ssa_obj.intensity_vec
#         time_vec_fixed = ssa_obj.time_vec_fixed
#         non_consider_time = ssa_obj.start_time

#         evaluating_inhibitor = ssa_obj.evaluating_inhibitor
#         evaluating_frap = ssa_obj.evaluating_frap
#         time_inhibit = ssa_obj.time_inhibit
#         truetime = ssa_obj.time
#         tstep = len(ssa_obj.time)

#         npoints = tstep #non_consider_time + tstep

#         rib_vec = []

#         solutions = []
        
#         pv = ssa_obj.probe
    
#         genelength = len(pv[0])-1
        



#         evf = int(evaluating_frap)
#         evi = int(evaluating_inhibitor)
#         try:
#             intime = float(time_inhibit)
#         except:
#             intime = 0


#         solutionssave = []
        
#         st = time.time() 
#         n_traj = n
#         force_python = False
        
#         try:
#             if force_python == True:
#                 st[0]
                
#             rib_vec = []
    
#             solutions = []            
#             solutionssave = []
#             N_rib = 200
#             all_results = np.zeros((n_traj, N_rib*len(time_vec_fixed)), dtype=np.int32)
#             all_ribtimes = np.zeros((n_traj,int(1.3*all_k[0]*truetime[-1])),dtype=np.float64)
#             result = np.zeros((len(time_vec_fixed)*N_rib), dtype=np.int32)
#             nribs = np.array([0],dtype=np.int32)
#             k = np.array(all_k)
#             seeds = np.random.randint(0, 0x7FFFFFF, n_traj)
#             all_frapresults = np.zeros((n_traj,N_rib*len(time_vec_fixed)),dtype=np.int32)
#             all_collisions = np.zeros((n_traj,int(1.3*all_k[0]*truetime[-1])),dtype=np.int32)
#             all_nribs = np.zeros((n_traj,1))
#             all_col_points = []
#             for i in range(n_traj):
#                 result = np.zeros((len(time_vec_fixed)*N_rib), dtype=np.int32)
#                 ribtimes = np.zeros((int(1.3*k[0]*truetime[-1])),dtype=np.float64)
#                 frapresult = np.zeros((len(time_vec_fixed)*N_rib),dtype=np.int32)
#                 coltimes = np.zeros((int(1.3*k[0]*truetime[-1])),dtype=np.int32)
#                 colpointsx = np.zeros(len(k[1:-1])*(int(1.3*k[0]*truetime[-1])),dtype=np.int32)
#                 colpointst = np.zeros(len(k[1:-1])*(int(1.3*k[0]*truetime[-1])),dtype=np.float64)
                
#                 nribs = np.array([0],dtype=np.int32)
                
#                 ssa_translation.run_SSA(result, ribtimes, coltimes, colpointsx,colpointst, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
#                 #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
#                 all_results[i, :] = result
#                 all_frapresults[i,:] = frapresult
#                 all_ribtimes[i,:] = ribtimes
#                 all_collisions[i,:] = coltimes
#                 all_nribs[i,:] = nribs
                
#                 endcolrec = np.where(colpointsx == 0)[0][0]
                
#                 colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
#                 all_col_points.append(colpoints.T)
                
                
    
#             for i in range(n_traj):
#                 soln = all_results[i, :].reshape((N_rib, len(time_vec_fixed)))
#                 validind = np.where(np.sum(soln,axis=1)!=0)[0]
#                 if np.max(validind) != N_rib-1:
#                     validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)
            
#                 so = soln[(validind,)]
                
#                 solutionssave.append(so)
#                 solutions.append(soln)
            
#             collisions = np.array([[]])
#             watched_ribs = []
#             for i in range(n_traj):
#                 totalrib = all_nribs[i]
            
#                 if totalrib > all_collisions.shape[1]:
#                     collisions = np.append(collisions, all_collisions[i][:])
#                     watched_ribs.append(int(all_collisions.shape[1]))
            
#                 else:
                   
#                     collisions = np.append(collisions, all_collisions[i][:int(totalrib[0])])
#                     watched_ribs.append(int(totalrib[0]))
            
#             sttime = time.time() - st

#         except:
            
#             print('C++ library failed, Using Python Implementation')
#             rib_vec = []
    
#             solutions = []            
#             solutionssave = []
#             N_rib = 200
#             collisions = np.array([[]])
#             all_results = np.zeros((n_traj, N_rib*len(time_vec_fixed)), dtype=np.int32)
#             all_col_points = []
#             watched_ribs = []
#             for i in range(n_traj):
                
#                 soln,all_ribtimes,Ncol,col_points = self.SSA(all_k, truetime, inhibit_time=time_inhibit+non_consider_time, FRAP=evaluating_frap, Inhibitor=evaluating_inhibitor)
#                 #soln = soln.reshape((1, (len(time_vec_fixed)*N_rib)))
                
#                 collisions = np.append(collisions,Ncol)
#                 watched_ribs.append(int(len(collisions)))
#                 validind = np.where(np.sum(soln,axis=1)!=0)[0]
#                 all_col_points.append(np.array(col_points))
#                 if np.max(validind) != N_rib-1:
#                     validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)
                
#                 so = soln[(validind,)]
               
#                 solutionssave.append(so)

#                 solutions.append(soln)
            
#                 result = soln.reshape((1, (len(time_vec_fixed)*N_rib)))
#                 all_results[i, :] = result
            
#             sttime = time.time() - st


#                 #rb = sparse.lil_matrix((len(time_vec_fixed),genelength),dtype=int)
#                 #for j in range(soln.shape[1]):

#                     #if len(np.where(soln[:,j]!=0)[0]) !=0:
#                     #print(np.where(soln[:,j]!=0)[0])


#                     #rb[j,np.where(soln[:,j]!=0)[0]] = 1


#                         #for value in soln[:,j][np.where(soln[:,j]!=0)[0]].astype(int):

#                             #rb[j, value-1] = 1

#                 #rib_vec.append(rb)

        



#         no_ribosomes = np.zeros((n_traj, (genelength+1)))
        
#         startindex = np.where(truetime >= non_consider_time)[0][0]
        
#         #all_results = all_results[:,startindex*N_rib:]

#         for i in range(len(solutions)):
#             for j in range(len(solutions[0][0][startindex:])):
#                 rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]
               
#                 no_ribosomes[i, rib_pos.astype(int)] += 1
#         no_ribosomes = no_ribosomes[:, 1:]

#         ribosome_means = np.mean(no_ribosomes, axis=0)
#         ribosome_density = ribosome_means/npoints

#         no_ribosomes_per_mrna = np.mean(no_ribosomes)
        
 

#         if pv.shape[0] <=1:
#             I = np.zeros((n_traj, len(time_vec_fixed[startindex:])))
         
            
#         else:
#             I = np.zeros((int(pv.shape[0]),n_traj, len(time_vec_fixed[startindex:])))
         

#         #I = np.zeros((1,tstep+1))
        
#         if evaluating_frap == False:
#             if pv.shape[0] <=1:
#                 for i in range(n_traj):
        
#                     traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
        
#                     I[i, :] = np.sum(pv[0][traj], axis=1)[startindex:].T
#             else:
#                 for j in range(pv.shape[0]):
#                     for i in range(n_traj):
            
#                         traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
            
#                         I[j,i, :] = np.sum(pv[j][traj], axis=1)[startindex:].T                
    
    
#             intensity_vec = I
        
#         else:
#             fraptime = time_inhibit
         
            
#             inds = np.where(truetime > fraptime)

#             inds2 = np.where(truetime  < fraptime+20)
#             inds = np.intersect1d(inds,inds2)
#             endfrap = inds[-1]-1
            

            
#             for i in range(n_traj):
    
#                 traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
                
#                 nribs = np.sum(solutionssave[i][:,endfrap]!=0)
             
#                 #ribloc = solutionssave[i][:,endfrap]
                
#                 #adj_pv = pv[solutionssave[i][:,inds[-1]][:nribs]]
#                 frap_app = 20

#                 revI = self.get_negative_intensity(traj,genelength,pv,truetime,fraptime+start_time,fraptime+start_time+frap_app)
                

#                 I[i, :] = np.sum(pv[traj], axis=1)[startindex:].T
                             
#                 I[i,inds[0]:inds[0]+20] = 0
#                 #I[i,endfrap-startindex:] = np.sum(pv[traj],axis=1)[endfrap-startindex:].T

#                 I[i,inds[0]+frap_app:len(revI)+inds[0]+frap_app] = I[i,inds[0]+frap_app:len(revI)+inds[0]+frap_app] + revI
                
                
      
                
                
    
    
#             intensity_vec = I




#         new_ssa_obj = ssa()
#         new_ssa_obj.no_ribosomes = np.vstack(( ssa_obj.no_ribosomes , no_ribosomes))
#         new_ssa_obj.n_traj = n_traj+ssa_obj.n_traj
#         new_ssa_obj.k = all_k
#         new_ssa_obj.no_rib_per_mrna = float(n_traj)/(n_traj+ssa_obj.n_traj) * no_ribosomes_per_mrna  + float(ssa_obj.n_traj)/(n_traj+ssa_obj.n_traj) * ssa_obj.no_rib_per_mrna
#         new_ssa_obj.rib_density = ribosome_density
#         new_ssa_obj.rib_means = ribosome_means
        
#         new_ssa_obj.rib_means = np.mean(np.vstack((ssa_obj.rib_means,ribosome_means)),0)
        
        
#         new_ssa_obj.rib_vec = rib_vec
#         new_ssa_obj.intensity_vec = np.vstack((ssa_obj.intensity_vec,intensity_vec))
#         new_ssa_obj.time_vec_fixed = time_vec_fixed
#         new_ssa_obj.time = truetime
#         new_ssa_obj.time_rec = truetime[startindex:]
#         new_ssa_obj.start_time = non_consider_time
#         new_ssa_obj.watched_ribs = ssa_obj.watched_ribs + watched_ribs
#         try:
#             new_ssa_obj.col_points = ssa_obj.col_points + all_col_points
#         except:
#             pass


#         new_ssa_obj.evaluating_inhibitor = evaluating_inhibitor
#         new_ssa_obj.evaluating_frap = evaluating_frap
#         new_ssa_obj.time_inhibit = time_inhibit
#         new_ssa_obj.solutions = ssa_obj.solutions + solutionssave
#         new_ssa_obj.solvetime = sttime
#         new_ssa_obj.collisions = np.hstack((ssa_obj.collisions,collisions))
        
        
#         try:
#             new_ssa_obj.ribtimes = np.hstack((ssa_obj.ribtimes, all_ribtimes[np.where(all_ribtimes > 0)]))
            
            
            
#         except:
#             pass


#         #solt = solutions.T

#         fragmented_trajectories = []
#         fragtimes = []
#         maxlen = 0
    
#         fragmentspertraj= []
#         for k in range(n_traj):
#             ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
#             changes = ind[1:] - ind[:-1]
#             addindexes = np.where(changes > 0)[0]
#             subindexes = np.where(changes < 0)[0]
            
#             sub = solutions[k][:,1:] - solutions[k][:,:-1]
#             neutralindexes = np.unique(np.where(sub < 0)[1])
#             neutralindexes = np.setxor1d(neutralindexes, subindexes)
            
#             for index in neutralindexes:
#                 pre = solutions[k][:,index]
#                 post = solutions[k][:,index+1]
#                 changecount = 0
#                 while len(np.where(post - pre < 0)[0]) > 0:
    
#                     post = np.append([genelength],post)
#                     pre = np.append(pre,0)
                    
#                     changecount+=1
                
#                 for i in range(changecount):
#                     addindexes = np.sort(np.append(addindexes,index))
#                     subindexes = np.sort(np.append(subindexes,index))
                    
#                 changes[index] = -changecount
#                 ind[index] += changecount
             
                
#             for index in np.where(np.abs(changes)>1)[0]:
#                 if changes[index] < 0:
#                     for i in range(np.abs(changes[index])-1):
#                         subindexes = np.sort(np.append(subindexes,index))
#                 else:
#                     for i in range(np.abs(changes[index])-1):
#                         addindexes = np.sort(np.append(addindexes,index))   
                
#             truefrags = len(subindexes)
     
                
        
           
#             if len(subindexes) < len(addindexes):
#                 subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(len(truetime)-1)).astype(int))
                
            
#             fragmentspertraj.append(len(subindexes))
            
#             for m in range(min(len(subindexes),len(addindexes))):
#                 traj = solutions[k][:, addindexes[m]:subindexes[m]+1]
#                 traj_ind = changes[addindexes[m]:subindexes[m]+1]
                
#                 startind = ind[addindexes[m]]
#                 minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
#                 fragment = np.array([])
            
                    
                
#                 iterind = startind
                
#                 if subindexes[m]-addindexes[m] > 0:
#                     if len(minusloc) > 1:
#                         if m <= truefrags:
#                             for n in range(len(minusloc)-1):
#                                 iterind = iterind + min(0,traj_ind[minusloc[n]])
#                                 fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                                
                                
                                
                  
                
                      
                            
#                             fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                            
#                         else:
#                             for n in range(len(minusloc)-1):

#                                 iterind = iterind + min(0,traj_ind[minusloc[n]])
                                
#                                 fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                  
                                
#                             fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
          
                        
                    
#                     else:

#                         fragment = solutions[k][startind][addindexes[m]:subindexes[m]+1].flatten()
                   
                
                    
#                     fragtimes.append(addindexes[m]+1)
                       
                    
#                     fragmented_trajectories.append(fragment)
#                     #if m <= truefrags:
#                         #kes.append(genelength/truetime[len(fragment)])
            
#                     if len(fragment) > maxlen:
#                         maxlen = len(fragment)
                    
    
#             fragarray = np.zeros((len(fragmented_trajectories), maxlen))
#             for i in range(len(fragmented_trajectories)):
#                 fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
                
                
#         fraglen_size = max(fragarray.shape[1],ssa_obj.fragments.shape[1])
        
#         if fragarray.shape[1] != fraglen_size:
#             fragarray =  np.hstack((fragarray, np.zeros((fragarray.shape[0],fraglen_size-fragarray.shape[1]))) )
#         if ssa_obj.fragments.shape[1] != fraglen_size:
#             ssa_obj.fragments =  np.hstack((ssa_obj.fragments, np.zeros((ssa_obj.fragments.shape[0],fraglen_size-ssa_obj.fragments.shape[1]))) )            

        
#         new_ssa_obj.fragments = np.vstack((ssa_obj.fragments,fragarray))
#         new_ssa_obj.fragtimes = ssa_obj.fragtimes+fragtimes
#         new_ssa_obj.frag_per_traj = fragmentspertraj
#         new_ssa_obj.full_frags = ssa_obj.full_frags + truefrags
#         new_ssa_obj.all_results = np.vstack((ssa_obj.all_results,all_results))
        
#         if pv.shape[0] > 1:
#             for i in range(pv.shape[0]):
#                 if i > 0:
#                     autocorr_vec2, mean_autocorr2, error_autocorr2, dwelltime2, ke_sim2  = self.get_autocorr(intensity_vec[i], truetime, 0, genelength)
#                     autocorr_vec = np.vstack((autocorr_vec,autocorr_vec2))
#                     mean_autocorr = np.vstack((mean_autocorr,mean_autocorr2))
#                     error_autocorr = np.vstack((error_autocorr,error_autocorr2))
#                     dwelltime.append(dwelltime2)
#                     ke_sim.append(ke_sim2)
#                 else:
#                     autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.get_autocorr(intensity_vec[i], truetime, 0, genelength)
#                     autocorr_vec_norm, mean_autocorr_norm, error_autocorr_norm, dwelltime, ke_sim = self.get_autocorr_norm(intensity_vec[i], truetime, 0, genelength)
#                     dwelltime = [dwelltime]
#                     ke_sim = [ke_sim]
            
#         else:
#             autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.get_autocorr(intensity_vec, truetime, 0, genelength)
#             autocorr_vec_norm, mean_autocorr_norm, error_autocorr_norm, dwelltime, ke_sim = self.get_autocorr_norm(intensity_vec, truetime, 0, genelength)
            
#             acov,nacov = self.get_all_autocovariances(intensity_vec,truetime,genelength )              
        
#         new_ssa_obj.autocorr_vec = autocorr_vec
#         new_ssa_obj.mean_autocorr = mean_autocorr
#         new_ssa_obj.error_autocorr = error_autocorr
#         new_ssa_obj.autocorr_vec_norm = autocorr_vec_norm
#         new_ssa_obj.mean_autocorr_norm = mean_autocorr_norm
#         new_ssa_obj.error_autocorr_norm = error_autocorr_norm
#         new_ssa_obj.dwelltime = dwelltime
        
#         new_ssa_obj.ke_sim = float(n_traj)/(n_traj+ssa_obj.n_traj) * ke_sim  + float(ssa_obj.n_traj)/(n_traj+ssa_obj.n_traj) * ssa_obj.ke_sim
#         new_ssa_obj.ke_true = float(genelength)/np.mean(   new_ssa_obj.ribtimes   )
#         new_ssa_obj.probe = ssa_obj.probe
        
        
        
#         new_ssa_obj.autocovariance_dict  = acov
#         new_ssa_obj.autocovariance_norm_dict = nacov

# #        try:
# #            probePosition = []
# #            for key in self.POI.tag_epitopes.keys():
# #                probePosition = probePosition + self.POI.tag_epitopes[key]
# #            probePosition = np.unique(probePosition).tolist()
# #        except:
# #            print('No POI found')
# #                #nt_seq = self.tag_full['T_flag'] + nt_seq
# #
# #
# #        nt_seq = self.POI.nt_seq
# #        genelength = int(len(nt_seq)/3)
# #
# #
# #
# #        pv = np.zeros((1, genelength)).astype(int).flatten()
# #
# #        for i in range(len(probePosition)):
# #            pv[probePosition[i]:] = i
# #
# #
# #
# #
# #
# #        npoints = len(time_vec_fixed)
# #        tstep = npoints-non_consider_time
# #        for i in range(n):
# #
# #            soln = self.SSA(all_k, time_vec_fixed, inhibit_time=time_inhibit+non_consider_time, FRAP=evaluating_frap, Inhibitor=evaluating_inhibitor)
# #
# #            rb = sparse.lil_matrix((len(time_vec_fixed), genelength), dtype=int)
# #            for j in range(soln.shape[1]):
# #
# #                #if len(np.where(soln[:,j]!=0)[0]) !=0:
# #                #print(np.where(soln[:,j]!=0)[0])
# #
# #
# #                #rb[j,np.where(soln[:,j]!=0)[0]] = 1
# #
# #
# #                    for value in soln[:, j][np.where(soln[:, j] != 0 )[0]].astype(int):
# #
# #                        rb[j, value-1] = 1
# #
# #            rib_vec.append(rb)
# #
# #
# #        no_ribosomes = np.zeros((len(rib_vec), genelength))
# #
# #
# #
# #        for i in range(len(rib_vec)):
# #            no_ribosomes[i] = np.sum(rib_vec[i].todense()[non_consider_time:], axis=0).flatten()
# #
# #        ribosome_means = np.mean(no_ribosomes, axis=0)
# #        ribosome_density = ribosome_means/npoints
# #
# #        no_ribosomes_per_mrna = np.mean(no_ribosomes)
# #
# #        intensity_vec = np.zeros((len(rib_vec), tstep+1))
# #
# #        I = np.zeros((1, tstep+1))
# #        for i in range(len(rib_vec)):
# #            for j in range(tstep):
# #                temp_output = rib_vec[i][non_consider_time + j, :].todense()
# #
# #                I[0, j] = np.sum(pv * temp_output.flatten().T)
# #            intensity_vec[i] = I
# #
# #
# #
# #        ssa_obj = ssa()
# #
# #        ssa_obj.n_traj = nRepetitions + n
# #        ssa_obj.k = all_k
# #        ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
# #        ssa_obj.rib_density = ribosome_density
# #        ssa_obj.rib_means = ribosome_means
# #        ssa_obj.rib_vec = rib_vec
# #        ssa_obj.intensity_vec = intensity_vec
# #        ssa_obj.time_vec_fixed = time_vec_fixed
# #        ssa_obj.start_time = non_consider_time
# #        ssa_obj.probe = probePosition
# #        ssa_obj.evaluating_inhibitor = evaluating_inhibitor
# #        ssa_obj.evaluating_frap = evaluating_frap
# #        ssa_obj.time_inhibit = time_inhibit
# #
# #
# #
# #        if evaluating_inhibitor == False:
# #            autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.get_autocorr(intensity_vec, time_vec_fixed, 0, genelength)
# #            ssa_obj.autocorr_vec = autocorr_vec
# #            ssa_obj.mean_autocorr = mean_autocorr
# #            ssa_obj.error_autocorr = error_autocorr
# #            ssa_obj.dwelltime = dwelltime
# #            ssa_obj.ke_sim = ke_sim

#         return new_ssa_obj

