# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 18:13:32 2020

@author: willi
"""


from . import SequenceManipMethods
from . import PropensityFactory, ProbeVectorFactory
PropensityFactory = PropensityFactory.PropensityFactory
ProbeVectorFactory = ProbeVectorFactory.ProbeVectorFactory

from . import CodonDictionaries
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
import numpy as np

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

                
                

    
    def visualize_probe(self, colors=None):
        probe = self.probe_loc
        fig,ax = plt.subplots(1)
        N = len(self.kelong)
        ncolors = probe.shape[0]
        
        cmap = cm.get_cmap('gist_rainbow')
        if colors == None:
            colors = cmap(np.linspace(.01,.95, ncolors))
        
        rectangle =  mpatches.Rectangle((0,.1), N ,.8,linewidth=1,edgecolor='k',facecolor='gray')

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