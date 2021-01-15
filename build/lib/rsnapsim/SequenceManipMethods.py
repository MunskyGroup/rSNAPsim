# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:35:43 2020

@author: willi
"""
from . import CodonDictionaries 
from . import FileParser
from . import poi as POI
import numpy as np
import re

class SequenceManipMethods():
    '''
    class that handles anything dealing with sequences
    '''
    def __init__(self,sequence=''):
        self.sequence = sequence
        self.codon_dicts = CodonDictionaries.CodonDictionaries()
        pass    

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
            nextstop = orf1_stops[np.where(orf1_stops > start)[0][0]]+3
            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    orfs['1'].append((start, nextstop))

                    laststop = nextstop

        laststop = 0
        for start in orf2_starts:
            nextstop = orf2_stops[np.where(orf2_stops > start)[0][0]]+3
            if (nextstop - start) > min_len:
                if nextstop != laststop:
                    orfs['2'].append((start, nextstop))
                    laststop = nextstop

        laststop = 0
        for start in orf3_starts:
            nextstop = orf3_stops[np.where(orf3_stops > start)[0][0]]+3

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
                
                protein = POI.poi()
                
                pro = self.nt2aa(seq[orfs[str(i+1)][j][0]:orfs[str(i+1)][j][1]])
                nt_seq = seq[orfs[str(i+1)][j][0]:orfs[str(i+1)][j][1]]
                # if pro[-1] == '*':
                #     pro = pro[:-1]
                #     nt_seq = nt_seq[:-3]
                
                
                
                
                
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
        
    def seq_to_protein_obj(self, sequence_str, min_codons=80):
        
        orfs = self.get_orfs(sequence_str,min_codons=min_codons)
        protein_strs,proteins, tagged_proteins  = self.get_proteins(orfs,sequence_str)
        
        return proteins
    

    def open_seq_file(self, seqfile,min_codons=80):
        
        '''
        Reads a sequence file, either a .txt file or a .gb genbank file

        *args*

            **seqfile**, sequence file either in txt, gb, gbk format
        '''

        fp = FileParser.FileParser()
        sequence_name = fp.get_name(seqfile)
        sequence_description = fp.get_description(seqfile)
        sequence_str= fp.get_sequence(seqfile).upper()
        
        
        orfs = self.get_orfs(sequence_str,min_codons=min_codons)
        protein_strs,proteins, tagged_proteins  = self.get_proteins(orfs,sequence_str)
        return protein_strs,proteins, tagged_proteins
        
                
        
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

