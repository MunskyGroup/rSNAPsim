# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 16:29:21 2020

@author: willi
"""

import re
try:
    from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord
except:
    pass
try:
    from Bio import SeqIO
    from Bio import Entrez
except:
    pass

class FileParser():
    '''
    Class to parse incoming files
    '''
    def __init__(self):
        self.file = None
        self.letter_dict = ['a','c','u','g']
        self.sub_dict = {'m':'a','w':'a','r':'g','y':'t','k':'g','s':'g','w':'a','h':'a','n':'a'}       
        self.test_seq = 'aucuguacguacguaucgaucguguacuggcaaaacguaguagcugagcaucaucuaug'                
        pass    

    def clean_seq(self,seq):
        '''
        clean the sequences to lowercase only a, u, g, c
        '''
        seq = seq.lower()
        
        for key in self.sub_dict.keys():            
            seq = seq.replace(key,self.sub_dict[key])
            
        seq = seq.replace('t','u')
        return seq

    
    def get_sequence(self,file):
        '''
        Return the sequence from several different file types
        
        Supported: .txt, .dna, .gb, .fasta
        '''
        
        extension = file.split('.')[-1] 
        if extension == 'dna':            
            try:
                seq_record = snapgene_file_to_seqrecord(file)
            except:
                print('To read .dna files please install snapegenereader: pip install snapgene_reader - https://github.com/IsaacLuo/SnapGeneFileReader' )
            
            sequence_str = str(seq_record.seq)  

        if extension == 'txt':
            sequence_str =  self.__get_seq_from_txt(file)       
                    

        if extension == 'gb':
            gb_record = SeqIO.read(open(file, "r"), "genbank")
            sequence_str = str(gb_record.seq)

        if extension == 'fasta':
            fastas = list(SeqIO.parse(file,'fasta'))
            if len(fastas) > 1:
                return 'Multiple line fastas not supported'
            else:
                sequence_str = (str(fastas[0].seq))
        cleaned_sequence_str = self.clean_seq(sequence_str)
        
        return cleaned_sequence_str
    
    
    def get_name(self,file):
        
        extension = file.split('.')[-1] 
        if extension == 'fasta':                  
            fastas = list(SeqIO.parse(file,'fasta'))
            if len(fastas) > 1:
                return 'Multiple line fastas not supported'
            else:
                name = (str(fastas[0].id))        
        
        if extension == 'gb': 
            gb_record = SeqIO.read(open(file, "r"), "genbank")
            name = str(gb_record.id)   
            
        if extension == 'txt':
            name = self.__get_name_from_text(file)
            
        if extension == 'dna':
            try:
                seq_record = snapgene_file_to_seqrecord(file)
            except:
                print('To read .dna files please install snapegenereader: pip install snapgene_reader - https://github.com/IsaacLuo/SnapGeneFileReader' )
            name = seq_record.name
            
        return name
    
    
    def get_description(self,file):
        extension = file.split('.')[-1] 
        
        if extension == 'fasta':                  
            fastas = list(SeqIO.parse(file,'fasta'))
            if len(fastas) > 1:
                return 'Multiple line fastas not supported'
            else:
                
                desc = str(fastas[0].description)

        if extension == 'gb':                  
            gb_record = SeqIO.read(open(file, "r"), "genbank")
            desc = str(gb_record.description)   
            
        if extension == 'dna':
            try:
                seq_record = snapgene_file_to_seqrecord(file)
            except:
                print('To read .dna files please install snapegenereader: pip install snapgene_reader - https://github.com/IsaacLuo/SnapGeneFileReader' )
            desc = seq_record.description
        if extension == 'txt':
            desc = '<unknown description>'

        return desc
    
    
            
    @classmethod
    def __get_seq_from_txt(cls,file):
        with open(file) as f:
            raw = f.readlines()
        
        raw = ''.join(raw)        
        onlychar = re.split(r'[^A-Za-z]', raw)
        validt = ['A', 'G', 'T', 'C']
        validu = ['A', 'G', 'U', 'C']
        
        sequence_str = ''
   
        for i in range(len(onlychar)):
            section = onlychar[i]
            if set(section.upper()) == set(validt):
                sequence_str += section.upper()
                
            
            elif set(section.upper()) == set(validu):
                sequence_str += section.upper()              
        return sequence_str
    
    @classmethod
    def __get_name_from_text(cls,file):
        with open(file) as f:
            raw = f.readlines()
            raw = ''.join(raw)
            
            validt = ['A', 'G', 'T', 'C']
            
            validu = ['A', 'G', 'U', 'C']
            onlychar = re.split(r'[^A-Za-z]', raw)

            namelen = 0

            for i in range(len(onlychar)):
                section = onlychar[i]
                if set(section.upper()) != set(validt):
                    if set(section.upper()) != set(validu):
                        if len(section)>namelen:
                            name = section
                            namelen = len(section)
                    
        return name
