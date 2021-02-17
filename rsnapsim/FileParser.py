# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 16:29:21 2020

@author: William Raymond
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

import os

class FileParser():
    '''
    Class to parse incoming files
    
    TODO: Extend this for psuedouridine /m6a / modification dictionary
    '''
    def __init__(self):
        self.file = None
        self.letter_dict = ['a','c','u','g']
        # IPUAC SUBSTITIONS - They need to  adjusted / extended for mrna modifications?
        self.sub_dict = {'m':'a','w':'a','r':'a','y':'t','k':'g','s':'g','w':'a','h':'a','n':'a','Ψ':'u'}    
        
        self.test_seq = 'aucuguacguacguaucgaucguguacuggcaaaacguaguagcugagcaucaucuaug'                
        pass    

    def clean_seq(self,seq):
        '''
        Return an mrna sequence of lowercase a,u,c,g from IPUAC substitutions
        
        .. warning:: this code will replace substitutive nucleotides with preferential order a, g , u , c. for example: N (any base) is set to A, W (T,U, or A) is set to A, S (C or G) is set to G
        
        

        Parameters
        ----------
        seq : str
            sequence string.

        Returns
        -------
        seq : str
            cleaned sequence str (only lowercase a,u,c,g).

        '''
        seq = seq.lower()
        
        for key in self.sub_dict.keys():            
            seq = seq.replace(key,self.sub_dict[key])
            
        seq = seq.replace('t','u')
        return seq

    
    def get_sequence(self,file):
        '''
        

        Parameters
        ----------
        file : path
            Path to the file to open.

        Returns
        -------
        str
            mRNA sequence string.

        '''
        self.__check_valid_file(file)
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
    
    def __check_valid_file(self,file):
        extension = file.split('.')[-1]

        if extension in ['fasta','gb','txt','dna']:
            return True
        else:
            raise Exception("Unrecognized File type, the sequence file must be a .txt, .dna, .gb, or .fasta")
    
    
    def get_name(self,file):
        '''
        attempt to find the transcript name from a file 

        Parameters
        ----------
        file : path
            path to file.

        Returns
        -------
        str
            transcript name or "unknown".

        '''
        self.__check_valid_file(file)
        
        name = 'unknown'
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
        '''
        attempt to find the transcript description from a file 

        Parameters
        ----------
        file : path
            path to file.

        Returns
        -------
        str
            transcript description or "unknown".

        '''
        self.__check_valid_file(file)
        
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
        name = ''
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
                            
        if name =='':
            
            name =  os.path.basename(file)[:-4]
                    
        return name
