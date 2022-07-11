# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:35:43 2020

@author: willi
"""
import re
import time
import os
import numpy as np
import itertools as it
from . import CodonDictionaries
from . import AuxCodonDicts
from . import FileParser
from . import poi as POI
from . import custom_errors as custom_err
from .core import SequenceCore
import warnings

try:
    from Bio import SeqIO
    from Bio import Entrez
except:

    print('BioPython is not installed, polling genbank will not be possible')
    pass


class SequenceManipMethods(SequenceCore):
    '''
    class that handles manipulation methods dealing with sequences
    '''
    def __init__(self, sequence=''):
        self.sequence = sequence
        #get the codon dictionaries
        self.codon_dicts = CodonDictionaries.CodonDictionaries()
        self.aux_dicts = AuxCodonDicts.AuxCodonDicts()
        self.previous_files = {}
        super().__init__()
        


    def open_seq_file(self, seqfile, min_codons=80, add_tag=True):

        '''
        Reads a sequence file, either a .txt file or a .gb genbank file

        *args*

            **seqfile**, sequence file either in txt, gb, gbk format
        '''

        fp = FileParser.FileParser()
        #TODO expose this to the user:
        #sequence_name = fp.get_name(seqfile)
        #sequence_description = fp.get_description(seqfile)
        
        sequence_str = fp.get_sequence(seqfile)
        if isinstance(sequence_str, str):
            sequence_str = sequence_str.upper()
            orfs = self.get_orfs(sequence_str, min_codons=min_codons)
    
            protein_strs, proteins, tagged_proteins = self.get_proteins(orfs,
                                                                        sequence_str,
                                                                        add_tag=add_tag)
        if isinstance(sequence_str, list):
            
            #multiple sequence fasta detected, split them up and handle them 
            #seperately, then merge them into the dictionaries later
            ps_list = []
            pro_list = []
            tagd_list = []
            for sequence in sequence_str:
                sequence = sequence.upper()
                orfs = self.get_orfs(sequence, min_codons=min_codons)
                
                protein_strs, proteins, tagged_proteins = self.get_proteins(orfs,
                                                                            sequence,
                                                                            add_tag=add_tag)
                ps_list.append(protein_strs)
                pro_list.append(proteins)
                tagd_list.append(tagged_proteins)
                
            # compress the dictionaries by keys
            flatten_dict = lambda dictionary: [[item for sublist in [y[x] for y in dictionary]
                                                for item in sublist] for x in ['1','2','3']]
            

            protein_strs = {i: flatten_dict(ps_list)[int(i)-1] for i in ['1','2','3']}
            proteins = {i: flatten_dict(pro_list)[int(i)-1] for i in ['1','2','3']}
            tagged_proteins = {i: flatten_dict(tagd_list)[int(i)-1] for i in ['1','2','3']}
                
        return protein_strs, proteins, tagged_proteins, sequence_str


    def get_proteins(self, orfs, seq, add_tag=True):
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
        proteins_strs = {'0':[], '+1':[], '+2':[], '-1':[]}
        protein_objs = {'0':[], '+1':[], '+2':[], '-1':[]}
        proteins_w_tags = {'0':[], '+1':[], '+2':[], '-1':[]}


        #tagged_proteins = {a:[] for a in cd.tag_dict.keys()}
        #tagged_protein_seq = {a:[] for a in cd.tag_dict.keys()}

        orf_keys = ['0','+1','+2','-1']
        for i in range(len(orfs)):
            
            for j in range(len(orfs[orf_keys[i]])):

                protein = POI.poi()

                pro = self.nt2aa(seq[orfs[orf_keys[i]][j][0]:orfs[orf_keys[i]][j][1]])
                nt_seq = seq[orfs[orf_keys[i]][j][0]:orfs[orf_keys[i]][j][1]]
                # if pro[-1] == '*':
                #     pro = pro[:-1]
                #     nt_seq = nt_seq[:-3]





                protein.aa_seq = pro
                protein.nt_seq = nt_seq

                proteins_strs[orf_keys[i]].append(pro)

                protein.gene_length = len(pro) #length of the gene
                protein.tag_length = 0   #length of the tags
                protein.total_length = len(pro)  #total length of the full amino acid sequence
                protein.source_seq = seq
                protein.UTR_5p = seq[:orfs[orf_keys[i]][j][0]]
                protein.UTR_3p =  seq[orfs[orf_keys[i]][j][1]:]
                protein.orf = orf_keys[i]
                protein.loc = (orfs[orf_keys[i]][j][0], orfs[orf_keys[i]][j][1])
                protein.tags = []

                protein_objs[orf_keys[i]].append(protein)




        for i in range(len(orfs)):
            for pr in protein_objs[orf_keys[i]]:

                tag_detected = False

                for tag in cd.tag_dict.keys():

                    if cd.tag_dict[tag] in pr.aa_seq:
                        tag_detected = True


                if tag_detected:
                    self.analyze_protein(pr)
                    pr.tag_added = False
                    proteins_w_tags[orf_keys[i]].append(pr)
                else:
                    if add_tag:
                        self.add_tag_to_protein(pr)
                        pr.tag_added = True
                    else:
                        pr.tag_added = False
                        self.analyze_protein(pr)



        return proteins_strs, protein_objs, proteins_w_tags


    def get_largest_poi(self,seqfile, min_codons=80, add_tag=True):
        '''
        Convenience file to get the largest poi if you know your file 
        has multiple orfs

        Parameters
        ----------
        seqfile : sequence file
            a file containting sequence data to get a mRNA sequence out of.
        min_codons : iny, optional
            minimum contingous codons to consider an ORF. The default is 80.
        add_tag : bool, optional
            Add a fluorescent tag if none is found. The default is True.

        Returns
        -------
        POI obj.

        '''
        fp = FileParser.FileParser()
        sequence_str = fp.get_sequence(seqfile).upper()
        orfs = self.get_orfs(sequence_str, min_codons=min_codons)
        protein_strs, proteins, tagged_proteins = self.get_proteins(orfs,
                                                                    sequence_str,
                                                                    add_tag=add_tag)
        
        sizes = [[len(y) for y in x] for x in protein_strs.values()]
        maxsize = max([item for sublist in sizes for item in sublist])
        
        for i in range(3):
            if maxsize in sizes[i]:
                frame = i
                pindex = sizes[i].index(maxsize)
                
        return proteins[str(int(frame+1))][pindex]
    
    def seq_to_protein_obj(self, nucleotide_sequence_str, min_codons=80, add_tag=True):
        '''

        Parameters
        ----------
        nucleotide_sequence_str : str
            nucleotide sequence string (.
        min_codons : TYPE, optional
            DESCRIPTION. The default is 80.
        add_tag : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        proteins : TYPE
            DESCRIPTION.

        '''
        nucleotide_sequence_str = self.clean_seq(nucleotide_sequence_str).upper()
        orfs = self.get_orfs(nucleotide_sequence_str, min_codons=min_codons)
        _, proteins, _ = self.get_proteins(
            orfs, nucleotide_sequence_str, add_tag=add_tag)

        return proteins
    

    def get_gb_file(self, accession_number, save_dir = '.'):
        '''
        A function to poll genbank given an accession number and pull the
        relevant gb file

        *args*

            **accession_number**, the accession number of the sequence to find.
            http://www.nslc.wustl.edu/elgin/genomics/bio4342/1archives/2006/AccReference.pdf

        *keyword args*

            **savetofile**, true or false to save the gb file in the same
            directory as sms for future use



        '''
        
        if not os.path.isdir(save_dir):
            msg = 'Specified save path does not exist, double check the path'\
            ' specified.'
            raise custom_err.PathDoesNotExistError(msg)
        
        # check if this file was already pulled
        if accession_number in self.previous_files.keys():
            if os.path.isfile(os.path.join(save_dir, self.previous_files[accession_number])):
                return os.path.join(save_dir, self.previous_files[accession_number])
             

        Entrez.email = "wsraymon@rams.colostate.edu"
        Entrez.tool = 'SingleMoleculeSimulator'
        er = False
        try:
            handle = Entrez.efetch(db="nucleotide", rettype="gb",
                                   retmode="text", id=accession_number)
            #using "gb" as an alias for "genbank"
            gb_record = SeqIO.read(handle, "genbank")
            handle.close()
        except:
            er = True
        time.sleep(2)
        if er == True:
            msg = 'Cannot find given ascession number for genbank, file re'\
                'quest failed.'
            raise custom_err.AscNumDoesNotExistError(msg)
            


        gb_rec = gb_record
        #gb_obj = gb_record

        #sequence_str = str(gb_record.seq)
        sequence_name = gb_record.name
        
        filename = os.path.join(save_dir, sequence_name+ '.gb')
        f = open(filename, 'w')


        f.write(gb_rec.format('gb'))

        f.close()
        
        self.previous_files[accession_number] = sequence_name + '.gb'
        return os.path.join(save_dir, sequence_name + '.gb')


    def add_tag_to_protein(self, POI, tag_type='T_Flag'):
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
        POI.tag_length = len(cd.tag_full[tag_type])

        self.analyze_protein(POI)


    def analyze_protein(self, poi_obj, epitope_loc='front',):
        cd = self.codon_dicts
        nt_seq = poi_obj.nt_seq
        aa_seq = poi_obj.aa_seq
        #self.POI.name = self.sequence_name
        total_length = len(poi_obj.aa_seq)

        '''
        for key in self.tagged_proteins:
            if protein in self.tagged_proteins[key]:
                self.POI.tag_types.append(key)
        '''
        poi_obj.tag_types = []
        for tag in cd.tag_dict.keys():
            if cd.tag_dict[tag] in aa_seq:
                poi_obj.tag_types.append(tag)

                #''.join(sms.poi[0].split('DYKDDDDK')

        poi_obj.tag_epitopes = {a:[] for a in poi_obj.tag_types}
        gs = poi_obj.aa_seq


        for i in range(len(poi_obj.tag_types)):

            try:
                nt_tag = cd.tag_full[poi_obj.tag_types[i]]
                aa_tag = self.nt2aa(nt_tag)
            except:
                epi = cd.tag_dict[poi_obj.tag_types[i]]
                firstep = poi_obj.aa_seq.find(epi)
                lastep = len(poi_obj.aa_seq) - poi_obj.aa_seq[::-1].find(epi[::-1])
                aa_tag = poi_obj.aa_seq[firstep:lastep]
                nt_tag = poi_obj.nt_seq[3*firstep:3*lastep]

            if epitope_loc == 'front':
                offset = 0
            if epitope_loc == 'middle':
                offset = int(len(cd.tag_dict[poi_obj.tag_types[i]])/2)
            if epitope_loc == 'back':
                offset = len(cd.tag_dict[poi_obj.tag_types[i]])

            poi_obj.tag_epitopes[poi_obj.tag_types[i]] = [
                m.start()+1+offset for m in re.finditer(
                    cd.tag_dict[poi_obj.tag_types[i]], poi_obj.aa_seq)]

            gs = gs.replace(aa_tag, '')
            
        poi_obj.gene_seq = gs
        poi_obj.gene_length = len(gs)
        poi_obj.total_length = total_length
        
        taglocs = np.array([x for x in poi_obj.tag_epitopes.values()])
        if taglocs.shape[0] > 0:
            tag_start,tag_stop = (int(np.min(taglocs)), int(np.max(taglocs)))
        
            poi_obj.tag_seq = aa_seq[tag_start:tag_stop]
            poi_obj.tag_length = len(aa_seq[tag_start:tag_stop])
        else:
            poi_obj.tag_seq = ''
            poi_obj.tag_length = 0

        codons = []
        for i in range(0, len(nt_seq), 3):
            codons.append(nt_seq[i:i+3])
        #POI.codons = codons

        #POI.codon_sensitivity, POI.CAI, POI.CAI_codons = self.codon_usage(POI.nt_seq)
        poi_obj.ki = .03
        poi_obj.ke = 10
        poi_obj.kt = 10
