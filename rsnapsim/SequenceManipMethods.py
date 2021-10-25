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
from . import FileParser
from . import poi as POI

try:
    from Bio import SeqIO
    from Bio import Entrez
except:

    print('BioPython is not installed, polling genbank will not be possible')
    pass




class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class PathDoesNotExistError(Error):
    """Exception raised for when trying to save a GB file to a directory 
    that doesnt exist

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
class AscNumDoesNotExistError(Error):
    """Exception raised for when trying to pull a gb from an ascession number
    that doesnt exist

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
        
class UnrecognizedFlagError(Error):
    """Exception raised for when trying to pull a gb from an ascession number
    that doesnt exist

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
                

class SequenceManipMethods():
    '''
    class that handles manipulation methods dealing with sequences
    '''
    def __init__(self, sequence=''):
        self.sequence = sequence
        self.codon_dicts = CodonDictionaries.CodonDictionaries()
        #get the codon dictionaries


    def optimize_ntseq(self, nt_seq, opt_dict=None):
        '''
        Optimizes a nucleotide sequence

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence string
        opt_dict : dictionary, optional
            a user defined dictionary to optimize over. The default is None.

        Returns
        -------
        opt_seq : str
            Optimized NT sequenced based on a given dictionary of rates

        '''

        if opt_dict is None:
            opt_dict = self.codon_dicts.human_codon_frequency_bias_nakamura

        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
        aa = [self.codon_dicts.aa_table[x] for x in seperated_codons]
        opt_seq = ''
        for i in range(0, len(aa)):
            ind = np.argmax([opt_dict[x] for x in self.codon_dicts.aa_table_r[aa[i]]])
            opt_codon = self.codon_dicts.aa_table_r[aa[i]][ind]
            opt_seq = opt_seq + opt_codon
        return opt_seq


    def deoptimize_ntseq(self, nt_seq, deopt_dict=None):

        '''
        Optimizes a nucleotide sequence

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence string
        deopt_dict : dictionary, optional
            a user defined dictionary to deoptimize over. The default is None.

        Returns
        -------
        deopt_seq : str
            Deoptimized NT sequenced based on a given dictionary of rates

        '''

        if deopt_dict is None:
            deopt_dict = self.codon_dicts.human_codon_frequency_bias_nakamura
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
        aa = [self.codon_dicts.aa_table[x] for x in seperated_codons]
        opt_seq = ''
        for i in range(0, len(aa)):
            ind = np.argmin([deopt_dict[x] for x in self.codon_dicts.aa_table_r[aa[i]]])
            opt_codon = self.codon_dicts.aa_table_r[aa[i]][ind]
            opt_seq = opt_seq + opt_codon
        return opt_seq


    def nt2aa(self, nt_seq):
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
            raise PathDoesNotExistError(msg)

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
            raise AscNumDoesNotExistError(msg)
            


        gb_rec = gb_record
        #gb_obj = gb_record

        #sequence_str = str(gb_record.seq)
        sequence_name = gb_record.name
        
        filename = os.path.join(save_dir, sequence_name+ '.gb')
        f = open(filename, 'w')


        f.write(gb_rec.format('gb'))

        f.close()


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
        allstarts = np.array([m.start() for m in re.finditer(
            '(?=A[TU]G((?:.{3})+?)[TU](?:AG|AA|GA))', nt_seq)])


        #allsegments = re.findall('(?=A[TU]G((?:.{3})+?)[TU](?:AG|AA|GA))',self.sequence_str)
        allstops = np.array(
            [m.start() for m in re.finditer('(?=[TU](?:AG|AA|GA))', nt_seq)])
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
    
    def get_codon_list(self, nt_seq):
        '''
        

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.

        Returns
        -------
        List of codons.

        '''
        #split codons by 3
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)]
        return seperated_codons
    
    def get_codon_count_dict(self, nt_seq):
        '''
        
        
        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.

        Returns
        -------
        dict of codon counts.

        '''
        #split codons by 3
        codons = nt_seq.upper()
        seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)]
        cdict = {}
        for i in range(len(seperated_codons)):
            if seperated_codons[i] in cdict.keys():
                cdict[seperated_codons[i]] +=1
            else:
                cdict[seperated_codons[i]] = 1
            
        return cdict

    def codon_usage(self, nt_seq, codon_dict=None):
        '''
        Analyzes codon useage from the nucleotide sequence

        *args*

            **nt_seq**,  nucleotide sequence as a string

        *returns*

            **codon_sensitivity**, a list of codon sensitivity for the nucleotide sequence

            **cai**, cai value

        '''

        if codon_dict == None:
            codon_dict = self.codon_dicts.human_codon_frequency_bias_nakamura

        codon_usage = np.zeros((1, 21))
        gene_len = len(nt_seq)/3
        aa_seq = self.nt2aa(nt_seq)

        for i in range(len(self.codon_dicts.aa_keys)-1):

            codon_usage[0, i] = len(
                re.findall(self.codon_dicts.aa_keys[i], aa_seq))
        codon_usage[0, 20] = len(re.findall('\*', aa_seq))
        codon_norm = codon_usage/gene_len
        codon_sensitivity = np.round(
            codon_norm*self.codon_dicts.sensitivity_fast_slow, 2)

        cai_codons = []
        for i in range(0, len(nt_seq), 3):
            synonmous_codons = self.codon_dicts.aa_table_r[
                self.codon_dicts.aa_table[nt_seq[i:i+3]]]

            max_freq = max([codon_dict[x] for x in synonmous_codons])

            cai_codons.append(codon_dict[nt_seq[i:i+3]] /max_freq)

        cai = self.geomean(cai_codons)

        return codon_sensitivity, cai, cai_codons



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
        proteins_strs = {'1':[], '2':[], '3':[]}
        protein_objs = {'1':[], '2':[], '3':[]}
        proteins_w_tags = {'1':[], '2':[], '3':[]}


        #tagged_proteins = {a:[] for a in cd.tag_dict.keys()}
        #tagged_protein_seq = {a:[] for a in cd.tag_dict.keys()}

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

                protein.gene_length = len(pro) #length of the gene
                protein.tag_length = 0   #length of the tags
                protein.total_length = len(pro)  #total length of the full amino acid sequence
                protein.source_seq = seq
                protein.orf = i
                protein.loc = (orfs[str(i+1)][j][0], orfs[str(i+1)][j][1]+3)
                protein.tags = []

                protein_objs[str(i+1)].append(protein)




        for i in range(len(orfs)):
            for pr in protein_objs[str(i+1)]:

                tag_detected = False

                for tag in cd.tag_dict.keys():

                    if cd.tag_dict[tag] in pr.aa_seq:
                        tag_detected = True


                if tag_detected:
                    self.analyze_protein(pr)
                    pr.tag_added = False
                    proteins_w_tags[str(i+1)].append(pr)
                else:
                    if add_tag:
                        self.add_tag_to_protein(pr)
                        pr.tag_added = True
                    else:
                        pr.tag_added = False
                        self.analyze_protein(pr)



        return proteins_strs, protein_objs, proteins_w_tags

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

    def seq_to_protein_obj(self, sequence_str, min_codons=80, add_tag=True):

        orfs = self.get_orfs(sequence_str, min_codons=min_codons)
        _, proteins, _ = self.get_proteins(
            orfs, sequence_str, add_tag=add_tag)

        return proteins


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
                
        return proteins[str(int(i))][pindex]


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


    def get_kmer_freq(self, nt_seq, kmer_length, substrings=False):
        '''
        return the K-mer frequency of a nuclotide sequence
    
        seq = 'AACGTACGTAGCTCATG...'
        
        kmer_dict with length 3 = 
        
        {'AAA':1, 'AAC': 2, 'AAG':0 ...}
        

        Parameters
        ----------
        nt_seq : str
            nuclotide sequence to get kmers from.
        kmer_length : int
            DESCRIPTION.
        substrings : bool, optional
            generate kmers for all lower length kmers as well, ie
            for k = 3, 2, and 1 instead of just 3. The default is False.

        Returns
        -------
        kmer_freq_vec : 1darray
            kmer frequency per k-mer key.
        kmer_ind : list of str
            list of keys for the k-mer frequency vector.

        '''
          
        nt_seq = nt_seq.upper()
        unique_char = list(set(nt_seq))

        combos =[x for x in it.product(unique_char, repeat=kmer_length)]   
        if substrings:
            for n in range(kmer_length-1, 0, -1):
                combos += [x for x in it.product(unique_char, repeat=n)] 
        kmer_ind = [''.join(y) for y in combos]

        kmer_freq_vec = np.zeros(len(combos)).astype(int)
        for i in range(len(nt_seq)-kmer_length):
            kmer_freq_vec[kmer_ind.index(nt_seq[i:i+kmer_length])] += 1
        if substrings:
            for n in range(kmer_length-1, 0, -1):
                for i in range(len(nt_seq)-n):
                    kmer_freq_vec[kmer_ind.index(nt_seq[i:i+n])] += 1            

        return kmer_freq_vec, kmer_ind
    
    
    def clean_nt_seq(self, nt_seq, upper=True, sub_dict=None, t_or_u='u',
                     random_sub=False):
        '''
        Return an mrna sequence of lowercase a,u,c,g from IPUAC substitutions

        .. warning:: this code will replace substitutive nucleotides with
        preferential order a, g , u , c. for example: N (any base) is set
        to A, W (T,U, or A) is set to A, S (C or G) is set to G



        Parameters
        ----------
        seq : str
            sequence string.

        Returns
        -------
        seq : str
            cleaned sequence str (only lowercase a,u,c,g).

        '''
        if sub_dict == None:
            sub_dict = self.codon_dicts.ipuac_nt_t

        seq = nt_seq.lower()

        if random_sub:
            new_str = []
            for i in range(len(seq)):
                if seq[i] not in ['a','t','g','c','u']:
                    idx = np.random.randint( len(sub_dict[seq[i]])   )
                    new_str += [sub_dict[seq[i]][idx] , ]
                else:
                    new_str += [seq[i], ]
            seq = ''.join(new_str)
        else:
            for key in sub_dict.keys():
                seq = seq.replace(key, sub_dict[key][0])
        
        if t_or_u in ['u','U']:
            seq = seq.replace('t', 'u')
        elif t_or_u in ['t','T']:
            seq = seq.replace('u', 't')
        else:
            msg = "Cannot recognize the flag for t_or_u, please use 'u' or 't'"\
                " to indicate which letter to use for the sequence."
            raise UnrecognizedFlagError(msg)
            
        if upper:
            seq = seq.upper()
        else:
            seq = seq.lower()
        return seq

    
    def get_kmer_freq_dict(self, nt_seq, kmer_length, substrings=False ):
        '''
        Get K-mer frequency as a dictionary, for example
        
        seq = 'AACGTACGTAGCTCATG...'
        
        kmer_dict with length 3 = 
        
        {'AAA':1, 'AAC': 2, 'AAG':0 ...}
        
        

        Parameters
        ----------
        nt_seq : str
            DESCRIPTION.
        kmer_length : int
            size kmers to generate.
        substrings : bool, optional
            generate kmers for all lower length kmers as well, ie
            for k = 3, 2, and 1 instead of just 3. The default is False.
        
        Returns
        -------
        dict
            dictionary of kmer frequencies.

        '''
        kmer_freq, kmer_key = self.get_kmer_freq(nt_seq, kmer_length, substrings=substrings)
        return  dict(zip(kmer_key, kmer_freq))

    def get_gc_content(self, nt_seq):
        '''
        Get the GC content of a nucleotide sequence

        Parameters
        ----------
        nt_seq : str
            nucleotide sequence.

        Returns
        -------
        GC content float (percentage).

        '''
        nt_seq = nt_seq.upper()
        return float((nt_seq.count('G') + nt_seq.count('C')))/len(nt_seq)
        

    def get_tag_loc(self, aa_seq, tag, epitope_loc='front'):
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
