# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 11:19:12 2021

@author: willi
"""

class AuxCodonDicts():
    '''
    This class contains auxillery codon dependency dictionaries,
    for use as rates / k elongation
    calculations.

    IBEN2015 tRNA Gene Copy Numbers (GCN):

        tRNA gene copy number variation in humans, James R Iben , Richard J Maraia
        Gene 2014 Feb 25; 536(2): 376-384

        https://pubmed.ncbi.nlm.nih.gov/24342656/

    hg19 tRNA Gene Copy Numbers (GCN):

        GtRNAdb 2.0: an expanded database of transfer RNA genes
        identified in complete and draft genomes
        Patricia P. Chan and Todd M. Lowe
        Nucleic Acids Research (online) 2015 Dec 15 - 2016 Jan4;
        44(Database issue) D184-D189

    Human codon frequency bias nakamura:

        https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606

    Human codon frequency bias HIVE:

        https://hive.biochemistry.gwu.edu/dna.cgi?cmd=tissue_codon_usage&id=586358&mode=cocoputs

        A new and updated resource for codon usage tables.
        Athey J. et al. BMC Bioinformatics Sept 02 2017 - 18, Article number: 391 (2017)
        https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1793-7

    '''
    def __init__(self):

        #https://www.pnas.org/doi/10.1073/pnas.1918145117
        self.Gobet2020_EPA_rates = {
            }
        
        self.Gobet2020_PA_rates = {
            }
        

        self.IBEN2015_tRNA_GCN_averages_by_anticodon = {
            'AGC': 36.23, 'CGC':4.68,
            'UGC':10.07, 'ACG': 6.87,
            'CCG':3.71, 'CCU': 5.48,
            'UCG':6.18, 'UCU':6.42,
            'AUU':1.71, 'GUU':31.70,
            'GUC':20.66, 'GCA':35.06,
            'CUG':27.85, 'UUG':13.73,
            'CUC':18.39, 'UUC':15.04,
            'CCC':9.27, 'GCC':21.87,
            'UCC':11.15, 'GUG':10.13,
            'AAU':13.44, 'GAU':.66,
            'UAU':4.52, 'AAG':11.20,
            'CAA':7.31, 'CAG':14.41,
            'UAA':7.71, 'UAG':2.96,
            'CUU':16.75, 'UUU':18.36,
            'CAU':22.07, 'GAA':12.08,
            'AGG':9.58, 'CGG':3.59,
            'UGG':5.81, 'UCA':3.32,
            'AGA':11.76, 'CGA':3.37,
            'GCU':8.87, 'UGA':5.17,
            'CUA':1.11, 'UUA':2.36,
            'AGU':10.15, 'CGU':4.71,
            'UGU':6.39, 'CCA':9.66,
            'AUA':1.21, 'GUA':16.09,
            'AAC':11.51, 'CAC':18.12,
            'UAC':5.46}


        #tmpdict = self.__convert_codon_to_anticodon(human_avtRNA_gene_number_anticodon_IBEN2015)
        #self.IBEN2015_tRNA_GCN_averages = self.__add_other_keys(tmpdict)

        self.IBEN2015_tRNA_GCN_stds_by_anticodon = {
            'AGC': 2.68, 'CGC':0.65,
            'UGC':1.02, 'ACG': .66,
            'CCG':.4, 'CCU': .59,
            'UCG':.39, 'UCU': .60,
            'AUU':.25, 'GUU':4.47,
            'GUC':3.72, 'GCA':2.37,
            'CUG':1.79, 'UUG':.64,
            'CUC':3.8, 'UUC':1.27,
            'CCC':1.15, 'GCC':5.23,
            'UCC': 2.87, 'GUG':1.34,
            'AAU':2.09, 'GAU':.25,
            'UAU':.62, 'AAG':1.35,
            'CAA':.54, 'CAG':3.12,
            'UAA':.64, 'UAG':.36,
            'CUU':1.38, 'UUU':.67,
            'CAU':1.83, 'GAA':.99,
            'AGG':1.13, 'CGG':.64,
            'UGG':.85, 'UCA':.55,
            'AGA':.91, 'CGA':.5,
            'GCU':1.04, 'UGA':.4,
            'CUA':.24, 'UUA':0.08,
            'AGU':.7, 'CGU':.48,
            'UGU':.42, 'CCA':.99,
            'AUA':.11, 'GUA':1.68,
            'AAC':.91, 'CAC':1.78,
            'UAC':.68}

        #tmpdict = self.__convert_codon_to_anticodon(human_stdtRNA_gene_number_anticodon_IBEN2015)
        #self.IBEN2015_tRNA_GCN_stds = self.__add_other_keys(tmpdict)


        self.hg19_tRNA_GCN_by_anticodon = {
            'AGC': 33, 'CGC':5,
            'UGC':11, 'ACG': 8, 'GGC':2,
            'CCG':4, 'CCU': 8,
            'UCG':6, 'UCU': 6,
            'AUU':2, 'GUU':36,
            'GUC':19, 'GCA':36,
            'CUG':21, 'UUG':10,
            'CUC':13, 'UUC':14,
            'CCC':13, 'GCC':15,
            'UCC': 11, 'GUG':10,
            'AAU':18, 'GAU':6,
            'UAU':5, 'AAG':13,
            'CAA':8, 'CAG':10,
            'UAA':6, 'UAG':4,
            'CUU':26, 'UUU':20,
            'CAU':11, 'GAA':15,
            'AGG':11, 'CGG':4,
            'UGG':9, 'UCA':2,
            'AGA':11, 'CGA':4,
            'GCU':9, 'UGA':5,
            'CUA':0, 'UUA':2,
            'AGU':10, 'CGU':6,
            'UGU':6, 'CCA':8,
            'AUA':1, 'GUA':16,
            'AAC':12, 'CAC':18, #duplicate for doubled up codons
            'UAC':7,}

        #tmpdict = self.__convert_codon_to_anticodon(self.hg19_tRNA_GCN_anticodon)
        #self.hg19_tRNA_GCN_codon = self.__add_other_keys(tmpdict)




        self.Human_codon_freq_nakumura = {}

        # Frequncy per thousand
        self.Human_codon_freq_nakumura = {
            'TTT': 17.6, 'TCT': 15.2, 'TAT': 12.2, 'TGT': 10.6, 'TTC': 20.3,
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

        self.Human_codon_freq_nakumura = self.__add_other_keys(self.Human_codon_freq_nakumura)



        self.Human_codon_freq_HIVE = {}

        # Frequncy per thousand
        self.Human_codon_freq_HIVE = {
            'TTT': 17.22, 'TCT': 16.96, 'TAT': 12.15, 'TGT': 10.47, 'TTC': 17.51,
            'TCC': 17.29, 'TAC': 13.47, 'TGC': 10.83, 'TTA': 8.82, 'TCA': 14.9,
            'TAA': .44, 'TGA': .8, 'TTG': 13.47, 'TCG':  4.05, 'TAG': 0.35,
            'TGG': 11.66, 'CTT': 14.7, 'CCT': 19.15, 'CAT': 11.9, 'CGT': 4.56,
            'CTC': 17.81, 'CCC': 18.86, 'CAC': 14.62, 'CGC': 8.71, 'CTA':  7.48,
            'CCA': 18.79, 'CAA': 14.17, 'CGA':  6.42, 'CTG': 36.01, 'CCG':  6.15,
            'CAG': 35.28, 'CGG': 10.62, 'ATT': 16.58, 'ACT': 14.3, 'AAT': 18.52,
            'AGT': 14.06, 'ATC': 18.68, 'ACC': 17.77, 'AAC': 18.27, 'AGC': 19.66,
            'ATA':  8.16, 'ACA': 16.57, 'AAA': 27.77, 'AGA': 13.4, 'ATG': 21.48,
            'ACG': 5.59, 'AAG': 31.79, 'AGG': 12.17,
            'GTT': 11.8, 'GCT': 18.89,
            'GAT': 24.09, 'GGT': 10.82, 'GTC': 13.44, 'GCC': 25.64, 'GAC': 24.25,
            'GGC': 19.73, 'GTA':  7.71, 'GCA': 17.09, 'GAA': 33.85, 'GGA': 17.19,
            'GTG': 25.75, 'GCG': 5.94, 'GAG': 39.4, 'GGG': 15.26}

        self.Human_codon_freq_HIVE = self.__add_other_keys(self.Human_codon_freq_HIVE)

    def __add_other_keys(self, dictionary):
        keys = list(dictionary.keys())

        for key in keys:
            if 'U' in key.upper():
                newkey = key.upper().replace('U', 'T')
                if newkey not in keys:
                    dictionary[newkey] = dictionary[key]
            if 'T' in key.upper():
                newkey = key.upper().replace('T', 'U')
                if newkey not in keys:
                    dictionary[newkey] = dictionary[key]
        return dictionary


    def __convert_codon_to_anticodon(self, dictionary):
        '''
        Convert a dictionary's keys to its complimentary anticodons
        (or anticodons back to codons)

        for example:

            {'GUU': 15, 'GTT': 15, 'GAA' : 10}

        becomes:

            {'AAC': 15, 'AAC': 15, 'UUC' : 10}

        Parameters
        ----------
        dictionary : dict
            dictionary with codons as the keys.

        Returns
        -------
        newdict : dict
            converted dictionary.

        '''
        newdict = {}
        for key in dictionary.keys():

            newkey = key[::-1]
            newkeys = [x for x in newkey]

            for i in range(3):
                if newkeys[i] == 'A':
                    newkeys[i] = 'U'

                elif newkeys[i] == 'U':
                    newkeys[i] = 'A'

                elif newkeys[i] == 'T':
                    newkeys[i] = 'A'

                elif newkeys[i] == 'G':
                    newkeys[i] = 'C'

                elif newkeys[i] == 'C':
                    newkeys[i] = 'G'
            codon = ''.join(newkeys)
            newdict[codon] = dictionary[key]
        return newdict

    def single_dict(self, dictionary, u_or_t='t'):
        '''
        Function that will return a dictionary with keys only by U or T,
        not both. For example say you have a dictionary with the following:

            {'GUU': 15, 'GTT': 15, 'GAA' : 10}

        this function will remove the duplicates and only use U or T

            {'GUU': 15, 'GAA' : 10}

        Parameters
        ----------
        dictionary : dict
            Dictionary of codons to value with both U and T.
        u_or_t : str, optional
            character to keep in the dictionary, ie all
            keys will be converted to t or u. The default is 't'.

        Returns
        -------
        new_dict : dict
            Singular dictionaries with codons that only use T or U for
            codons.

        '''
        new_dict = {}
        if u_or_t.lower() == 't':
            switchchar = 'u'
        else:
            switchchar = 't'
        for key in dictionary.keys():
            if switchchar.upper() in key:
                newkey = key.replace(switchchar.upper(), u_or_t.upper())

                new_dict[newkey] = dictionary[key]


        return new_dict

'''
    def load_from_bed(self, bed_file, anticodon=True):
        #TODO
        return 1

    def load_freq_from_ccds(self, ccds_fa):
        #TODO
        return 1
'''