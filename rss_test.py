# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 14:24:28 2020

@author: willi
"""

import unittest
import numpy as np
import FileParser

#Default Test example using strings
class TestStringMethods(unittest.TestCase):

    def test_upper(self):    #test for a Equal Condition
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self): #Test for a true and a False
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):  #test for a specific error
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)
            
            
class TestFileParser(unittest.TestCase):
    def setUp(self):
        self.parser = FileParser.FileParser()


    # Txt parsing
    
    def test_get_seq_txt(self):
        seq = self.parser.get_sequence('example.txt')
        self.assertEqual(seq,'AAUGAUCUAGUCGUGUGACUUACUGGGGAUCGGUCAGUGUCGUUGGGCAUGU'.lower())
        
    def test_get_name_txt(self):
        seq = self.parser.get_name('example.txt')
        self.assertEqual(seq,'sequencename')

    def test_get_desc_txt(self):       
        seq = self.parser.get_description('example.txt')
        self.assertEqual(seq,'<unknown description>')
     
    #genbank file parsing        
        
    def test_get_seq_gb(self):
        seq = self.parser.get_sequence('example.gb')
        self.assertEqual(seq[:40],'uguaacgaacggugcaauagugauccacacccaacgccug'.lower())
        
    def test_get_name_gb(self):       
        seq = self.parser.get_name('example.gb')
        self.assertEqual(seq,'NC_005816.1')

    def test_get_desc_gb(self):       
        seq = self.parser.get_description('example.gb')
        self.assertEqual(seq,'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence')
       
    #DNA file parsing    
        
    def test_get_seq_dna(self):
        seq = self.parser.get_sequence('example.dna')
        self.assertEqual(seq[:40],'uggaagggcuaauucacucccaaagaagacaagauauccu'.lower())
        
    def test_get_name_dna(self):       
        seq = self.parser.get_name('example.dna')
        self.assertEqual(seq,'<unknown name>')

    def test_get_desc_dna(self):       
        seq = self.parser.get_description('example.dna')
        self.assertEqual(seq,'<unknown description>')

    #fasta file parsing   

    def test_get_seq_fasta(self):
        seq = self.parser.get_sequence('example_fasta.fasta')
        self.assertEqual(seq,'AAUAUGAUGUAGCUGUCGUAGCUAGUC'.lower())
        
    def test_get_name_fasta(self):       
        seq = self.parser.get_name('example_fasta.fasta')
        self.assertEqual(seq,'MCHU')

    def test_get_desc_fasta(self):       
        seq = self.parser.get_description('example_fasta.fasta')
        self.assertEqual(seq,'MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken')



class TestSequenceMethods(unittest.TestCase):
    def setUp(self):
        from rss import SequenceManipMethods
        self.smm = SequenceManipMethods('')
        self.test_str = 'AAUGAUCUAGUCGUGUGACUUACUGGGGAUCGGUCAGUGUCGUUGGGCAUG'

        self.test_seq = 'AUGGACUACAAGGACGACGACGACAAAGGUGACUACAAAGAUGAUG'\
        'ACGAUAAAGGCGACUAUAAGGACGAUGACGACAAGGGCGGAAACUCACUGAUCAAGGAAAACAUGCGGA'\
        'UGAAGGUGGUGAUGGAGGGCUCCGUGAAUGGUCACCAGUUCAAGUGCACCGGAGAGGGAGAGGGAAACCC'\
        'GUACAUGGGAACUCAGACCAUGCGCAUUAAGGUCAUCGAAGGAGGUCCGCUGCCGUUCGCUUUCGAUAUC'\
        'CUGGCCACUUCGUUCGGAGGAGGGUCGCGCACGUUCAUCAAGUACCCGAAGGGAAUCCCGGACUUCUUUA'\
        'AGCAGUCAUUCCCGGAAGGAUUCACUUGGGAACGGGUGACCCGGUAUGAAGAUGGAGGUGUGGUGACUGU'\
        'CAUGCAAGAUACUUCGCUGGAGGAUGGGUGCCUCGUGUACCACGUCCAAGUCCGCGGAGUGAAUUUCCCGU'\
        'CCAACGGACCAGUGAUGCAGAAAAAGACGAAGGGUUGGGAACCUAAUACUGAAAUGAUGUACCCCGCAGAC'\
        'GGAGGGCUGAGGGGCUACACCCACAUGGCGCUGAAGGUCGACGGAGGAGAUUACAAGGAUGACGACGAUAA'\
        'GCAACAAGAUUACAAAGACGAUGAUGACAAGGGCCAGCAGGGCGACUACAAGGACGACGACGACAAGCAGC'\
        'AGGACUACAAAGAUGACGAUGAUAAAGGAGGAGGACAUCUGUCCUGUUCGUUCGUGACCACCUACAGAUCA'\
        'AAGAAAACCGUGGGAAACAUCAAGAUGCCGGGCAUUCAUGCCGUCGACCACCGCCUGGAGCGGCUCGAAGA'\
        'AUCAGACAAUGAGAUGUUCGUCGUGCAAAGAGAACAUGCCGUGGCCAAGUUCGCGGGACUGGGAGGCGGUG'\
        'GAGGCGAUUACAAAGACGAUGAUGACAAGGGUGACUAUAAAGACGACGAUGACAAAGGGGAUUACAAGGAU'\
        'GAUGAUGAUAAGGGAGGCGGUGGAUCAGGUGGAGGAGGUUCACUGCAGGAUGAUGAUAUCGCCGCGCUCGU'\
        'CGUCGACAACGGCUCCGGCAUGUGCAAGGCCGGCUUCGCGGGCGACGAUGCCCCCCGGGCCGUCUUCCCCU'\
        'CCAUCGUGGGGCGCCCCAGGCACCAGGGCGUGAUGGUGGGCAUGGGUCAGAAGGAUUCCUAUGUGGGCGAC'\
        'GAGGCCCAGAGCAAGAGAGGCAUCCUCACCCUGAAGUACCCCAUCGAGCACGGCAUCGUCACCAACUGGGA'\
        'CGACAUGGAGAAAAUCUGGCACCACACCUUCUACAAUGAGCUGCGUGUGGCUCCCGAGGAGCACCCCGUGC'\
        'UGCUGACCGAGGCCCCCCUGAACCCCAAGGCCAACCGCGAGAAGAUGACCCAGAUCAUGUUUGAGACCUUC'\
        'AACACCCCAGCCAUGUACGUUGCUAUCCAGGCUGUGCUAUCCCUGUACGCCUCUGGCCGUACCACUGGCAU'\
        'CGUGAUGGACUCCGGUGACGGGGUCACCCACACUGUGCCCAUCUACGAGGGGUAUGCCCUCCCCCAUGCCA'\
        'UCCUGCGUCUGGACCUGGCUGGCCGGGACCUGACUGACUACCUCAUGAAGAUCCUCACCGAGCGCGGCUAC'\
        'AGCUUCACCACCACGGCCGAGCGGGAAAUCGUGCGUGACAUUAAGGAGAAGCUGUGCUACGUCGCCCUGGA'\
        'CUUCGAGCAAGAGAUGGCCACGGCUGCUUCCAGCUCCUCCCUGGAGAAGAGCUACGAGCUGCCUGACGGCC'\
        'AGGUCAUCACCAUUGGCAAUGAGCGGUUCCGCUGCCCUGAGGCACUCUUCCAGCCUUCCUUCCUGGGCAUG'\
        'GAGUCCUGUGGCAUCCACGAAACUACCUUCAACUCCAUCAUGAAGUGUGACGUGGACAUCCGCAAAGACCU'\
        'GUACGCCAACACAGUGCUGUCUGGCGGCACCACCAUGUACCCUGGCAUUGCCGACAGGAUGCAGAAGGAGA'\
        'UCACUGCCCUGGCACCCAGCACAAUGAAGAUCAAGAUCAUUGCUCCUCCUGAGCGCAAGUACUCCGUGUGG'\
        'AUCGGCGGCUCCAUCCUGGCCUCGCUGUCCACCUUCCAGCAGAUGUGGAUCAGCAAGCAGGAGUAUGACGA'\
        'GUCCGGCCCCUCCAUCGUCCACCGCAAAUGCUUCUAG'

    
    def test_nt2aa(self):
        seq = self.smm.nt2aa(self.test_str)
        self.assertEqual(seq,'NDLVV*LTGDRSVSLGM'.upper())

    def test_nt2aa_lower(self):
        seq = self.smm.nt2aa(self.test_str.lower())
        self.assertEqual(seq,'NDLVV*LTGDRSVSLGM'.upper())
        
    def test_get_orfs(self):
        orfs = self.smm.get_orfs(self.test_seq)
        self.assertEqual(orfs['1'][0],(0,2133) )
        
    def test_get_protein(self):
        orfs = self.smm.get_orfs(self.test_seq)
        proteins_strs, proteins, proteins_w_tags = self.smm.get_proteins(orfs,self.test_seq)
        self.assertEqual(1,len(proteins['1']))
        self.assertEqual(1,len(proteins['2']))
        self.assertEqual(0,len(proteins['3']))
        
    def test_protein_obj_tagging(self):
        orfs = self.smm.get_orfs(self.test_seq)
        proteins_strs, proteins, proteins_w_tags = self.smm.get_proteins(orfs,self.test_seq)
        
        poi = proteins['1'][0]
        poi2 = proteins['2'][0]

        self.assertEqual(poi.tag_types[0],'T_Flag')
        self.assertEqual([2, 11, 20, 196, 206, 218, 228, 300, 309, 318], poi.tag_epitopes['T_Flag'])
    
        
        self.assertEqual(poi.tag_added,False)
        self.assertEqual(poi2.tag_added,True)
        
    def test_protein_obj_lengths(self):
        orfs = self.smm.get_orfs(self.test_seq)
        proteins_strs, proteins, proteins_w_tags = self.smm.get_proteins(orfs,self.test_seq)
        
        poi = proteins['1'][0]
        poi2 = proteins['2'][0]
        
        self.assertEqual(poi.tag_length, len(self.smm.nt2aa(self.test_seq)[:337] ) )
        self.assertEqual(poi.gene_length, len(self.smm.nt2aa(self.test_seq)[337:] ) )
        self.assertEqual(poi.total_length, len(self.smm.nt2aa(self.test_seq) ) )        
        
        
    def test_protein_obj_seqs(self):
        orfs = self.smm.get_orfs(self.test_seq)
        proteins_strs, proteins, proteins_w_tags = self.smm.get_proteins(orfs,self.test_seq)
        
        poi = proteins['1'][0]
        poi2 = proteins['2'][0]        
                
        self.assertEqual(poi.aa_seq, self.smm.nt2aa(self.test_seq) )     
        self.assertEqual(poi.nt_seq, self.test_seq )  
        self.assertEqual(poi.gene_seq, self.smm.nt2aa(self.test_seq)[337:] )
        self.assertEqual(poi.tag_seq, self.smm.nt2aa(self.test_seq)[:337] )
        

if __name__ == '__main__':
    #unittest.main(exit=False)
    
    runner = unittest.TextTestRunner()
    file_parser_suite = unittest.TestSuite()
    
    
    
    file_parser_suite.addTest(TestFileParser('test_get_seq_txt' ))
    file_parser_suite.addTest(TestFileParser('test_get_name_txt' ))
    file_parser_suite.addTest(TestFileParser('test_get_seq_gb' ))
    file_parser_suite.addTest(TestFileParser('test_get_name_gb' ))
    file_parser_suite.addTest(TestFileParser('test_get_seq_dna' ))
    file_parser_suite.addTest(TestFileParser('test_get_name_dna' ))
    
    file_parser_suite.addTest(TestFileParser('test_get_desc_fasta' ))
    file_parser_suite.addTest(TestFileParser('test_get_name_fasta' ))
    file_parser_suite.addTest(TestFileParser('test_get_seq_fasta' ))
    file_parser_suite.addTest(TestFileParser('test_get_desc_dna' ))
    file_parser_suite.addTest(TestFileParser('test_get_desc_gb' ))
    file_parser_suite.addTest(TestFileParser('test_get_desc_txt' ))    
    
    sequence_methods_suite = unittest.TestSuite()
    sequence_methods_suite.addTest(TestSequenceMethods('test_nt2aa' ))
    sequence_methods_suite.addTest(TestSequenceMethods('test_nt2aa_lower' ))
    sequence_methods_suite.addTest(TestSequenceMethods('test_get_orfs' ))
    sequence_methods_suite.addTest(TestSequenceMethods('test_get_protein' ))
    sequence_methods_suite.addTest(TestSequenceMethods('test_protein_obj_tagging' ))
    sequence_methods_suite.addTest(TestSequenceMethods('test_protein_obj_lengths' ))
    sequence_methods_suite.addTest(TestSequenceMethods('test_protein_obj_seqs' ))
    
    
    print('testing file parser...')
    runner.run(file_parser_suite)
    
    print('testing sequence manipulation methods...')
    runner.run(sequence_methods_suite)
    
    
    