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
    
    
    
    runner.run(file_parser_suite)
    
    
    