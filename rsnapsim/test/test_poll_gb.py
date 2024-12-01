# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 12:54:52 2021

@author: willi
"""
import os
cwd = os.getcwd()
os.chdir('../../..')

import rsnapsim as rss
from rsnapsim import seqmanip

import numpy as np
import time
import matplotlib.pyplot as plt

os.chdir(cwd)

Ascession_Number = 'MN908947.3'


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
        
#test catching wrong ascession



class AscNumDoesNotExistError(Error):
    """Exception raised for when trying to pull a gb from an ascession number
    that doesnt exist

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        
import unittest

class TestGBpulling(unittest.TestCase):

    def test_pull_file(self):
        if os.path.isfile(os.path.join('.', 'content', 'MN908947.gb')):
            os.remove(os.path.join('.', 'content', 'MN908947.gb'))
                               
                               
        rss.seqmanip.get_gb_file(Ascession_Number, '.\content')
        
        self.assertTrue(os.path.isfile(os.path.join('.',
                                                         'content',
                                                         'MN908947.gb')))

    def test_wrong_path(self):
        with self.assertRaises(rss.SequenceManipMethods.PathDoesNotExistError) as context:
            rss.seqmanip.get_gb_file('aaa', '.\content')
            msg = 'Specified save path does not exist, double check the path'\
            ' specified.'
            self.assertEqual(
                context.exception.msg,
                msg)
            
    def test_wrong_asc(self):
        with self.assertRaises(rss.SequenceManipMethods.AscNumDoesNotExistError) as context:
            rss.seqmanip.get_gb_file('113', '.\contents')
            msg = 'Cannot find given ascession number for genbank, file re'\
                'quest failed.'
            
            print('')
            print(context.exception)
            print('')
            self.assertEqual(
                context.exception.msg,
                msg)

if __name__ == '__main__':
    unittest.main()