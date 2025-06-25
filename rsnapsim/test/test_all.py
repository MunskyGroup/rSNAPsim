# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:49:23 2025

@author: willi
"""

# This should be replaced with a dedicated list for security reasons
import os
import sys
test_files = ['test_poi.py',
              'test_seqmanip.py',
              'test_file_parser.py',
              'test_default_models.py']

for f in test_files:
    if f != 'test_all.py':
        print('running %s'%f)
        os.system('python %s'%f)