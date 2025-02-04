# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:49:23 2025

@author: willi
"""

# This should be replaced with a dedicated list for security reasons
import os
test_files = [f for f in os.listdir('.') if f[:5] == 'test_' and f[-3:] == '.py']

for f in test_files:
    if f != 'test_all.py':
        print('running %s'%f)
        os.system(f)