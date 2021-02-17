# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:48:33 2020

@author: William Raymond
"""

import rsnapsim as rSNAPsim
import time
import os
import sys
import platform
'''
Method to tag solution objects with some metadata
'''
class GenericMetaData():
    def __init__(self):
        self.id = ''
        self.created_at = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime( time.time() ))
        self.user = os.path.expanduser("~")
        self.rss_version = rSNAPsim.__version__
        self.platform = platform.platform()
        self.python_version = sys.version
        
    def get(self):
        '''
        generate and return a metadata dictionary for a solver object

        Returns
        -------
        dict
            a dictionary of metadata such as solution id, time ran, user, platform and rss version.

        '''
        return self.__dict__
        