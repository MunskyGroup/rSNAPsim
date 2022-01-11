# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:48:33 2020

@author: William Raymond
"""

import time
import sys
import platform
import rsnapsim as rSNAPsim

class GenericMetaData():
    '''
    Provides a metatdata dictionary to the user about rsnapsim files or
    models being generated.

    Attributes
    ----------

    id
    created_at
    rss_version
    platform
    python_version

    '''
    def __init__(self):
        self.id = ''
        self.created_at = time.strftime('%Y-%m-%d %H:%M:%S',
                                        time.localtime(time.time()))
        #self.user = os.path.expanduser("~")
        self.rss_version = rSNAPsim.__version__
        self.platform = platform.platform()
        self.python_version = sys.version

    def get(self):
        '''
        generate and return a metadata dictionary for a solver object

        Returns
        -------
        dict
            a dictionary of metadata such as solution id, time ran, user,
            platform and rss version.

        '''
        return self.__dict__

    def get_as_str(self):
        metadata = self.get()
        # convert to a metadata string
        meta_str = '######################################################'
        meta_str += 'model ID : ' + metadata['id'] + '\n'
        meta_str += 'files created at : ' + metadata['created_at']+ '\n'
        meta_str += 'platform : ' + metadata['platform']+ '\n'
        meta_str += 'python version : ' + metadata['python_version']+ '\n'
        meta_str += 'rsnapsim version : ' + metadata['rss_version']+ '\n'
        meta_str += '#####################################################'
        return meta_str
