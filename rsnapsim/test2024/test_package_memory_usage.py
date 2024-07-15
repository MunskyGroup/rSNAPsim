# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 12:09:54 2021

@author: willi
"""

from guppy import hpy; 
import os

heap = hpy()
heap_status1 = heap.heap()
heap.setref()
heap_status2 = heap.heap() #reference memory usage

cwd = os.getcwd()
os.chdir('../../..')

import rsnapsim as rss

os.chdir(cwd)
heap_status3 = heap.heap()
print("Memory Usage of module: rsnapsim -  ", (heap_status3.size - heap_status2.size)/1e6, " mb")
