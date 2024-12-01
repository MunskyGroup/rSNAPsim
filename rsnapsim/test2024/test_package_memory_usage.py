# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 12:09:54 2021

@author: willi
"""



#heap = hpy()
#heap_status1 = heap.heap()
#heap.setref()
#heap_status2 = heap.heap() #reference memory usage

#heap_status3 = heap.heap()
#print("Memory Usage of module: rsnapsim -  ", (heap_status3.size - heap_status2.size)/1e6, " mb")



from collections import Counter
import linecache
import os
import tracemalloc

def display_top(snapshot, key_type='lineno', limit=3):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        print("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))

tracemalloc.start()

cwd = os.getcwd()
os.chdir('../../..')

import rsnapsim as rss

os.chdir(cwd)


snapshot = tracemalloc.take_snapshot()
display_top(snapshot)