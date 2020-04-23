# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:27:39 2020

@author: willi
"""

import rSNAPsim as rss
import matplotlib.pyplot as plt
import time
import numpy as np

r = rss.rSNAPsim()
k_bin,k_lens = r.get_binned_k(np.ones(110)*2,50)



r.open_seq_file('gene_files/HUMINSR.gb')
r.run_default()
pv,pl = r.get_probvec()
bpv,bpl = r.get_binned_probe_vec(pl,100)
st = time.time()
ssa_1 = r.ssa_binned(bins=100)
print('time to run 100 bin resolution simulation')
print(time.time()-st)
print(' ')
st = time.time()
ssa_2 = r.ssa_solver()

print('time to run full resolution 100 traj')
print(time.time()-st)

plt.plot(ssa_1.intensity_vec.T, 'g',alpha=.1)
plt.plot(ssa_2.intensity_vec.T, 'm',alpha=.1)
plt.xlabel('time')
plt.ylabel('intensity')
plt.legend(['100 bin','full'])

plt.figure()
plt.plot(ssa_1.autocorr_vec,'g',alpha=.1)
plt.plot(ssa_2.autocorr_vec,'m',alpha=.1)
plt.xlabel('tau')
plt.ylabel('ACOV')
plt.legend(['100 bin','full'])