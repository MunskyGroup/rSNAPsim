# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 12:43:13 2021

@author: willi
"""

import os
os.chdir('..')
os.chdir('..')
print(os.getcwd())
import rsnapsim as rss
import numpy as np
os.chdir('rsnapsim')
os.chdir('test')


import matplotlib.pyplot as plt
import time


poi_strs, poi_objs, tagged_pois,seq = rss.seqmanip.open_seq_file('../gene_files/Bactin_withTags.txt')
bactin = tagged_pois['1'][0]

rss.solver.protein = bactin #pass the protein object
t = np.linspace(0,1000,1001)
solution = rss.solver.solve_ssa(bactin.kelong, t,  n_traj= 300, ki=.013, kt = 10, low_memory=False,record_stats=True)
solver = rss.solver

X = solution.ribosome_locations

bactin_mrna_blank_mw = rss.diffcalc.calculate_rna_strand_base_mw(bactin.nt_seq, n_loops=10, )
aa_mw_vec = rss.diffcalc.calculate_single_rib_mw(bactin.aa_seq,bactin.probe_loc,)

mw_vec = rss.diffcalc.mw_over_time(X, bactin.nt_seq, bactin.probe_loc, bactin_mrna_blank_mw)
df_vec = rss.diffcalc.calculate_diffusion_constant(mw_vec)


def brownian_motion(df_vec,t ):
    traj = np.zeros( (2,) + df_vec.shape ) 
    print(traj.shape)
    current_x = traj[:,:,0]
    m = 2
    current_t = t[0]
    for i in range(0, traj.shape[-1]):
        
        x = brownian_step(current_x, df_vec[:,i], current_t, t[i] )
        current_t = t[i]
        traj[:,:,i] = x
        current_x = x
    return traj

def brownian_step(current_x, df_vec, t0, tf):
    dt = np.float(tf-t0)
    s = np.sqrt (4.0* df_vec* dt) * np.random.randn(*current_x.shape)
    return current_x + s



traj = brownian_motion(df_vec, t)
plt.plot(traj[0,:,:].T,traj[1,:,:].T,'r',alpha=.1)
print(np.mean(np.abs(traj[:,:,1:] - traj[:,:,:-1])))


solution = rss.solver.solve_ssa(bactin.kelong, t, n_traj= 300,  ki=.083, kt = 10, low_memory=False,record_stats=True)
solver = rss.solver

X = solution.ribosome_locations

bactin_mrna_blank_mw = rss.diffcalc.calculate_rna_strand_base_mw(bactin.nt_seq, n_loops=10, )
aa_mw_vec = rss.diffcalc.calculate_single_rib_mw(bactin.aa_seq,bactin.probe_loc,)

mw_vec = rss.diffcalc.mw_over_time(X, bactin.nt_seq, bactin.probe_loc, bactin_mrna_blank_mw)
df_vec = rss.diffcalc.calculate_diffusion_constant(mw_vec)


traj = brownian_motion(df_vec, t)
plt.plot(traj[0,:,:].T,traj[1,:,:].T,'b',alpha=.1)
print(np.mean(np.abs(traj[:,:,1:] - traj[:,:,:-1])))


solution = rss.solver.solve_ssa(bactin.kelong, t, n_traj= 300, ki=.183, kt = 10, low_memory=False,record_stats=True)
solver = rss.solver

X = solution.ribosome_locations

bactin_mrna_blank_mw = rss.diffcalc.calculate_rna_strand_base_mw(bactin.nt_seq, n_loops=10, )
aa_mw_vec = rss.diffcalc.calculate_single_rib_mw(bactin.aa_seq,bactin.probe_loc,)

mw_vec = rss.diffcalc.mw_over_time(X, bactin.nt_seq, bactin.probe_loc, bactin_mrna_blank_mw)
df_vec = rss.diffcalc.calculate_diffusion_constant(mw_vec)


traj = brownian_motion(df_vec, t)
print(np.mean(np.abs(traj[:,:,1:] - traj[:,:,:-1])))
plt.plot(traj[0,:,:].T,traj[1,:,:].T,'g',alpha=.1)

plt.legend()

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='r', lw=4),
                Line2D([0], [0], color='b', lw=4),
                Line2D([0], [0], color='g', lw=4)]

plt.legend(custom_lines, ['.013', '.83', '.183'])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Brownian motion vs Ki')
'''
def brownian_motion_simulation ( m = 2, n = 1001, d = 10.0, t = 1.0 ):


#
#  Set the time step.
#
  dt = t / float ( n - 1 )
#
#  Compute the individual steps.
#
  x = np.zeros ( [ m, n ] )

  for j in range ( 1, n ):
#
#  S is the stepsize
#
    s = np.sqrt ( 2.0 * m * d * dt ) * np.random.randn ( 1 )
#
#  Direction is random.
#
    if ( m == 1 ):
      dx = s * np.ones ( 1 );
    else:
      dx = np.random.randn ( m )
      norm_dx = np.sqrt ( np.sum ( dx ** 2 ) )
      for i in range ( 0, m ):
        dx[i] = s * dx[i] / norm_dx
#
#  Each position is the sum of the previous steps.
#
    x[0:m,j] = x[0:m,j-1] + dx[0:m]

  return x
'''