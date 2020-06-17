# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 11:23:46 2020

@author: willi
"""

import numpy as np
import os



from rss import rSNAPsim
from rss import ProbeVectorFactory as pvf
from rss import PropensityFactory as pff
from rss import TranslationSolvers as tss
from rss import IntensityAnalyses as ia
import matplotlib.pyplot as plt
import time


rsim = rSNAPsim()
rsim.open_seq_file('./gene_files/H2B_2tags.txt')

poi = rsim.proteins['1'][0]  #protein object

solver = tss()  #solver class
# solver.colors = 2


poi.tag_epitopes['T_Flag'] = [10,20,30,40,50,60,70]
poi.tag_epitopes['T_Hemagglutinin'] = [300,330,340,350]


solver.protein=poi





t = np.linspace(0,500,501)

taus,means,var =  solver.solve_ballistic_model(.033,10)

#with recording and low memory
sttime = time.time()
ssa_soln = solver.solve_ssa([.033] + poi.kelong + [10],t,n_traj=100)
solvetime = time.time()-sttime
print(ssa_soln.intensity_vec.shape)
plt.plot(np.mean(ssa_soln.intensity_vec[0],axis=1),color='seagreen',alpha=.4)
plt.plot(np.mean(ssa_soln.intensity_vec[1],axis=1),color='violet',alpha=.4)

mu_I1 = means[0]*len( [10,20,30,40,50,60,70])
mu_I2 = means[1]*4
plt.plot([0,500],[mu_I1,mu_I1],color='k')
plt.plot([0,500],[mu_I2,mu_I2],color='k')

plt.xlabel('time')
plt.ylabel('intensity')
print("Low memory, no recording: solved in %f seconds" % solvetime)