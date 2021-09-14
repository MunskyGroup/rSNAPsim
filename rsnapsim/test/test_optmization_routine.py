# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 16:13:17 2021

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
Bactin_obj = tagged_pois['1'][0]

rss.solver.protein = Bactin_obj #pass the protein object
t = np.linspace(0,1000,1001)
solution = rss.solver.solve_ssa(Bactin_obj.kelong, t, ki=.033, kt = 10, low_memory=False,record_stats=True)
solver = rss.solver


ki = [.08]
ke = [ 5]
colors = ['r']
i = 0
      
Bactin_obj.ke_mu = ke[i] #set the protein object average kelong
solver.default_conditions['burnin'] = 1000 #lets burnin for 1000s
solver.t = np.linspace(0,2000,2001)
solver.n_traj = 30
      
t = np.linspace(0,2000,2001)
ssa_soln = solver.solve_ssa( Bactin_obj.kelong,t, kt=10, ki= ki[0],n_traj=30)  
      
plt.plot(ssa_soln.intensity_vec[0],color=colors[i],alpha=.4)
plt.title('Intensity Trajectories')
plt.ylabel('Time')      
      
    
sim_data = rss.optimizer.IntensityData()      #create a new data object
sim_data.add_data(ssa_soln.time,ssa_soln.intensity_vec) #add some intensity data to it
sim_data.get_stats()  #generate stats such as intensity mu and var

      
sim_acov, sim_acov_err = rss.inta.get_autocov(sim_data.intensity_vec,norm='ind')
sim_acc, sim_acc_err = rss.inta.get_autocorr(sim_acov)
sim_data.acorr = sim_acc
sim_data.acorr_err = sim_acc_err
sim_data.histogram = np.histogram(sim_data.intensity_vec,bins=30)[0]
sim_data.histogram_bins = np.histogram(sim_data.intensity_vec,bins=30)[1]

opt = rss.optimizer.TranslationOptimization()  #Optimization object
opt.data_obj = sim_data
opt.parnames = ['ki','ke']
true_par = [.08,5]

def model_fun(parameters):
    Bactin_obj.ke_mu = parameters[1]
    soln = rss.solver.solve_ssa(Bactin_obj.kelong,
                                t,
                                ki=parameters[0],
                                kt = 10,
                                low_memory=True,
                                record_stats=False,
                                n_traj=30)
    return soln



opt.opts['bounds'] = ([0.01,.17],[0.1,12])
opt.initial_params = np.array([.033,10])
opt.params = np.array([.033,10])
opt.args['LL_acorr'] = (200,'ind','G0')
opt.args['LL_I_distb'] = (1,)
opt.run_optimization('LL_I_distb','MH', model_fun, stepsize=[.1,.1],
                     disp=True,
                     mut_rate=.5,
                     logspace=True,
                     niter=500,)
                     


chain = opt.chain
intensity = opt.intensity_fun(chain.bestpar)
fit_acorr,fit_acorr_error = opt.autocorrelation_fun(intensity)
fit_acorr = np.mean(fit_acorr[0],axis=1)

plt.plot(fit_acorr,'#1cfff7',lw=3)
plt.plot(fit_acorr-fit_acorr_error[0],'--', color = '#1cfff7',label='_nolegend_')
plt.plot(fit_acorr+fit_acorr_error[0],'--', color = '#1cfff7',label='_nolegend_')


plt.plot(np.mean(opt.data_obj.acorr[0],axis=1),'r',lw=3)
plt.plot(np.mean(opt.data_obj.acorr[0],axis=1) -opt.data_obj.acorr_err[0] ,'r--',lw=1,label='_nolegend_')
plt.plot(np.mean(opt.data_obj.acorr[0],axis=1) +opt.data_obj.acorr_err[0] ,'r--',lw=1,label='_nolegend_')

plt.plot([200,200],[-.3,1],'k--')
plt.legend(['Model','Data'])
plt.ylabel('tau')
plt.xlabel('Autocorrelation')

plt.figure()
int_hist = opt.intensity_distribution(intensity,bins = sim_data.histogram_bins)
counts1, bins = np.histogram(intensity,density=True,bins=20)
plt.hist(bins[:-1], bins, weights=counts1)

counts2, bins = np.histogram(sim_data.intensity_vec,density=True,bins=20)
counts3, bins3 = np.histogram(sim_data.intensity_vec,bins=20)
plt.hist(bins3[:-1], bins3, weights=counts2,fill=False,edgecolor='r')
plt.ylabel('Probability')
plt.xlabel('Intensity')

plt.figure()
plt.plot(intensity[0],'#1cfff7',alpha=.3)
plt.plot(sim_data.intensity_vec[0],'r',alpha=.3)


plt.xlabel('Time')
plt.ylabel('Intensity')

np.dot(counts3,np.log(counts1))
#np.sum(np.log(counts3*counts1))


