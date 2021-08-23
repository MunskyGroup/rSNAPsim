# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:34:04 2020

@author: willi
"""
import numpy as np
from . import IntensityAnalyses
from . import PropensityFactory, ProbeVectorFactory
import scipy as sci
import copy
import time
import pandas as pd

import pickle
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import multivariate_normal,chi2
from matplotlib.patches import Ellipse

class TranslationOptimization():
    '''
    optimization class adjusted from Micheal May's GenericOpt
    
    
    Usage: one must set a solver_obj containing a model, and a data_obj containing any data to fit
    
    ::
        opt = rss.optimizer.TranslationOptimization()  #Optimization object
        opt.solver_obj = solver
        opt.data_obj = sim_data
        opt.parnames = ['ki','ke']
        true_par = [.08,5]    
        
        Then you set the options  and optimization function you would like to fit
        all arguments for these are stored in a dictionary and then passed to run_optimization
        
        opt.opts['bounds'] = ([0.01,.17],[0.1,12])
        opt.initial_params = np.array([.033,10])
        opt.params = np.array([.033,10])
        opt.args['LL_acorr'] = (200,'ind','G0')
        opt.args['LL_I_distb'] = (1,)
        opt.run_optimization('LL_I_distb','MH',stepsize=[1,1],disp=True,mut_rate=.99,logspace=True,niter=500)
                    
    '''        
    def __init__(self):
        
        self.solver_obj = None
        self.data_obj = None
        self.int_a = IntensityAnalyses.IntensityAnalyses()
        
        self.parnames = []
        self.params = np.array([])
        self.initial_params = np.array([])
        
        self.pff = PropensityFactory.PropensityFactory()
        self.pvf = ProbeVectorFactory.ProbeVectorFactory()
        
        
        self._int_a = IntensityAnalyses.IntensityAnalyses()
        self.methods={'met_hast':self.methast,\
                      'methast':self.methast,\
                      'MH':self.methast,\
                      'lin_opt':self.lin_opt,\
                      'linopt':self.lin_opt,\
                      'genetic':self.genetic,\
                      'sim_anneal':self.sim_anneal,\
                      'simanneal':self.sim_anneal}
            
        self.objective_funs={'sse':self.sse,\
                   'SSE':self.sse,\
                   'chi':self.chisq,\
                   'chisq':self.chisq,\
                   'LL_acorr':self.get_loglikelihood_autocorr,\
                   'I_mu_sse':self.get_intensity_sse,\
                   'LL_I_distb':self.get_loglikelihood_intensity_distribution,\
                   'LL_acorr_ode':self.get_loglikelihood_autocorr_ode}
            
        self.opts={'bounds':([0,10],)*len(self.params)}  
        
        self.args = {'LL_acorr': (20,'raw','G0'), 'I_mu_sse':(), 'combined_objective':[],'LL_I_distb':(30,True)}
        
        self._obj_weights = self.args.copy()
        for key in self._obj_weights.keys():
            self._obj_weights[key] = 1
        
     
        
        self.last_run = None
        self.chain = None
        
    def print_args(objfun):
        if objfun == 'LL_acorr':
            print('LL_acorr: Log Likelihood Autocorrelation')
            print('timepoints: default 20 | how many points of decorrelation to compare the data to the model.  ')
            print('normalization: raw, ind, global | normalization for the autocorrelation based on individual trajectories or global dataset means / var')
            print('Normalization point: G0, G1, interp, max | which point to use as 1, g0 meaning the first point in the autocorrelation funciton, and g1 being the second point. Interpolation will calculate g0 from g1,g2,g3 to ignore shot noise. Max will use whatever is the maximum of the autocorrelation')
        if objfun == 'I_mu_sse':
            print("I_mu_sse: intensity mean sum of squared errors")
        if objfun == 'LL_I_distb':
            print("LL_I_distb: Log Likelihood Intensity Distribution")      
            print("number of bins: Number of bins to bin model data and exp data over")
            print("Normalization value: Number to divide all data by, for example to convert epitopes to UMP")
        
    
    def add_bounds(self,lower_bounds, upper_bounds):
        '''
        Adds bounds to the optimizer object

        Parameters
        ----------
        lower_bounds : List or ndarray
            lower bound of each parameter.
        upper_bounds : List or ndarray
            upper bound of each parameter.

        Returns
        -------
        None.

        '''
        if isinstance(lower_bounds,np.ndarray):
            lower_bounds = lower_bounds.tolist()
        if isinstance(upper_bounds,np.ndarray):
            upper_bounds = upper_bounds.tolist() 
            
        self.opts['bounds'] = tuple(list(x) for x in zip(lower_bounds,upper_bounds))
        
        
    def intensity_fun(self,x):
        self.solver_obj._poi.ke_mu = x[1]
        self.solver_obj._poi.ki = x[0]
        ssa_soln = self.solver_obj.solve_ssa_set_conditions()
        return ssa_soln.intensity_vec
    
    def autocovariance_fun(self,intensity,norm='ind'):
        acov,err_acov = self._int_a.get_autocov(intensity,norm=norm)        
        return acov,err_acov

    def autocorrelation_fun(self,intensity,norm='ind',g0='G0'):
        acov,err_acov = self._int_a.get_autocov(intensity,norm=norm)
        acorr,err_acorr = self._int_a.get_autocorr(acov,g0=g0)
        return acorr,err_acorr    

    def intensity_distribution(self,intensity,bins=None,density=True,norm=1):
        
       
        int_dist = np.histogram(intensity/norm, bins=bins, density=density,) 
        int_dist_bins = int_dist[1]
        int_dist_heights = int_dist[0]
        return int_dist_heights,int_dist_bins
    
    def analytical_autocorrelation(self,x,bins=None,bin_method='intellegent'):
        
        self.solver_obj._poi.ke_mu = x[1]
        self.solver_obj._poi.ki = x[0]       
        
        if not isinstance(bins,None):
            
            
            inds = self.pff.intellegent_bin(np.atleast_2d(self.solver_obj._poi.probe_loc),100)
            bpl,bpv = self.pvf.bin_probe_vecs(self.solver_obj._poi.probe_loc,inds)
            k_bin = self.pff.bin_k(self.solver_obj._poi.kelong, inds)
            x0 = np.zeros((k_bin.shape[0],1))
            t = self.solver_obj.t
            ode_soln_bin = self.solver_obj.solve_ode(k_bin,   t, x0, self.solver_obj._poi.ki, bpl,corr=True)
        
        
    def genetic(self,objfun,**kwargs):
        
        return sci.optimize.differential_evolution(objfun,self.opts['bounds'],**kwargs)
    
    def sim_anneal(self,objfun,**kwargs):
        
        minimizer_kwargs = {"bounds":self.opts['bounds'] }
        return sci.optimize.basinhopping(objfun,self.initial_params,minimizer_kwargs =minimizer_kwargs,**kwargs)
    
    
    def methast(self, optfun, optfun_type, niter=1000, burnin=100, stepsize=None,mut_rate=.3, disp=False, logspace=True, proposal = None):
        '''

        Parameters
        ----------
        optfun : function
            objective function to optimize.
        optfun_type : str
            the name of the objective function for __.args to access its arguments
        niter : int, optional
            Number of iterations. The default is 1000.
        burnin : int, optional
            Number of iterations to burn (disregard from the final chain). The default is 100.
        stepsize : list, ndarray, optional
            The stepsize of the proposal distribution. The default is None.
        mut_rate : float 0-1, optional
            mutation probability. The default is .3.
        disp : boolean, optional
            return a display of the current progress. The default is False.
        logspace : boolean, optional
            use logspace to search. The default is True.
        proposal : function, optional
            function to draw proposal distributions from. The default is None.

        Returns
        -------
        OptimizeResult
            scipy result object.

        '''
        if logspace:
            bounds = tuple([np.log10(x).tolist() for x in self.opts['bounds']])
            initial_par = np.log10(self.initial_params)
        else:
            bounds = self.opts['bounds']
            initial_par = self.initial_params
        
        if stepsize is None:
            stepsize = (initial_par/10).tolist()
            
        if proposal is None:
            proposal = np.random.normal
            
        
        def evolvepars(p,stepsize):
            new=np.copy(p)
            for j in range(len(new)):
                    if np.random.rand() < mut_rate:
                        while 1:
                            new[j]=p[j]+proposal(0,stepsize[j])
                    #print p
                    #print new
                            if (bounds[j][0]<new[j]<bounds[j][1]):
                                break
            return new
                        
        objective_fun_args = self.args[optfun_type]
        

        if logspace:
            oldpars=initial_par
            f_old=f_best=optfun(10**oldpars, *objective_fun_args)             
            bestpars=initial_par
        else:
            oldpars=initial_par
            f_old=f_best=optfun(oldpars, *objective_fun_args) 
            bestpars=initial_par
            
        if np.isnan(f_old):
            f_old = np.inf
        
        
    
        if disp:
            print('Burning in....')
        
        for i in range(-burnin,niter):
            if disp:
                if i > 0 and i%(niter/10) == 0:
                    self.__mh_print_report(i,bestpars,f_best,niter)
            
            newpars=evolvepars(oldpars,stepsize)

            if logspace:
                tmp_pars = copy.deepcopy(newpars)
                f_new=optfun(10**tmp_pars,*objective_fun_args)
            else:
                f_new=optfun(newpars,*objective_fun_args)
            #print("newpars = "+str(newpars)+ "     oldpars = "+str(oldpars))
            #print("fnew="+str(f_new)+", fold="+str(f_old))
            
            if np.isnan(f_new):
                f_new = np.inf
     
            if f_new<f_old or np.random.rand()<np.exp(f_old-f_new):

                #raw_input("Press Enter to continue...")
                f_old=f_new
                oldpars=newpars
                self.__update_mh_chain(newpars,f_new)
                if f_new<f_best:
                    f_best=f_new
                    bestpars=newpars
                    
        result = sci.optimize.OptimizeResult()
        
        if logspace:
            result.x = 10**bestpars
        else:
            result.x = bestpars
        
        result.fun = f_best
        result.success = True
        result.nit = niter
        result.logspace = logspace
        return result
    
    def __mh_print_report(self,iteration,bestpars,fbest,niter):
        print('current iteration: %d out of %d | best_parameters: %s | best evaulation: %f' % (iteration,niter,''.join(str(bestpars.tolist())),fbest )  )
    
    def __update_mh_chain(self,pars,funeval):
       
        self.chain.accepted_parchain = np.vstack( (self.chain.accepted_parchain, pars) )
        self.chain.accepted_evalchain = np.vstack((self.chain.accepted_evalchain,np.sum(funeval)))
        #self.chain.accepted_objfunchain  = np.vstack((self.chain.accepted_objfunchain,np.atleast_1d(funeval)))
        self.chain.accepted = self.chain.accepted + 1
        

            
    def lin_opt(self):
        return 1
    
    
    def combined_objective(self,x,objfun_list,intensity_fun):
        obj_sum = 0
        intensity = intensity_fun(x)
        obj_fun_evals = np.zeros(len(objfun_list))
        k = 0
        for i in range(len(objfun_list)):
            objective = objfun_list[i]
            obj_args = self.args[objective]
            
            obj_fun_evals[i] = self._obj_weights[objective]*self.objective_funs[objective](intensity,*obj_args)
         
            #obj_sum += self._obj_weights[objective]*self.objective_funs[objective](intensity,*obj_args)  
        obj_sum = np.sum(obj_fun_evals)
        self.__update_chain(x, obj_fun_evals )
        
        return obj_sum
    
    def run_optimization(self, objective_fun_list, method ,model = None, data = None, intensity_fun = None,**kwargs):
        if isinstance(objective_fun_list,str):
            objective_fun_list = [objective_fun_list]
        
        if model == None:
            model = self.solver_obj
        if data == None:
            data = self.data_obj
        if intensity_fun == None:
            intensity_fun = self.intensity_fun
                        
        obj_fun = self.combined_objective
                
        method_fun = self.methods[method]
        
        self.chain = OptChain()
        
        self.chain.parchain = self.initial_params
        self.chain.initial_params = self.initial_params
        self.chain.iterations = 0
        self.chain.parnames = self.parnames
        
        self.chain.bestpar = None
        self.chain.besteval = None
        self.chain.opt_method = method
        

        self.chain.evalchain  = np.array([11110])
        self.chain.objfunchain = np.zeros((1,len(objective_fun_list)))
        
        starttime = time.time()
        
        args = (objective_fun_list,intensity_fun)
        
        if method in ['met_haste','methaste','MH']:
            self.args['combined_objective'] = [objective_fun_list,intensity_fun]
            self.chain.accepted = 0
            self.chain.accepted_evalchain = np.array([11110])
            self.chain.accepted_objfunchain =  np.zeros((1,len(objective_fun_list)))
            self.chain.accepted_parchain = self.initial_params
            result = method_fun(obj_fun,'combined_objective', **kwargs)
            self.chain.accepted_evalchain = self.chain.accepted_evalchain[1:]
            #self.chain.accepted_objfunchain =  self.chain.accepted_objfunchain[1:,:]
            self.chain.accepted_parchain = self.chain.accepted_parchain[1:,:]     
            
        else:
            kwargs['args'] = args
            result =  method_fun(obj_fun,**kwargs)
            
        self.chain.runtime = time.time()-starttime
        
        #result =  method_fun(obj_fun,**kwargs)
        
        self.chain.bestpar = result.x
        self.chain.besteval = result.fun
        
        self.chain.parchain = self.chain.parchain[1:,:]
        self.chain.objfunchain = self.chain.objfunchain[1:,:]
        self.chain.evalchain  = self.chain.evalchain[1:]
        self.chain.objective_fun_list = objective_fun_list
        self.chain.objective_args = [self.args[x] for x in objective_fun_list ] 
        self.chain.opt_args = kwargs

    def __update_chain(self,pars, funeval):   

  
        self.chain.parchain = np.vstack( (self.chain.parchain, pars) )

        self.chain.evalchain = np.vstack((self.chain.evalchain,np.sum(funeval)))
        self.chain.objfunchain  = np.vstack((self.chain.objfunchain,np.atleast_1d(funeval)))
        self.chain.iterations = self.chain.iterations + 1
        
        
    def __intensity_generator(self,pars,objective_fun_list):
        intensity = self.intensity_fun(pars)
        
        return intensity
    
    def get_loglikelihood_autocorr(self, intensity, n_points,norm,g0):
        
        model_acorr,model_acorr_err = self.autocorrelation_fun(intensity,norm=norm,g0=g0)

        total_n_spots = self.solver_obj.n_traj
        data_autocorrelation = self.data_obj.acorr
        data_acc_err = self.data_obj.acorr_err
        LL = self.loglikelihood_acc(model_acorr[:,:n_points,:],  data_autocorrelation[:,:n_points,:], data_acc_err[:,:n_points], total_n_spots)
    
        
        
        return LL

    def get_loglikelihood_autocorr_ode(self, intensity, n_points,norm,g0):
        
        model_acorr,model_acorr_err = self.autocorrelation_fun(intensity,norm=norm,g0=g0)

        total_n_spots = self.solver_obj.n_traj
        data_autocorrelation = self.data_obj.acorr
        data_acc_err = self.data_obj.acorr_err
        LL = self.loglikelihood_acc(model_acorr[:,:n_points,:],  data_autocorrelation[:,:n_points,:], data_acc_err[:,:n_points], total_n_spots)
    
       
        
        return LL

    def get_intensity_sse(self,intensity):
        
        return (np.mean(intensity) - np.mean(self.data_obj.I_mu))**2
    
    def get_loglikelihood_intensity_distribution(self,intensity,norm):
        
        dist_sim_data = self.intensity_distribution(intensity,bins=self.data_obj.histogram_bins,density=True,norm=norm)[0]
        dist_sim_data[dist_sim_data==0] = 1e-7
        
        LL = -np.dot(self.data_obj.histogram,np.log(dist_sim_data))
        
        if LL == -np.inf:
            LL = np.inf
        return LL
        
    @staticmethod
    def sse(model_data,data):        
        return np.sum((model_data-data)**2)
    
    @staticmethod
    def chisq(model_data,model_var, data):
        return np.sum(((model_data-data)**2)/model_var)
    
    
    @staticmethod
    def loglikelihood_acc( model, data,data_err, nspots):
        '''
        Parameters
        ----------
        model : ndarray
            autocorrelation model generated data array.
        data : ndarray
            autocorrelation data array.
        data_err : ndarray
            error array (SEM, STD) of the data.
        nspots : int
            number of spots.

        Returns
        -------
        Loglikelihood (float)
            returns the logliklihood from the formula:
                
                log L(G|M) = Const. - 1/Nt * sum( (G_D(t_i) - G_M(t_i))^2 / sigma_G_D(t_i) )

        '''
        
        d = np.mean(data,axis=-1)
        m = np.mean(model,axis=-1)
       

    
        return (nspots/2) * np.sum(( d[:,2:] - m[:,2:])**2/ data_err[:,2:])
    
    @staticmethod
    def loglikelihood_distb( model_intensity, data_intensity,nbins=30,norm=1):
        '''
    
        Parameters
        ----------
        model_intensity : ndarray
            array of model simulated intensity over time.
        data_intensity : ndarray
            array of data intensity over time.
        nbins : int, optional
            Number of bins. The default is 20.

        Returns
        -------
        float
            LogLikelihood of the intensity distributions.

        '''
        dist_sim_data = np.histogram(model_intensity/norm, bins=nbins, density=True)[1]
        hist_exp_data = np.histogram(data_intensity, bins=nbins)
        return -np.dot(hist_exp_data,np.log(dist_sim_data))
    

class IntensityData():
    def __init__(self):
      
        self.ragged = False
        self.head = None
        self.ssa_obj = None    
        
    
    def add_data(self, t, intensity_vec):
        '''

        Parameters
        ----------
        intensity_vec : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.intensity_vec = intensity_vec
        if isinstance(intensity_vec,np.ndarray):
            self.ragged=False
            
        else:
            self.ragged=True
        self.times = t
        
         
    def get_stats(self):
        if self.ragged == False:
            self.I_mu = np.mean(self.intensity_vec,axis=2)
            self.I_var = np.mean(self.intensity_vec,axis=2)
        else:
            spot_no = len(self.intensity_vec)
            def autocorr(x):
                result = np.correlate(x, x, mode='full')
                return result[result.size // 2:]        
            # Calculating the autocovariance
            n_lags = 300
            n_selected_Particles = np.max(spot_no)
            ac_array = np.zeros((n_selected_Particles , n_lags ))
            norm_ac_array = np.zeros((n_selected_Particles , n_lags ))
            counter = 0
            # Calculating Autocorrelation.
            for i in range(1,n_selected_Particles+1):
                intensity_green = self.intensity_vec[i][0]
                temp_ac = autocorr(intensity_green)
                size_ac_temp  = len(temp_ac)
                ac_array[counter, 0:size_ac_temp] = autocorr(intensity_green)
                norm_ac_array[counter, 0:size_ac_temp] = autocorr(intensity_green)/ float(ac_array[counter, :] .max())
                counter += 1
            # Plotting mean autocovariance
            lag_time = [*range(0, self.max_lag_output, 1)]
            lag_time = [element * self.sampling_rate for element in lag_time]
            mean_ac_data = norm_ac_array.mean(0)
            std_ac_data = norm_ac_array.std(0)
            # normalized autocovariance, removing shot noise
            
            if self.remove_shotnoise ==1:
                #mean_ac_data_norm = mean_ac_data[1:-1]/mean_ac_data[1]
                return lag_time[1:self.max_lag_output], mean_ac_data[1:self.max_lag_output], std_ac_data[1:self.max_lag_output]
            else:
                #mean_ac_data_norm = mean_ac_data[1:-1]/mean_ac_data[0]
                return lag_time[0:self.max_lag_output], mean_ac_data[0:self.max_lag_output], std_ac_data[0:self.max_lag_output]
            
        
    def load_data(self, file):
        extension = file.split('.')[-1]
        if extension in ['xls', 'xlsx']:
            df = pd.read_excel(file)
            
            self.head = df.head()
            n_cells = max(df['Cell_#'])
            max_spots =  max(df['Spot_#'])
            total_n_spots = 0
            data_struct = []
            for i in range(0,n_cells):
                i_g = []
                i_r = []
                intensitys = []
                times = []
                for j in range(0,max_spots):
                    temp_len_spot =  df[(df['Cell_#'] == i) & (df['Spot_#'] ==j) ]
                    if len(temp_len_spot) >0:
                        time = df[(df['Cell_#'] ==i) & (df['Spot_#'] ==j)]['Time_(sec)'].values
                        intensity_red = df[(df['Cell_#'] ==i) & (df['Spot_#'] ==j)]['R_GaussianInt'].values
                        intensity_green = df[(df['Cell_#'] ==i) & (df['Spot_#'] ==j)]['G_GaussianInt'].values
            
                        i_g.append(intensity_green)
                        i_r.append(intensity_red)
                        
                        intensity = np.vstack((intensity_red,intensity_green))
                        intensitys.append(intensity)
                        times.append(time)
                data_struct.append([i_g,i_r,intensitys,times ])
                
            time_arrays = []
            intensity_arrays = []
            cells = []
            for n in range(len(data_struct)):
                cells = cells + [n for x in range( len(data_struct[n][3])  )]
                time_arrays = time_arrays + data_struct[n][3]
            
            for n in range(len(data_struct)):
                intensity_arrays = intensity_arrays + data_struct[n][2]
            
            len_t = len(time_arrays[0])
            self.ragged = False
            for n in range(len(time_arrays)):
                if len(time_arrays[n]) != len_t:
                    self.ragged = True
                    
            self.intensity_vec = intensitys
            self.times = time_arrays
            self.cells = cells
            
            
    def __guess_headers(self,dataframe):
        columns = dataframe.columns.tolist()
        
        cell_header = None; spot_header = None; time_header = None;
        R_header=None; G_header = None; B_header=None;
    
        kwargs = {}
        for col in columns:
            if 'spot' in col.lower():
                spot_header = col
                kwargs['spot_header'] = spot_header
            if np.sum([x in col.lower() for x in ['green','g_']]) >0:
                green_header = col
                kwargs['green_header'] = green_header
            if np.sum([x in col.lower() for x in ['red','r_']]) >0:
                red_header = col
                kwargs['red_header'] = red_header
            if np.sum([x in col.lower() for x in ['blue','b_']]) >0:
                blue_header = col  
                kwargs['blue_header'] = blue_header
            if 'cell' in col.lower():
                cell_header = col
                kwargs['cell_header'] = cell_header
            if np.sum([x in col.lower() for x in ['time','sec']]) >0:
                time_header = col   
                kwargs['time_header'] = time_header
        
        return kwargs
        
    
    def load_dataframe(self,dataframe):
        kwargs = self.__guess_headers(dataframe)
        print(kwargs)
        self.__load_dataframe(dataframe, **kwargs)
        
            
    def __load_dataframe(self, dataframe, cell_header = 'Cell_No', spot_header = 'Spot_No',time_header = 'Time_sec', red_header='Red_Int', green_header = 'Green_int', blue_header='Blue_int' ):
            df = dataframe
            
            self.head = df.head()
            n_cells = max(df[cell_header])+1
            max_spots =  max(df[spot_header])+1
            total_n_spots = 0
            data_struct = []
            

            for i in range(0,n_cells):
                i_g = []
                i_r = []
                i_b = []
                intensitys = []
                times = []
                for j in range(0,max_spots):
                    temp_len_spot =  df[(df[cell_header] == i) & (df[spot_header] ==j) ]
    
                    if len(temp_len_spot) >0:
                        time = df[(df[cell_header] ==i) & (df[spot_header] ==j)][time_header].values
                        intensity_red = df[(df[cell_header] ==i) & (df[spot_header] ==j)][red_header].values
                        intensity_green = df[(df[cell_header] ==i) & (df[spot_header] ==j)][green_header].values
                        intensity_blue = df[(df[cell_header] ==i) & (df[spot_header] ==j)][blue_header].values
                        i_g.append(intensity_green)
                        i_r.append(intensity_red)
                        i_b.append(intensity_blue)
                        
                        intensity = np.vstack((intensity_red,intensity_green,intensity_blue))
                        intensitys.append(intensity)
                        times.append(time)
                data_struct.append([i_r,i_g,i_b,intensitys,times ])
         
            time_arrays = []
            intensity_arrays = []
            cells = []
            for n in range(len(data_struct)):
                cells = cells + [n for x in range( len(data_struct[n][4])  )]
                time_arrays = time_arrays + data_struct[n][4]

            len_t = len(time_arrays[0])
            self.ragged = False
            for n in range(len(time_arrays)):
                if len(time_arrays[n]) != len_t:
                    self.ragged = True
                    
            self.intensity_vec = intensitys
            self.data_struct = data_struct
            self.times = time_arrays
            self.cells = cells               


class OptChain():
    '''
    Container class for parameter optimization chains
    '''
    def __init__(self):
        self.parchain = np.array([])
        self.iterations = 0
        self.parnames = np.array([])
        self.evalchain  = np.array([])
        self.objfunchain = np.array([])
        self.bestpar = None
        self.besteval = None
        self.opt_method = None
        self.opt_args = None
        self.objective_args = None
        self.objective_fun_list  = None
        self.intensity_fun = None
        self.logspace = False
        
    def report(self):
        print('=====================')
        print('Optimizer: %s ran for %d iterations ' % (self.opt_method,self.iterations))
        print('Optimizer arguments: ' + str(self.opt_args))
        print('Objective function: ' + str(self.objective_fun_list))
        print('Objective arguments: ' + str(self.objective_args))
        print('_____________________')
        print('Best Parameter Set: %s, feval: %d'%  (''.join(str(self.bestpar.tolist())),self.besteval ) )
        print('=====================')
        
    def check_parameter_convergence(self):
        
        def get_acc2(data, trunc=False):
            '''
            Get autocorrelation function
    
            *NOT* multi-tau
            '''
            N = len(data)
            fvi = np.fft.fft(data, n=2*N)
            acf = fvi*np.conjugate(fvi)
            acf = np.fft.ifft(acf)
            acf = np.real(acf[:N])/float(N)
            if trunc:
                acf[acf < 0]=0
                for i in range(1, len(acf)):
                    if acf[i] > acf[i-1]:
                        acf[i] = acf[i-1]
            return acf    
        
        for i in range(0,self.parchain.shape[1] ):
            acc = get_acc2((self.parchain[:,i] - np.mean(self.parchain[:,i] ))  / np.var(self.parchain[:,i] ) )       
            if i == 0:
                self.par_acc = acc
            else:
                self.par_acc = np.vstack( (self.par_acc, acc))
        self.par_acc = self.par_acc.T
            
    def check_objfun_convergence(self):
        
        def get_acc2(data, trunc=False):
            '''
            Get autocorrelation function
    
            *NOT* multi-tau
            '''
            N = len(data)
            fvi = np.fft.fft(data, n=2*N)
            acf = fvi*np.conjugate(fvi)
            acf = np.fft.ifft(acf)
            acf = np.real(acf[:N])/float(N)
            if trunc:
                acf[acf < 0]=0
                for i in range(1, len(acf)):
                    if acf[i] > acf[i-1]:
                        acf[i] = acf[i-1]
            return acf    
        
        for i in range(0,self.objfunchain.shape[1] ):
            acc = get_acc2((self.objfunchain[:,i]  -   np.mean(self.objfunchain[:,i]  )) / np.var(self.objfunchain[:,i] )  )       
            if i == 0:
                self.objfun_acc = acc
            else:
                self.objfun_acc = np.vstack( (self.objfun_acc, acc))
                
        self.objfun_acc = self.objfun_acc.T
        
    
    def __clear_invalid_values(self):
        evalchain = self.evalchain
        parchain = self.parchain
        trimmed_evals = evalchain[~np.isnan(evalchain)]        
        final_evals = trimmed_evals[np.isfinite(trimmed_evals)]
        trimmed_parchains = parchain[~np.isnan(evalchain)[:,0]]
        final_parchains = trimmed_parchains[np.isfinite(trimmed_evals)]
        
        
        return final_evals, final_parchains
        
        
    def parplot(self,ellipse=True,logspace=False):
        n_par = len(self.bestpar)
       
        fig, ax = plt.subplots( n_par,n_par,dpi=300)
        
        plotnum = np.arange(n_par**2).reshape((n_par,n_par)) + 1
        triangle_inds = np.tril(plotnum)[np.where(np.tril(plotnum) !=0 )]
        nplots = len(triangle_inds)
        
        
        if np.sum(np.isnan(self.evalchain))>0 or np.sum(~np.isfinite(self.evalchain)) >0 :
            print('Warning: NaN or Infinite values detected within function evaluation chain, these parameter sets are left out of the plot ')
        eval_chain,par_chain = self.__clear_invalid_values()
         
        
        viridis = cm.get_cmap('viridis', int(np.ceil(np.max(eval_chain))))
        colors = eval_chain
        covariances = np.cov(par_chain.T)
        used_pairs = []
        for i in range(n_par-1,-1,-1):
            for j in range(n_par-1,-1,-1):
                if set([i,j]) not in used_pairs:
                    if i != j:
                        a = ax[i][j].scatter(par_chain[:,i],par_chain[:,j], marker='.',c= colors)
                        if ellipse == True:
                            self.__get_ellipse([np.mean(par_chain[:,i]),np.mean(par_chain[:,j])], covariances,ax=ax[i][j] )
                                         
                    else:
                        b = ax[i][j].hist(par_chain[:,i],bins=40,density=True)
                        
           
                    used_pairs.append(set([i,j]))
                else: 
                    ax[i][j].axis('off')
                
                if len(self.parnames) >0:
                    if i == n_par-1:
                        ax[i][j].set_xlabel(self.parnames[i])
                    if j == 0:
                        ax[i][j].set_ylabel(self.parnames[j])
                        
                if logspace:
                    if i != j:
                        ax[i][j].set_yscale('log')
                    ax[i][j].set_xscale('log')
 
        fig.tight_layout()
        fig.colorbar(a, ax=ax)
        fig.show()
        
    def __get_ellipse(self,mu,cov,ax=None,crosshairs=False,*args,**kwargs):
        ''' Command to plot the fit instance 
        '''
        cmap= cm.viridis
        cmap = cm.coolwarm
        ci = .95
        vals,vecs = np.linalg.eig(cov)
        theta = float(abs((360/(2*np.pi))*np.arctan(vecs[1,0]/vecs[0,0])))
        # able to change CI now.
        scale = chi2.ppf(ci,2)
        w = np.sqrt(vals[0]*scale)*2
        h = np.sqrt(vals[1]*scale)*2
        

        e = Ellipse(xy = tuple(mu) ,width = w, height = h,angle=90-theta,**kwargs )

        ax.add_artist(e)
    
        e.set_clip_box(ax.bbox)
        # e.set_alpha(.75)
        e.set_edgecolor(('red'))
        e.set_linestyle(('-'))
        e.set_facecolor(('none'))
        return ax
    
    def save(self,filename):
        ext = filename.split('.')[-1]
        if ext == 'txt':
            x=1
        if ext in ['p','pickle']:
            pickle.dump(self,open('filename','wb'))
        if ext == 'csv':
            x=1
        if ext == 'npz':
            x=1
            
    def load(self,filename):
        ext = filename.split('.')[-1]
        if ext == 'txt':
            x=1
        if ext in ['p','pickle']:
            tmp_obj = pickle.load(open('filename','wb'))
            for item in tmp_obj.__dict__.keys():
                self.__dict__[item] = tmp_obj.__dict__[item]
            
        if ext == 'csv':
            x=1
        if ext == 'npz':
            x=1
                    

class suite():
    '''
    a coterie of constructs
    '''

    def __init__(self):
        self.pois = []
        self.discernable = False
        self.models = []
        self.combo_list = []
        