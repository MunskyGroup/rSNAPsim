# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:49:42 2020

@author: William Raymond
"""

import numpy as np
class IntensityAnalysesRagged():
    def __init__(self):
        pass
    
    @staticmethod
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
    
    def get_g0(self,covariance, mode = 'interp'):
        '''
        
        '''
        if mode.lower() in ['interp','inter','extrapolate','interpolate']:
            X = [1,2,3,4]
            G0 = np.zeros((covariance.shape[0], covariance.shape[2]))
            for i in range(covariance.shape[0]):
                for n in range(covariance.shape[0]):
                    V = covariance[i,X,n]
                    print(V.shape)
                    G0[i,n] = np.interp(0,X,V)
           
            
        if mode.lower() in ['g1','1']:
            G0 = covariance[:,1,:]
            
        if mode.lower() in ['g0','0']:
            G0 = covariance[:,0,:]

        if mode.lower() in ['max','maximum']:
            G0 = np.max(covariance,axis=1)
            
        return G0
    
    def get_autocorr(self,autocov,g0='G0'):
        '''
        normalize the autocovariance over g0
        '''
        autocorr = np.copy(autocov)
        n_traj = autocorr.shape[-1]
        autocov_error =  1.0/np.sqrt(n_traj)*np.std(autocov,ddof=1,axis=2)
        
        g0 = self.get_g0(autocov,g0)
        for n in range(autocov.shape[0]):
            autocorr[n] = autocorr[n]/g0[n]
        
        err_autocorr =  1.0/np.sqrt(n_traj)*np.std(autocorr,ddof=1,axis=2)
        return autocorr,err_autocorr
                

    def get_autocov(self,intensity_vec,max_lag,norm='raw'):
        
        colors, n_traj = self.__get_ivec_dims(intensity_vec)
        
        autocorr_vec = np.zeros((colors, max_lag, n_traj    ) )
        autocorr_err = np.zeros((colors, max_lag, n_traj  ) )
        

        
        
        for n in range(colors):
            if norm in [ 'Individual','I','individual','ind']:
                for i in range(n_traj):
                    ivec = intensity_vec[i][n]
                    print(ivec.shape)
                    trace_len = min(max_lag, len(ivec) )
                    autocorr_vec[n,:trace_len,i] = self.get_acc2( (ivec-np.mean(ivec))/np.var(ivec)  )[:trace_len]
                    
            elif norm in ['global','Global','g','G']:
                means = []
                var = []
                for i in range(n_traj):
                    ivec = intensity_vec[i][n]
                    means.append( np.mean(ivec) )
                    var.append(np.var(ivec))         
                    
                global_mean = np.mean(means)
                global_var = np.var(var)
                
                for i in range(n_traj):
                    ivec = intensity_vec[i][n]
                    trace_len = min(max_lag, len(ivec) )
                    autocorr_vec[n,:trace_len,i] = self.get_acc2((ivec-global_mean)/global_var )[:trace_len]
            elif norm in ['raw','Raw']:
                for i in range(n_traj):
                    ivec = intensity_vec[i][n]
                    trace_len = min(max_lag, len(ivec) )
                    autocorr_vec[n,:trace_len,i] = self.get_acc2(ivec)[:trace_len]  
                    
            # elif norm in ['max_ind','Max','maximum','Maximum','max']:
            #     for i in range(n_traj):
            #         ivec = intensity_vec[i][n]
            #         trace_len = min(max_lag, len(ivec) )
            #         autocorr_vec[n,:trace_len,i] = self.get_acc2(ivec)[:trace_len]
            #         autocorr_vec[n,:trace_len,i] = autocorr_vec[n,:trace_len,i]/np.max(autocorr_vec[n,:trace_len,i])
                    
            else:
                print('unrecognized normalization, please use individual, global, or none')
                return

        autocorr_err =  1.0/np.sqrt(n_traj)*np.std(autocorr_vec,ddof=1,axis=2)
                
        return autocorr_vec, autocorr_err
    
    
    def get_intensity_histogram(self,intensity_vec,n_bins,scale_factor,time_slice=1):
        colors, n_traj = self.__get_ivec_dims(intensity_vec)
        hist_list = []
        hist_edges_list = []
        for n in range(colors):
            for i in range(n_traj):
                if i == 0:
                    ivec = intensity_vec[i][n][::time_slice]
                    
                else:
                    tmp_vec = intensity_vec[i][n][::time_slice]
                    
                    ivec = np.hstack((ivec, tmp_vec))
                    
                    
            exp_int_dist_ump = np.divide(ivec,scale_factor)
            exp_int_dist_ump[exp_int_dist_ump<0] = 0
            exp_int_dist_ump[np.isnan(exp_int_dist_ump)] = 0
            
            temp_hist = np.histogram(exp_int_dist_ump, bins=n_bins)
            hist_list.append(temp_hist[1])
            hist_edges_list.append(temp_hist[0])
            
        hist_data = np.array(hist_list)
        hist_bins = np.array(hist_edges_list)   
        
        return hist_bins, hist_data
            
    
    def __get_ivec_dims(self,ivec):
        colors = ivec[0].shape[0]
        n_traj = len(ivec)           
        return colors,n_traj
     