# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:50:56 2020

@author: William Raymond
"""
import numpy as np
from . import custom_errors as custom_err


class IntensityAnalyses():
    def __init__(self):
        pass

    @staticmethod
    def get_acc2(data, trunc=False):
        '''
        return autocorrelation from a data vector via Fourier transform

        Parameters
        ----------
        data : ndarray
            numpy array of data to get the correlation from.
        trunc : bool, optional
            Remove non zero entries and only return decreasing entries.
            The default is False.

        Returns
        -------
        autocorrelation : ndarray
            autocorrelation function.

        '''

        N = len(data)
        fvi = np.fft.fft(data, n=2*N)
        acf = fvi*np.conjugate(fvi)
        acf = np.fft.ifft(acf)
        acf = np.real(acf[:N])/float(N)
        if trunc:
            acf[acf < 0] = 0
            for i in range(1, len(acf)):
                if acf[i] > acf[i-1]:
                    acf[i] = acf[i-1]
        return acf



    # def get_frags_ssa_obj(self,ssa_obj, total_length = None):

    #     fragmented_trajectories = []
    #     fragtimes = []
    #     maxlen = 0

    #     fragmentspertraj= []
    #     n_traj = ssa_obj.n_traj

    #     solutions = ssa_obj.solutions
    #     if total_length == None:
    #         total_length = 0
    #         for k in range(n_traj):
    #             if total_length < np.max(solutions[k]):
    #                 total_length = np.max(solutions[k])

    #     nsteps = solutions[0].shape[1]
    #     for k in range(n_traj):
    #         ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
    #         changes = ind[1:] - ind[:-1]
    #         addindexes = np.where(changes > 0)[0]
    #         subindexes = np.where(changes < 0)[0]

    #         sub = solutions[k][:,1:] - solutions[k][:,:-1]
    #         neutralindexes = np.unique(np.where(sub < 0)[1])
    #         neutralindexes = np.setxor1d(neutralindexes, subindexes)

    #         for index in neutralindexes:
    #             pre = solutions[k][:,index]
    #             post = solutions[k][:,index+1]
    #             changecount = 0
    #             while len(np.where(post - pre < 0)[0]) > 0:

    #                 post = np.append([total_length],post)
    #                 pre = np.append(pre,0)

    #                 changecount+=1

    #             for i in range(changecount):
    #                 addindexes = np.sort(np.append(addindexes,index))
    #                 subindexes = np.sort(np.append(subindexes,index))

    #             changes[index] = -changecount
    #             ind[index] += changecount


    #         for index in np.where(np.abs(changes)>1)[0]:
    #             if changes[index] < 0:
    #                 for i in range(np.abs(changes[index])-1):
    #                     subindexes = np.sort(np.append(subindexes,index))
    #             else:
    #                 for i in range(np.abs(changes[index])-1):
    #                     addindexes = np.sort(np.append(addindexes,index))

    #         truefrags = len(subindexes)




    #         if len(subindexes) < len(addindexes):
    #             subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(nsteps-1)).astype(int))


    #         fragmentspertraj.append(len(subindexes))

    #         for m in range(min(len(subindexes),len(addindexes))):
    #             traj = solutions[k][:, addindexes[m]:subindexes[m]+1]
    #             traj_ind = changes[addindexes[m]:subindexes[m]+1]

    #             startind = ind[addindexes[m]]
    #             minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
    #             fragment = np.array([])



    #             iterind = startind

    #             if subindexes[m]-addindexes[m] > 0:
    #                 if len(minusloc) > 1:
    #                     if m <= truefrags:
    #                         for n in range(len(minusloc)-1):
    #                             iterind = iterind + min(0,traj_ind[minusloc[n]])
    #                             fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten())







    #                         fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())

    #                     else:
    #                         for n in range(len(minusloc)-1):

    #                             iterind = iterind + min(0,traj_ind[minusloc[n]])

    #                             fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten())


    #                         fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())



    #                 else:

    #                     fragment = solutions[k][startind][addindexes[m]:subindexes[m]+1].flatten()



    #                 fragtimes.append(addindexes[m]+1)


    #                 fragmented_trajectories.append(fragment)
    #                 #if m <= truefrags:
    #                     #kes.append(genelength/truetime[len(fragment)])

    #                 if len(fragment) > maxlen:
    #                     maxlen = len(fragment)


    #         fragarray = np.zeros((len(fragmented_trajectories), maxlen))
    #         for i in range(len(fragmented_trajectories)):
    #             fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]


    #     ssa_obj.fragments = fragarray
    #     ssa_obj.fragtimes = fragtimes
    #     ssa_obj.frag_per_traj = fragmentspertraj
    #     ssa_obj.full_frags = truefrags


    ##TODO Convert excitation of the intensity vector (realistic probe)

    def get_g0(self, covariance, mode='interp'):
        '''
        return the normalization point for autocorrelations, g0 delay
    
        Parameters
        ----------
        correlation : ndarray
            numpy array of a cross or autocovariance.
        mode : string, optional
            the type of G0 shot noise to return,

            * 'Interp' - will interpolate the g0 position from the G1, G2, and G3 points.
            
            * 'g1' - (second point) g1 will be returned
            
            * 'g0' - g0 will be returned (first point)
            
            * 'max' -maximum of the correlation will be returned
            
            The default is 'interp'.

        Returns
        -------
        G0 : float
            point to normalize correlation over.

        '''
        if mode.lower() in ['interp', 'inter', 'extrapolate', 'interpolate']:
            X = [1, 2, 3, 4]
            V = covariance[:, X, :]
            G0 = np.interp(0, X, V)

        if mode.lower() in ['g1', '1']:
            G0 = covariance[:, 1, :]

        if mode.lower() in ['g0', '0']:
            G0 = covariance[:, 0, :]

        if mode.lower() in ['max', 'maximum']:
            G0 = np.max(covariance, axis=1)
        return G0

    def normalize_cc(self, correlation, mode='global_max'):
        '''
        normalize a cross correlation

        Parameters
        ----------
        correlation : ndarray
            numpy array of the cross correlation to .
        mode : str, optional
            The type of normalization. The default is 'max'.

            * 'max' - normalize by maximum point
            * 'middle' - maximize by the middle point

        Returns
        -------
        norm_cor : ndarray
            normalized cross correlation.

        '''
        if mode.lower() in ['global_max', 'global_maximum', 'gmax']:
            norm_cor = correlation/np.max(correlation)

        if mode.lower() in ['indiv_max', 'indiv_maximum', 'imax',
                            'individual_maximum', 'individual_max']:
            norm_cor = correlation/np.max(correlation, 1)[:, np.newaxis, :]

        if mode.lower() in ['indiv_center', 'indiv_middle',
                            'individual_middle',
                            'individual_center']:
            centerpoint = int((correlation.shape[1]+1)/2)-1
            norm_cor = correlation/(
                correlation[:, centerpoint][:, np.newaxis, :])

        if mode.lower() in ['global_center', 'global_middle']:
            centerpoint = int((correlation.shape[1]+1)/2)-1
            norm_cor = correlation/(np.mean(correlation[:, centerpoint]))
        if mode.lower() in ['raw']:
            norm_cor = correlation

        return norm_cor


    def get_crosscorr(self, intensity_vecs, norm='indiv', g0='indiv_center', scale_fix=False):
        '''
        return a cross correlation tensor from an intensity tensor of (ncolor, ntime, ntraj)
    
        Given signals :math:`Signals X(t), Y(t)`:
            
        .. math:: 
            ACOV(t, \tau) = cov(X_{t}, Y_{\tau}) = E{X_{t}, Y_{\tau}} 
        
        Parameters
        ----------
        intensity_vecs : ndarray
            (ncolor, ntime, ntraj) intensity tensor.
        norm : str, optional
            Normalization to apply,
            
            * global - subtract by the global intensity mean for the correlation
            
            * indiv - subtract by each trajectory's mean for the correlation
            
            * raw - do no mean subtraction
            
            The default is 'indiv'.
        g0 : str, optional
            point to normalize the correlation by after correlating
            
            * global_middle - divide the correlation by the average of the center point of all trajectories
            
            * indiv_middle - divide the correlation trajectories by their individual center point
            
            * global_max - divide the correlation trajectories by the global max point
            
            * indiv_max - divide the correlation trajectories by their individual max point
            
            * raw - do not divide by a point.

        Returns
        -------
        cross_corr : ndarray
            cross correlation tensor,  size (ncolor**2, 2*ntime-1, ntraj).
        err_crosscorr : ndarray
            cross correlation standard error estimation.
        inds : list
            indices describing which colors were correlated with which.

        '''

        i = 0
        k = 0
        inds = []

        if isinstance(intensity_vecs,np.ndarray):
            ncolors = intensity_vecs.shape[0]
            time_pts = intensity_vecs.shape[1]
            traj = intensity_vecs.shape[2]
            cross_corr = np.zeros((ncolors**2, time_pts*2-1, traj))
            
        # check for jagged arrays
        if isinstance(intensity_vecs, list):
            if not self.__check_jagged_array(intensity_vecs):
               intensity_vec = np.array(intensity_vecs)
            else:
                #jagged array detected, calculate the acov and return with zero padding
                cross_corr, err_crosscorr, inds = self.__calculate_ccov_jagged_array(intensity_vecs, norm=norm, scale_fix=scale_fix)
                return cross_corr, err_crosscorr, inds
            

        if scale_fix:
            scale_array = np.hstack( [np.linspace(1, time_pts, time_pts)[1:] / time_pts, [1], np.linspace(1, time_pts, time_pts)[::-1][1:] / time_pts] )
        else:
            scale_array = 1
            
        for n in range(intensity_vecs.shape[0]):
            for m in range(intensity_vecs.shape[0]):
                iv1 = intensity_vecs[n].T
                iv2 = intensity_vecs[m].T
                inds.append((n, m))

                if norm in ['global', 'Global', 'g', 'G']:
                    global_mean1 = np.mean(iv1)
                    global_mean2 = np.mean(iv2)

                # slen = np.correlate(iv1[i]-np.mean(iv1[i]),iv2[i]-np.mean(iv2[i]),'full').shape[0]
                # crosscorr_vec = np.zeros((iv1.shape[0],slen))

                for i in range(traj):
                    if norm in ['Individual', 'I', 'individual', 'ind', 'indiv']:
                        cross_corr[k, :, i] = np.correlate(
                            iv1[i, :]-np.mean(iv1[i, :]),
                            iv2[i, :]-np.mean(iv2[i, :]), 'full')/time_pts/scale_array

                    elif norm in ['global', 'Global', 'g', 'G']:
                        cross_corr[k, :, i] = np.correlate(
                            iv1[i, :]-global_mean1, iv2[i, :]-global_mean2, 'full')/time_pts/scale_array

                    elif norm in ['raw', 'Raw', None, 'none', 'None']:
                        cross_corr[k, :, i] = np.correlate(
                            iv1[i, :], iv2[i, :], 'full')/time_pts/scale_array

                    else:
                        msg = 'unrecognized normalization, please use '\
                              'individual, global, or none for norm arguement'
                        raise custom_err.UnrecognizedNormalizationError(msg)
                        

                k += 1
        if g0 != None:
            cross_corr = self.normalize_cc(cross_corr, mode=g0)

        err_crosscorr = 1.0/np.sqrt(traj)*np.std(cross_corr, ddof=1, axis=2)

        return cross_corr, err_crosscorr, inds

    # def get_crosscorr2(self, iv1,iv2):
    #     '''
    #     returns the autocorrelations
    #     '''

    #     i = 0
    #     slen = np.correlate(iv1[i]-np.mean(iv1[i]),iv2[i]-np.mean(iv2[i]),'full').shape[0]
    #     crosscorr_vec = np.zeros((iv1.shape[0],slen))

    #     for i in range(iv1.shape[0]):
    #         crosscorr_vec[i,:] = np.correlate(iv1[i]-np.mean(iv1[i]),iv2[i]-np.mean(iv2[i]),'full')/len(iv1)

    #     normalized_autocorr = crosscorr_vec.T/ crosscorr_vec[:,len(iv1[i])-1]
    #     mean_autocorr = np.mean(normalized_autocorr, axis=1)



    #     return crosscorr_vec, mean_autocorr


    '''
    def get_autocorr_norm(self, intensity_vec, time_vec,
                          geneLength, normalization='Individual'):


        autocorr_vec = np.zeros((intensity_vec.shape))

        if normalization in ['Individual', 'I', 'individual', 'ind']:
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i, :] = self.get_acc2(intensity_vec[i] - np.mean(intensity_vec[i]))
        elif normalization in ['global', 'Global', 'g', 'G']:
            global_mean = np.mean(intensity_vec)
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i, :] = self.get_acc2(
                    intensity_vec[i] - global_mean)

        else:
            print('unrecognized normalization, using indivdual means')
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i, :] = self.get_acc2(
                    intensity_vec[i] - np.mean(intensity_vec[i]))

        normalized_autocorr = autocorr_vec.T / autocorr_vec[:, 0]
        mean_autocorr = np.mean(normalized_autocorr, axis=1)

        error_autocorr = np.std(normalized_autocorr, axis=1)/np.sqrt(intensity_vec.shape[0])

        dwelltime = None

        try:
            dwelltime = time_vec[np.where(mean_autocorr < .01)[0][0]]

        except:
            try:
                dwelltime = time_vec[np.where(mean_autocorr < .05)[0][0]]

            except:
                dwelltime = 1


        try:
            zeroind = np.where(mean_autocorr < 0)[0][0]
            length = int(.3*len(mean_autocorr))
            zeromean = np.mean(mean_autocorr[-length:])
            zeromean2 = np.mean(mean_autocorr[zeroind:])

            normalized_autocorr = normalized_autocorr - zeromean2
            mean_autocorr = np.mean(normalized_autocorr, axis=1)

            error_autocorr = np.std(normalized_autocorr, axis=1)/np.sqrt(intensity_vec.shape[0])

        except:
            pass



        ke_exp = np.round(geneLength/dwelltime, 1)

        return normalized_autocorr, mean_autocorr, error_autocorr, dwelltime, ke_exp
    '''


    def get_autocov(self, intensity_vec, norm='global', scale_fix=False):
        '''
        Return the autocovariance of an intensity vector as defined by:
        
        .. math:: 
            
            ACOV(t, \tau) = cov(X_{t}, X_{\tau}) = E{X_{t}, X_{\tau}} 
        
        There are also several normalization options:
            
        Raw - perform no normalization to each intensity trajectory
        
        Global - subtract the global intensity mean and divide by global
        variance
        
        .. math:: 
            
            \mu_I = E(Intensity) \\
            \sigma_I^2 = V(Intensity) \\
            X_{normalized} = (X(t) - \mu_I) / \sigma_I^2            

        Individual - subtact each trajectory by its mean and divide by
        its variance

        .. math:: 
            
            X_{normalized} = (X_{i}(t) - E(X_{i}(t))) / V(X_{i}(t))   

        Parameters
        ----------
        intensity_vec : ndarray
            intensity tensor of shape (ncolor, ntimes, ntraj).
        norm : str, optional
            normalization to use, 'raw' for no normalization, 'global' to 
            normalize by the gobal intensity moments, 
            'individual' to normalize each trajectory by
            its individual moments. The default is 'global'.
        scale_fix : bool, optional
            rescale the autocovariance function to ignore the zero padding, 
            that is divide all lags by the percentage of that total length. Default is 'False'
            Example: G(t=N) / ((signal_length - N) / signal_length )

        Returns
        -------
        autocorr_vec : ndarray
            returns autorcovariance array of size (ncolor, ntime-1, ntraj).
        autocorr_err : ndarray
            returns autorcovariance SEM array of size (ncolor, ntime-1, ntraj).

        '''
        
        if isinstance(intensity_vec,np.ndarray):
            autocorr_vec = np.zeros((intensity_vec.shape))
            autocorr_err = np.zeros((intensity_vec.shape))
            colors = intensity_vec.shape[0]
            n_traj = intensity_vec.shape[2]
            t = intensity_vec.shape[1]
            
        # check for jagged arrays
        if isinstance(intensity_vec, list):
            if not self.__check_jagged_array(intensity_vec):
               intensity_vec = np.array(intensity_vec)
            else:
                #jagged array detected, calculate the acov and return with zero padding
                autocorr_vec, autocorr_err = self.__calculate_acov_jagged_array(intensity_vec, norm=norm, scale_fix=scale_fix)
                return autocorr_vec, autocorr_err 
            

        if scale_fix:
            scale_array = np.linspace(1, t, t)[::-1] / t
        else:
            scale_array = 1
        
        for n in range(colors):
            if norm in ['Individual', 'I', 'individual', 'ind','i']:
                for i in range(intensity_vec.shape[2]):
                    ivec = intensity_vec[n, :, i]
                    autocorr_vec[n, :, i] = self.get_acc2(
                        (ivec - np.mean(ivec))/np.var(ivec))/scale_array

            elif norm in ['global', 'Global', 'g', 'G']:
                global_mean = np.mean(intensity_vec[n])
                global_var = np.var(intensity_vec[n])
                for i in range(intensity_vec.shape[2]):
                    autocorr_vec[n, :, i] = self.get_acc2(
                        (intensity_vec[n, :, i]-global_mean)/global_var )/scale_array
            elif norm in ['raw', 'Raw']:
                for i in range(intensity_vec.shape[2]):
                    autocorr_vec[n, :, i] = self.get_acc2(
                        intensity_vec[n, :, i]) /scale_array
            else:
                print('unrecognized normalization,'/
                      ' please use individual, global, or none')
                return

        autocorr_err = 1.0/np.sqrt(n_traj)*np.std(
            autocorr_vec, ddof=1, axis=2)

        return autocorr_vec, autocorr_err


    def get_autocorr(self, autocov, norm_type='interp', norm = 'individual'):
        '''
        Given an autocovariance tensor, normalize to the autocorrelation
        
        .. math:: 
            
            ACORR(X(t)) = ACOV(X(t)) / Normalization Constant
            
        where Normalization constant is defined as the delay to divide all 
        autocrrelations by:
            
            * G0 - the first delay without shot noise correction
            
            * G1 - the second delay of the autocorrelation
            
            * interp - interpolated G0, take G1-4 and calculate the G0 without shot noise
            
            
        norm = global will normalize the autocovariance by the global average normalization_constant
        norm = individual will normalize each trajectory by its own normalization_constant
            
        Parameters
        ----------
        autocov : ndarray
            autocovariance tensor of shape (Ncolor, Ntimes, Ntrajectories).
        norm_type : str, optional
            Delay to normalize by, G0, G1 or interp for interpolated G0. The default is 'interp'.
        norm : str, optional
            globally normalize autocovariance or individually normalize.
            Normalize each trajectory by its own g0 or by the global g0. 
            The default is 'individual'

        Returns
        -------
        autocorr : ndarray
            autocorrelation tensor of shape (Ncolor, Ntimes, Ntrajectories).
        err_autocorr : ndarray
            SEM autocorrelation tensor of shape (Ncolor, Ntimes, Ntrajectories).

        '''
        autocorr = np.copy(autocov)
        n_traj = autocorr.shape[-1]

        if norm_type.lower() in ['individual','indiv','i']:
            g0 = self.get_g0(autocov, norm)
            for n in range(autocov.shape[0]):
                autocorr[n] = autocorr[n]/g0[n]
        elif norm_type.lower() in ['global','g']:
            g0 = self.get_g0(autocov, norm)
            g0_mean = np.mean(g0)
            for n in range(autocov.shape[0]):
                autocorr[n] = autocorr[n]/g0_mean     
                
        else: 
        
            msg = 'unrecognized normalization, please use '\
                  'individual, or global for norm arguement'
            raise custom_err.UnrecognizedNormalizationError(msg)                
        err_autocorr =  1.0/np.sqrt(n_traj)*np.std(autocorr, ddof=1, axis=2)
        return autocorr, err_autocorr

    @staticmethod
    def standardize_signal(signal, axis=0, norm='global'):
        '''
        Perform standardization to set signal mean to 0 with unit variance of 
        1. This can be applied globally or individually to each trajectory
        with the flag norm ='global' or norm='indiv'

        Parameters
        ----------
        signal : ndarray
            intensity or signal array.
        axis : int, optional
            axis to apply the standardization over. The default is 0.
        norm : str, optional
            apply this using global mean ('global') and var or individual trajectory
            mean and var ('indiv'). The default is 'global'.

        Raises
        ------
        UnrecognizedNormalizationError
            Throws an error if the norm is not global or individual.

        Returns
        -------
        ndarray
            standardized array over axis desired.

        '''
        if norm in ['global', 'Global', 'g', 'G']:
            return (signal - np.mean(signal)) / np.std(signal)
        if norm in ['Individual', 'I', 'individual', 'ind','indiv','i']:
            return (signal - np.mean(signal, axis=axis)) / np.std(signal,axis=axis)
        else:
            msg = 'unrecognized normalization, please use '\
                  'individual, global, or none for norm arguement'
            raise custom_err.UnrecognizedNormalizationError(msg)    
            
    #sets individual 0-1
    @staticmethod
    def minmax_signal(signal, axis=0, norm='global'):
        '''
        normalize a singal by its min and max, leaving trajectories from 0 to 1
        This can be applied globally or individually to each trajectory
        with the flag norm ='global' or norm='indiv'


        .. code-block::
            
            #global normalization
            S_95= np.quantile(0.95)
            S_normalized = (S - np.min(S)) / (np.max(S) - np.min(S))
            
            #individual normalization
            S_95= np.quantile(0.95)
            S_normalized = (S - np.min(S,axis=axis)) / (np.max(S,axis=axis) - np.min(S,axis=axis))
        

        Parameters
        ----------
        signal : ndarray
            intensity or signal array.
        axis : int, optional
            axis to apply the standardization over. The default is 0.
        norm : str, optional
            apply this using global mean ('global') and var or individual trajectory
            mean and var ('indiv'). The default is 'global'.

        Raises
        ------
        UnrecognizedNormalizationError
            Throws an error if the norm is not global or individual.

        Returns
        -------
        ndarray
            minmaxed array over axis desired.

        '''
        if norm in ['global', 'Global', 'g', 'G']:
            return (signal -  np.min(signal) ) / ( np.max(signal)  - np.min(signal) )
        if norm in ['Individual', 'I', 'individual', 'ind','indiv','i']:
            return (signal -  np.min(signal,axis=axis) ) / ( np.max(signal,axis=axis)  - np.min(signal,axis=axis) )
        else:
            msg = 'unrecognized normalization, please use '\
                  'individual, global, or none for norm arguement'
            raise custom_err.UnrecognizedNormalizationError(msg)    
            
            
    def __check_jagged_array(self,list_of_lists):
        _,t_sizes,_ = self.__check_jagged_array_sizes(list_of_lists)
        if len(set(t_sizes)) ==1:
            return False
        else:
            return True

    def __check_jagged_array_sizes(self,list_of_lists):
        color_size = len(list_of_lists)
        t_sizes = []
        n_traj = len(list_of_lists[0])
        for i in range(color_size):
            for j in range(n_traj):
                t_sizes.append(len(list_of_lists[i][j]))
        return color_size, t_sizes, n_traj
    
    
    def __calculate_ccov_jagged_array(self, intensity_vecs, norm='global', g0='indiv_center', scale_fix=False):
        
        colors, t_sizes, n_traj = self.__check_jagged_array_sizes(intensity_vecs)
        
        cross_corr = np.zeros([colors**2, np.max(t_sizes)*2-1, n_traj])
        center = int(np.ceil((np.max(t_sizes)*2-1 )/2) )
        total_size = np.max(t_sizes)*2-1 
        scale_array = 1
        k = 0
        inds = []
        
        
        for n in range(colors):
            for m in range(colors):
                inds.append((n, m))
                
                

                if norm in ['global', 'Global', 'g', 'G']:
                    global_mean1 = np.mean([np.mean(x) for x in  intensity_vecs[n]])
                    global_mean2 = np.mean([np.mean(x) for x in  intensity_vecs[m]])

                # slen = np.correlate(iv1[i]-np.mean(iv1[i]),iv2[i]-np.mean(iv2[i]),'full').shape[0]
                # crosscorr_vec = np.zeros((iv1.shape[0],slen))

                for i in range(n_traj):
                    
                    iv1 = intensity_vecs[n][i]
                    iv2 = intensity_vecs[m][i]
                
                   
                    time_pts = min(len(iv1), len(iv2)) #sweep the smaller sized array
                    if scale_fix:
                        scale_array = np.hstack( [np.linspace(1, time_pts, time_pts)[1:] / time_pts, [1], np.linspace(1, time_pts, time_pts)[::-1][1:] / time_pts] )
                    else:
                        scale_array = 1
                    
                    if norm in ['Individual', 'I', 'individual', 'ind', 'indiv']:
                        corr = np.correlate(
                            iv1-np.mean(iv1),
                            iv2-np.mean(iv2), 'full')/time_pts/scale_array
                        corr_len = len(corr)
                        start = int((total_size -  corr_len)/2)
                        stop = total_size-start
                        cross_corr[k,start:stop  , i] = corr

                    elif norm in ['global', 'Global', 'g', 'G']:
                        cross_corr[k, center-time_pts+1:center+time_pts, i] = np.correlate(
                            iv1-global_mean1, iv2-global_mean2, 'full')/time_pts/scale_array

                    elif norm in ['raw', 'Raw', None, 'none', 'None']:
                        cross_corr[k, center-time_pts+1:center+time_pts, i] = np.correlate(
                            iv1,iv2, 'full')/time_pts/scale_array

                    else:
                        msg = 'unrecognized normalization, please use '\
                              'individual, global, or none for norm arguement'
                        raise custom_err.UnrecognizedNormalizationError(msg)
                        

                k += 1
        if g0 != None:
            cross_corr = self.normalize_cc(cross_corr, mode=g0)

        err_crosscorr = 1.0/np.sqrt(n_traj)*np.std(cross_corr, ddof=1, axis=2)

        return cross_corr, err_crosscorr, inds        
        

    
    def __calculate_acov_jagged_array(self,list_of_lists, norm='global', scale_fix=False,):
        
        colors, t_sizes, n_traj = self.__check_jagged_array_sizes(list_of_lists)
        
        acov_array = np.zeros([colors, np.max(t_sizes), n_traj])
        scale_array = 1
        
        
        for n in range(colors):
            if norm in ['Individual', 'I', 'individual', 'ind','i']:
                for i in range(n_traj):
                    if scale_fix:
                        scale_array = np.linspace(1, t_sizes[n*n_traj + i], t_sizes[n*n_traj + i])[::-1] / t_sizes[n*n_traj + i]
                    ivec = np.array(list_of_lists[n][i]).flatten()
                    acov_array[n, :t_sizes[n*n_traj + i], i] = self.get_acc2(
                        (ivec - np.mean(ivec))/np.var(ivec))/scale_array

            elif norm in ['global', 'Global', 'g', 'G']:
                
                global_mean = np.mean([np.mean(x) for x in list_of_lists[n]])
                global_var = np.mean([np.var(x) for x in list_of_lists[n]])
                
                for i in range(n_traj):
                    if scale_fix:
                        scale_array = np.linspace(1, t_sizes[n*n_traj + i], t_sizes[n*n_traj + i])[::-1] / t_sizes[n*n_traj + i]
                                        
                    ivec = np.array(list_of_lists[n][i]).flatten()
                    acov_array[n, :t_sizes[n*n_traj + i], i] = self.get_acc2(
                        (ivec-global_mean)/global_var )/scale_array
                    
            elif norm in ['raw', 'Raw']:
                for i in range(n_traj):
                    if scale_fix:
                        scale_array = np.linspace(1, t_sizes[n*n_traj + i], t_sizes[n*n_traj + i])[::-1] / t_sizes[n*n_traj + i]
                    
                    ivec = np.array(list_of_lists[n][i]).flatten()
                    acov_array[n, :t_sizes[n*n_traj + i], i] = self.get_acc2(
                        ivec) /scale_array
            else:
                print('unrecognized normalization,'/
                      ' please use individual, global, or none')
                return

        acov_error = 1.0/np.sqrt(n_traj)*np.std(
            acov_array, ddof=1, axis=2)

        return acov_array, acov_error
        


    @staticmethod
    def minmax_quantile_signal(signal, axis=0, norm='global', quantile=.95, max_outlier=1.5):
        '''
        normalize a singal by its min and max from a quantile,
        leaving trajectories from 0 to 1 at that quantile. Set all outliers
        over a maximum to max_outlier.
        This can be applied globally or individually to each trajectory
        with the flag norm ='global' or norm='indiv'
        
        
        .. code-block::
            
            #global normalization
            S_95= np.quantile(0.95)
            S_normalized = (S - np.min(S)) / (np.quantile(0.95) - np.min(S))
            
            #individual normalization
            S_95= np.quantile(0.95)
            S_normalized = (S - np.min(S,axis=axis)) / (np.quantile(0.95) - np.min(S,axis=axis))
        

        Parameters
        ----------
        signal : ndarray
            intensity or signal array.
        axis : int, optional
            axis to apply the standardization over. The default is 0.
        norm : str, optional
            apply this using global mean ('global') and var or individual trajectory
            mean and var ('indiv'). The default is 'global'.
        quantile : float (0-1)
            the % quantile to set to the new 1
        max_outlier : float
            the maximum value to use for any outlier after normalization

        Raises
        ------
        UnrecognizedNormalizationError
            Throws an error if the norm is not global or individual.

        Returns
        -------
        ndarray
            minmaxed per quantile array over axis desired.

        '''
        if norm in ['global', 'Global', 'g', 'G']:
            
            max_95 = np.quantile(signal, quantile)
            sig = (signal -  np.min(signal) ) / ( max_95  - np.min(signal) )
            return np.minimum(sig,max_outlier)
        
        if norm in ['Individual', 'I', 'individual', 'ind','indiv','i']:
            
            max_95 = np.quantile(signal, quantile,axis=axis)
            sig = (signal -  np.min(signal,axis=axis) ) / ( max_95 - np.min(signal,axis=axis) )
            return np.minimum(sig,max_outlier)
        
        else:
            msg = 'unrecognized normalization, please use '\
                  'individual, global, or none for norm arguement'
            raise custom_err.UnrecognizedNormalizationError(msg)    
            
            
    '''

    def get_autocorr2(self, intensity_vec, time_vec, totalSimulationTime,
                      geneLength, normalization='Individual'):


        autocorr_vec = np.zeros((intensity_vec.shape))


        if normalization in ['Individual', 'I', 'individual', 'ind']:
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i, :] = self.get_acc2(
                    intensity_vec[i]-np.mean(intensity_vec[i]))
        elif normalization in ['global', 'Global', 'g', 'G']:
            global_mean = np.mean(intensity_vec)
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i, :] = self.get_acc2(
                    intensity_vec[i]-global_mean)

        else:
            print('unrecognized normalization, using indivdual means')
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i, :] = self.get_acc2(
                    intensity_vec[i]-np.mean(intensity_vec[i]))
        autocorr = autocorr_vec.T
        mean_autocorr = np.mean(autocorr, axis=1)
        error_autocorr = np.std(
            autocorr, axis=1)/np.sqrt(intensity_vec.shape[0])


        dwelltime = None

        try:
            dwelltime = time_vec[np.where(mean_autocorr < .01)[0][0]]

        except:
            try:
                dwelltime = time_vec[np.where(mean_autocorr < .05)[0][0]]
            except:
                dwelltime = 1

        ke_exp = np.round(geneLength/dwelltime, 1)

        return autocorr, mean_autocorr, error_autocorr, dwelltime, ke_exp
    '''