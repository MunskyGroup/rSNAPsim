# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:50:56 2020

@author: William Raymond
"""
import numpy as np
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
            Remove non zero entries and only return decreasing entries. The default is False.

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
            acf[acf < 0]=0
            for i in range(1, len(acf)):
                if acf[i] > acf[i-1]:
                    acf[i] = acf[i-1]
        return acf    
        
    
    def get_fragments(self,position_tensor, total_length = None):
        '''
        Get individual ribosome trajectories for kymograph generation
        
        .. warning:: This is not perfect process, if the time resolution of a simulation is too coarse this will not be able to rebuild ribosome trajectories. This will fail if ribosomes are moving fast enough between timesteps to pass previously recorded positions.

        Parameters
        ----------
        position_tensor : ndarray
            Ribosome position tensor:  Nribosomes x Ntime, accessed from SSA_Soln.solutions if a full statistical model was run.
        total_length : int, optional
            Length of the transcript. If left blank it will use the maximum detected ribosome position

        Returns
        -------
        fragtimes : ndarray
            N_ribosomes x N_longest_time. The individual trajectory times for the positions fragments
        fragarray : ndarray
            N_ribosomes x N_position. The individual ribosome trajectory positions 

        '''
        if total_length == None:
            total_length = np.max(position_tensor)
            
        fragmented_trajectories = []
        fragtimes = []
        maxlen = 0
    
        fragmentspertraj= []
        solutions = [position_tensor]
        nsteps = position_tensor.shape[0]
        k = 0


        ind = np.zeros(position_tensor.shape[1]).astype(int)
        for i in range(position_tensor.shape[1]):
          tmp = np.where(position_tensor[:,i] == 0)[0]
          if len(tmp) > 0:
            ind[i] = np.where(position_tensor[:,i] == 0)[0][0]
          else:
            ind[i] = position_tensor.shape[0]-1
        
        changes = ind[1:] - ind[:-1]
        addindexes = np.where(changes > 0)[0]
        subindexes = np.where(changes < 0)[0]
        
        sub = solutions[k][:,1:] - solutions[k][:,:-1]
        neutralindexes = np.unique(np.where(sub < 0)[1])
        neutralindexes = np.setxor1d(neutralindexes, subindexes)
        
        for index in neutralindexes:
            pre = solutions[k][:,index]
            post = solutions[k][:,index+1]
            changecount = 0
            while len(np.where(post - pre < 0)[0]) > 0:

                post = np.append([total_length],post)
                pre = np.append(pre,0)
                
                changecount+=1
            
            for i in range(changecount):
                addindexes = np.sort(np.append(addindexes,index))
                subindexes = np.sort(np.append(subindexes,index))
                
            changes[index] = -changecount
            ind[index] += changecount
         
            
        for index in np.where(np.abs(changes)>1)[0]:
            if changes[index] < 0:
                for i in range(np.abs(changes[index])-1):
                    subindexes = np.sort(np.append(subindexes,index))
            else:
                for i in range(np.abs(changes[index])-1):
                    addindexes = np.sort(np.append(addindexes,index))   
            
        truefrags = len(subindexes)
 
            
    
       
        if len(subindexes) < len(addindexes):
            subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(nsteps-1)).astype(int))
            
        
        fragmentspertraj.append(len(subindexes))
        
        for m in range(min(len(subindexes),len(addindexes))):
            traj = solutions[k][:, addindexes[m]:subindexes[m]+1]
            traj_ind = changes[addindexes[m]:subindexes[m]+1]
            
            startind = ind[addindexes[m]]
            minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
            fragment = np.array([])
        
                
            
            iterind = startind
            
            if subindexes[m]-addindexes[m] > 0:
                if len(minusloc) > 1:
                    if m <= truefrags:
                        for n in range(len(minusloc)-1):
                            iterind = iterind + int(min(0,traj_ind[minusloc[n]]))

                            fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                            
                            
                            
              
            
                  
                        
                        fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                        
                    else:
                        for n in range(len(minusloc)-1):

                            iterind = iterind + min(0,traj_ind[minusloc[n]])
                            
                            fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
              
                            
                        fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
      
                    
                
                else:

                    fragment = solutions[k][startind][addindexes[m]:subindexes[m]+1].flatten()
               
            
                
                fragtimes.append(addindexes[m]+1)
                   
                
                fragmented_trajectories.append(fragment)
                #if m <= truefrags:
                    #kes.append(genelength/truetime[len(fragment)])
        
                if len(fragment) > maxlen:
                    maxlen = len(fragment)
                

        fragarray = np.zeros((len(fragmented_trajectories), maxlen))
        for i in range(len(fragmented_trajectories)):
            fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
        
        return fragtimes,fragarray
            
            
        
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
        
        
    def get_g0(self,correlation, mode = 'interp'):
        '''
        

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
        if mode.lower() in ['interp','inter','extrapolate','interpolate']:
            X = [1,2,3,4]
            V = correlation[:,X,:]
            G0 = np.interp(0,X,V)      
            
        if mode.lower() in ['g1','1']:
            G0 = correlation[:,1,:]
            
        if mode.lower() in ['g0','0']:
            G0 = correlation[:,0,:]
            
        if mode.lower() in ['max','maximum']:
            G0 = np.max(correlation,axis=1)
        return G0
    
    def normalize_cc(self, correlation,mode='global_max'):
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
        if mode.lower() in ['global_max','global_maximum','gmax']:
            norm_cor = correlation/np.max(correlation)
        
        if mode.lower() in ['indiv_max','indiv_maximum','imax','individual_maximum','individual_max']:
            norm_cor = correlation/np.max(correlation,1)[:,np.newaxis,:]
            
        if mode.lower() in ['indiv_center','indiv_middle','individual_middle','individual_center']:
            centerpoint = int((correlation.shape[1]+1)/2)-1
            norm_cor = correlation/( correlation[:,centerpoint][:,np.newaxis,:] )
            
        if mode.lower() in ['global_center','global_middle']:
            centerpoint = int((correlation.shape[1]+1)/2)-1
            norm_cor = correlation/(np.mean(correlation[:,centerpoint]))   
        if mode.lower() in ['raw']:
            norm_cor = correlation
            
        return norm_cor 
        
    
    def get_crosscorr(self,intensity_vecs, norm='indiv', g0='indiv_center'):
        '''
        return a cross correlation tensor from an intensity tensor of (ncolor, ntime, ntraj)

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
        err_crosscorr : TYPE
            cross correlation standard error estimation.
            
            .. math:: \sigma_{cc} = \frac{\sigma}{\sqrt(n_{trajectories})}
            
        inds : list
            indices describing which colors were correlated with which.

        '''
        
        ncolors = intensity_vecs.shape[0]
        time_pts = intensity_vecs.shape[1]
        traj = intensity_vecs.shape[2]
        cross_corr = np.zeros( (ncolors**2,time_pts*2-1,traj ))
        
        i = 0
        k = 0
        inds = []

        
        for n in range(intensity_vecs.shape[0]):
            for m in range(intensity_vecs.shape[0]):
                iv1 = intensity_vecs[n].T
                iv2 = intensity_vecs[m].T
                inds.append((n,m))
                
                if norm in ['global','Global','g','G']:
                    global_mean1 = np.mean(iv1)
                    global_mean2 = np.mean(iv2)
                
                # slen = np.correlate(iv1[i]-np.mean(iv1[i]),iv2[i]-np.mean(iv2[i]),'full').shape[0]
                # crosscorr_vec = np.zeros((iv1.shape[0],slen))
                
                for i in range(traj):
                    if norm in [ 'Individual','I','individual','ind','indiv']:       
                        cross_corr[k,:,i] = np.correlate(iv1[i,:]-np.mean(iv1[i,:]),iv2[i,:]-np.mean(iv2[i,:]),'full')/time_pts
                    
                    elif norm in ['global','Global','g','G']:
                        cross_corr[k,:,i] = np.correlate(iv1[i,:]-global_mean1,iv2[i,:]-global_mean2,'full')/time_pts

                    elif norm in ['raw','Raw', None,'none','None']:
                        cross_corr[k,:,i] = np.correlate(iv1[i,:],iv2[i,:],'full')/time_pts

                    else:
                        raise(ValueError, "unrecognized normalization, please use individual, global, or none for norm arguement")
                    
                        return
                
                k +=1
        if g0 != None:
            cross_corr = self.normalize_cc(cross_corr,mode=g0)
                              
        err_crosscorr =  1.0/np.sqrt(traj)*np.std(cross_corr,ddof=1,axis=2)
        
        return cross_corr  , err_crosscorr,inds

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
        

    def get_autocorr_norm(self, intensity_vec, time_vec, totalSimulationTime, geneLength,normalization= 'Individual'):
        '''
        returns the autocorrelations
        '''

        autocorr_vec = np.zeros((intensity_vec.shape))


        if normalization in [ 'Individual','I','individual','ind']:
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-np.mean(intensity_vec[i]))
        elif normalization in ['global','Global','g','G']:
            global_mean = np.mean(intensity_vec)
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-global_mean)
            
        else:   
            print('unrecognized normalization, using indivdual means')
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-np.mean(intensity_vec[i]))
                
                


        normalized_autocorr = autocorr_vec.T/ autocorr_vec[:,0]
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
            zeroind = np.where(mean_autocorr<0)[0][0]
            length = int(.3*len(mean_autocorr))
            zeromean = np.mean(mean_autocorr[-length:])
            zeromean2 = np.mean(mean_autocorr[zeroind:])
    
            normalized_autocorr = normalized_autocorr-zeromean2
            mean_autocorr = np.mean(normalized_autocorr, axis=1)
            
            error_autocorr = np.std(normalized_autocorr, axis=1)/np.sqrt(intensity_vec.shape[0])
     
        except:
            pass
        


        ke_exp = np.round(geneLength/dwelltime ,1)

        return normalized_autocorr, mean_autocorr, error_autocorr, dwelltime, ke_exp
    
    
    
    def get_autocov(self,intensity_vec,norm='raw'):
        autocorr_vec = np.zeros((intensity_vec.shape))
        autocorr_err = np.zeros((intensity_vec.shape))
        colors = intensity_vec.shape[0]
        n_traj = intensity_vec.shape[2]
        
        for n in range(colors):
            if norm in [ 'Individual','I','individual','ind']:
                for i in range(intensity_vec.shape[2]):
                    ivec = intensity_vec[n,:,i]
                    autocorr_vec[n,:,i] = self.get_acc2( (ivec-np.mean(ivec))/np.var(ivec)  )
                    
            elif norm in ['global','Global','g','G']:
                global_mean = np.mean(intensity_vec[n])
                global_var = np.var(intensity_vec[n])
                for i in range(intensity_vec.shape[2]):
                    autocorr_vec[n,:,i] = self.get_acc2((intensity_vec[n,:,i]-global_mean)/global_var )
            elif norm in ['raw','Raw']:
                for i in range(intensity_vec.shape[2]):
                    autocorr_vec[n,:,i] = self.get_acc2(intensity_vec[n,:,i])     
            else:
                print('unrecognized normalization, please use individual, global, or none')
                return

        autocorr_err =  1.0/np.sqrt(n_traj)*np.std(autocorr_vec,ddof=1,axis=2)
                
        return autocorr_vec, autocorr_err
    
    
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
                
                
                
    def get_autocorr2(self, intensity_vec, time_vec, totalSimulationTime, geneLength, normalization='Individual'):
        '''
        returns the autocorrelations
        '''

        autocorr_vec = np.zeros((intensity_vec.shape))
        
        
        if normalization in [ 'Individual','I','individual','ind']:
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-np.mean(intensity_vec[i]))
        elif normalization in ['global','Global','g','G']:
            global_mean = np.mean(intensity_vec)
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-global_mean)
            
        else:   
            print('unrecognized normalization, using indivdual means')
            for i in range(intensity_vec.shape[0]):
                autocorr_vec[i,:] = self.get_acc2(intensity_vec[i]-np.mean(intensity_vec[i]))
        autocorr = autocorr_vec.T
        mean_autocorr = np.mean(autocorr, axis=1)
        error_autocorr = np.std(autocorr, axis=1)/np.sqrt(intensity_vec.shape[0])
        
 
        dwelltime = None

        try:
            dwelltime = time_vec[np.where(mean_autocorr < .01)[0][0]]

        except:
            try:
                dwelltime = time_vec[np.where(mean_autocorr < .05)[0][0]]
            except:
                dwelltime = 1




        ke_exp = np.round(geneLength/dwelltime ,1)

        return autocorr, mean_autocorr, error_autocorr, dwelltime, ke_exp  


