# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 19:21:12 2021

@author: willi
"""

import warnings
import numpy as np
from . import SequenceManipMethods
smm = SequenceManipMethods.SequenceManipMethods
from scipy import stats

class RibosomeAnalyses():

    '''
    A class containing functions for analyzing ribosomal position tensors
    of the form (Ntraj, Nt, Recording_index)
    '''

    def __init__(self):
        self.aa_link = .35 #nm on average
        self.nt_link = .33 #nm on average
        pass
    
    
    
    
    def circular_probe_positions(self, position_tensor, probe_vec, utr_3p_length=100,
                                 utr_5p_length=0,
                                 utr_tag='MS2', utr_fluorophore='mCherry', n_utr_probes=24):
        
        construct_length = 3*probe_vec.shape[-1] + utr_3p_length #nucleotides
        
        base_radius = construct_length*self.nt_link/(2*np.pi) #nanometers
        center = 0
        def angle_to_ribosome_location(r,angle):
            x,y = r*np.cos(angle), r*np.sin(angle)
            return x,y
        
        def ribosome_location_to_angle(codon_location,construct_length):
            angle = codon_location*3 / construct_length * 2*np.pi
            return angle
        
        angle_tensor = ribosome_location_to_angle(position_tensor, construct_length)
        for i in range(np.atleast_2d(probe_vec).shape[0]):
            if i ==0:
                
                intensity_tensor = probe_vec[0][ np.maximum(position_tensor -1,0) ]
                intensity_tensor = np.expand_dims(intensity_tensor,0)
                shape_dim = intensity_tensor.shape
            else:
                intensity_tensor = np.concatenate((intensity_tensor,
                                                   np.expand_dims(probe_vec[i][np.maximum(position_tensor -1,0) ],0 )) )
        utr_angle = ribosome_location_to_angle(probe_vec.shape[-1],construct_length)
        utr_intensity = np.array((*angle_to_ribosome_location(base_radius,utr_angle  ), n_utr_probes))
    
    
        #intensity_tensor Ncolor, Ntraj, Nt, Nrib
        #angle_tensor Ntraj, Nt, Nrib
        
        x,y = angle_to_ribosome_location(base_radius, angle_tensor)
        return x,y,intensity_tensor, utr_intensity
        
        
        
    def linear_probe_positions(self, position_tensor, probe_vec,
                                 utr_tag='MS2', utr_fluorophore='mCherry',
                                 utr_3p_length=100,
                                 n_utr_probes=24):
        
        construct_length = 3*probe_vec.shape[-1] + utr_3p_length #nucleotides
        
        for i in range(np.atleast_2d(probe_vec).shape[0]):
            if i ==0:
                
                intensity_tensor = probe_vec[0][ np.maximum(position_tensor -1,0) ]
                intensity_tensor = np.expand_dims(intensity_tensor,0)
                shape_dim = intensity_tensor.shape
            else:
                intensity_tensor = np.concatenate((intensity_tensor,
                                                   np.expand_dims(probe_vec[i][np.maximum(position_tensor -1,0) ],0 )) )

        x = 3*position_tensor*self.nt_link
        utr_intensity = np.array( (3*probe_vec.shape[-1]*self.nt_link, n_utr_probes))
        return x, intensity_tensor, utr_intensity
    

    def spring_probe_positions(self, position_tensor, probe_vec, utr_3p_length=100,
                                 utr_5p_length=0,
                                 utr_tag='MS2', utr_fluorophore='mCherry', n_utr_probes=24,
                                 stretched_length=10, compressed_length=1):
        
        construct_length = 3*probe_vec.shape[-1] + utr_3p_length #nucleotides
        
        base_radius = construct_length*self.nt_link/(2*np.pi) #nanometers
        center = 0
        def angle_tensor_to_ribosome_tensor(r,angle):
            r = np.repeat(np.expand_dims(r,-1), angle.shape[-1],axis=-1)
            x,y = np.multiply(r,np.cos(angle)), np.multiply(r,np.sin(angle))
            return x,y
        
        def angle_to_ribosome_location(r,angle):
            x,y = r*np.cos(angle), r*np.sin(angle)
            return x,y
              
        def ribosome_location_to_angle(codon_location,construct_length):
            angle = codon_location*3 / construct_length * 2*np.pi
            return angle
        
        def ribtensor_to_angletensor(position_tensor, stretched_tensor, construct_length):
            return position_tensor*3 / (construct_length + np.expand_dims(stretched_tensor,-1)*stretched_length)*2*np.pi
        
        stretched_tensor = np.sum((position_tensor !=0),axis=-1)
        angle_tensor = ribtensor_to_angletensor(position_tensor, stretched_tensor, construct_length)
        
        for i in range(np.atleast_2d(probe_vec).shape[0]):
            if i ==0:
                
                intensity_tensor = probe_vec[0][ np.maximum(position_tensor -1,0) ]
                intensity_tensor = np.expand_dims(intensity_tensor,0)
                shape_dim = intensity_tensor.shape
            else:
                intensity_tensor = np.concatenate((intensity_tensor,
                                                   np.expand_dims(probe_vec[i][np.maximum(position_tensor -1,0) ],0 )) )
        utr_angle = ribosome_location_to_angle(probe_vec.shape[-1],construct_length)
        utr_intensity = np.array((*angle_to_ribosome_location(base_radius,utr_angle  ), n_utr_probes))
    
    
        #intensity_tensor Ncolor, Ntraj, Nt, Nrib
        #angle_tensor Ntraj, Nt, Nrib
        radius_tensor = (construct_length + stretched_tensor*stretched_length)*self.nt_link/(2*np.pi) #nanometers
        x,y = angle_tensor_to_ribosome_tensor(radius_tensor, angle_tensor)
        return x,y,intensity_tensor, utr_intensity
            
    def WLC_positions_2D(self, position_tensor, probe_vec, dt,
                                 utr_tag='MS2', utr_fluorophore='mCherry', 
                                 stretched_length=5, compressed_length=1,
                                 utr_3p_length = 100, stretched_codons=5):
        
        construct_length = 3*probe_vec.shape[-1] + utr_3p_length
        
        pL = self.nt_link* compressed_length
        pL_stretched = self.nt_link*stretched_length
        stretched_region = stretched_codons*3

        ptheta = lambda theta: np.sqrt(pL/(2*np.pi))*np.exp(-pL*theta**2 /2 )
        ptheta_stiff = lambda theta: np.sqrt(pL_stretched/(2*np.pi))*np.exp(-pL_stretched*theta**2 /2 )
        angles = np.linspace(-6.28, 6.28, 5000);
        cdf = ptheta(np.linspace(-6.28,6.28,5000))
        ccdf = np.cumsum(cdf)
        cdf_sum = ccdf/np.max(ccdf)

        cdf = ptheta_stiff(np.linspace(-6.28,6.28,5000))
        ccdf = np.cumsum(cdf)
        cdf_sum_stiff = ccdf/np.max(ccdf)
        get_angle = lambda: angles[np.where(np.random.rand() < cdf_sum)[0][0]]
        get_angle_stretched = lambda: angles[np.where(np.random.rand() < cdf_sum_stiff)[0][0]]
        
        #generate an inital chain for each trajectory
        def init_chains(position_tensor, L):
            x,y = np.zeros([position_tensor.shape[0], L]), np.zeros([position_tensor.shape[0], L])
            angles = np.zeros([position_tensor.shape[0], L-1])
            for j in range(position_tensor.shape[0]): 
                xs,ys = 0,0
                occupied = position_tensor[j,0,0][position_tensor[j,0,0]!=0]*3
                
                x[j,0], y[j,0] = xs, ys
                for i in range(1,L):
                    if np.any(np.abs(i - occupied) < stretched_region ):
                        angle = get_angle_stretched()
                    else:
                        angle = get_angle()
                    xs += pL*np.cos(angle)
                    ys += pL*np.sin(angle)
                    x[j,i], y[j,i] = xs, ys
                    angles[j,i-1] = angle
                    
                
            return x, y,angles

        def get_angle_2pi_rad(pt1,pt2):
            a  = np.arctan2( pt1[1] - pt2[1], pt1[0] - pt2[0])
            if a < 0:
                a = a+ np.pi*2
            return a
    
        def update_chain( chainx, chainy, dt, nt, L, achain):
            
            for j in range(position_tensor.shape[0]): 
                xs,ys = chainx[j][0],chainy[j][0]
                occupied = position_tensor[j,nt,:][position_tensor[j,nt,:]!=0]*3
                
                in_occupied_region = False
                for i in range(1,L):
                    
                    p_ptx, p_pty = chainx[j][i-1], chainy[j][i-1]
                    old_ptx, old_pty = chainx[j][i], chainy[j][i]
    
                    #current_angle = np.arctan((old_pty - p_pty )/(old_ptx - p_ptx ))
                    current_angle = achain[j][i-1]
                    
    
                    if np.any(np.abs(i - occupied) < stretched_region ):
                    #if 0:
                        
                        if not in_occupied_region:
                            # inside first link of occupied region
                            dangle =  get_angle_stretched()*dt
                            angle = dangle+  current_angle
                            occupied_angle  = angle
                            new_ptx, new_pty = p_ptx + pL*np.cos(angle) , p_pty + pL*np.sin(angle)
                            in_occupied_region = True            
                        else:
                            # inside 2+ link of occupied region
                            dangle =  get_angle_stretched()*dt
                            angle =   occupied_angle + dangle
                            new_ptx, new_pty = p_ptx + pL*np.cos(angle) , p_pty + pL*np.sin(angle)
                            in_occupied_region = True

                    else:
                        dangle = get_angle()*dt
                        angle = current_angle + dangle
                        in_occupied_region = False
                    
                        new_ptx, new_pty = p_ptx + pL*np.cos(angle) , p_pty + pL*np.sin(angle)
                    
                    dx, dy = new_ptx - old_ptx, new_pty - old_pty
    
                    chainx[j,i:], chainy[j,i:] = chainx[j,i:] + dx,  chainy[j,i:] + dy
                    achain[j,i-1] = angle
                    achain[j,i:] += dangle
                    
            return chainx, chainy, achain        
        
        xchain,ychain = np.zeros([position_tensor.shape[0],position_tensor.shape[1],construct_length ]), np.zeros([position_tensor.shape[0],position_tensor.shape[1], construct_length])
        
        
        xi,yi,ai  = init_chains(position_tensor, construct_length)
        xchain[:, 0, :] = xi
        ychain[:, 0, :] = yi
        
        for i in range(1, position_tensor.shape[1]):
            xi, yi, ai = update_chain(xi, yi, dt,i, construct_length,ai)
            xchain[:, i, :] = xi
            ychain[:, i, :] = yi        
        
        return xchain,ychain


            
    def gaussian_jump_angles(self, position_tensor, mu_angle, var_angle):
        x=1


    def get_fragments(self, position_tensor, total_length=None):
        #TODO
        '''
        ..TODO refactor this code!!
        
        Get individual ribosome trajectories for kymograph generation

        .. warning:: This is not perfect process, if the time resolution of a
        simulation is too coarse this will not be able to rebuild ribosome
        trajectories. This will fail if ribosomes are moving fast enough 
        between timesteps to pass previously recorded positions.

        Parameters
        ----------
        position_tensor : ndarray
            Ribosome position tensor:  Nribosomes x Ntime, accessed from
            SSA_Soln.solutions if a full statistical model was run.
        total_length : int, optional
            Length of the transcript. If left blank it will use the maximum
            detected ribosome position

        Returns
        -------
        fragtimes : ndarray
            N_ribosomes x N_longest_time. The individual trajectory times
            for the positions fragments
        fragarray : ndarray
            N_ribosomes x N_position. The individual ribosome trajectory
            positions

        '''
        if total_length == None:
            total_length = np.max(position_tensor)

        fragmented_trajectories = []
        fragtimes = []
        maxlen = 0

        fragmentspertraj = []
        solutions = [position_tensor]
        nsteps = position_tensor.shape[0]
        k = 0


        ind = np.zeros(position_tensor.shape[1]).astype(int)
        for i in range(position_tensor.shape[1]):
          tmp = np.where(position_tensor[:, i] == 0)[0]
          if len(tmp) > 0:
            ind[i] = np.where(position_tensor[:, i] == 0)[0][0]
          else:
            ind[i] = position_tensor.shape[0]-1

        changes = ind[1:] - ind[:-1]
        addindexes = np.where(changes > 0)[0]
        subindexes = np.where(changes < 0)[0]

        sub = solutions[k][:, 1:] - solutions[k][:, :-1]
        neutralindexes = np.unique(np.where(sub < 0)[1])
        neutralindexes = np.setxor1d(neutralindexes, subindexes)

        for index in neutralindexes:
            pre = solutions[k][:, index]
            post = solutions[k][:, index+1]
            changecount = 0
            while len(np.where(post - pre < 0)[0]) > 0:

                post = np.append([total_length], post)
                pre = np.append(pre, 0)

                changecount += 1

            for i in range(changecount):
                addindexes = np.sort(np.append(addindexes, index))
                subindexes = np.sort(np.append(subindexes, index))

            changes[index] = -changecount
            ind[index] += changecount


        for index in np.where(np.abs(changes) > 1)[0]:
            if changes[index] < 0:
                for i in range(np.abs(changes[index])-1):
                    subindexes = np.sort(np.append(subindexes, index))
            else:
                for i in range(np.abs(changes[index])-1):
                    addindexes = np.sort(np.append(addindexes, index))

        truefrags = len(subindexes)




        if len(subindexes) < len(addindexes):
            subindexes = np.append(
                subindexes, (
                    np.ones((
                        len(addindexes) - len(subindexes)))*(nsteps-1)).astype(int))


        fragmentspertraj.append(len(subindexes))

        for m in range(min(len(subindexes), len(addindexes))):
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
                            iterind = iterind + int(min(0, traj_ind[minusloc[n]]))

                            fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten())

                        fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())

                    else:
                        for n in range(len(minusloc)-1):

                            iterind = iterind + min(0, traj_ind[minusloc[n]])

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

        return fragtimes, fragarray

