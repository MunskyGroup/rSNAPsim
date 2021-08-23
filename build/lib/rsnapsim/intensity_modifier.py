# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 09:24:35 2020

@author: willi
"""

import pandas as pd
import matplotlib.pyplot as plt
import time
import numpy as np





def modify_intensity(intensity_vector, SNR, noise_type='AGWN'):
    '''
    given an intensity vector and signal to noise ratio this will add noise 
    according to SNR = mu_signal**2 / var_noise**2
    
    
    awgn = additive white gaussian noise

    
    '''
    if noise_type.lower() in ['agwn']:
        av_intensity = np.mean(intensity_vector)
        var_intensity = np.var(intensity_vector)
        
        sig_intensity = 10 * np.log10(av_intensity)
        noise_sig_db = sig_intensity - SNR
        noise_par =  10 ** (noise_sig_db / 10)
        
        
        intensity_with_modifications = intensity_vector + noise_par*np.random.randn(*intensity_vector.shape)

    if noise_type.lower() in ['poisson']:
        
        var_intensity = np.var(intensity_vector)
        
        var_poission = var_intensity / SNR
        
        intensity_with_modifications = intensity_vector + np.random.poission(var_poission, size = intensity_vector.shape) 

    if noise_type.lower() in ['pink']:
        ivec = np.copy(intensity_vector)
        av_intensity = np.mean(ivec)
        sig_intensity = 10 * np.log10(av_intensity)
        noise_sig_db = sig_intensity - SNR
        noise_par =  10 ** (noise_sig_db / 10)      
       
        for i in range(ivec.shape[0]):
            
            
            ivec[i,:] = ivec[i,:] + voss(ivec.shape[1])*noise_par
            
        intensity_with_modifications = ivec


    return intensity_with_modifications




def slice_intensity(intensity_vector, total_frames, framerate):
    
    '''
    Given an intensity vector, this will slice the vector into the total frames and framerate desired
    
    Note this will throw away some time points at the end to do the correct reshaping
    
    '''
    remove = intensity_vector.shape[1]%int((framerate*total_frames))
    
    intensity_vector = intensity_vector[:,:-remove]
    s0 = intensity_vector.shape[0]*intensity_vector.shape[1] / int((framerate*total_frames))

    intensity_vector = intensity_vector.reshape(int(s0), int((framerate*total_frames)) )
    
    
    #intensity_with_modifications = intensity_vector[0,::framerate]
    
    
    return intensity_vector



def voss(nrows,ncols=16):
    """Generates standardized pink noise using the Voss-McCartney algorithm.
    
    https://www.dsprelated.com/showarticle/908.php
    
    nrows: number of values to generate
    rcols: number of random sources to add
    
    returns: NumPy array of standardized pink noise
    """
    array = np.empty((nrows, ncols))
    array.fill(np.nan)
    array[0, :] = np.random.random(ncols)
    array[:, 0] = np.random.random(nrows)
    
    # the total number of changes is nrows
    n = nrows
    cols = np.random.geometric(0.5, n)
    cols[cols >= ncols] = 0
    rows = np.random.randint(nrows, size=n)
    array[rows, cols] = np.random.random(n)

    df = pd.DataFrame(array)
    df.fillna(method='ffill', axis=0, inplace=True)
    total = df.sum(axis=1)

    signal = total.values - np.mean(total.values)
    signal =  signal/ np.std(signal)
    return signal




def make_training_dataset(dataset, framerate, total_frames, SNR, noise_type = 'agwn' ):
    ivec = slice_intensity(dataset,total_frames,framerate)
    ivec_with_noise = modify_intensity(ivec, SNR, noise_type = noise_type)
    
    return ivec_with_noise
    
def make_testing_dataset(dataset,framerate,total_frames,SNR, noise_type = 'agwn' ):
    ivec = dataset[:,::framerate][:,:total_frames] 
    ivec_with_noise = modify_intensity(ivec, SNR, noise_type = noise_type)
    
    return ivec_with_noise


###############################
# How to set up the datasets
###############################

# Loop over framerate and SNR and save accuracy.

framerate = 1
total_frames = 2000
SNR = .2

A = np.loadtxt('C:/Users/willi/Documents/GitHub/rSNAPsim/large_no_noise/trainingA_1s.txt' ,delimiter = ',')
A_training = make_training_dataset(A, framerate,total_frames,SNR)

B = np.loadtxt('C:/Users/willi/Documents/GitHub/rSNAPsim/large_no_noise/trainingB_1s.txt' ,delimiter = ',')
B_training = make_training_dataset(B, framerate,total_frames,SNR)

C = np.loadtxt('C:/Users/willi/Documents/GitHub/rSNAPsim/large_no_noise/trainingC_1s.txt' ,delimiter = ',')
C_training = make_training_dataset(C,framerate,total_frames,SNR)

# testing = np.loadtxt('./mixture_framerate1s.txt' ,delimiter = ',')
# testing = make_testing_dataset(testing, framerate,total_frames,SNR)

#plt.scatter(A_training[:1,0:-10],A_training[:1,10:] ,alpha=.3 ); plt.scatter(B_training[:1,0:-10],B_training[:1,10:] ,alpha=.3 );plt.scatter(C_training[:1,0:-10],C_training[:1,10:] ,alpha=.3 )



