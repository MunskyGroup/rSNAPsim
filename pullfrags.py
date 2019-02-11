# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 14:44:49 2019

@author: wsraymon
"""

'''
test fragment pulling script
'''


import numpy as np
import rSNAPsim


case1 = np.zeros((3,100))
case2 = np.zeros((3,100))
case3 = np.zeros((4,100))

case1[0,0:26] = np.linspace(0,25,26)
case1[1,10:26] = np.linspace(0,15,16)
case1[0,26:26+10] = np.linspace(16,25,10)


case2[0,0:26] = np.linspace(0,25,26)
case2[1,10:26] = np.linspace(0,15,16)
case2[0,26:26+10] = np.linspace(16,25,10)
case2[1,26:36] =  np.linspace(1,10,10)
case2[0,36:51] = np.linspace(11,25,15)



solutions = [case1,case2]


fragmented_trajectories = []
fragtimes = []
maxlen = 0
kes = []

k = 1
ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
changes = ind[1:] - ind[:-1]
addindexes = np.where(changes > 0)[0]
subindexes = np.where(changes < 0)[0]

if len(subindexes) <= len(addindexes):
    for m in range(len(subindexes)):
        traj = solutions[k][:, addindexes[m]:subindexes[m]+1]
    
        traj_ind = changes[addindexes[m]:subindexes[m]+1]

        startind = ind[addindexes[m]]
        minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
        fragment = np.array([])


        if subindexes[m]-addindexes[m] > 0:
            if len(minusloc) > 1:
                for n in range(len(minusloc)-1):
                    print(fragment)
                    potential_frag = traj[startind-n, minusloc[n]+1:minusloc[n+1]].flatten()
                    
                    discontinuity = np.where(potential_frag[1:]-potential_frag[:-1]<0)[0]
                    #if len(discontinuity) !=0:
                        
                    
                    fragment = np.append(fragment, traj[startind-n, minusloc[n]+1:minusloc[n+1]+1].flatten())
                    
                    
                    
  

                fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                print(fragment)
                
                #print(traj[0, minusloc[-1]:].flatten())
                
            else:
                fragment = solutions[k][startind][addindexes[m]:subindexes[m]].flatten()
                print([addindexes[m]])
            
            fragtimes.append(addindexes[m])
               
            
            fragmented_trajectories.append(fragment)
            
            #kes.append(genelength/truetime[len(fragment)])

            if len(fragment) > maxlen:
                maxlen = len(fragment)

else:
    for m in range(len(addindexes)):
        traj = solutions[k][:, addindexes[m]:subindexes[m]]
        traj_ind = changes[addindexes[m]:subindexes[m]]
     
        startind = ind[addindexes[m]]
        minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
        fragment = np.array([])


        if subindexes[m]-addindexes[m] > 0:
            if len(minusloc) > 1:
                for n in range(len(minusloc)-1):
                    fragment = np.append(fragment, traj[startind-n, minusloc[n]:minusloc[n+1]].flatten())

                fragment = np.append(fragment, traj[0, minusloc[-1]:].flatten())
            else:
                fragment = solutions[k][startind][addindexes[m]:subindexes[m]].flatten()
            fragmented_trajectories.append(fragment)
            fragtimes.append(addindexes[m])
         
            #kes.append(25/truetime[len(fragment)])

            if len(fragment) > maxlen:
                maxlen = len(fragment)


fragarray = np.zeros((len(fragmented_trajectories), maxlen))
for i in range(len(fragmented_trajectories)):
    fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]