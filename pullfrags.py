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




case3 = np.loadtxt('fragerror2.txt')


solutions = [case1,case2,case3]

fragmented_trajectories = []
fragtimes = []
maxlen = 0
kes = []

k = 2
ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
changes = ind[1:] - ind[:-1]
addindexes = np.where(changes > 0)[0]
subindexes = np.where(changes < 0)[0]
genelength = 465
sub = solutions[k][:,1:] - solutions[k][:,:-1]
neutralindexes = np.unique(np.where(sub < 0)[1])
neutralindexes = np.setxor1d(neutralindexes, subindexes)


for index in neutralindexes:
    pre = solutions[k][:,index]
    post = solutions[k][:,index+1]
    changecount = 0
    while len(np.where(post - pre < 0)[0]) > 0:

        post = np.append([genelength],post)
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
    subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(1000-1)).astype(int))
    



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
                    iterind = iterind + min(0,traj_ind[minusloc[n]])
                    fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                    
                    
                    
                    
      
    
          
                
                fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                
            else:
                print('above m')
                print(addindexes[m])
                print(subindexes[m])
                for n in range(len(minusloc)-1):
                    

                    iterind = iterind + min(0,traj_ind[minusloc[n]])
                    
                   
                    fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                    
                    
                fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
                
    
            #print(traj[0, minusloc[-1]:].flatten())
            
        else:

            fragment = solutions[k][startind][addindexes[m]:subindexes[m]].flatten()
           
        
        fragtimes.append(addindexes[m])
           
        
        fragmented_trajectories.append(fragment)
        
        #kes.append(genelength/truetime[len(fragment)])

        if len(fragment) > maxlen:
            maxlen = len(fragment)
'''
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

'''
fragarray = np.zeros((len(fragmented_trajectories), maxlen))
for i in range(len(fragmented_trajectories)):
    fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]