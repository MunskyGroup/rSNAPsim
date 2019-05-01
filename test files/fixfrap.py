# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 16:52:37 2019

@author: wsraymon
"""

import numpy as np

test_ivec = np.load('test_ivec.npy')
test_results = np.load('test_ribs.npy')
pv = np.load('pv.npy')

ti = 500
stop_frap = ti+20
tvec = np.linspace(0,999,1000)

startindex = np.where(tvec >= ti)[0][0]
stopfrap = np.where(tvec >= stop_frap)[0][0]

traj = test_results[0, :].reshape((200, len(tvec)))
frapribs = traj[:,startindex:stopfrap]


genelength = 405




def pullfrags(solution,genelength):
    solutions = [solution]
    fragmented_trajectories = []
    fragtimes = []
    endfragtimes = []
    maxlen = 0
    
    fragmentspertraj= []
    for k in range(1):
        ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
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
            subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(len(tvec)-1)).astype(int))
            
        
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
                            iterind = iterind + min(0,traj_ind[minusloc[n]])
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
                if addindexes[m]+1  + len(fragment) > len(tvec):
                    endfragtimes.append(len(tvec))
                else:
                    endfragtimes.append(addindexes[m]+1  + len(fragment))
                   
                
                fragmented_trajectories.append(fragment)
                #if m <= truefrags:
                    #kes.append(genelength/truetime[len(fragment)])
        
                if len(fragment) > maxlen:
                    maxlen = len(fragment)
                
    
        fragarray = np.zeros((len(fragmented_trajectories), maxlen))
        for i in range(len(fragmented_trajectories)):
            fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
            
            
    return fragarray, fragtimes, fragmentspertraj,endfragtimes







fragarray,fragtimes,fragmentspertraj,endfragtimes = pullfrags(traj,405)
affected_frags = []
fragindexes = []
for i in range(len(fragtimes)):
   if  np.sum([fragtimes[i]> np.array([startindex, stop_frap]), endfragtimes[i] > np.array([startindex, stop_frap])]) in [1,2,3]:
       affected_frags.append(i)
       fragindexes.append([fragtimes[i],endfragtimes[i]])
    

#affected_frags = np.intersect1d(np.where(np.array(fragtimes) >=  startindex), np.where(np.array(fragtimes)<= stop_frap))
af = fragarray[affected_frags]
findexes = np.array(fragindexes)
frange = findexes[:,1]-stop_frap
afterfrapribs = findexes[np.where(frange > 0 )]
relevantfrags = np.array(affected_frags)[np.where(frange > 0 )]
stopfrapindex = stop_frap - afterfrapribs[:,0]

rfrags = fragarray[relevantfrags]
np.diag(rfrags[:,stopfrapindex])
laglen = afterfrapribs[:,1] - stop_frap
posistions_at_end_of_FRAP = np.diag(rfrags[:,stopfrapindex])

lenfrap = stop_frap - startindex
offset = pv[posistions_at_end_of_FRAP.astype(int)]
trailing_intensity = np.zeros((max(laglen)+ lenfrap))


for j in len(fragindexes):
    frag = fragindexes[j]
    fragpos = af[j]
    fragstart = frag[0]
    fragend = frag[1]
    for i in range(fragend - fragstart):
        if fragstart + i > startindex:
            trailing_intensity[i] -= pv[fragpos[i]]
            
    
    


for i in range(len(laglen)):
    trailing_intensity[:laglen[i]] -= offset[i]


