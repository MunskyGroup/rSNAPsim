# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 18:08:58 2020

@author: willi
"""

import numpy as np
from . import PropensityFactory, ProbeVectorFactory
PropensityFactory = PropensityFactory.PropensityFactory
ProbeVectorFactory = ProbeVectorFactory.ProbeVectorFactory
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
import matplotlib.patches as mpatches
import time


path_to_cpp = ''
path_to_gen = ''
path_to_trna = ''

import os
for root, dirs, files in os.walk(".", topdown=False):
   for branch in dirs:
       if 'ssa_cpp' in branch:
           path_to_cpp = os.path.join(root, branch)
       if 'generalized_cpp' in branch:
           path_to_gen = os.path.join(root, branch)
       if 'trna_ssa' in branch:
           path_to_trna = os.path.join(root, branch)
if path_to_cpp != '':
    try:
        cwd = os.getcwd()
        os.chdir(path_to_cpp)
        
        import ssa_translation
        import ssa_translation_lowmem
        os.chdir(cwd)
    except:
        os.chdir(cwd)
    

if path_to_gen != '':
    try:
        cwd = os.getcwd()
        
        os.chdir(path_to_gen)
        print('importing C++ models')
        import ssa_translation_generic
        import ssa_translation_generic_lowmem
        print('c++ models loaded successfully')
        os.chdir(cwd)
    except:
        os.chdir(cwd)

if path_to_trna !='':
    try:
        cwd = os.getcwd()
        
        os.chdir(path_to_gen)
        print('importing C++ tRNA models')
        import ssa_trna
        print('c++ models loaded successfully')
        os.chdir(cwd)
    except:
        os.chdir(cwd)   
        
        
    
try: 
    import ssa_translation_lowmem
except:
    pass

class ModelBuilder():
    '''
    Container class for the solvers
    '''
    def __init__(self, time=None,xi=None):
        self.k = None
        self.k_bind = None
        self.k_term = None
        self.multiframe = True
        self.additional_rxns = {}
        self.probe_locations = None
        self.colors = 1
        
        self.default_conditions = {'low_mem':True,
                                   'perturb':[0,0,0],
                                   'leaky_probes':False,
                                   'bins':None,
                                   'k_probe':0,
                                   'footprint':9,
                                   'burnin':500,
                                   'record_stats':False}
        
        self._poi = None
        self.k_enters = np.array([])
        self.k_stops = np.array([])
        self.k_jumps = np.array([])
        self.k_pauses = np.array([])
        self.frames_used = 1
        self.probe = None
        
        
    @property    
    def n_extra_rxns(self):
        return len(self.k_enters) + len(self.k_stops) + len(self.k_jumps) + len(self.k_pauses) 
        
    def add_k(self, poi, kelong_mean = 10):
        self._poi = poi
        ks = PropensityFactory().get_k_3_frame(poi.nt_seq, kelong_mean)        
        self.k = ks

        L = len(poi.aa_seq)
        self.k[L-1] = 0
        self.k[L*2-1  -1] = 0
        self.k[L*3-2 -1] = 0
        self.k_stops = np.array([[L,0,10],[L-1,1,10],[L-1,2,10]],dtype=np.float64)
        
    
    def add_jumps(self, locs, frames, dest_locs, dest_frames, rates):
        if isinstance(locs,int):
            locs = [locs]
        if isinstance(frames,int):
            frames = [frames]
        if isinstance(dest_locs,int):
            dest_locs = [dest_locs]
        if isinstance(dest_frames,int):
            dest_frames = [dest_frames]
        if isinstance(rates,float):
            rates = [rates]
        self.__check_frames(frames)
        self.k_jumps = np.array((locs,frames,dest_locs,dest_frames,rates )).T.astype(np.float64)
        
    
    def add_enters(self,locs, frames, rates):
        if isinstance(locs,int):
            locs = [locs]
        if isinstance(frames,int):
            frames = [frames]
        if isinstance(rates,float):
            rates = [rates]        
        self.__check_frames(frames)
        self.k_enters = np.array((locs,frames,rates )).T.astype(np.float64)
    
    def add_stops(self,locs, frames, rates):
        if isinstance(locs,int):
            locs = [locs]
        if isinstance(frames,int):
            frames = [frames]
        if not isinstance(rates,list):
            rates = [rates]        
        self.__check_frames(frames)

        self.k_stops = np.array((locs,frames,rates )).T.astype(np.float64)
        
    
    def clear_model(self):
        self._poi = None
        self.k_enters = np.array([])
        self.k_stops = np.array([])
        self.k_jumps = np.array([])
        self.k_pauses = np.array([])
        self.frames_used = 1        
        self.probe = None
    
    def add_pauses(self,locs, frames, rates):
        if isinstance(locs,int):
            locs = [locs]
        if isinstance(frames,int):
            frames = [frames]
        if isinstance(rates,float):
            rates = [rates]        
            
        self.__check_frames(frames)
        self.k_pauses = np.array((locs,frames,rates )).T.astype(np.float64)
    
    def __check_frames(self,frames):
        for frame in frames:
            if frame > self.frames_used:
                self.frames_used = frame
        
    def visualize_transcript(self):
        
        probe = self._poi.probe_loc
        #self._poi.generate_3frame_tags()
        
        keys = []
        for dics in self._poi.multiframe_epitopes:
            for key in dics.keys():
                if key not in keys:
                    keys.append(key)
        
        
        fig,ax = plt.subplots(1,1,dpi=300)
        N = len(self._poi.aa_seq)
        ncolors = len(keys)
        xs = [.1, 2.1, 4.1][::-1]
        fr = ["+0","+1","+2"]
        
        
        opt = {'head_width': .1, 'head_length': .1,
        'length_includes_head': True}
        
        
        cmap = cm.get_cmap('viridis')
        colors = cmap(np.linspace(.01,.95, ncolors))
        
        go_color = '#2A9D8F'
        stop_color = '#E76F51'
        jump_color = '#843b62'
        pause_color = '#e9c46a'
        
        
        
        if ncolors <4:
            colors = ['#64ff21','#08ffff','#fb00ff']
        
        for i in range(3):
            

            
            rectangle =  mpatches.Rectangle((0,xs[i]), N ,.7,linewidth=1,edgecolor='k',facecolor='darkgray')
    
            ax.add_patch(rectangle)
            
            color = []
            color_dict = self._poi.multiframe_epitopes[i]
            location = []
            for n in range(len(keys)):
                if keys[n] in color_dict.keys():
                
                    location.append(color_dict[keys[n]])
                    color.append( [n for x in color_dict[keys[n]]])
                
          
            
            colorlabels = ['Color %d'% j for j in range(ncolors)    ]
            
        
            color = [j for i in color for j in i]
            location = [j for i in location for j in i]
            
            for c in range(ncolors):
                ax.plot([-10,-10],[-4,-4],color = colors[c]  )  #fix the legend colors
            
            
            for c,loc in zip(color,location):
                ax.plot([loc,loc],[xs[i] ,xs[i]+.8],color = 'k',lw=2 )
                ax.plot([loc,loc],[xs[i] ,xs[i]+.8],color = colors[c]  )
                
            ax.set_ylim([-.1,8])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.set_xlabel('codon')
            
            ax.axes.get_yaxis().set_visible(False)
            if i == 0:
                leg1 = ax.axes.legend(colorlabels,loc=9)
                ax.text(-20,7.5,'Transcript Name: %s' % self._poi.name)
                ax.text(-20,6.7,'Total Length: %d aa' % self._poi.total_length)
                ax.text(-20,5.9,'Seq: %s ...' % self._poi.aa_seq[:10])
                

                
            ax.text(-N*.07 ,xs[i] + .4,fr[i])
                
        for j in range(self.k_enters.shape[0]):
            yloc = xs[int(self.k_enters[j][1])]
            loc = int(self.k_enters[j][0])
            ax.plot([loc,loc],[yloc, yloc+.8],color = 'k',lw=3  )
            ax.plot([loc,loc],[yloc, yloc+.8],color = go_color,lw=2  )

            ax.annotate("",
                    xy=(loc, yloc+.85), xycoords='data',
                    xytext=(loc-1, yloc+.85), textcoords='data',
                    arrowprops=dict(arrowstyle="-|>",lw=1.5,facecolor='k',edgecolor='k',
                                    connectionstyle="arc3"),
                    )            
            
            ax.annotate("",
                    xy=(loc, yloc+.85), xycoords='data',
                    xytext=(loc-1, yloc+.85), textcoords='data',
                    arrowprops=dict(arrowstyle="-|>",lw=1,facecolor=go_color,edgecolor=go_color,
                                    connectionstyle="arc3"),
                    )
            

        for j in range(self.k_stops.shape[0]):
            yloc = xs[int(self.k_stops[j][1])]
            loc = int(self.k_stops[j][0])
            ax.plot([loc,loc],[yloc, yloc+.8],color = 'k',lw=3  )
            ax.plot([loc,loc],[yloc, yloc+.8],color = stop_color,lw=2  )
            
        for j in range(self.k_jumps.shape[0]):
            
            yloc1 = xs[int(self.k_jumps[j][1])]
            loc1 = int(self.k_jumps[j][0])
            yloc2 = xs[int(self.k_jumps[j][3])]
            loc2 = int(self.k_jumps[j][2])   

            ax.plot([loc1,loc1],[yloc1, yloc1+.8],color = 'k',lw=3  )
            ax.plot([loc2,loc2],[yloc2, yloc2+.8],color = 'k',lw=3  )
            
            ax.plot([loc1,loc1],[yloc1, yloc1+.8],color = jump_color,lw=2  )
            ax.plot([loc2,loc2],[yloc2, yloc2+.8],color = jump_color,lw=2  )
            
            if yloc2 == yloc1:
                ax.plot([loc1,loc1],[yloc1, yloc1-.3],color=plt.rcParams['axes.edgecolor'])
                ax.plot([loc1,loc2],[yloc1-.3, yloc1-.3],color=plt.rcParams['axes.edgecolor'])
                
                ax.annotate("",
                        xy=(loc2, yloc2), xycoords='data',
                        xytext=(loc2, yloc1-.3), textcoords='data',
                        arrowprops=dict(arrowstyle="-|>",lw=2,
                                        connectionstyle="arc3"),
                        )
                
            else:
                dx = loc2 - loc1
            
                dy = yloc2- yloc1 
                
                if yloc1 > yloc2:
                    dx = loc2 - loc1
            
                    dy = yloc2+.8- yloc1 
                
                   # ax.arrow(loc1,yloc1,dx,dy,color = 'k',lw=2)#,  head_width=5, head_length=.2)
                    ax.annotate("",
                        xy=(loc2, yloc2+.7), xycoords='data',
                        xytext=(loc1, yloc1), textcoords='data',
                        arrowprops=dict(arrowstyle="-|>",lw=2,
                                        connectionstyle="arc3"),
                        )
                else:
                    dx = loc2 - loc1
            
                    dy = yloc2- yloc1+.8
                    #ax.arrow(loc1,yloc1,dx,dy,color = 'k',lw=1)
                    
                    ax.annotate("",
                        xy=(loc2, yloc2), xycoords='data',
                        xytext=(loc1, yloc1+.8), textcoords='data',
                        
                        
                        arrowprops=dict(arrowstyle="-|>",lw=2,
                                        connectionstyle="arc3"),
                        )

        for j in range(self.k_pauses.shape[0]):
            yloc = xs[int(self.k_pauses[j][1])]
            loc = int(self.k_pauses[j][0])
            ax.plot([loc,loc],[yloc, yloc+.8],color = 'k',lw=3  )
            ax.plot([loc,loc],[yloc, yloc+.8],color = pause_color,lw=2  )         
            
            
            
            
        custom_lines = [Line2D([0], [0], color=go_color, lw=2),
                Line2D([0], [0], color=stop_color, lw=2),
                Line2D([0], [0], color=pause_color, lw=2), 
                Line2D([0], [0], color=jump_color, lw=2)]
        
        ax.legend(custom_lines, ['Enter', 'Stop', 'Pause','Jump'])
        ax.add_artist(leg1)
        #ax.set_facecolor('darkgray')
        fig.show()    
        
        

    def __generate_additional_ks(self,enters,pauses,jumps,stops,L):
        
        def frame_check_1(L,arr):        
            return (L- arr[:,1]+1)*(arr[:,1]>0) + L*(arr[:,1]>1)
        
        def frame_check_3(L,arr):        
            return (L- arr[:,3]+1)*(arr[:,3]>0) + L*(arr[:,3]>1)            
                    
        def gen_ks_1_loc(L,arr):
            arr[:,0] = arr[:,0]+frame_check_1(L,arr)
            arr[:,1] = arr[:,2]    
            arr = arr[:,0:2]
            max_arr = np.max( arr[:,0])     
            return arr,max_arr
        
        def gen_ks_3_loc(L,arr):
            arr[:,0] = arr[:,0]+ frame_check_1(L,arr)     
            arr[:,1] = arr[:,2]+ frame_check_3(L,arr)
            arr[:,2] = arr[:,4]
            arr = arr[:,0:3]
            max_arr = max([np.max( arr[:,0]),np.max( arr[:,1])])
            return arr,max_arr
    
        max_enter = 0
        max_pause = 0
        max_stop = 0
        max_jump = 0
        k_jumps = np.copy(jumps)
        k_pauses = np.copy(pauses)
        k_stops = np.copy(stops)
        k_enters = np.copy(enters)
        if len(k_enters) != 0:
            k_enters,max_enter = gen_ks_1_loc(L,k_enters)
    
        if len(k_pauses) != 0:
            k_pauses,max_pause = gen_ks_1_loc(L,k_pauses)
    
        if len(k_stops) != 0:
            k_stops,max_stop = gen_ks_1_loc(L,k_stops)
        
        if len(k_jumps) != 0:
            k_jumps,max_jump = gen_ks_3_loc(L,k_jumps)
            
        max_loc = max(max_jump,max_stop,max_pause,max_enter)
        
        if max_loc <=L: 
            frames_used = 0
        if max_loc > L:
            frames_used = 1
        if max_loc > 2*L-1 :
            frames_used = 2
        
        return k_enters, k_pauses, k_stops, k_jumps, frames_used            

    def generate_probe(self):
        keys = []
        for frame in self._poi.multiframe_epitopes:
            for key in frame.keys():
                if key not in keys:
                    keys.append(key)
                    
        ncolors = len(keys)
        
        N = len(''.join(self._poi.multiframe_aa_seq))
        Ns = [0,len(self._poi.multiframe_aa_seq[0]),len(self._poi.multiframe_aa_seq[1])+len(self._poi.multiframe_aa_seq[0])   ]
        probe = np.zeros((ncolors,N) )
        
        for j in range(3):
            dic = self._poi.multiframe_epitopes[j]
            for i in range(len(keys)):
                if keys[i] in dic.keys():
                    epis = dic[keys[i]]
                    
                    probe[i,Ns[j] + np.array(epis)] = 1
        return probe

    def run_ssa(self,t, n_traj=10,low_mem = True ):
        
        try: 
            self.probe.shape
        except:
            self.probe = self.generate_probe()
        
        
        probe = self.probe.astype(int).copy(order='C')
        t_array = t
        t0 = t[0]
        ncolors = probe.shape[0]

        N_rib = 200
        #result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
        #kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
        n_trajectories = n_traj
        start = time.time()
        
        lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))
        
        all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
        
        if low_mem == False:
            all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
            all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
            all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
            nribs = np.array([0],dtype=np.int32)
            all_ribs = np.zeros((n_trajectories,1))
            
        else:
            all_results = np.zeros((n_trajectories,ncolors,len(t_array)),dtype=np.int32)
        
        seeds = np.random.randint(0,0x7FFFFFF,n_trajectories, dtype= np.int32)
        
        k_enters,k_pauses,k_stops,k_jumps,frames_used = self.__generate_additional_ks(self.k_enters,self.k_pauses,self.k_jumps,self.k_stops,self._poi.total_length)
        
        #k_add = np.hstack((self.k_enters.flatten(),self.k_pauses.flatten(),self.k_stops.flatten(),self.k_jumps.flatten() ))
        k_add = np.hstack((k_enters.flatten(),k_pauses.flatten(),k_stops.flatten(),k_jumps.flatten() ))
        
        

        max_ribs = 200
        
        

        n_enters = self.k_enters.shape[0]
        n_pauses = self.k_pauses.shape[0]
        n_stops = self.k_stops.shape[0]
        n_jumps = self.k_jumps.shape[0]
        

                
        if low_mem == False:
            kelong = np.array(self.k).astype(np.float64)
            for i in range(n_trajectories):
                result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
                frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
                
                ribtimes = np.zeros((400),dtype=np.float64)
                coltimes = np.zeros((400),dtype=np.int32)
                ssa_translation_lowmem.run_SSA_generic(result,ribtimes,coltimes, kelong,frapresult,t_array, np.array([0,0,0],dtype=np.float64), seeds[i],nribs, k_add.flatten()  ,n_enters,n_pauses,n_stops,n_jumps)
                all_results[i,:] = result
                all_frapresults[i,:] = frapresult
                all_coltimes[i,:] = coltimes
                all_ribtimes[i,:] = ribtimes
                all_ribs[i,:] = nribs[0]
            
        else:
            
            
        
            kelong = np.array(self.k).astype(np.float64)
            
            for i in range(n_trajectories):
                result = np.zeros((ncolors,len(t_array)),dtype=np.int32)    
                frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
                
                

                ssa_translation_lowmem.run_SSA_generic_lowmem(result,kelong ,
                                                              frapresult,
                                                              t_array, 
                                                              np.array([0,0,0],dtype=np.float64),
                                                              int(seeds[i]),
                                                              k_add.flatten()  ,
                                                              int(n_enters),
                                                              int(n_pauses),
                                                              int(n_stops),
                                                              int(n_jumps), 
                                                              probe.astype(np.int32),
                                                              int(ncolors) ,
                                                              int(max_ribs) )
                all_results[i,:,:] = result
                #all_frapresults[i,:] = frapresult
        

      
        #traj = all_results[0,:].reshape((N_rib,len(t_array))).T
        return all_results
    
    
    