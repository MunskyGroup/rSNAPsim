# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 10:28:04 2019

@author: wsraymon
"""

#adjusting simulated cell

import rSNAPsim
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

import matplotlib.patches as mpatches
import matplotlib.animation as animation
from matplotlib.collections import PatchCollection

from matplotlib.colors import LinearSegmentedColormap
import scipy


def normal_pdf(x, mean, var):
    return np.exp(-(x - mean)**2 / (2*var))

from matplotlib.colors import Normalize


class GenericSSA():

    '''
    Generic SSA solver - used for the simulated cell animations
    '''

    def __init__(self,type='linear'):

        self.time_variant = False
        self.xi=np.array([])
        self.ti= None
        self.tf=None
        self.S=np.array([])
        self.type=type
        self.ptimes=100
        self.params={}
        if type=='linear':
            #self.fast_rxn = 0.5
            self.W0=np.array([])
            self.W1=np.array([])
        if type == 'nonlinear':
            #self.fast_rxn = 0.5
            self.P=lambda x,t:None



    def gettvec(self):
        return np.linspace(self.ti,self.tf,self.ptimes)

    def updateparams(self):
        for k, v in self.params.iteritems():
            setattr(self, k, v)

    def _run_trajectory(self):#!!!!!!!!!!!!!!!!!renamed run to solve(big deal)
        x=self.xi
        t=self.ti
        __n=len(x)
        self.time=self.gettvec()
        data=np.zeros((len(self.xi),self.ptimes))
        ip=0

        if self.type=='linear':
            if self.time_variant == False:
                while t<self.tf:
                    rate=np.atleast_2d(np.dot(self.W1,x))+self.W0
                    rate=np.cumsum(rate)
                    with np.errstate(divide='ignore', invalid='ignore'):
                        t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()
                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            data[:,ip]=x.reshape(__n)
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i
                            break
                    x=x+np.atleast_2d(self.S[:,event]).T

            else:


                x = np.concatenate((np.array([0]),self.xi.flatten()))

                #x = self.xi
                __n=len(x)
                self.time=self.gettvec()
                data=np.zeros((len(self.xi),self.ptimes))
                a,b = self.S.shape
                S = np.vstack((np.zeros(b),self.S))
                S = np.hstack((np.zeros((a+1,1)),S))
                while t<self.tf:
                    __n=len(x)
                    self.time=self.gettvec()
                    data=np.zeros((len(self.xi),self.ptimes))
                    a,b = self.S.shape
                    S = np.vstack((np.zeros(b),self.S))
                    S = np.hstack((np.zeros((a+1,1)),S))
                    while t<self.tf:
                        trate=self.get_P(x[1:],t)

                        rate = np.concatenate((np.array([self.fast_rxn]),trate))
                        rate=np.cumsum(rate)


                        t=(t-np.log(np.random.rand(1))/rate[-1])
                        ro=rate[-1]*np.random.rand()

                        while t>self.time[ip]:
                            if t>self.tf:
                                b = len(self.time[ip:])
                                fill = np.repeat(x[1:],b)
                                data[:,ip:]=fill.reshape(__n-1,b)
                                return data
                            else:
                                #data[:,ip]=x.reshape(__n)
                                data[:,ip]=x[1:]
                                ip=ip+1
                        for i in range(len(rate)):
                            if rate[i]>=ro:
                                event=i

                                break

                        x=x+S[:,event].ravel()
                    '''
                    rate=np.atleast_2d(np.dot(self.W1(t),x))+self.W0(t)

                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    print(t)
                    ro=rate[-1]*np.random.rand()
                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            data[:,ip]=x.reshape(__n)
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i
                            break
                    x=x+np.atleast_2d(self.S[:,event]).T

                    '''



        elif self.type=='nonlinear':
            if self.time_variant == True:  #if time variant use fast reaction
                x = np.concatenate((np.array([0]),self.xi.flatten()))

                #x = self.xi
                __n=len(x)
                self.time = self.gettvec()
                data = np.zeros((len(self.xi), self.ptimes))
                a, b = self.S.shape
                S = np.vstack((np.zeros(b), self.S))
                S = np.hstack((np.zeros((a+1, 1)), S))
                while t < self.tf:
                    trate=self.get_P(x[1:],t)
                    rate = np.concatenate((np.array([self.fast_rxn]),trate))
                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()

                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x[1:],b)
                            data[:,ip:]=fill.reshape(__n-1,b)
                            return data
                        else:
                            #data[:,ip]=x.reshape(__n)
                            data[:,ip]=x[1:]
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i

                            break

                    x=x+S[:,event].ravel()

            else:   #if not time variant ignore fast reaction

                x = self.xi.flatten()

                #x = self.xi
                __n=len(x)
                self.time=self.gettvec()

                while t<self.tf:
                    rate=self.get_P(x,t)

                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()

                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            #data[:,ip]=x.reshape(__n)
                            data[:,ip]=x
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i

                            break

                    x=x+self.S[:,event].ravel()



        else:
            'Error'
        self.data=data
        return data

    def _solve(self,n):
        __data=np.zeros((len(self.xi),self.ptimes,n))

        for i in range(n):
            __d=self._run_trajectory()
            __data[:,:,i]=__d
        self.data = __data

        return __data

    def setpar(self,key,val):
        self.params[key]=val

    def get_dist(self,specID=0):
        '''
        build distribution (non-normalized and pdf)
        of rna for the model)
        '''
        n_specs, n_times, n_traj = self.data.shape
        n_specs = int(n_specs)
        n_times = int(n_times)
        n_traj = int(n_traj)
        specID = int(specID)
        max_rna = int(np.max(self.data[specID,:,:]))
        self.pdf = np.zeros((n_times,max_rna+1))
        self.fdist = np.zeros((n_times,max_rna+1))
        for i in range(n_times):
            ind = int(i)
            for j in range(n_traj):
                jnd = int(j)
                self.fdist[ind,int(self.data[int(specID),ind,jnd])] +=1
            self.pdf[ind,:] = self.fdist[ind,:] / np.sum(self.fdist[ind,:])

    def get_traj(self,specID=0,ntraj='all'):
        n_specs, n_times, n_traj = self.data.shape
        n_specs = int(n_specs)
        n_times = int(n_times)




        if isinstance(specID,int):
            if ntraj == 'all':
                return self.data[specID],ntraj
            else:
                try:
                    ntraj = ntraj.flatten().astype(int).tolist()
                except:
                    ntraj = int(ntraj)
                    pass

                return self.data[specID][:,ntraj],ntraj
        else:
            if specID == 'all':
                if ntraj == 'all':
                    return self.data,ntraj
                else:

                    try:
                        ntraj = ntraj.flatten().astype(int).tolist()
                    except:
                        pass

                    return self.data,ntraj

    def get_means(self,specID=0):
        '''
        get the first moment.
        '''
        n_specs, n_times, n_traj = self.data.shape
        max_rna = int(np.max(self.data[int(specID),:,:])+1)
        self.means = np.zeros(n_times)
        for i in range(n_times):
            self.means[i] = np.sum(np.arange(max_rna)*self.pdf[i,:])

    def get_variances(self,specID=0):
        '''
        get the second moment.
        '''
        self.get_dist()
        n_specs, n_times, n_traj = self.data.shape
        max_rna = int(np.max(self.data[specID,:,:])+1)
        self.variances = np.zeros(n_times)
        self.covariances = np.zeros((n_specs,n_specs,n_times))
        for i in range(n_times):
            self.variances[i] = np.sum((np.arange(max_rna)**2)*self.pdf[i,:])-(np.sum(np.arange(max_rna)*self.pdf[i,:])**2)
            self.covariances[:,:,i] = np.cov(self.data[:,i,:])

    def return_all_var(self):
        all_members = self.__dict__.keys()

        return [ (item, self.__dict__[item]) for item in all_members if not item.startswith("_")]

    def return_names(self):
        all_members = self.__dict__.keys()
        return [ item for item in all_members if not item.startswith("_")]




def simulate_cell(diffusion_constant, kon, koff, kRNA, kdecay, ti=0, tf=1000, tstep=1000, cell_radius=50, imagesize=5, dpi=90, filename='simulated_test.gif', ssa_obj=None, fcolor='#00FF00', rnacolor='#FF0000'):
    '''
    [DNA] ==kRNA==> [RNA] <==koff== [RNA*] ==translation simulation==> [Protein]===> null
                    // ||             /\
                    || `'=====kon====='`
                    ||
                    \/
                    null

    '''
    print('simulating RNA creation....')
    t = np.linspace(ti, tf, tstep)


    dna_s = np.array([[ 0,  0],
                      [ 1, -1]])

    dna_w1 = np.array([[kRNA, 0],
                       [0, 0]],dtype=float)


    dna_w0 = np.array([[0], [0]])


    dna_si = GenericSSA(type='linear' )
    dna_si.W1 = dna_w1
    dna_si.W0 = dna_w0
    dna_si.S = dna_s

    dna_si.ti = t[0]
    dna_si.tf = t[-1]
    dna_si.n = 1
    xi = np.zeros((2, 1))
    xi[0] = 1
    dna_si.xi = xi
    dna_si.ptimes = len(t)

    dna_si.time_variant = False
    dna_si._solve(1)
    rna_creation_data = dna_si.data




    stoich = np.array([[  0,    0,  1],
                       [  -1,  1, -1],
                       [  1, -1, 0]])

    propensity = np.array([[0, kon, 0],
                          [0, 0,koff],
                          [0,kdecay, 0]], dtype=float)

    w0 = np.array([[0],[0],[0]])

    solver_instance = GenericSSA(type='linear' )
    solver_instance.W1 = propensity
    solver_instance.W0 = w0
    solver_instance.S = stoich

    solver_instance.ti = t[0]
    solver_instance.tf = t[-1]
    solver_instance.n = 1
    xi = np.zeros((3,1))
    xi[1] = 1
    solver_instance.xi = xi
    solver_instance.ptimes = len(t)

    solver_instance.time_variant = False




    print('simulating RNA activation....')



    R = cell_radius
    squarelen = float(R/np.sqrt(2))
    
    theta = np.linspace(0,2*np.pi,1000)
    rnuc = 1.*R-2.*R/3
    
    ncircx,ncircy = rnuc*np.cos(theta), rnuc*np.sin(theta)
    
    
    n_RNA_t = np.zeros((len(t),int(np.max(rna_creation_data[1]))))

    nRNA = 0
    nparticles = (int(np.max(rna_creation_data[1])))
    for i in range(len(t)):

        while nRNA != rna_creation_data[1][i]:
            data = solver_instance._solve(1)


            rnaonoff = data[2] + 1 - data[0]



            n_RNA_t[i:, nRNA] = rnaonoff[:-i].flatten()
            nRNA += 1


    rna_particles = n_RNA_t.T
    rna_exist = np.where(rna_particles >0,1,0)
    rnaex = data



    print('simulating RNA motion....')
    rna_locations = np.empty((nparticles, len(t), 2))

    dt = t[-1]/len(t)

    delta = diffusion_constant



    def linecirc(m, b, xc, yc, r):

        if np.isinf(m) == False:
            a = 1+m**2
            e = 2*(m*(b-yc)-xc)
            c = yc**2+xc**2 + b**2-2*yc*b-r**2
            x = np.roots([a, e, c])

            if np.isreal(x).all() == False:
                x = [np.nan, np.nan]
                y = [np.nan, np.nan]
            else:
                y = [b + m*x[0], b+m*x[1]]

        elif abs(xc-b) > r:
            x = [np.nan, np.nan]
        else:
            x = [b, b]
            step = np.sqrt(r**2-(b-xc)**2)
            y = [yc + step, yc-step]

        return [x[0], y[0]], [x[1], y[1]]


    def dist(x1, y1, x2, y2):
        return np.sqrt((x1-x2)**2+(y1-y2)**2)


    for i in range(nparticles):
        x = np.empty((2,len(t) - np.where(rna_exist[i] != 0 )[0][0]  ))
        centers = np.zeros(x.shape)
        #x[:,0] = ncircx[int(np.random.uniform(0,1000))]
        x0 = [ncircx[int(np.random.uniform(0,1000))],  ncircy[int(np.random.uniform(0,1000))]]
        #x0 = [  ((R+squarelen/4) - (R-squarelen/4))*np.random.random() + (R-squarelen/4),((R+squarelen/4) - (R-squarelen/4))*np.random.random() + (R-squarelen/4) ]


        #x0 = x0 - np.array([0, 0])
        
        x[:,0] =x0
        r = scipy.stats.norm.rvs(size=np.array(x0).shape + (len(t) - np.where(rna_exist[i] !=0 )[0][0],), scale=delta*np.sqrt(dt))



        out = np.empty(r.shape)

        np.cumsum(r, axis=-1, out=out)
        out += np.expand_dims(np.array(x0), axis=-1)

        #out = np.array([[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36],
                        #[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36]])

        centers = np.zeros(out.shape)
        dists = np.zeros((x.shape[1], 1)).flatten()

        incirc = np.hypot(out.T[:, 0]-centers.T[:, 0], out.T[:, 1]-centers.T[:, 1])
        dists[np.where(out[0] != 0)] = incirc[np.where(out[0] != 0)]

        while len(np.where(dists>R)[0]) != 0:   #trajectory left the cell
            out = out.T
            left_cell = np.where(dists > R)[0][0]



            pts = [[out[left_cell][0], out[left_cell][1]], [out[left_cell-1][0], out[left_cell-1][1]]]

            p = np.polyfit([out[left_cell][0], out[left_cell-1][0]], [out[left_cell][1], out[left_cell-1][1]], 1)
            m = p[0]
            b = p[1]

            intercepts = linecirc(m, b, 0, 0, R)
            if dist(*tuple(intercepts[0])+tuple(pts[0])) > dist(*tuple(intercepts[1])+tuple(pts[0])):
                inter = np.array(intercepts[1])
            else:
                inter = np.array(intercepts[0])

            a = out[left_cell] - inter


            out[left_cell-1:] = out[left_cell-1:] - 2*(np.dot(inter, a)/np.linalg.norm(inter)**2)*inter




            dists = np.zeros((x.shape[1], 1)).flatten()
            out = out.T
            incirc = np.hypot(out.T[:, 0]-centers.T[:, 0], out.T[:, 1]-centers.T[:, 1])
            dists[np.where(out[0] != 0)] = incirc[np.where(out[0] != 0)]



        data = ((out.T).T*rna_exist[i][np.where(rna_exist[i] != 0)[0][0]:].T).T
        data[np.where(rna_exist[i] != 0)[0][-1]- np.where(rna_exist[i] != 0)[0][0]+1 :] = -R-1


        rna_locations[i, np.where(rna_exist[i] != 0)[0][0]:, :] =  data
        rna_locations[i, :np.where(rna_exist[i] != 0)[0][0], :] =  -R-1


    print(nparticles)
    rna_loc_compressed = rna_locations[np.where(np.sum(np.sum(rna_locations+R, axis=1), axis=1) > 0)]

    if ssa_obj == None:
        print('no ssa data given')
        print('simulating translation....')
        print(int(rna_loc_compressed.shape[0]))
        ssa_obj = ssa_solver(n_traj=int(rna_loc_compressed.shape[0]),tf=tf,tstep=tstep)

        ivec = ssa_obj.intensity_vec/np.max(ssa_obj.intensity_vec)
        ivec = ivec.T  #get the intensity vec for the "fluorescence"


    else:
        print('Translation data given')
        print('Given ' + str(ssa_obj.n_traj) + ' Needed '+str(int(rna_loc_compressed.shape[0])) )
        if ssa_obj.n_traj  < int(rna_loc_compressed.shape[0]):
            print('simulating ' + str(int(rna_loc_compressed.shape[0]) - ssa_obj.n_traj) + ' additional trajectories....')
            ssa_obj =ssa_solver_append(ssa_obj, n=int(rna_loc_compressed.shape[0]) - ssa_obj.n_traj)
            ivec = ssa_obj.intensity_vec/np.max(ssa_obj.intensity_vec)
            ivec = ivec.T  #get the intensity vec for the "fluorescence"

        else:
            ivec = ssa_obj.intensity_vec[0:int(rna_loc_compressed.shape[0])]/np.max(ssa_obj.intensity_vec[0:int(rna_loc_compressed.shape[0])])
            ivec = ivec[0:int(rna_loc_compressed.shape[0])].T  #get the intensity vec for the "fluorescence"



            
    def customcmap(color1,color2):
    
        ccmap = LinearSegmentedColormap.from_list('custom', [(0,color1),(.75,color1),(1,color2)], N=256)
        
        
        return ccmap

    
    import time
    print(time.time())
    print('making movie...')
    #simulate brownian motion
    
    def update_line(num, xpos, ypos, line):
        if num !=0:
            
            for child in ax.get_children():  #remove the previous patch collection (green spots)
                
                if isinstance(child, mpl.image.AxesImage):
                    
                    child.remove()
        if num == 1:
            print(time.time() - st )*len(t)
            
       
        radi = 5*ivec[inds[num]]+.0001   
        xmin, xmax, ymin, ymax = (-R-10, R+10, -R-10, R+10)
        n_bins = 300
        xx = np.linspace(xmin, xmax, n_bins)
        yy = np.linspace(ymin, ymax, n_bins)
        ncount = 0
        cmap = plt.cm.viridis
        cmap = customcmap('#b619ff','#fbff28')
        cmap2 = plt.cm.OrRd
        
        for x1, y1, r in zip(xpos[inds[num]],ypos[inds[num]], radi):   #make circle objects of radius based on ivec
 
            

            
            '''
            if ncount ==0:
                weightsx = normal_pdf(xx, x1, r)
                weightsy = normal_pdf(yy, y1, r)                
                weights = np.array(np.meshgrid(weightsx, weightsy)).prod(0)
                greys = np.empty(weights.shape + (3,), dtype=np.uint8)
                greys.fill(0)  
                ncount+=1 
            '''
            
                
            if r > .001:
                
                if x1 !=-R-1 and y1 !=-R-1:
                    
                    weightsx = normal_pdf(xx, x1, r)
                    weightsy = normal_pdf(yy, -y1, r)
                    
                    weights = np.array(np.meshgrid(weightsx, weightsy)).prod(0)
                 
        
        
                
                    vmax = np.abs(weights).max()
                    vmin = -vmax
                    
                    
                    alphas = Normalize(0, 1, clip=True)(np.abs(weights))
                    
                    colors = Normalize(vmin, vmax)(weights)
                    colors = cmap(colors)
                    
                    colors[..., -1] = alphas
        
                
                
                    ax.imshow(colors, extent=(xmin, xmax, ymin, ymax),zorder=4) 
                    
                    
            weightsx = normal_pdf(xx, x1, .5)
            weightsy = normal_pdf(yy, -y1, .5)
            
            weights = np.array(np.meshgrid(weightsx, weightsy)).prod(0)
         


        
            vmax = np.abs(weights).max()
            vmin = -vmax
            
            
            alphas = Normalize(0, 1, clip=True)(np.abs(weights))
            
            colors = Normalize(vmin, vmax)(weights)
            colors = cmap2(colors)
            
            colors[..., -1] = alphas
    
        
        
            ax.imshow(colors, extent=(xmin, xmax, ymin, ymax),zorder=3)   
        
        
        theta = np.linspace(0,2*np.pi,1000)
        rnuc = R
    
        ncircx,ncircy = rnuc*np.cos(theta), rnuc*np.sin(theta)
        
        n_bins = 140
        xx = np.linspace(xmin, xmax, n_bins)
        yy = np.linspace(ymin, ymax, n_bins)        
        
        
        weights = np.array(np.meshgrid(xx, yy)).prod(0)
        greys = np.empty(weights.shape + (4,), dtype=np.float32)
        greys.fill(70)  
         
        
        greys =  (1.5*np.random.poisson(greys)).astype(int)
        
        
        alphas = Normalize(0, 1, clip=True)(np.abs(greys))
        alphas[:,:,:] = .5
        
        
        greys = Normalize(0, 256)(greys)
        greys[..., -1] = .4

        
        ax.imshow(greys, extent=(xmin, xmax, ymin, ymax),zorder=5)  
        
        
        '''
        line.set_data(xpos[inds[num]], ypos[inds[num]])
        line.set_linewidth(0)
        line.set_marker('o')
        line.set_markersize(1)
        line.set_color(rnacolor)
        line.set      
        '''
        
        
               
            
        
        plt.xlabel(str(inds[num])) 
        
        
        return line,
    
    '''
    def update_line(num, xpos,ypos, line):  #function for the FuncAnimation
        if num !=0:

            for child in ax.get_children():  #remove the previous patch collection (green spots)

                if isinstance(child, PatchCollection):
                    child.remove()
                if isinstance(child, mpatches.Ellipse):
                    child.remove()

        patches = []
        radi = 3*ivec[inds[num]]   #create a max radius of 3 for intensity vecs


        for x1, y1, r in zip(xpos[inds[num]],ypos[inds[num]], radi):   #make circle objects of radius based on ivec
            circle = mpatches.Circle((x1, y1), r, color=fcolor)
            patches.append(circle)
            #fig.gca().add_artist(circle)

        line.set_data(xpos[inds[num]], ypos[inds[num]])
        line.set_linewidth(0)
        line.set_marker('o')
        line.set_markersize(3)
        line.set_color(rnacolor)
        line.set
        p = PatchCollection(patches, zorder=3, facecolors=(fcolor,))  #create a patch collection to add to axis
        ax.add_collection(p)  #adds the circles to axis


        p = mpatches.Circle((0,0), radius=R, color='black')  #add the black circle
        ax.add_patch(p)


        whitep = mpatches.Ellipse((-R, -R), width=7, height=7, color='white', zorder=5)  #add the black circle
        ax.add_patch(whitep)

        plt.xlabel(str(inds[num]))  #update time label


        return line,
    
    '''
    
    xpos = rna_loc_compressed.T[0]
    ypos = rna_loc_compressed.T[1]

    for i in range(int(rna_loc_compressed.shape[0])):
        start = np.where(xpos[:,i] !=-R-1)[0][0]

     
        ivec[:,i][start:] = ivec[:,i][:ivec[:,i].shape[0]-start]
        ivec[:,i][:start] = 0
            
    
    filetype='.gif'

    if filetype == '.gif':

        Writer = animation.writers['pillow']
    if filetype == '.html':
        Writer = animation.writers['html']
    if filetype == '.mov':

        Writer = animation.writers['FFMpeg']

    writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

    fig1 = plt.figure(figsize=(imagesize, imagesize),dpi=dpi)  #make figure
    fig1.tight_layout()

    ax= fig1.add_subplot('111')
    plt.yticks([])
    plt.xticks([])
    p = mpatches.Rectangle((-R-10,-R-10),width=(R+10)*2,height=(R+10)*2,color='black',zorder=1)
    ax.add_patch(p)
    p = mpatches.Circle((0, 0), radius=R+1, color='white',zorder=1)
    ax.add_patch(p)
    p = mpatches.Circle((0, 0), radius=R, color='black',zorder=1)  #add the black circle
    ax.add_patch(p)
    
    
    
    p = mpatches.Circle((0, 0), radius=R*1.-2.*R/3, color='#014abf',zorder=1)  #add the black circle
    ax.add_patch(p)
    plt.gca().set_aspect('equal', adjustable='box')
    
    
    xmin, xmax, ymin, ymax = (-R-10, R+10, -R-10, R+10)
    n_bins = 300
    xx = np.linspace(xmin, xmax, n_bins)
    yy = np.linspace(ymin, ymax, n_bins)
    weightsx = normal_pdf(xx, 0, 20)
    weightsy = normal_pdf(yy, 0, 20)
    
    weights = np.array(np.meshgrid(weightsx, weightsy)).prod(0)
 



    vmax = np.abs(weights).max()
    vmin = -vmax
    
    
    alphas = Normalize(0, 1, clip=True)(np.abs(weights))
    
    colors = Normalize(vmin, vmax)(weights)
    cmap = plt.cm.YlGnBu
    colors = cmap(colors)
    
    colors[..., -1] = alphas



    ax.imshow(colors, extent=(xmin, xmax, ymin, ymax),zorder=2)   
            
            
    p = mpatches.Circle((-R-1, -R-1), radius=8, color='black',zorder=5)  #add the black circle
    ax.add_patch(p)
    plt.gca().set_aspect('equal', adjustable='box')
    
    
    
    

    l, = plt.plot([], [], 'r-',zorder=2)
    plt.xlim(-R-10, R+10)
    plt.ylim(-R-10, R+10)
    plt.xlabel('0')
    plt.title('Simulated Cell')

    inds = np.linspace(0, len(t)-1, len(t)).astype(int)
    st = time.time()
    itercount = 0
    #creates the animation
    
    '''
    for num in range(0,100):    
        radi = 3*ivec[inds[num]]   +.0001
        xmin, xmax, ymin, ymax = (R-10, R+10, R+10, R-10)
        n_bins = 100
    
        for x1, y1, r in zip(xpos[inds[num]],ypos[inds[num]], radi):   #make circle objects of radius based on ivec
    
            
            xx = np.linspace(xmin, xmax, n_bins)
            yy = np.linspace(ymin, ymax, n_bins)
            
            weightsx = normal_pdf(xx, x1, r)
            
            weightsy = normal_pdf(yy, y1, r)
            
            weights = np.array(np.meshgrid(weightsx, weightsy)).prod(0)
            if i == 0:
                greys = np.empty(weights.shape + (3,), dtype=np.uint8)
                greys.fill(0)  
                i+=1
                ax.imshow(greys)
        
            vmax = np.abs(weights).max()
            vmin = -vmax
            cmap = plt.cm.viridis
            
            alphas = Normalize(0, 1, clip=True)(np.abs(weights))
            
            colors = Normalize(vmin, vmax)(weights)
            colors = cmap(colors)
            
            # Now set the alpha channel to the one we created above
            colors[..., -1] = alphas
    
        # Note that the absolute values may be slightly different
        
       
            
            # Create the figure and image
            # Note that the absolute values may be slightly different
        

            ax.imshow(colors, extent=(xmin, xmax, ymin, ymax))     
            
            plt.xlabel(str(inds[num]))     
    '''
    
    line_ani = animation.FuncAnimation(fig1, update_line, tstep, fargs=(xpos,ypos, l),
                                       interval=10, blit=True,repeat=False)

    line_ani.save((filename + filetype), writer=writer)  #save the animation


    print(time.time())
    #return solver_instance,n_RNA_t,rna_creation_data,data,rna_locations
    return rna_locations, rna_loc_compressed, rna_particles, rna_creation_data, rna_exist, rnaonoff, rnaex,ivec





Writer = animation.writers['pillow']

sms = rSNAPsim.rSNAPsim()
sms.open_seq_file("gene_files/KDM5B.txt")
sms.get_orfs(sms.sequence_str, min_codons = 80)
sms.get_temporal_proteins()
sms.analyze_poi(sms.pois[0],sms.pois_seq[0])

#ssa_obj = sms.ssa_solver(n_traj=300)
ssa_obj = rSNAPsim.ssa()
ssa_obj.load_from_txt('simcelltestobj.txt')

rna_locations, rna_loc_compressed, rna_particles, rna_creation_data, rna_exist, rnaonoff, rnaex,ivec = simulate_cell(.3,3,1,.03,.01,ssa_obj=ssa_obj,tf=500,tstep=500)




def plot_gauss(centers,rs,weight):

    xmin, xmax, ymin, ymax = (0, 100, 0, 100)
    n_bins = 100
    xx = np.linspace(xmin, xmax, n_bins)
    yy = np.linspace(ymin, ymax, n_bins)
    fig, ax = plt.subplots()
  
        
    i = 0
    for center in centers:
        weightsx = normal_pdf(xx, center[0], 10)
        weightsy = normal_pdf(yy, ymax - center[1], 10)
        
        weights = np.array(np.meshgrid(weightsx, weightsy)).prod(0)
        if i == 0:
            greys = np.empty(weights.shape + (3,), dtype=np.uint8)
            greys.fill(0)  
            i+=1
            ax.imshow(greys)
    
        vmax = np.abs(weights).max()
        vmin = -vmax
        cmap = plt.cm.viridis
        
        alphas = Normalize(0, 1, clip=True)(np.abs(weights))
        
        colors = Normalize(vmin, vmax)(weights)
        colors = cmap(colors)
        
        # Now set the alpha channel to the one we created above
        colors[..., -1] = alphas
        
        # Create the figure and image
        # Note that the absolute values may be slightly different
        
       
        ax.imshow(colors, extent=(xmin, xmax, ymin, ymax))
        
        
