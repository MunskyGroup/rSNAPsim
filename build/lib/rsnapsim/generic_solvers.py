# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:52:28 2020

@author: willi
"""
import numpy as np
from scipy.sparse.linalg import expm
from scipy.integrate import odeint
#from scipy.integrate import ode
import scipy as sp
import scipy.optimize
#import random as random
import collections as c
#import matplotlib.pyplot as plt
import copy as cp
#from multiprocessing import process as Process
import matplotlib.pyplot as plt

'''
Generic solvers written by Micheal May, Zach Fox and Will Raymond ~2018
'''


class GenericSSA:
    def __init__(self,type='linear'):
        self.xi=np.array([])
        self.ti=None
        self.tf=None
        self.S=np.array([])
        self.type=type
        self.ptimes=100
        self.params={}
        if type=='linear':
            self.W0=np.array([])   
            self.W1=np.array([])
        if type == 'nonlinear':
            self.P=lambda x,t:None 
    def gettvec(self):
        return np.linspace(self.ti,self.tf,self.ptimes)
        
    def updateparams(self):
        for k, v in self.params.iteritems():
            setattr(self, k, v)

    def _run_trajectory(self):#!!!!!!!!!!!!!!!!!renamed run to solve(big deal)
        x=self.xi
        t=0
        __n=len(x)
        self.time=self.gettvec()
        data=np.zeros((len(self.xi),self.ptimes))
        ip=0
        if self.type=='linear':
            while t<self.tf:
                rate=np.atleast_2d(np.dot(self.W1,x))+self.W0
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
                        data[:,ip]=x.reshape(__n)
                        ip=ip+1
                for i in range(len(rate)):
                    if rate[i]>=ro:
                        event=i
                        break
                x=x+np.atleast_2d(self.S[:,event]).T
        elif self.type=='nonlinear':
            pass
        else:
            'Error'
        return data

    def _solve(self,n):
        soln=np.zeros((len(self.xi),self.ptimes,n))
        for i in range(n):
            d=self._run_trajectory()
            soln[:,:,i]=d
        self.soln = soln
        return soln
        
    def getW0(self):
        pass
    
    def getW1(self):
        pass
    
    def setpar(self,key,val):
        self.params[key]=val

class GenericFSP:
    def __init__(self,ti=[],tf = [],xi=[],A=np.array([]),ptimes=100):
        self.ti=ti
        self.tf=tf
        self.xi=xi
        self.A=A
        self.errors=np.empty(ptimes)
        self.ptimes=ptimes
        self.params={}
        
    def updateparams(self):
        for k, v in self.params.iteritems():
            setattr(self, k, v)
        
    def gettvec(self):
        '''
        get vector of times 
        '''
        return np.linspace(0,self.tf-self.ti,self.ptimes)

    def _solve(self):
        '''
        solve the FSP at each time. 
        '''
        tvec=self.gettvec()
        self.model_solns = np.zeros((len(self.xi),len(tvec)))
        for i in range(len(tvec)):
            mexp = expm(self.A*tvec[i])
            self.soln = mexp.dot(self.pi)
#           self.errors[i]=1.0-np.sum(self.model_solns[:,i])
        return self.soln

class GenericODE:
    def __init__(self,ti=None,tf=None,ODE=None,ptimes=50,xi=None):
        self.ti=ti
        self.tf=tf
        self.ODE=ODE
        self.ptimes=ptimes
        self.xi=xi
        self.tvec=None
        self.soln=None
        self.parameters={}
    
    def updateparams(self):
        try:
            for k, v in self.params.items():
                setattr(self, k, v)
        except:
             for k, v in self.params.iteritems():
                setattr(self, k, v)           
            
    def gettvec(self):
        '''
        get vector of times
        '''
        return np.linspace(self.ti,self.tf,self.ptimes)

    def _solve(self):
        '''
        solve the ODEs at each time
        '''
        self.tvec=self.gettvec()
        try: 
            solution=odeint(self.ODE,self.xi,self.tvec,Dfun=self.Jacobian)
        except: 
            solution=odeint(self.ODE,self.xi,self.tvec)
        self.soln = solution.T
        return solution.T
        
    def _expm_solve(self):
        self.tvec = self.gettvec()
        
        
        
class GenericStats():
    '''
    a class of static functions to analyze data
    '''
    def __init__(self,model=None):
        if model is not None:
            self.model=model
            self.type=model.__class__.__name__
            self.simdata={}
            self.simdata['data']=self.model.solve()
            self.simdata['time']=self.model.gettvec()
            
            self.means=self.get_means(**self.simdata)
            self.acc=self.get_acc(self.simdata)
            if self.__class__.__bases__ is 'GenericSSA' or self.__class__.__name__ is 'GenericSSA':
                self.var=self.get_var(self.simdata['data'])

    @staticmethod
    def get_var(data,time):
        '''
        a function that returns statistical variance of data
        '''
        if type(data) is not np.ndarray :
            print("Error: data is not a numpy array")
            return
        s=data.shape
        if len(s) is 1:
            s=s+1
        var=np.zeros((s[0],s[1]))
        for j in range(len(time)):
            for i in range(s[0]):
                var[i,j]=np.var(data[i,j,:])
        return var
        
    @staticmethod
    def get_means(data,time):
        '''
        a function that returns the statistical mean of data
        '''
        if type(data) is not np.ndarray :
            print("Error: data is not a numpy array")
            return
        __s=data.shape
        means=np.zeros((__s[0],__s[1]))
        for j in range(len(time)):
            for i in range(__s[0]):
                means[i,j]=np.mean(data[i,j,:])
        return means
    
    @staticmethod 
    def get_acc2(data,trunc=True): 
        N = len(data)    
        fvi = np.fft.fft(data,n=2*N)     
        acf = fvi*np.conjugate(fvi)
        acf = np.fft.ifft(acf)
        acf = np.real(acf[:N])/float(N)
        if trunc:
            acf[acf<0]=0
            for i in range(1,len(acf)): 
                if acf[i] > acf[i-1]:
                    acf[i] = acf[i-1]
        return acf
        
    @staticmethod 
    def get_fastacc2(data,trunc=True): 
        N = len(data)    
        fvi = np.fft.fft(data,n=2*N)     
        acf = fvi*np.conjugate(fvi)
        acf = np.fft.ifft(acf)
        acf = np.real(acf[:N])/float(N)
        if trunc:
            acf[acf<0]=0
            for i in range(1,len(acf)): 
                if acf[i] > acf[i-1]:
                    acf[i] = acf[i-1]
        return acf
    @staticmethod
    def get_acc(data,one_sided=True,trunc=True):
        '''
        a function that return the statistical autocorreletion of data
        '''
        half = len(data)/2
        z = np.correlate(data[:half+1],data,mode='valid') / float(half)
        if trunc:
            z[z<0]=0
            for i in range(1,len(z)): 
                if z[i] > z[i-1]:
                    z[i] = z[i-1]
        return z 

    @staticmethod
    def get_xcc(data1,data2=None,axis=1,norm=False,trunc=False):
        '''
        a function that returns the statistical crosscorrelation(or autocorrelation)of data
        '''
        if data2 is None:
            data2=data1
        if type(data1) is not np.ndarray:
            return "Type error, data must be a numpy array"
        corr = np.correlate(data1,data2,mode='full')
        if norm and corr.max()>0:
            corr = corr/corr.max()
        lc = len(corr)/2
        if trunc:
            corr[corr<0]=0
            tmp = corr[lc:]
            for i in range(1,lc): 
                if tmp[i] > tmp[i-1]:
                    tmp[i] = tmp[i-1]
            corr[lc:] = tmp
        return corr 

    @staticmethod
    def get_kld(p,q):
        '''
        returns the kullback-leibler divergence between the distributions 
        p and q. 
        '''
        support = (p>0)*(q>0)
        #if np.max(self.target)<np.min(self.predicted) or np.min(self.target)>np.max(self.predicted):
        #    kld = np.inf
        kld = np.sum( p[support] * (np.log(p[support]) - np.log(q[support])))
        return kld    
    
    @staticmethod
    def normalize(data):
        return data/sum(data)
        
    @staticmethod
    def get_curvature(y,h):
        k=np.zeros(len(y)-2)
        for i in range(len(y)-2):
            k[i]=((y[i+2]-2*y[i+1]+y[i])/(h**2.0))/(1+((y[i+2]-y[i])/(2*h))**2)**(1.5)
        return k
        
    @staticmethod
    def get_sse(data1,data2,axis=None):
        '''calculates sum of squared error of two numpy array datasets'''
        return np.linalg.norm(data1-data2)
        
class GenericOpt:
    def __init__(self,model):
        self.model=model
        self.modeltype=model.__class__.__name__
        self.activeobjfun='chisq'
        self.activemethod='MH'
        self.activefun='solve'
        self._mdt={'met_hast':self.methast,\
                      'methast':self.methast,\
                      'MH':self.methast,\
                      'lin_opt':self.lin_opt,\
                      'linopt':self.lin_opt,\
                      'genetic':self.genetic,\
                      'sim_anneal':self.sim_anneal,\
                      'simanneal':self.sim_anneal}
        self._odt={'sse':self.sse,\
                   'SSE':self.sse,\
                   'chi':self.chisq,\
                   'chisq':self.chisq}
        self.opts={'bounds':([0,10],)*len(self.model.params)}
        self.data=[]
        self.time=[]
        self.cache=c.deque([])
        self.bestpars=[]
    def fit(self,*args):
        #if method is not None:
        #    return self.methods[method](self,args)
        if self.activemethod is not None and self.activeobjfun is not None:
            print("Running: " +self.activemethod +" with "+self.activeobjfun)
            return self._mdt[self.activemethod](*args)
        else:
            print("Please pick an active solver and objective function")
            
    def getsimdata(self):
        return self.model.solve()
    
    def sse(self):
        simdata=self.getsimdata()
        return np.sum((self.data-simdata)**2)

    def chisq(self):
        simdata=self.getsimdata()
        return np.sum(((self.data-simdata)**2)/self.model.variance)
    
    def kld(self,x):
        pass#return np.stats.kld(1,2,3)
    def methast(self,N=1000,burn=100,stepsize=5):
        def evolvepars(p,stepsize):
            new=cp.copy(p)
            for j in range(len(new)):
                while 1 :
                    new[j]=p[j]+np.random.normal(0,stepsize)
                    #print p
                    #print new
                    if (self.opts['bounds'][j][0]<new[j]<self.opts['bounds'][j][1]):
                        break
            return new
                        
        optfun=self.getoptfun()
        oldpars=self.bestpars=self.getpars()
        f_old=f_best=optfun(oldpars)   
        for i in range(-burn,N):
            newpars=evolvepars(oldpars,stepsize)
            f_new=optfun(newpars)
            #print("newpars = "+str(newpars)+ "     oldpars = "+str(oldpars))
            #print("fnew="+str(f_new)+", fold="+str(f_old))
            if f_new<f_old or np.random.rand()<np.exp(f_old-f_new):
                #raw_input("Press Enter to continue...")
                f_old=f_new
                oldpars=newpars
                print("Accpeted "+str(i))
                if i>0:
                    self.cacheitem(newpars)    
                if f_new<f_best:
                    f_best=f_new
                    self.bestpars=newpars
                

        return self.bestpars
    
    def genetic(self):
        f=self.getoptfun()
        logfun = lambda x: f(x)
        return sp.optimize.differential_evolution(logfun,**self.opts)
    def lin_opt(self):
        pass
    def sim_anneal(self):
        return sp.optimize.minimize(self.getoptfun(),self.getpars(),method='Anneal',options={'lower':[0,0]})
    
    def getoptfun(self):
        def fun(pars):
            self.setpars(pars)
            #print self.model.params
            return self._odt[self.activeobjfun]()
        return fun
        
    def getpars(self):
        i=0
        pars=np.zeros(len(self.model.params))
        for key,val in self.model.params.iteritems():
                pars[i]=val
                i+=1
        return pars
    def clearcache(self):
        self.cache=c.deque([])
    def cachepars(self):
        self.cache.append((self.getpars()))
    def cacheitem(self,item):
        self.cache.append(item)
    def setpars(self,parval):
        i=0
        for key,val in self.model.params.iteritems():
                self.model.params[key]=parval[i]
                i+=1
    def parplot(self,s=2,alpha=.5):
        N=len(self.model.params)
        D=zip(*self.cache)
        #print D

        for i in range(N):
            for j in range(N):
                plt.subplot(N,N,j+N*i+1)
                if i is not j:                    
                    plt.scatter(D[i],D[j],linewidth=0.0,s=s,alpha=alpha)
                else:
                    plt.hist(D[i])
        
class GenericModel(GenericSSA,GenericODE):
    '''
    define a generic model to run SSA,ODE~may add FSP analysis here one day
    '''
    def __init__(self,ti=None,tf=None,ODE=None,ptimes=50,xi=None):
        self.ti=ti
        self.tf=tf
        self.ODE=ODE
        self.ptimes=ptimes
        self.xi=xi
        self.tvec=None
        self.soln=None
        self.parameters={}
        self.type=None
    def updateparams(self):
        for k, v in self.params.iteritems():
            setattr(self, k, v)
            
    def gettvec(self):
        '''
        get vector of times
        '''
        return np.linspace(self.ti,self.tf,self.ptimes)

    def _solve(self):
        '''
        solve the ODEs at each time
        '''
        self.tvec=self.gettvec()
        try: 
            solution=odeint(self.ODE,self.xi,self.tvec,Dfun=self.Jacobian)
        except: 
            solution=odeint(self.ODE,self.xi,self.tvec)
        self.soln = solution.T
        return solution.T
    
