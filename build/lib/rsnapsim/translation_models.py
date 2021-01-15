import numpy as np
#import stochasticSimulationAlgorithm as ssa
import time
from scipy.sparse import spdiags
from scipy.linalg import solve_lyapunov
#import generic_solvers
#reload(generic_solvers)
from .generic_solvers import *


import scipy as sp
from . import expv

from scipy.sparse.linalg import expm

class TranslateModels:
    def getW1(self):
        '''
        build a propensity vector that works with my SSA. 
        '''
        W1 = np.diag(self.ke)
        f_rxn = np.zeros(self.N)
        return np.vstack((f_rxn,W1))

    def getW0(self):
        W0=np.atleast_2d(np.zeros(self.N+1)).T
        W0[0,0]=self.kb
        return W0

    def getselongation(self):
        __stoich_matrix_elongation = np.diag(np.ones(self.N),k=0)+np.diag(-1*np.ones(self.N-1),k=1)
        __f_rxn = np.atleast_2d(np.zeros(self.N))
        __f_rxn[-1,-1] = -1
        return np.hstack((__stoich_matrix_elongation,__f_rxn.T))

    def map_to_fluorescence(self,soln=None):
        self.probe_design_vector = self.get_fluor_vec()
     
        if soln is None:
            self.fluorescence = np.dot(self.probe_design_vector,self.soln)
        else:
            print(soln.shape)

            self.fluorescence = np.dot(np.dot(self.probe_design_vector,soln),self.probe_design_vector.T)
        return self.fluorescence 

    def get_fluor_vec(self,binary=None):
        if binary is not None:
            self.binary = binary
        if self.binary is None:
            n=len(self.getxi())
            return np.array(range(n))/float(n)
        else:
            return np.array([np.cumsum(self.binary,axis=1).dot(self.fi)])[0]
    
    def map_to_fluorescence2(self,color_1_bin,color_2_bin,soln):
        self.probe_design_vector_1 = self.get_fluor_vec(color_1_bin)
        self.probe_design_vector_2 = self.get_fluor_vec(color_2_bin)
        
        self.probe_design_vector = np.vstack((self.probe_design_vector_1,self.probe_design_vector_2))
        self.fluorescence = np.dot(np.dot(self.probe_design_vector,soln),self.probe_design_vector.T)
        
        return self.fluorescence   

    def map_to_fluorescence3(self,soln=None):
        '''
        Get the fluorescence at steady state
        '''
        self.probe_design_vector =  self.get_fluor_vec()
        print(self.probe_design_vector)
        #if self.soln is None:
           # print("no solution saved to memory, run self.solve() wtf")
        if soln is None:
            self.fluorescence = np.dot(self.probe_design_vector,self.soln)
        else:
            self.fluorescence = np.dot(self.probe_design_vector,soln)
        return self.fluorescence 

class TranslateODE(TranslateModels,GenericODE):
    def __init__(self,ptimes=None,tf=None,ke=None,kb=None,fi=None,N=None):
        GenericODE.__init__(self)
        if N: 
            self.ptimes=ptimes
            self.ti = 0 
            self.tf = tf
            self.params={'ke':ke,\
                         'kb':kb,\
                         'fi':fi}
            self.N=N
        else:
            self.ptimes=100
            self.ti = 0; self.tf =10.0;
            self.params={'ke':5.0,\
                         'kb':5.0,\
                         'fi':0.1}
            self.N=100
        #self.W1 = self.getW1()
        #self.W0 = self.getW0()   
        #self.S=self.getselongation()
        #self.ODE=lambda x,t:np.ravel(np.dot(self.S,(np.dot(self.W1,np.atleast_2d(x).T)))+np.dot(self.S,self.W0))
        #self.build_probe()
        #self.summary = {}
        #self.soln=None
        #self.xi=self.getxi()
    
    def solve(self):
        '''
        wrap generic run to reinitialize variables
        '''
        self.xi=np.zeros(self.N)
        self.W1=self.getW1()
        self.W0=self.getW0()
        self.S=self.getselongation()
        
        if self.N > 100:
            self.soln = self.expv_solve_means()
        
        #self.build_probe()
        else:
            self.ODE=lambda x,t:np.ravel(np.dot(np.dot(self.S,self.W1),np.atleast_2d(x).T)+np.dot(self.S,self.W0))
            solution=self._solve()
            self.soln=solution
            
        return self.map_to_fluorescence()
    
    
    def expv_solve_means(self):
        '''
        expv solver for systems that are greater than 100 states 
        
        This is specific for the Translational models as it can be written as A*X
        '''
        
        # reshape and add Kbind to the propensity and stoich, for x' = A*x
        
        S_1 = np.zeros((self.S.shape[0]+1,self.S.shape[1]))
        S_1[1:,:] = self.S
        W1_1 = np.zeros((self.W1.shape[0],self.W1.shape[1] + 1))
        
        W1_1[:,1:] = self.W1
        W1_1[0,0] = self.kb
        W1_1[0,1] = 0

        phi1 = sp.sparse.csc_matrix(np.dot(S_1,W1_1)) 
                
        x0 = np.zeros((self.N + 1,1))
        x0[0] = 1
        t = self.gettvec()
        xv =np.copy(x0) 
        #intensity = np.zeros(len(t))
        soln = np.zeros((self.N+1,len(t) ))
        
        
        
        for i in range(1,len(t)):
    
            xv,m,v = expv(t[i]-t[i-1],phi1,xv,tol = 1e-6,m=30) 

            soln[:,i] = xv
            #intensity[i] = np.dot(bpv,soln[1:,i])
        print(soln.shape)
        return soln[1:,:]
                    

    def build_probe(self,__placements=None):
        '''
        build a default probe on the beginning and end of RNA
        '''
        self.updateparams()
        if __placements is None:
            __placements = np.zeros(self.N)
            __placements[3:10] = np.ones(7)
            #__placements[-10:-3] = np.ones(7)
            self.probe_design_vector = np.cumsum(__placements)
            #self.probe_design_vector = __placements
        else:
            self.probe_design_vector = np.cumsum(__placements)
        
    def build_W1(self):
        self.updateparams()
        __main_diag = [-1]*(self.N)
        W1=np.multiply(spdiags(np.array([__main_diag,[1]*self.N]),np.array([0,-1]),self.N,self.N),self.ke)
        return W1.toarray()

    def build_W0(self):
        self.updateparams()
        W0= np.atleast_2d(np.zeros(self.N)).T
        W0[0]=self.kb
        return W0

    def getfc(self):
        return np.dot(range(self.N),self.fi)
    
    # def map_to_fluorescence(self,soln=None):
    #     try:
    #         self.soln.shape
    #     except:
    #         print("no solution saved to memory, run self.solve()")
    #     if soln is None:
    #         self.fluorescence = np.dot(self.probe_design_vector,self.soln)
    #     else:
    #         self.fluorescence = np.dot(self.probe_design_vector,soln)
    #     return self.fluorescence 
        
    def get_flu_acc(self):
        fdata=self.map_to_fluorescence()
        return self.get_xcc(fdata)

    def getxi(self):
        self.updateparams()
        return np.zeros(self.N)

class TranslateVars(TranslateODE):    
    def __init__(self):
        
        GenericODE.__init__(self)
        TranslateODE.__init__(self)
        self.ptimes=100
        self.ti = 0; self.tf =10.0;

        self.N=100        
        self.kb=.33
        self.ke = np.ones(100)
        self.kb = .033
        self.fi = .1
        

        self.S = self.getselongation()
        self.W0 = self.getW0()
        self.W1 = self.getW1()
        self.ODE = self.variance_ODE
        self.xi=np.zeros(self.N+self.N**2)
    
    def vsolve(self,steady_state=False):  
        '''
        wrap solve function to have changing parameters
        if steady_state is true, solve the covariances by solving 
        the Lyapunov function. 
        '''
        self.S = self.getselongation()
        self.W0 = self.getW0()
        self.W1 = self.getW1()
        self.ODE = self.variance_ODE
        self.xi=np.zeros(self.N+self.N**2)
        #self.build_probe()
        self._solve()
        
    def variance_ODE(self,x,t):
        '''
        build variance/covariance ODEs.
        '''
        # dSIG = S*W1*SIG + SIG*W1'*S' + S*diag((W1*MU + W0)')*S'         
        __MU = x[:self.N]
        __SIG = np.reshape(x[self.N:],(self.N,self.N))
        # Build diagonal matrix
        inner = np.dot(self.W1,np.array([__MU]).T)+self.W0
        self.diag_mat = np.diag(np.ravel(inner),k=0)
        # Compute RHS
        RHS_vars = np.ravel(np.dot(np.dot(self.S,self.W1),__SIG) + np.dot(np.dot(__SIG,self.W1.T),self.S.T) + np.dot(np.dot(self.S,self.diag_mat),self.S.T))
        RHS_means = np.ravel(np.dot(self.S,(np.dot(self.W1,np.atleast_2d(__MU).T)))+np.dot(self.S,self.W0))
        return np.concatenate((RHS_means,RHS_vars))    

class TranslateCorrs(TranslateVars):
    '''
    get auto and cross correlations
    '''    
    def __init__(self):
        GenericODE.__init__(self)
     
        #### Default Parameters ####
#        self.pnum=100
#        self.ti = 0; self.tf =10.0;
#        self.N = 100; self.ke = 5  
#        self.kb = 5
#        self.params = {}
#        self.params['ke'] = self.ke 
#        self.params['kb'] = self.kb
#        self.binary = None
#        # Get fundamental matrices.        
#        self.S = self.getselongation()
#        self.W0 = self.getW0()
#        self.W1 = self.getW1()
#        self.mu_ss = self.get_mean_SS()
#        self.var_ss = self.get_var_SS()
#        
#        # solve for covariance steady state using ode
#        self.xi= np.ravel(self.var_ss)
#        self.ODE = self.correlation_ODE

    def get_mean_SS(self):       
        '''
        get the steady state mean values for the system
        '''
        M1 = -1*np.linalg.inv(np.dot(self.S,self.W1))
        return np.dot(M1,np.dot(self.S,self.W0))
   
    def get_var_SS(self):
        '''
        get the steady state variance for the system
        by solving the lyapunov fn. 
        '''
        A = np.dot(self.S,self.W1)
        inner = np.dot(self.W1,self.mu_ss)+self.W0
        diag_mat = np.diag(np.ravel(inner),k=0)
        Q = np.dot(np.dot(self.S,diag_mat),self.S.T)
        return solve_lyapunov(A,-Q)
         
    def csolve(self):
        '''
        update parameters and solve correlations
        '''
        self.get_autonomous_matrix()
        # Get steady state covariance as IC. 
        #self.vsolve() 
        #self.xi= np.ravel(self.soln[self.N:,-1])
        self.mu_ss = self.get_mean_SS()
        self.var_ss = self.get_var_SS()
        #self.vsolve()
        self.xi= np.ravel(self.var_ss)
        self.ODE = self.correlation_ODE
#        self._solve()
        self.faster_solve()

    def get_autonomous_matrix(self):
        '''
        Gets the phi matrix
        '''
        # Update fundamental matrices if parameters have been changed
        self.S = self.getselongation()
        self.W0 = self.getW0()
        self.W1 = self.getW1()
        self.phi = np.dot(self.W1.T,self.S.T) 
        self.phi = -np.diag(self.ke)+np.diag(self.ke,-1)[:-1,:-1]

    def correlation_ODE(self,x,t):
        '''
        Build variance/covariance ODEs. 
        '''
        _tmp = np.dot(self.phi,x.reshape(self.N,self.N)) 
        return np.ravel(_tmp.reshape(self.N**2,1) )

    def faster_solve(self):
        '''
        a faster way to solve the system of ODEs using 
        expokit.
        '''
        phi = sp.sparse.csc_matrix(self.phi)
        self.tvec = self.gettvec()
        #self.soln = np.zeros((self.N,self.N,len(self.tvec)))
        
        self.soln = sp.sparse.csc_matrix( (self.N*self.N,len(self.tvec)) , dtype=float )
        xi = self.xi.reshape(self.N,self.N)
        xtmp = np.copy(xi)
        self.intensity = np.zeros(len(self.tvec))
 
        self.intensity[0] = self.map_to_fluorescence(xi)
        #print(self.phi)
        for i in range(1,len(self.tvec)):
            for j in range(self.N):
                xv = xtmp[:,j] 
                try:
                    xv,m,v = expv(self.tvec[i]-self.tvec[i-1],phi,xv,tol = 1e-6,m=30)
                except:
                    xv = np.zeros(xv.shape)
                #xv[np.abs(xv)<1e-8]=0
                xtmp[:,j] = xv

            self.intensity[i] = self.map_to_fluorescence(xtmp)
            if self.intensity[i]<1e-6:
                break
                self.soln[:,j,i] = xv 
        #self.soln =self.soln.reshape(self.N**2,len(self.tvec))
        return self.intensity 

class TranslateCoarseCorrs(TranslateCorrs):
    '''
    get auto and cross correlations for Coarse Grained model
    '''    
    def __init__(self,model,N):
        self.N=N
        self.binsize=int(np.ceil(model.N/float(N)))
        self.loadfullmodel(model)
        
    def loadprobe(self,probedesign):
        self.Nreal=len(probedesign)
        binary=np.zeros(self.N)
        for i in range(self.N):
            binary[i]=np.sum(probedesign[i*self.binsize:(i+1)*self.binsize])
        return binary
        
    def loadfullmodel(self,model):
        GenericODE.__init__(self)
        self.ptimes=model.ptimes
        self.ti = model.ti
        self.tf = model.tf
        self.Nreal = model.N 
        self.ke = model.ke/self.N
        self.kb = model.kb
        self.fi = model.fi
        self.params = {}
        for key,val in model.params.iteritems():
            if key is 'ke':
                self.params['ke'] = float(model.ke)*self.binsize
            else:
                self.params[key] = val
        self.binary = self.loadprobe(model.binary)
        self.build_probe(self.binary)
        
        
       
        # Get fundamental matrices.        
        self.S = self.getselongation()
        self.W0 = self.getW0()
        self.W1 = self.getW1()
        self.mu_ss = self.get_mean_SS()
        self.var_ss = self.get_var_SS()
        
        # solve for covariance steady state using ode
        self.xi= np.ravel(self.var_ss)
        self.ODE = self.correlation_ODE

        #_tmp = np.dot(np.dot(self.W1.T,self.S.T),x.reshape(self.N,self.N))
        #return np.ravel(_tmp.reshape(self.N**2,1) )
       # Ci = np.reshape(cc.soln[:,i],(cc.N,cc.N))
