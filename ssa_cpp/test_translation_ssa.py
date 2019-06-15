import numpy as np
import ssa_translation
import matplotlib.pyplot as plt
import time

# load the elongation 
kelong = np.loadtxt('elongationrates.txt')
kbind = kelong[0]
kcompl = kelong[-1]
kelong = kelong[1:-1]

t_array = np.array([0,100,500],dtype=np.float64)
t0 = 15
t_array = np.linspace(0,1000,1000,dtype=np.float64)
N_rib = 200
result = np.zeros((len(t_array)*N_rib),dtype=np.int32  )
#kelong = np.array([3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5],dtype=np.float64)
n_trajectories = 10
start = time.time()
all_results = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)
lenfrap = len(np.intersect1d(np.where(t_array>0)[0],np.where(t_array<20)[0]))

all_frapresults = np.zeros((n_trajectories,N_rib*len(t_array)),dtype=np.int32)




all_ribtimes = np.zeros((n_trajectories,400),dtype=np.float64)
all_coltimes = np.zeros((n_trajectories,400),dtype=np.int32)
nribs = np.array([0],dtype=np.int32)
all_ribs = np.zeros((n_trajectories,1))
seeds = np.random.randint(0,0x7FFFFFF,n_trajectories)
all_col_points = []
for i in range(n_trajectories):
    result = np.zeros((len(t_array)*N_rib),dtype=np.int32)    
    frapresult = np.zeros((len(t_array)*N_rib),dtype=np.int32)
    
    ribtimes = np.zeros((400),dtype=np.float64)
    coltimes = np.zeros((400),dtype=np.int32)
    
    colpointsx = np.zeros(len(kelong)*400,dtype=np.int32)
    colpointst = np.zeros(len(kelong)*400,dtype=np.float64)
    
    ssa_translation.run_SSA(result,ribtimes,coltimes,colpointsx,colpointst, kelong,frapresult,t_array,.03,kcompl, 1,0,300, seeds[i],nribs)
    all_results[i,:] = result
    all_frapresults[i,:] = frapresult
    all_coltimes[i,:] = coltimes
    all_ribtimes[i,:] = ribtimes
    all_ribs[i,:] = nribs[0]

    endcolrec = np.where(colpointsx == 0)[0][0]
    
    colpoints = np.vstack((colpointsx[:endcolrec],colpointst[:endcolrec]))
    all_col_points.append(colpoints.T)
print('time for {0} trajectories {1}'.format(n_trajectories,time.time()-start))
#plt.hist(result[result>0])
#plt.show()
#traj = result.reshape((N_rib,len(t_array))).T
##print('The result is \n {0}'.format(result.reshape((N_rib,len(t_array))).T))
#plt.plot(traj[-1,:])
#plt.show()

# map to fluorescence.
ntimes = len(t_array)
intensity_vec = np.zeros(ntimes)
pv = np.loadtxt('probe_design.txt')
tstart = 0
I = np.zeros((n_trajectories,ntimes-tstart))
for i in range(n_trajectories):
    traj = all_results[i,:].reshape((N_rib,len(t_array))).T
    for j in range(tstart,ntimes):
        temp_output = traj[j,:]
        I[i,j] = np.sum(pv[temp_output[temp_output>0]-1])

# Plotting
all_traj = np.loadtxt('ivec_1000t')
#all_traj = ssa_obj_01.intensity_vec
f,ax = plt.subplots(2,1)
ax[0].plot(I[0:50,-500:].T)
ax[1].plot(all_traj[0:50,-500:].T)
f2,ax2 = plt.subplots(1,2)
ax2[0].hist(I[:,-500:].ravel())
ax2[1].hist(all_traj[:100,-500:].ravel())
plt.show()

