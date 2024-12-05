# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 12:06:39 2024

@author: willi
"""

import sys, os
#print(os.listdir('.'))
print(os.getcwd())
sys.path.append('.')
cwd = os.getcwd()
os.chdir('../..')
import rsnapsim as rsnp
os.chdir(cwd)
import numpy as np

import multiprocessing
import time
import inspect

example_mRNA = '''ATGGCGAACCTTGGCTGCTGGATGCTGGTTCTCTTTGTGGCCACATGGAGTGACCTGGGC
                    CTCTGCAAGAAGCGCCCGAAGCCTGGAGGATGGAACACTGGGGGCAGCCGATACCCGGGG
                    CAGGGCAGCCCTGGAGGCAACCGCTACCCACCTCAGGGCGGTGGTGGCTGGGGGCAGCCT
                    CATGGTGGTGGCTGGGGGCAGCCTCATGGTGGTGGCTGGGGGCAGCCCCATGGTGGTGGC
                    TGGGGACAGCCTCATGGTGGTGGCTGGGGTCAAGGAGGTGGCACCCACAGTCAGTGGAAC
                    AAGCCGAGTAAGCCAAAAACCAACATGAAGCACATGGCTGGTGCTGCAGCAGCTGGGGCA
                    GTGGTGGGGGGCCTTGGCGGCTACATGCTGGGAAGTGCCATGAGCAGGCCCATCATACAT
                    TTCGGCAGTGACTATGAGGACCGTTACTATCGTGAAAACATGCACCGTTACCCCAACCAA
                    GTGTACTACAGGCCCATGGATGAGTACAGCAACCAGAACAACTTTGTGCACGACTGCGTC
                    AATATCACAATCAAGCAGCACACGGTCACCACAACCACCAAGGGGGAGAACTTCACCGAG
                    ACCGACGTTAAGATGATGGAGCGCGTGGTTGAGCAGATGTGTATCACCCAGTACGAGAGG
                    GAATCTCAGGCCTATTACCAGAGAGGATCGAGCATGGTCCTCTTCTCCTCTCCACCTGTG
                    ATCCTCCTGATCTCTTTCCTCATCTTCCTGATAGTGGGATGA'''.replace('\n','').replace(' ','')



# Make the mRNA object
poi = rsnp.seqmanip.seq_to_protein_obj(example_mRNA) # convert a given sequence to a protein of interest object
mRNA = poi['0'][0]                              # pull out the main open reading frame
mRNA_length = len(mRNA.kelong)                  # get the length of the mRNA
mRNA.generate_3frame_tags()                     # call to generate all open reading frames

# manually adding epitopes for two colors
mRNA.multiframe_epitopes[0] = {'T_Flag': [1, 10, 19, 195, 205, 217, 227, 299, 308, 317], 'T_HA':[400,410,420,430]}

#mRNA.visualize_mrna_strand() # plot open reading frame 1

model = rsnp.tasep_model(mRNA, 'dropoff') # model object




# Make the kelong mat (manually adding an extra location that is equal to zero, so particles dont run over the simulation)
kelong_mat = np.zeros([3, mRNA_length+1])
kelong_mat[0, :-1] = rsnp.propf.get_k(mRNA.nt_seq, .033, 10, 10)[1:-1]
kelong_mat[1, :-2] = rsnp.propf.get_k(mRNA.nt_seq[1:-2], .1, 10, 10)[1:-1]
kelong_mat[2, :-2] = rsnp.propf.get_k(mRNA.nt_seq[2:-1], .1, 10, 10)[1:-1]
kelong_mat[0, -2] = 0

model._kelong_mat = kelong_mat #override the current kelongation mat

footprint = 9
parameters = [0.03, 10, .4, 0]

init = lambda k,t,p,ke,o,l,pr,s,r,nr: (1-np.any(l[0:0+footprint]))*k[0]
model.add_lattice_reaction(init, parameters, rxn_name = 'initiation', exclusion=1, frame=0, loc=0, dexist=1,)


leave = lambda k,t,p,ke,o,l,pr,s,r,nr: l[590]*k[1] #(lattice location 590 = 1) * parameter
model.add_lattice_reaction(leave, parameters, rxn_name='termination', exclusion=0, frame=0, loc=590, dexist=-1,)


# ribosome has a random chance to disconnect

random_leave = lambda k,t,p,ke,o,l,pr,s,r,nr: [k[3] for i in range(nr)]
model.add_ribosome_reaction(random_leave, parameters, rxn_name='drop_off', exclusion=0, dexist=-1, ) #random chance to disconnect


# DEFAULT STEPPING OF ELONGATION USING THE ELONGATION MATRIX
#default stepping
elong_step = lambda k,t,p,ke,o,l,pr,s,r,nr: [ (ke[p[i,2], p[i,3]])*(1 - sum(l[p[i,3]+1:p[i,3]+footprint])) for i in range(nr)]
model.add_ribosome_reaction(elong_step, parameters, rxn_name='elongation', exclusion=1, dloc=1) 


# finally specify which reactions are ribosome specific
model._ribosome_reactions = [2,3]
model._constant_reactions = [0,1,]
model._lattice_arr0 = np.zeros([model._length+1], dtype=int)


model.compile_model_c()

class CustomSSASoln:
    def __init__(self, mRNA_model, rib_array, state_array, resource_array, t, burnin, n_traj, solve_time):
        self.ribosome_array = rib_array
        self.state_array = state_array
        self.resource_array = resource_array
        self.t = t
        self.burnin = burnin
        self.L = mRNA_model._length
        self.kelong_mat = mRNA_model._kelong_mat
        self.model_id = mRNA_model.model_id
        self.n_traj = n_traj
        self.n_colors = mRNA_model._n_colors
        self.solve_time = solve_time
        
    @property
    def lattice_arr(self):
        lattice_arr = np.zeros([self.n_traj, len(self.t), self.L ] ,dtype=int)
        for i in range(self.n_traj):
            for t in range(len(self.t)):
                rtraj = self.ribosome_array[i,t,:,3][self.ribosome_array[i,t,:,1]==1]
                if len(rtraj) > 0:
                    lattice_arr[i,t,rtraj] = 1
        return lattice_arr
    
    @property
    def intensity_arr(self):
        arr = np.zeros([self.n_traj, len(self.t), self.n_colors ] ,dtype=int)
        for i in range(self.n_traj):
            for t in range(len(self.t)):
                traj = np.sum(self.ribosome_array[i,t,:,4:4+self.n_colors],axis=0)
                arr[i,t,:] = traj
        return arr
    
    @property
    def I(self):
        return self.intensity_arr






class SSASolver:
    def __init__(self, NUMBER_OF_CORES=None):
        if NUMBER_OF_CORES == None:
            self.NUMBER_OF_CORES = int(multiprocessing.cpu_count()/2) #claim half the cores
        else:
            self.NUMBER_OF_CORES = 1
        

    def solve_ssa(self, mRNA_model, t, n_traj=1, burnin=0, seed=None, parallel=False, cplus=False):
        st = time.time()
        if seed == None:
            seeds = np.random.randint(0x7FFFFF, size=n_traj)
        
        if not cplus:
            if parallel:
                solns = Parallel(n_jobs = self.NUMBER_OF_CORES, prefer="threads")(delayed(self.__run)(*args, **kwargs) for args, kwargs in ([[(model, t), {'burnin':burnin, 'seed':seeds[i]}] for i in range(0,n_traj)]))
    
            else:
                solns = []
                for i in range(n_traj):
                    print(seeds[i])
                    solns.append(self.__run(mRNA_model, t, burnin=burnin, seed=seeds[i]))
            
            rib_array = np.array([solns[i][0] for i in range(len(solns))])
            resource_array = np.array([solns[i][2] for i in range(len(solns))])
            state_array = np.array([solns[i][1] for i in range(len(solns))])
            solve_time = time.time() - st
            soln = CustomSSASoln(mRNA_model, rib_array, state_array, resource_array, t, burnin, n_traj, solve_time) 
            
            
        if cplus:
            
            mRNA_model.cmodel

            parameters = model._parameters
            pars = []
            npars = []
            for i in range(len(parameters)):
                if isinstance(parameters[i], list):
                    npars.append(len(parameters[i]))
                    for j in range(len(parameters[i])):
                        pars.append(parameters[i][j])
                else:
                    npars.append(1)
                    pars.append(parameters[i])
                
            
            if len(model._state_arr0) == 0:
                temp_state_arr0 = np.zeros([1], dtype=np.int32)
            else:
                temp_state_arr0 = model._state_arr0
                
            if len(model._resource_arr0) == 0:
                temp_resource_arr0 = np.zeros([1], dtype=np.int32)
            else:
                temp_resource_arr0 = model._resource_arr0
            n_states = 0
            n_resources = 0 
            
            st = time.time()
            def convert_c_pa(pa):
                return pa.swapaxes(1,0).reshape([1,len(t), model._rib_arr0.shape[1], model._rib_arr0.shape[0]]).swapaxes(-1,-2)
            
            
            ribosome_array = np.zeros([n_traj, len(t), model._rib_arr0.shape[0], model._rib_arr0.shape[1] ])
            state_array = np.zeros([n_traj, len(t), len(model._state_arr0)+1])
            resource_array = np.zeros([n_traj, len(t), len(model._resource_arr0)+1])
            
            for i in range(n_traj):
                seed = seeds[i]

                pa = np.zeros([ model._rib_arr0.shape[0]*model._rib_arr0.shape[1], len(t)], dtype=np.int32, order='C')
                sa = np.zeros([len(t), len(model._state_arr0)+1], dtype=np.int32)
                ra = np.zeros([len(t), len(model._resource_arr0)+1], dtype=np.int32)
                
                mRNA_model.cmodel.run_ssa_cpp(pa, sa, ra, model._rib_arr0, temp_state_arr0, temp_resource_arr0,
                                           model._rxn_mat.astype(np.int32), np.array(model._constant_reactions + model._ribosome_reactions),
                                           model._kelong_mat,
                                           model._probe_mat.astype(np.int32), np.array(pars),
                                           t,
                                           model._rib_arr0.shape[0], model._n_states, model._n_resources, len(model._constant_reactions),
                                           len(model._ribosome_reactions),
                                           burnin, seed, )
                
                ribosome_array[i] = convert_c_pa(pa)
                state_array[i] = np.copy(sa)
                resource_array[i] = np.copy(ra)
                
                #print(pa)
            solve_time = time.time() - st
            soln = CustomSSASoln(model, ribosome_array, state_array, resource_array, t, burnin, n_traj, solve_time)         

        return soln
    

        
        
    def __initalize_trajs(self, mRNA_model):
        rib_arr, lattice_arr, state_arr, resource_arr = mRNA_model._x0()
        NR = np.sum(lattice_arr)
        occupied = rib_arr[:,3]
        return np.copy(rib_arr), np.copy(lattice_arr), np.copy(state_arr), np.copy(resource_arr), np.copy(occupied), NR
        
    def __constants(self, mRNA_model):
        particle_size = mRNA_model.particle_size
        rxn_mat = mRNA_model._rxn_mat.astype(int)
        n_colors = int(np.max(mRNA_model._probe_mat))
        n_states = mRNA_model._n_states
        n_resources = mRNA_model._n_resources
        n_rxns = int(rxn_mat.shape[0])
        n_ribosome_reactions = np.sum(rxn_mat[:,0] == 0)
        max_rib = mRNA_model._max_particles
        n_constant_reactions = len(rxn_mat) - n_ribosome_reactions
        L = kelong_mat.shape[1]
        
        return particle_size, max_rib, rxn_mat, n_colors, n_states, n_rxns, n_resources, n_constant_reactions, n_ribosome_reactions, L
        




        
    def __run(self, mRNA_model, t, burnin=0, seed=None):
        
        # initalize constants, flags, and initial state of the simulation
        footprint, max_rib, rxn_mat, n_colors, n_states, n_rxns, n_resources, n_constant_reactions, n_ribosome_reactions, L = self.__constants(mRNA_model)
        reaction_taken, dexist, rib_id, ribosome_moved, tindex = [0,]*5
        rib_arr, lattice_arr, state_arr, resource_arr, occupied, NR = self.__initalize_trajs(mRNA_model)
        kelong_mat = mRNA_model._kelong_mat
        probe_mat = mRNA_model._probe_mat.astype(int)
        #rint(state_arr)
        if seed != None:
            np.random.seed(seed)
            
        # initalize propensities
        constant_props = [mRNA_model._propensities[i] for i in mRNA_model._constant_reactions]
        constant_parameters = [mRNA_model._parameters[i] for i in mRNA_model._constant_reactions]

        def get_constant_props(*props):
            return [props[i](constant_parameters[i],tc, rib_arr, kelong_mat, occupied, lattice_arr, probe_mat, state_arr, resource_arr, NR) for i in range(len(props))]
            
        ribosome_props = [mRNA_model._propensities[i] for i in mRNA_model._ribosome_reactions]
        ribosome_parameters = [mRNA_model._parameters[i] for i in mRNA_model._ribosome_reactions]
        
        def get_ribosome_props(*props):
            rates = [props[i](ribosome_parameters[i],tc, rib_arr, kelong_mat, occupied, lattice_arr, probe_mat, state_arr, resource_arr, NR) for i in range(len(props))]
            return [item for sublist in rates for item in sublist]
        
        probe_function = mRNA_model._probe_function
        probe_fun = inspect.getsourcelines(probe_function)[0][0].split(':')[-1].replace('\n','').replace(' ','')
        if probe_fun == '1':
            use_probe_fun = False
        else:
            use_probe_fun = True
        
        probe_parameters = mRNA_model._probe_parameters
        
        reaction_ids = mRNA_model._constant_reactions + mRNA_model._ribosome_reactions 
        
        # Initalize arrays to store the trajectory
        #intensity_array = np.zeros([len(t), n_colors], dtype=int)
        ribosome_array = np.zeros([len(t), max_rib, 4+n_colors+n_constant_reactions+n_ribosome_reactions], dtype=int)
        state_array = np.zeros([len(t), n_states], dtype=int)
        resource_array = np.zeros([len(t), n_resources], dtype=int)
        
        
        tc = 0-burnin # current time
        tf = t[-1] # final time point
        
        while tc < tf:

            # get propensities and where
            rates = get_constant_props(*constant_props) + get_ribosome_props(*ribosome_props)
            if sum(n < 0 for n in rates):
                msg = 'Negative Rate detected, double check model design.'\
                    ''
                raise custom_err.NegativeRateError(msg)

            # select propensity

            rate_sum = np.cumsum(rates)

            tc = (tc-np.log(np.random.rand())/rate_sum[-1]) # Update the time point randomly
            ro = rate_sum[-1]*np.random.rand()  #draw random number for reaction

            # record
            while tc >= t[tindex]:


                ribosome_array[tindex] = rib_arr
                state_array[tindex] = state_arr
                resource_array[tindex] = resource_arr
                #intensity_array[tindex] = np.sum(rib_arr[:,4:4+n_colors],axis=0)
                tindex += 1
                if tindex == len(t):
                    return ribosome_array, state_array, resource_array
                    break


            for i in range(len(rates)): #pick the next reaction rate
                if rate_sum[i] >= ro:
                    event = i #which raction happened
                    rid = i
                    break
            if event>=n_constant_reactions: #if its a ribosome reaction, on which ribosome did it occur

                rib_ind = (event-n_constant_reactions)%NR  #edit event to match reaction matrix
                event = int((event-n_constant_reactions)/NR) + n_constant_reactions

            event = reaction_ids.index(event) # map the event to the reactions (reaction matrix and propensities may not match)
            # do the reaction
            ribosome_moved = 0
            if rxn_mat[event][0] == 2: # lattice reaction
                # find the matching ribosome
                fr = rxn_mat[event][2]
                loc = rxn_mat[event][3]


                if rxn_mat[event][4] == 1: #ribosome arriving
                    rib_ind = NR
                    rib_id += 1
                    rib_arr[NR,0] = rib_id
                    rib_arr[NR,1] = 1
                    rib_arr[NR,2] += rxn_mat[event][2]
                    rib_arr[NR,3] += rxn_mat[event][3]
                    rib_arr[NR,4+n_colors+event] += 1
                    NR += 1
                    ribosome_moved = 1



                elif rxn_mat[event][4] == -1: #ribosome leaving
                    rib_ind = np.where(rib_arr[:,3] == loc)[0][0]
                    rib_arr[rib_ind,[1,2,3]] = rib_arr[rib_ind, [1,2,3]] + rxn_mat[event,4:7]
                    NR -= 1
                    rib_arr[rib_ind] = 0
                    rib_arr[rib_ind:-1] = rib_arr[rib_ind+1:]
                    rib_arr[-1] = 0
                    lattice_arr[:] = 0
                    lattice_arr[occupied[:NR]] = 1

                else: # ribosome moving
                    rib_ind = np.where(rib_arr[:,3] == loc)[0][0]
                    #rib_arr[rib_ind,[1,2,3]] = rib_arr[rib_ind, [1,2,3]] + rxn_mat[event,4:7]
                    rib_arr[rib_ind, 1:(n_colors+4+1)] = rib_arr[rib_ind, 1:(n_colors+4+1)] + rxn_mat[event,4:n_colors+8]
                    rib_arr[NR,4+n_colors+event] += 1
                #change states
                if n_states > 0:
                    state_arr = state_arr +  rxn_mat[event][n_colors+7:n_colors+7+n_states]
                    if sum(n < 0 for n in state_arr):
                        msg = 'Negative state detected, check model design.'\
                            ''
                        raise custom_err.StatesError(msg)
                    #print('***')
                    #print(rxn_mat[event][n_colors+7:n_colors+7+n_states]  )

                #change resources
                if n_resources > 0:
                    resource_arr = resource_arr +  rxn_mat[event][n_colors+7+n_states:]
                    if sum(n < 0 for n in resource_arr):
                        msg = 'Negative resource detected, check model design.'\
                            ''
                        raise custom_err.NegativeResourcesError(msg)


                if rxn_mat[event][6] !=0:
                    ribosome_moved = 1

            if rxn_mat[event][0] == 0: #ribosome reaction
                # change movement and colors
                # dexist, dframe, dloc, dprobe1... dprobeN
                rib_arr[rib_ind, 1:(n_colors+4)] = rib_arr[rib_ind, 1:(n_colors+4)] + rxn_mat[event,4:n_colors+7]

                #change states
                if n_states > 0:
                    state_arr = state_arr +  rxn_mat[event][n_colors+7:n_colors+7+n_states]
                    if sum(n < 0 for n in state_arr):
                        msg = 'Negative state detected, check model design.'\
                            ''
                        raise custom_err.StatesError(msg)

                #change resources
                if n_resources > 0:
                    resource_arr = resource_arr +  rxn_mat[event][n_colors+7+n_states:]
                    if sum(n < 0 for n in resource_arr):
                        msg = 'Negative resource detected, check model design.'\
                            ''
                        raise custom_err.NegativeResourcesError(msg)

                ribosome_moved = 0
                if rxn_mat[event][6] !=0:
                    ribosome_moved = 1

                if rxn_mat[event][4] == -1: #ribosome left
                    NR -= 1
                    rib_arr[rib_ind] = 0
                    rib_arr[rib_ind] = 0
                    rib_arr[rib_ind:-1] = rib_arr[rib_ind+1:]
                    rib_arr[-1] = 0
                    lattice_arr[:] = 0
                    lattice_arr[occupied[:NR]] = 1

                rib_arr[rib_ind,4+n_colors+event] += 1


            if rxn_mat[event][0] == 1: # state reaction
                state_arr = state_arr +  rxn_mat[event][n_colors+7:n_colors+7+n_states]
                if sum(n < 0 for n in state_arr):
                    msg = 'Negative state detected, check model design.'\
                        ''
                    raise custom_err.StatesError(msg)

            if rxn_mat[event][0] == 3: # resource reaction
                resource_arr = resource_arr +  rxn_mat[event][n_colors+7+n_states:]
                if sum(n < 0 for n in resource_arr):
                    msg = 'Negative resource detected, check model design.'\
                        ''
                    raise custom_err.NegativeResourcesError(msg)

            if rxn_mat[event][0] == 4: # probe reaction
                state_mat[rxn_mat[event][1]] += rxn_mat[event][2]
                state_mat[rxn_mat[event][3]] += rxn_mat[event][4]
                state_mat[rxn_mat[event][5]] += rxn_mat[event][6]


            # Check probes
            if ribosome_moved:
                # update occupied vector
                occupied = rib_arr[:,3]

                # update the lattice vector
                lattice_arr[:] = 0
                lattice_arr[occupied[:NR]] = 1

                # check for probes
                pr = probe_mat[rib_arr[rib_ind,2], rib_arr[rib_ind,3]]
                if pr != 0:
                    if use_probe_fun:
                        if probe_function(ribosome_parameters[i],tc, rib_arr, kelong_mat, occupied, lattice_arr, probe_mat, state_arr, resource_arr, NR):
                            rib_arr[rib_ind, 3+pr] +=1 #add a probe if it passed a location and passes the probe function
                    else:
                        rib_arr[rib_ind, 3+pr] +=1
                ribosome_moved = 0


            #print('-----')
            #print(rib_arr)
            #print(rates)
            #print(event)
            #print(rid)
            #print(state_arr)






solver = SSASolver()
t = np.linspace(0,500,5001)
for i in range(1):
    soln = solver.solve_ssa(model, t)

#model.compile_model_c()



solver = SSASolver()
t = np.linspace(0,500,5001)
solnc = solver.solve_ssa(model, t, cplus=True)


