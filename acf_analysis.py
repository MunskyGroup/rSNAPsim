import numpy as np
import matplotlib.pyplot as plt 
import translation_models as models
import time

import rSNAPsim


def get_acf_analytical(kelong,probe_placements,tf,ptimes):
    '''
    Compute the autocorrelation function for the specified times, with the specified parameters. 
    '''
    m = models.TranslateCorrs()
    m.N = len(kelong)-1
    m.tf = tf
    m.ptimes = ptimes
    m.ke = kelong[1:]
    #m.ke = 13.567*np.ones(kelong[1:].shape[0])
   # m.ke[0] = 0.0
    #m.kb = kelong[0] 
    m.kb = 0.0033 
    m.fi = 1
    m.ti = 0

    # Solve correlations
    print("*****SOLVING MOMENT EQUATIONS*****")
    m.binary = probe_placements
    start = time.time()
    m.csolve()
    solve_time = time.time()-start
    print("Time to solve: %f" %solve_time)
    print("Done.")
    mean_I = m.map_to_fluorescence3(m.mu_ss)
    var_I = m.map_to_fluorescence(m.var_ss)   
    print(mean_I)
    print(var_I)
    return  m.tvec,np.ravel((m.intensity)/var_I),m
#    return m


if __name__ == '__main__':
    
    
    sms = rSNAPsim.rSNAPsim()
    sms.open_seq_file('gene_files/Bactin_withTags.txt')
    
    sms.get_orfs(sms.sequence_str, min_codons = 80)
    sms.get_temporal_proteins()
    sms.analyze_poi(sms.pois[0],sms.pois_seq[0])
    rates = sms.get_k(sms.POI.nt_seq,.03,.1)
    print(rates)
    
    kelong = rates
    

    #    mu = np.mean(kelong)
#    kelong = np.ones(101)*10
#    kelong = kelong+5*np.random.randn(101)
#    kelong[kelong<0] = 10
#    print(kelong)
    probes = np.loadtxt('bactin_binary.txt').ravel()
#    probes = np.concatenate((probes[1:],[0]))
    probes = probes[1:]
    print(probes)
#    probes = np.zeros(100)
#    probes[[0,1,2,3,4,5,6,7]] = 1
#    m = get_acf_analytical(kelong,probes,10,20)
    #tau,I = get_acf_analytical(kelong,probes,80,81)
    tau,I,m = get_acf_analytical(kelong,probes,1200,1201)
    data = np.vstack((tau,I))
    np.savetxt('bactin_acf_data_01.txt',data)
    np.save('bactin_intensity_01',I)
    # load the autocorrelation from the SSA. 
    plt.plot(tau,I)

#    acf_ssa = np.loadtxt('h2b_ssa_acf.txt')
#    plt.plot(acf_ssa[:80])

    plt.show()
    
