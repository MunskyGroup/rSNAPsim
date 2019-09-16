# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:12:17 2019

@author: wsraymon
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 11:58:41 2019

@author: William
"""

'''
Generate supplemental figure for luis
'''

import rSNAPsim
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})




def plot_sequence(ax,gene_length,total_length,epitopes_pos,tag_name,gene_name,codon):
    ax.cla()
    ax.plot([0,total_length],[0,0],color='k')
    ax.plot([0,total_length],[1,1],color='k')
    ax.plot([0,0],[0,1],color='k')
    ax.plot([total_length,total_length],[0,1],color='k')
    ax.axis([-10,total_length+10,-1,2])
    #ax.plot([tag_length,tag_length],[0,1],color='k',linewidth=1)
    tag_length = total_length - gene_length
    for i in range(len(epitopes_pos)):
        ax.plot([epitopes_pos[i],epitopes_pos[i]],[0,1],color='black',linewidth=2)
    ax.fill_between([0,0,tag_length,tag_length],[0,1,1,0],color='#CCCCCC')

    ticks = np.linspace(0,total_length,10).astype(int)
    ax.set_xticks(ticks)
    ax.set_xlabel('Codon Position')
    ax.get_yaxis().set_visible(False)
    ax.set_title((codon + ' location - H2B w/ T-flag'))
    
    ax.text(0,1.2,('grey: ' + tag_name), fontdict = dict(fontsize=10) )
    #ax.text(tag_length,1.3,('POI: ' + gene_name),fontdict = dict(fontsize=13))


changes = ['CTC']
percent = [1]#,.1,.03,.01,.001]
limit = [250]



files = ['gene_files/KDM5B_withTags.txt','gene_files/Bactin_withTags.txt','gene_files/H2B_withTags.txt']
rates = ['natural','slow','fast']
names = ['KDM5B','B-actin','H2B']
ki = [.03,.0603,.0612]
ke = [10.75,12,10.25]


font = {'family' : 'arial',
        
        'size'   : 25}

plt.rc('font', **font)

sms = rSNAPsim.rSNAPsim()

'''

k = 0
for f in files:
    sms = rSNAPsim.rSNAPsim()

    #sms.strGeneCopy['ACC'] = 3
    sms.open_seq_file(f)
    sms.get_orfs(sms.sequence_str, min_codons = 80)
    sms.get_temporal_proteins()
    sms.analyze_poi(sms.pois[0],sms.pois_seq[0])
    

    for i in range(len(rates)): 
        print(rates[i])
        constructk = sms.get_k_construct(sms.POI.nt_seq,ki[i],ke[i],codon_types=rates[i])
        testssa = sms.ssa_solver(n_traj=1,tf=1000,tstep=1000,force_python=True,all_k=constructk)
        plt.figure(figsize=(8,8))
        plt.tight_layout()
        plt.rc('font', **font)
        fig = plt.figure(figsize=(8,8))
        sms.kymograph(testssa,0,color='white',bg_intense=False,show_intense=False,lw=.8 )
        print(('figpics/'+names[k] + '_' +str(rates[i])+str(i)+'ky.svg'))
        plt.savefig(('figpics/'+names[k] + '_' +str(rates[i])+str(i)+'ky.svg'),dpi=300)
    k+=1
        
       
'''       
    


'''
testssa = sms.ssa_solver(n_traj=1,tf=120000,tstep=120000)
plt.hist(testssa.collisions[:testssa.full_frags], bins = np.linspace(0,limit[0],70).astype(int), density=True,align='left',facecolor='grey' )
plt.xlim(0,limit[0])

plt.xlabel('No. Collisions')
plt.ylabel('Probability')
plt.savefig(('figpics/'+codon+str(per).replace('.','')+'hst.svg'),dpi=300)
'''


sms.open_seq_file('gene_files/KDM5B_withTags.txt')
sms.get_orfs(sms.sequence_str, min_codons = 80)
sms.get_temporal_proteins()
sms.analyze_poi(sms.pois[0],sms.pois_seq[0])


for codon in changes:
    ovalue = sms.strGeneCopy[codon]
    codoninds = np.where(np.array(sms.POI.codons) == codon)[0]
    fig = plt.figure(figsize=(10,10))
    
    ax = fig.add_subplot('111')
    plot_sequence(ax,sms.POI.gene_length,sms.POI.total_length,codoninds,'T_FLAG','design3',codon)
    plt.savefig(('figpics/'+codon+'loc.svg'),dpi=300)
    for per in percent: 
        sms.strGeneCopy[codon] = ovalue*per
        testssa = sms.ssa_solver(n_traj=1,tf=1000,tstep=1000,force_python=True,k_initiation=.022,k_elong_mean=10.6)
        plt.figure(figsize=(8,7))
        plt.tight_layout()
        plt.rc('font', **font)
        fig = plt.figure(figsize=(13,8))
        sms.kymograph(testssa,0,color='white',bg_intense=False,show_intense=True,lw=.8 )
        
        plt.savefig(('figpics/'+codon+str(per).replace('.','')+'ky.svg'),dpi=300)
        plt.figure(figsize=(8,7))
        plt.tight_layout()
        
        
        
        '''
        testssa = sms.ssa_solver(n_traj=1,tf=120000,tstep=120000)
        plt.hist(testssa.collisions[:testssa.full_frags], bins = np.linspace(0,limit[0],70).astype(int), density=True,align='left',facecolor='grey' )
        plt.xlim(0,limit[0])
        
        plt.xlabel('No. Collisions')
        plt.ylabel('Probability')
        plt.savefig(('figpics/'+codon+str(per).replace('.','')+'hst.svg'),dpi=300)
        '''
   
        
'''     
    
print('big figure')
    
sms.strGeneCopy[codon] = ovalue
        

percents = np.logspace(-1,-3,20)
ncols2 = []
ovalue = sms.strGeneCopy['CTC']
i = 0
for per in percents:
    print(i)
    
    sms.strGeneCopy[codon] = ovalue*per
    testssa = sms.ssa_solver(n_traj=1,tf=400000,tstep=400000)
    ncols2.append(np.mean(testssa.collisions[:testssa.full_frags]))
    print(testssa.full_frags)
    i+=1
    
fig = plt.figure()
ax = plt.gca()
ax.scatter(percents ,ncols2 , c='grey', )

ax.set_xscale('log')
ax.set_xlabel('Percent tRNA depletion')
ax.set_ylabel('Average Collisions')
ax.set_xlim([9.9e-5,1.1])
plt.savefig(('figpics/'+codon+str(per).replace('.','')+'log.svg'),dpi=300)

'''