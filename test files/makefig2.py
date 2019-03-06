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
    ax.set_title((codon + ' location - Bactin w/ T-flag'))
    
    ax.text(0,1.2,('grey: ' + tag_name), fontdict = dict(fontsize=10) )
    #ax.text(tag_length,1.3,('POI: ' + gene_name),fontdict = dict(fontsize=13))


changes = ['CTA','CGG','GAT']
percent = [1,.01,.005]

sms = rSNAPsim.rSNAPsim()








#sms.strGeneCopy['ACC'] = 3
sms.open_seq_file("gene_files/Bactin_withTags.txt")
sms.get_orfs(sms.sequence_str, min_codons = 80)
sms.get_temporal_proteins()
sms.analyze_poi(sms.pois[0],sms.pois_seq[0])

font = {'family' : 'arial',
        
        'size'   : 25}

plt.rc('font', **font)


for codon in changes:
    ovalue = sms.strGeneCopy[codon]
    codoninds = np.where(np.array(sms.POI.codons) == codon)[0]
    fig = plt.figure(figsize=(10,4))
    
    ax = fig.add_subplot('111')
    plot_sequence(ax,sms.POI.gene_length,sms.POI.total_length,codoninds,'T_FLAG','Bactin',codon)
    plt.savefig(('figpics/'+codon+'loc.svg'),dpi=300)
    for per in percent: 
        sms.strGeneCopy[codon] = ovalue*per
        testssa = sms.ssa_solver(n_traj=1,tf=10000,tstep=10000)
        plt.figure(figsize=(8,7))
        plt.tight_layout()
        plt.rc('font', **font)
        sms.kymograph(testssa,0,color='white',bg_intense=False,show_intense=False,lw=.8)
        plt.savefig(('figpics/'+codon+str(per).replace('.','')+'ky.svg'),dpi=300)
        plt.figure(figsize=(8,7))
        plt.tight_layout()

    
        plt.hist(testssa.collisions, bins = int(np.max(testssa.collisions))/5+1, density=True,align='left',facecolor='grey' )
        plt.savefig(('figpics/'+codon+str(per).replace('.','')+'hst.svg'),dpi=300)
        
    

        
    sms.strGeneCopy[codon] = ovalue
        