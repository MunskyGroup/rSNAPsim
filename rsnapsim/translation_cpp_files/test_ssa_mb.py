
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 15:43:07 2021

@author: willi
"""
import numpy as np
import rsnapsim as rsim
h2b_gene_file = '../gene_files/.txt'

poi_seqs, pois, tagged_pois, seq = rsim.seqmanip.open_seq_file(h2b_gene_file)
h2b_construct = pois['1'][0]
h2b_construct.generate_3frame_tags()



model = rsim.model_builder
model.add_k(h2b_construct)

model._poi.multiframe_epitopes[0]['T_Flag'] = [310, 320, 330,340,350]
model._poi.multiframe_epitopes[1]['T_New'] = [ 240,250,230,260,270]
model._poi.multiframe_epitopes[0]['T_Hemagglutinin']= [100,110,120,130,140, ]
model._poi.multiframe_epitopes[2]['T_New2']=  [310, 320, 330,340,350]


model.add_enters([5,50], [0,1], [.08,.03])


model.add_jumps([299], [1], [300], [2], [2])
model.visualize_transcript()

t_array = np.linspace(0,1000,1001,dtype=np.float64)
len(h2b_construct.kelong)
#r = model.run_ssa(t_array,n_traj=10,low_mem=False)