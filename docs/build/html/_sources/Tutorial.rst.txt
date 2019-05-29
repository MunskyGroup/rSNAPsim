===============
rSNAPsim tutorial 
===============

.. contents::
	:depth: 5
	

Importing rSNAPsim commandline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your Python enviroment has rSNAPsim.py's directory within the enviroment paths, simply import directly.

::

	import rSNAPsim
	
	
Otherwise when importing ssit, the current working directory will have to be where ssit.py is on your system.

::

	import os    			   #import os
	os.getcwd()					  #print your current working directory
	os.chdir('path to directory containing rSNAPsim.py')  #how to change your directory 
	
	import rSNAPsim
	
	
or you can add rSNAPsim to your current sys path

::

	import sys
	sys.path.append('path to directory containing rSNAPsim.py')
	import rSNAPsim
	
	

Getting started
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now with rSNAPsim imported, we can start a rSNAPsim instance and load a sequence file (.txt or .gb)

::


	model = rSNAPsim.rSNAPsim()
	model.open_seq_file('HUMINSR.gb')
	
	

The sequence file was successfully read and stored in sequence_str

::

	>>> model.sequence_str
	'GGGGGGCTGCGCGGCCGGGTCGGTGCG...'
	


From here we can analzye this sequence using default settings or with specified settings 


::

	>>> model.run_default()
	>>> model.POI.aa_seq	
	'MGTGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHLQILLMFKTRPEDFRDLSFPKLIMITDYLLLFRVYGLESLKDLFPNLTVIRGSRLFFNYALVIFEMVHLKELGLYNLMNITRGSVRIEKNNELCYLATIDWSRILDSVEDNHIVLNKDDNEECGDICPGTAKGKTNCPATVINGQFVERCWTHSHCQKVCPTICK...'
	
	#Lets tag the protien with a T_Flag
	>>> model.analyze_poi(model.tag_dict['T_Flag']+model.POI.aa_seq, model.tag_full['T_Flag']+model.POI.nt_seq)
	
	

SSA trajectories can now be run and saved from the stored protien of intrest

::
	
	>>> ssa_traj = sms1.ssa_solver()
	>>> ssa_traj.ivec
	array([[ 0.,  0.,  0., ..., 24., 25., 27.],
       [ 0.,  0.,  0., ..., 50., 52., 53.],
       [ 0.,  0.,  0., ..., 40., 40., 40.],
       ...,
       [ 0.,  0.,  0., ..., 33., 33., 33.],
       [ 0.,  0.,  0., ..., 30., 30., 30.],
       [ 0.,  0.,  0., ..., 30., 30., 30.]])

   
