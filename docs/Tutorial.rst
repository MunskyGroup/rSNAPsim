===============
SMS tutorial 
===============

.. contents::
	:depth: 5
	

Importing SMS commandline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your Python enviroment has sms.py's directory within the enviroment paths, simply import directly.

::

	import sms
	
	
Otherwise when importing ssit, the current working directory will have to be where ssit.py is on your system.

::

	import os    					  #import os
	os.getcwd()					  #print your current working directory
	os.chdir('path to directory containing sms.py')  #how to change your directory 
	
	import sms
	
	
or you can add sms to your current sys path

::

	import sys
	sys.path.append('path to directory containing sms.py')
	import sms
	
	

Getting started
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now with sms imported, we can start a sms instance and load a sequence file (.txt or .gb)

::


	sms1 = sms.sms()
	sms1.open_seq_file('HUMINSR.gb')
	
	

The sequence file was successfully read and stored in sequence_str

::

	>>> sms1.sequence_str
	'GGGGGGCTGCGCGGCCGGGTCGGTGCG...'
	


From here we can analzye this sequence using default settings or with specified settings 


::

	>>> sms1.run_default()
	>>> sms1.POI.aa_seq	
	'MGTGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHLQILLMFKTRPEDFRDLSFPKLIMITDYLLLFRVYGLESLKDLFPNLTVIRGSRLFFNYALVIFEMVHLKELGLYNLMNITRGSVRIEKNNELCYLATIDWSRILDSVEDNHIVLNKDDNEECGDICPGTAKGKTNCPATVINGQFVERCWTHSHCQKVCPTICK...'
	
	#Lets tag the protien with a T_Flag
	>>> sms1.analyze_poi(sms1.tag_dict['T_Flag']+sms1.POI.aa_seq, sms1.tag_full['T_Flag']+sms1.POI.nt_seq)
	
	

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

   
