==========================
rSNAPsim Module Tutorial 
==========================

.. contents::
	:depth: 5
	

Importing rSNAPsim commandline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your Python enviroment has rSNAPsim.py's directory within the enviroment paths, simply import directly.

::

	import rSNAPsim
	
	
Otherwise when importing rSNAPsim, the current working directory will have to be where rSNAPsim.py is on your system.

::

	import os    					  #import os
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


	rss = rSNAPsim.rSNAPsim()
	rss.open_seq_file('HUMINSR.gb')
	
Alternatively you can pull the human insulin resceptor protein INSR from genebank

::

	>>> rss.get_gb_file('M10051')
	
	>>> rss.gb_obj 	
	SeqRecord(seq=Seq('GGGGGGCTGCGCGGCCGGGTCGGTGCGCACACGAGAAGGACGCGCGGCCCCCAG...AAA', IUPACAmbiguousDNA()), id='M10051.1', name='HUMINSR', description='Human insulin receptor mRNA, complete cds', dbxrefs=[])	
	


The sequence file was successfully read and stored in sequence_str in both cases

::

	>>> rss.sequence_str
	'GGGGGGCTGCGCGGCCGGGTCGGTGCG...'
	
	>>> rss.sequence_name
	'HUMINSR'
	


The steps to set up a protein for simulation are as follows:

::
   
	#analyze the sequence given for open reading frames and fluorescent tags
	>>> rss.get_orfs(rss.sequence_str)
	
	>>> rss.orfs
	{'1': [(78, 4284)], '2': [(667, 979)], '3': []}
	
	>>> rss.starts[0]   #all the start codons in frame 1
	array([  78,  138,  249,  330,  384,  510,  546, 1098, 1542, 1728, 1875,
       2670, 2745, 3243, 3369, 3444, 3453, 3543, 3552, 3576, 3633, 3675,
       3744, 3792, 3885, 3942, 3948, 3972, 4092, 4107, 4227, 4383])
	
Now the rss instance contains the open reading frames, from the dictionary we can see that 
there is a protein from codon 78 to codon 4284 in frame 1 and a smaller potential protein from 667 to 979 in the 2nd frame

Next step is to translate those proteins and check for fluorescent tags

::

	>>> rss.get_temporal_proteins()
	
	>>> rss.pois	
	['MDYKDDDDKGDYKDDDDKGDYKDDDDKGGNSLIKENMRMKVVMEGSVNGHQFKCTGEGEGNPYMGTQT\
	MRIKVIEGGPLPFAFDILATSFGGGSRTFIKYPKGIPDFFKQSFPEGFTWERVTRYEDGGVVTVMQDTSL\
	EDGCLVYHVQVRGVNFPSNGPVMQKKTKGWEPNTEMMYPADGGLRGYTHMALKVDGGDYKDDDDKQQDYK\
	DDDDKGQQGDYKDDDDKQQDYKDDDDKGGGHLSCSFVTTYRSKKTVGNIKMPGIHAVDHRLERLEESDNE\
	MFVVQREHAVAKFAGLGGGGGDYKDDDDKGDYKDDDDKGDYKDDDDKGGGGSGGGGSLQMGTGGRRGAAA\
	APLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHLQILLMFKTRPEDFRDLSFPKL\
	IMITDYLLLFRVYGLESLKDLFPNLTVIRGSRLFFNYALVIFEMVHLKELGLYNLMNITRGSVRIEKNNE\
	LCYLATIDWSRILDSVEDNHIVLNKDDNEECGDICPGTAKGKTNCPATVINGQFVERCWTHSHCQKVCPT\
	ICKSHGCTAEGLCCHSECLGNCSQPDDPTKCVACRNFYLDGRCVETCPPPYYHFQDWRCVNFSFCQDLHH\
	KCKNSRRQGCHQYVIHNNKCIPECPSGYTMNSSNLLCTPCLGPCPKVCHLLEGEKTIDSVTSAQELRGCT\
	VINGSLIINIRGGNNLAAELEANLGLIEEISGYLKIRRSYALVSLSFFRKLRLIRGETLEIGNYSFYALD\
	NQNLRQLWDWSKHNLTTTQGKLFFHYNPKLCLSEIHKMEEVSGTKGRQERNDIALKTNGDKASCENELLK\
	FSYIRTSFDKILLRWEPYWPPDFRDLLGFMLFYKEAPYQNVTEFDGQDACGSNSWTVVDIDPPLRSNDPK\
	SQNHPGWLMRGLKPWTQYAIFVKTLVTFSDERRTYGAKSDIIYVQTDATNPSVPLDPISVSNSSSQIILK\
	WKPPSDPNGNITHYLVFWERQAEDSELFELDYCLKGLKLPSRTWSPPFESEDSQKHNQSEYEDSAGECCS\
	CPKTDSQILKELEESSFRKTFEDYLHNVVFVPRKTSSGTGAEDPRPSRKRRSLGDVGNVTVAVPTVAAFP\
	NTSSTSVPTSPEEHRPFEKVVNKESLVISGLRHFTGYRIELQACNQDTPEERCSVAAYVSARTMPEAKAD\
	DIVGPVTHEIFENNVVHLMWQEPKEPNGLIVLYEVSYRRYGDEELHLCVSRKHFALERGCRLRGLSPGNY\
	SVRIRATSLAGNGSWTEPTYFYVTDYLDVPSNIAKIIIGPLIFVFLFSVVIGSIYLFLRKRQPDGPLGPL\
	YASSNPEYLSASDVFPCSVYVPDEWEVSREKITLLRELGQGSFGMVYEGNARDIIKGEAETRVAVKTVNE\
	SASLRERIEFLNEASVMKGFTCHHVVRLLGVVSKGQPTLVVMELMAHGDLKSYLRSLRPEAENNPGRPPP\
	TLQEMIQMAAEIADGMAYLNAKKFVHRDLAARNCMVAHDFTVKIGDFGMTRDIYETDYYRKGGKGLLPVR\
	WMAPESLKDGVFTTSSDMWSFGVVLWEITSLAEQPYQGLSNEQVLKFVMDGGYLDQPDNCPERVTDLMRM\
	CWQFNPKMRPTFLEIVNLLKDDLHPSFPEVSFFHSEENKAPESEELEMEFEDMENVPLDRSSHCQREEAG\
	sGRDGGSSLGFKRSYEEHIPYTHMNGGKKNGRILTLPRSNPS']
	
The stored protein from the genbank file format was automatically tagged with a 10X FLAG tag (DYDDDDK) since no tag was autodetected

Finally we can analyze this protein




::

	>>> rss.analyze_poi(rss.pois[0],rss.pois_seq[0]) 
	
	>>> rss.POI.gene_length
	1382
	
	>>> rss.POI.tag_epitopes
	{'T_Flag': [2, 11, 20, 196, 206, 218, 228, 300, 309, 318]}
	
Now the POI object is filled with the information needed to run SSA simulations!

For conveinece all these functions are built into ``run_default()``
 
::

	>>> rss.run_default()
	
	>>> rss.POI.aa_seq	
	'MGTGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHLQILLMFKTRPEDFRDLSFPKLIMITDYLLLFRVYGLESLKDLFPNLTVIRGSRLFFNYALVIFEMVHLKELGLYNLMNITRGSVRIEKNNELCYLATIDWSRILDSVEDNHIVLNKDDNEECGDICPGTAKGKTNCPATVINGQFVERCWTHSHCQKVCPTICK...'


Running SSA trajectories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
SSA trajectories can now be run and saved from the stored protien of intrest

If the users just runs a blank solver with no arguments, the rSNAPsim automatically uses the stored POI object to populate the 
rates and propensities with default settings of 

- k_elong_mean = 10
- k_initiation = .03
- n_traj = 100
- tf = 1000
- ti = 0


::
	
	>>> ssa_traj = rss.ssa_solver()
	>>> ssa_traj.ivec
	array([[ 0.,  0.,  0., ..., 24., 25., 27.],
       [ 0.,  0.,  0., ..., 50., 52., 53.],
       [ 0.,  0.,  0., ..., 40., 40., 40.],
       ...,
       [ 0.,  0.,  0., ..., 33., 33., 33.],
       [ 0.,  0.,  0., ..., 30., 30., 30.],
       [ 0.,  0.,  0., ..., 30., 30., 30.]])
	   
The ssa_solver returns an object containing all the data from the simulation


ssa_obj attributes
===================

===========================  ========================================================================
Attribute                    Description
===========================  ========================================================================
``ssa_traj.all_results``     all the results from the simulation (ribosome positions per time per trajectory per max ribosomes)
``ssa_traj.solvetime``       the time it took to solve the simulation
``ssa_traj.mean_autocorr``   mean autocorrelation of fluorescence 
``ssa_traj.error_autocorr``  error of the autocorrelation
``ssa_traj.autocorr_vec``    the autocorrelation vectors
``ssa_traj.solutions``       list of solutions for each trajectory (ribosomes pos per time)
``ssa_traj.rib_means``       list of mean ribosomes per posistion
``ssa_traj.rib_density``     list of ribosome probability per codon
``ssa_traj.k``               the rates used for the simulation
``ssa_traj.fragimes``        the start time of each ribosome trajectory
``ssa_traj.fragments``       each seperate ribosome trajectory
``ssa_traj.full_frags``      the number of fully completed ribosomes (finished translation)
``ssa_traj.ke_true``         the true k_elongation rate recorded by the simulation
``ssa_traj.ke_sim``          the simulated ke from FCS autocorrelation
``ssa_traj.dwelltime``       the dwell time of each ribosome
``ssa_traj.time``            the time vector for the simulation
``ssa_traj.collisions``      recorded collisions 
===========================  ========================================================================

additionally you can choose to save the ssa object to a txt or json file. 

::

	ssa_traj.save_txt('data.txt')
   
Data Structures
~~~~~~~~~~~~~~~~~~~~~~~

The rSNAPsim contains a few different dictionaries used for the codon dependancies, these can be manipulated by the user as well

rSNAPsim data
================

===========================  ========================================================================
Attribute                    Description
===========================  ========================================================================
``rss.tag_dict``			 dictionary of tag epitopes for fluorescent tags (SunTag, Flag, and Hemaglglutinin)
``rss.tag_full``             dictionary of the full tag sequences
``rss.aa_keys``			     list of one letter amino acid characters
``rss.aa_table``			 codon to amino acid dictionary
``rss.aa_table_r``			 amino acid to its possible codons dictionary
``rss.strGeneCopy``			 GeneCopy ratios of every codon in the human genome
``rss.strGeneCopy_fast``	 Fastest codon value for each amino acid
``rss.strGeneCopy_slow``	 Slowest codon value for each amino acid
===========================  ========================================================================


Plotting and Kymographs
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

	import matplotlib.pyplot as plt
	import numpy as np
	plt.figure()
	plt.plot(ssa_traj.intensity_vec.T, alpha=.3,color='green')  #plot all trajectories
	
	plt.plot(np.mean(ssa_traj.intensity_vec.T,axis=1),color='purple') #plot the mean
	


.. figure:: tutorial_plot1.png


To plot a kymograph the kymograph function can be used

::

	#with no bg intensity
	rss.kymograph(ssa_traj,0,color='white',bg_intense=False)
	
	#with bg intensity colormap
	rss.kymograph(ssa_traj,0,color='purple')  
	
	
.. figure:: tutorial_plot2.png


.. figure:: tutorial_plot3.png


Autocorrelation with error


::
	
	#Steady State trajectories
	ssa_traj = rss.ssa_solver(n_traj = 250, start_time = 1000, tf = 2000, tstep = 2000)
	
    #just the mean trajectory
	
	plt.plot(ssa_traj.mean_autocorr,color='b')
	plt.plot(ssa_traj.mean_autocorr- ssa_traj.error_autocorr,color='b',ls='--')
	plt.plot(ssa_traj.mean_autocorr+ ssa_traj.error_autocorr,color='b',ls='--')
	plt.xlabel('time')
	plt.ylabel('normalized autocorrelation')
	plt.plot([0,1000],[.01,.01],color='red',ls=':')
	
	#With all the trajectories
	
	normalized_autocorr = ssa_traj.autocorr_vec.T/ ssa_traj.autocorr_vec[:,0]
	plt.plot(normalized_autocorr,alpha=.1,color='b')		
	plt.plot(ssa_traj.mean_autocorr,color='b')
	plt.plot(ssa_traj.mean_autocorr- ssa_traj.error_autocorr,color='b',ls='--')
	plt.plot(ssa_traj.mean_autocorr+ ssa_traj.error_autocorr,color='b',ls='--')
	plt.xlabel('time')
	plt.ylabel('normalized autocorrelation')
	plt.plot([0,1000],[.01,.01],color='red',ls=':')

.. figure:: tutorial_plot4.png

.. figure:: tutorial_plot5.png

	




