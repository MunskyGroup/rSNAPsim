===============
Installation 
===============

.. contents::
	:depth: 5
	
Installing Single Molecule Simulator 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dependencies 
=============

* Python 2.7 or Python 3.5+
* `SciPy <https://www.scipy.org/>`_
	- NumPy
	- SymPy
	- Matplotlib
* `BioPython <https://biopython.org/>`_

	
Recommended Enviroment
=========================

`Anaconda 2.7 or 3.5+ <https://conda.io/docs/user-guide/install/download.html>`_

Manually Compiling the Stochastic Simulation for C++
=====================================================

Download :download:`translation_ssa <translation_ssa_v01.zip>` and uncompress it.

Download `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ and place it in the translation_ssa folder

Windows:

* open a new command prompt 
* change directory to the translation_ssa 
* run the following command: "python setup.py build_ext"
	* if this doesnt work, try adding --inplace at the end of the command

Mac:

* open a new terminal window
* change directory to the translation_ssa 
* run the following command: "python setup.py build_ext"
	* if this doesnt work, try adding --inplace at the end of the command

then add ssa_translation to your python path

::

	import sys
	sys.path   #prints current packages
	sys.path.append('directory of translation_ssa folder')
	
From here the C++ library for under sms.ssa_solver or sms.ssa

Download
===============
[Githublink here]

Adding SMS to System path directly
=====================================

Python enviroment variables can be accessed and adjusted by opening a new Python console and appending sys.path

::

	import sys
	sys.path   #prints current packages
	sys.path.append('directory of downloaded SMS')
	

Adding SMS to Spyder IDE
============================

Tools > PYTHONPATH Manager > Add path > select directory to SMS
	
.. figure:: install_spyder_tools.png

.. figure:: install_spyder_manager.png

.. figure:: install_spyder_path.png
	
	
Pip Install
===============

[Eventual pip installation including C++ library]
	
	
	

