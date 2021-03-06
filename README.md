Python 2.7 or 3.5+ Version of the Single Molecule Translation Simulator (MatLab) by Dr. Luis Aguilera 

[Computational Design and Interpretation of Single-RNA Translation Experiments](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6816579/)

Translated by Will Raymond - 2018/2019

------------------------------
### **rSNAPsim** - **R**NA **S**equence to **NA**scent **P**rotein **Sim**ulation
-------------------------------
### Project Goal

Provide a Python module that takes nucleotide sequence as an input and does the following: 
   * Choose a file or pull a file from GeneBank
   * Analyzes the sequence and identifies proteins 
   * Detects or adds fluorescent tags
   * Simulates translation trajectories and converts to intensity vectors of A.U. under various conditions
      * Constructs with Rare codons only or Common codons, FRAP or Harringtonite assays
   * Provides analyses of the trajectories 
   * Allows the user to save or export the data
   * Commandline / GUI implementations

---------------------------------
### Documentation

Tutorials, Module Documentation, Installiation and more [LINK TO MUNSKY GROUP WEBSITE]

Dependencies: 
 * [NumPy](https://www.numpy.org/) 
 * [SciPy](https://www.scipy.org/)
 * [BioPython](https://biopython.org/)
 * [matplotlib](https://matplotlib.org/)
 * [pandas](https://pandas.pydata.org/)
 * [SnapGene Reader](https://github.com/IsaacLuo/SnapGeneFileReader)
 
-----------------------------------
### Instillation 

#### Pip install

upon release the current release will be added to PyPI for pip installation

##### General

open up a Terminal or Command Prompt and run 
```pip install rSNAPsim ```

##### Anaconda

Open up the Anaconda Prompt (anaconda console) and run 
```conda pip install rSNAPsim ```

##### Compilation of the C++ 

* Download the latest stable [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) release
* Place the eigen-version into the ssa_cpp folder and rename it to "eigen"
* open a command prompt/terminal or anaconda prompt and cd to this directory
* run "python setup.py build_ext --inplace"

 if it says Cannot open include file: 'eigen/Eigen/Dense': and exit status 2 you have eigen in the wrong directory

 

-------------------------------------

#### Future work

- Adding the two color visualization and probe editing to the GUI
- Addting I(t) vs I(t+tau) visualization plot
