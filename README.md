Python 2.7 or 3.5+ Version of the Single Molecule Translation Simulator (MatLab) by Dr. Luis Aguilera 

[PAPER LINK HERE]

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

-------------------------------------

#### Future work

- Adding the two color visualization and probe editing to the GUI
- Addting I(t) vs I(t+tau) visualization plot
