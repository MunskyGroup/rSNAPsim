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

#### Within a conda enviroment:

```conda install eigen ```

```pip install -i https://test.pypi.org/simple/ rsnapsim-ssa-cpp ```

```pip install -i https://test.pypi.org/simple/ rsnapsim```

##### Compilation of the C++ 

The c++ model should attempt to compile when you pip install the ssa-cpp module, however in the event that it cannot here are some common errors:

* cannot include eigen3/Eigen/Dense
  * This means eigen was not installed correctly from the conda installiation, you may have to manually download [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) and pass the argument to the setup.py command. ```python setup.py build_ext --inplace -I[PATH TO EIGEN FOLDER]```
  
* gcc not found

-------------------------------------

#### Future work

- Example notebooks of all functions
