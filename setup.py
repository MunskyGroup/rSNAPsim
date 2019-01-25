# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 11:44:48 2018

@author: willi
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rSNAPsim", 
    
    #potential names:   CodonDependentTranslation - CDTsim (cdt taken)
    #                   mRNAsim
    #                   seqsim
    #                   seq2sim - NT sequence to simulation
    #                   seqsim
    #                   SSAseq
     
    
    
    version="0.0.1",
    author="Dr. Luis Aguilera, William Raymond, Dr. Brian Munsky",
    author_email="wsraymond@rams.colostate.edu",
    description="A package for mRNA sequence translation stochastic simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="githublink",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.5+ or 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)