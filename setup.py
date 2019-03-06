# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 11:44:48 2018

@author: william
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rSNAPsim", 
         
    
    
    version="0.0.1a",
    author="Dr. Luis Aguilera, William Raymond, Dr. Brian Munsky",
    author_email="wsraymon@rams.colostate.edu",
    description="A package for mRNA sequence translation stochastic simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    #url="githublink",
    packages=setuptools.find_packages(),
    install_requires = ['BioPython','NumPy','SciPy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)