# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 11:44:48 2018

@author: william
"""

import setuptools


import re
VERSIONFILE="rsnapsim/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rsnapsim", 
    version=verstr,
    author="Dr. Luis Aguilera, William Raymond, Dr. Brian Munsky",
    author_email="wsraymon@rams.colostate.edu",
    description="A package for mRNA sequence translation stochastic simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    #url="githublink",
    packages=setuptools.find_packages(exclude=('ssa_cpp.*','test.*','models.*','trna_ssa.*','ssa_cpp')),
    include_package_data=True,
    package_data = {"rsnapsim":['*.h', '*.cpp','*.pyx']},
    install_requires = ['BioPython','numpy','scipy','cython','matplotlib','pandas', 'SnapGeneFileReader','dna_features_viewer'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)