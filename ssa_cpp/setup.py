try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
    
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import os

# Call cythonize in advance so a single module can be compiled from a single Cython
# file along with other C++ files.
import numpy as np
sources = ['ssa_translation.pyx','ssa_translation_c_w.cpp']
cythonize('*.pyx', language='c++')
setup(name='SSA',ext_modules=[Extension('ssa_translation', sources, language='c++',include_dirs = ['/usr/local/include','/Library/Developer/CommandLineTools/usr/bin','/anaconda/lib/python3.6/site-packages/numpy/core/include',np.get_include(),os.getcwd(),('C:\\Users\\' + os.environ['USERNAME'])])],cmdclass = {'build_ext': build_ext})
