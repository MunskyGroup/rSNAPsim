try:
    from setuptools import find_packages
    from setuptools import setup
    from setuptools import Extension
    
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
    
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import sysconfig
import os

# Call cythonize in advance so a single module can be compiled from a single Cython
# file along with other C++ files.
import numpy as np


with open('README.md','r') as f:
    long_description = f.read()

#sources = ['ssa_translation.pyx','ssa_translation_c_w.cpp','ssa_translation_lowmem.pyx','ssa_translation_c_w_lowmem.cpp','ssa_translation_lowmem_leaky.pyx','ssa_translation_c_w_lowmem_leaky.cpp','ssa_translation_lowmem_nostats.pyx','ssa_translation_c_w_lowmem_nostats.cpp','ssa_translation_lowmem_leaky_nostats.pyx','ssa_translation_c_w_lowmem_leaky_nostats.cpp','ssa_translation_lowmem_bursting.pyx','ssa_translation_c_w_lowmem_bursting.cpp']
sources = ['ssa_translation_c_w_full.cpp','ssa_translation_lowmem.pyx','ssa_translation_c_w_lowmem.cpp',  'ssa_translation_generic_lowmem_c_w.cpp', 'ssa_translation_generic_c_w.cpp',  'ssa_trna_lowmem_c_w.cpp', 'ssa_trna_c_w.cpp']

cythonize('*.pyx', language='c++')


if not sysconfig.get_config_var('LIBS'):
    libs = sysconfig.get_config_var('LIBS')
else:
    libs = []

#try:
#
#    setup(name='SSA',ext_modules=[Extension('ssa_translation', sources, language='c++',include_dirs = ['/usr/local/include','/Library/Developer/CommandLineTools/usr/bin','/anaconda/lib/python3.6/site-packages/numpy/core/include','/anaconda/lib/python2.7/site-packages/numpy/core/include','/anaconda/lib/python3.6/site-packages/numpy/core/include',np.get_include(),'.',os.getcwd()])],cmdclass = {'build_ext': build_ext})
#except:
#    setup(name='SSA',ext_modules=cythonize([Extension('ssa_translation', sources,language='c++',extra_compile_args=[    "-stdlib=libc++"],    include_dirs = ['usr/include','/usr/local/include','/Library/Developer/CommandLineTools/usr/bin','/anaconda/lib/python3.6/site-packages/numpy/core/include',np.get_include(),'.',os.getcwd()])]),cmdclass = {'build_ext': build_ext})
#

class DependencyError(Exception):
    pass

try:
    packages=find_packages()
except:
    pass

try:
    eca = sysconfig.get_config_var('CPPFLAGS').split()
except:
    eca = []
    
# Try to find eigen directory if possible
include_list = [sysconfig.get_paths()['include'],np.get_include(),'.', os.getcwd()]

env_location = sysconfig.get_config_vars()['prefix']
print(env_location)
potential_eigens = []
for root, dirs, files in os.walk(env_location):
    for folder in dirs:
        if 'eigen3' in folder:
			
            potential_eigens.append( os.path.dirname (os.path.join(root, folder)))

print(potential_eigens)
final_dirs = []
for directory in potential_eigens:
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file == 'Dense':
                final_dirs.append(directory)

print('Located eigen3/Eigen/Dense:')
print(final_dirs)    

if len(final_dirs) > 0:
    include_list.append(final_dirs[0])
else:
    raise DependencyError("Cannot find Eigen installed on enviroment, please conda install eigen or provide a path to eigen with the setup command: python setup.py build_ext --inplace -I[path to eigen, no space after I, no brackets]")
     
setup(name='rsnapsim_ssa_cpp',
      ext_modules=[Extension('ssa_translation_lowmem', 
                sources, language='c++',
                include_dirs = include_list ,
                library_dirs = libs,
                extra_compile_args= eca)]
    ,cmdclass = {'build_ext': build_ext}
    ,author='William Raymond'
	,author_email='wsraymon@rams.colostate.edu'
    ,description= 'mRNA translation Stochastic Simulation Algorithms (SSA) for the rSNAPsim module.'
    ,version = "0.0.16"
    ,long_description = long_description
    ,long_description_content_type='text/markdown'
    ,url = 'https://github.com/MunskyGroup/rSNAPsim'
    ,install_requires = ['numpy',"Cython>=0.26.0"]
    ,packages=packages
    
    )


