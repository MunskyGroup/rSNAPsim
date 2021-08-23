try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
    
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import os
import sysconfig


# Call cythonize in advance so a single module can be compiled from a single Cython
# file along with other C++ files.
import numpy as np
sources = ['ssa_translation_generic.pyx','ssa_translation_generic_c_w.cpp','ssa_translation_generic_lowmem.pyx','ssa_translation_generic_lowmem_c_w.cpp']
cythonize('*.pyx', language='c++')



#setup(name='SSA',ext_modules=[Extension('ssa_translation_generic', sources, language='c++',include_dirs = ['/usr/local/include','/Library/Developer/CommandLineTools/usr/bin','/anaconda/lib/python3.6/site-packages/numpy/core/include','/anaconda/lib/python2.7/site-packages/numpy/core/include',np.get_include(),os.getcwd(),('C:\\Users\\' + os.environ['USERNAME'])])],cmdclass = {'build_ext': build_ext})

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


try:
    eca = sysconfig.get_config_var('CPPFLAGS').split()
except:
    eca = []

setup(name='SSA',
      ext_modules=[Extension('ssa_translation_generic', 
                sources, language='c++',
                include_dirs = [sysconfig.get_paths()['include'],
                                np.get_include(),
                                '.',
                                os.getcwd()],
                library_dirs = libs,
                extra_compile_args= eca)
    
    ,Extension('ssa_translation_generic_lowmem', 
                sources, language='c++',
                include_dirs = [sysconfig.get_paths()['include'],
                                np.get_include(),
                                '.',
                                os.getcwd()],
                library_dirs = libs,
                extra_compile_args= eca)]    ,cmdclass = {'build_ext': build_ext})
