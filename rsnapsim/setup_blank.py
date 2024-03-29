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

#sources_go_here
#cythonize_goes_here

if not sysconfig.get_config_var('LIBS'):
    libs = sysconfig.get_config_var('LIBS')
else:
    libs = []


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

#eigen_path_goes_here
#model_name_goes_here
#setup_goes_here
      ext_modules=[Extension(model_name, 
                sources, language='c++',
                include_dirs = include_list ,
                library_dirs = libs,
                extra_compile_args= eca)]
    ,cmdclass = {'build_ext': build_ext}
    ,author='rsnapsim'
    ,description= 'auto_generated_model_cpp_lib'
    ,version = "0.0a"
    ,long_description = 'This file was autogenerated for an rSNAPsim c++ model '
    ,long_description_content_type='text/markdown'
    ,url = ''
    ,install_requires = ['numpy',"Cython>=0.26.0"]
    ,packages=packages
    
    )


