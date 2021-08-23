# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 18:27:45 2021

@author: William Scott Raymond
"""

import os
import shutil
import time
import sys
import re
import subprocess
import platform
import hashlib
import numpy as np

class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class EigenMissingError(Error):
    """Exception raised for when an eigen instillation cannot be found

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class ExistenceError(Error):
    """Exception raised for when a requesting making a model that already
    exists without overwite == True

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class ModelNameError(Error):
    """Exception raised for when a requesting making a model that already
    exists without overwite == True

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message



class MisMatchedBrackets(Error):
    """Exception raised for when a requesting making a model that already
    exists without overwite == True

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class MisMatchedParenthesis(Error):
    """Exception raised for when a requesting making a model that already
    exists without overwite == True

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class GenericMetaData():
    '''
    Class that generates some generic metadata and returns it in dictionary

    meta data currently returned:
        * user
        * id
        * datetime this function was called
        * python version
        * platform this is being run on

    '''

    def __init__(self):
        self.id = ''
        self.created_at = time.strftime('%Y-%m-%d %H:%M:%S',
                                        time.localtime(time.time()))
        self.user = os.path.expanduser("~")
        self.platform = platform.platform()
        self.python_version = sys.version

    def get(self):
        '''
        generate and return a metadata dictionary for a solver object

        Returns
        -------
        dict
            a dictionary of metadata such as solution id, time ran, user,
            platform and rss version.

        '''
        return self.__dict__

class ModelFactory():
    '''
    Model factory class, this class parses, makes, and compiles custom
    Nascent Chain Tracking (NCT) TASEP models for the rsnapsim package

    the goal of this class is to provide the user a quick way to edit only
    the propensity function of a blank TASEP model such that they dont have to
    write and recompile a full c++/cython model themselves
    '''
    def __init__(self,):
        self.file_path = os.path.dirname(os.path.realpath(__file__))
        self.test_prop = '''

                            if (X_states(0,0) == 1){ // if off turn on rate
                                wn(0) = parameters[0];
                                forward_rate_matrix(0,49) = 0;  // turn on the pause secondary structure is present
                            }
                            if (X_states(0,1) == 1){ // if on turn off rate
                                if (X_full.block(0,50,1,10).sum() == 0){
                                    wn(1) = parameters[1];
                                }
                            }


                            int loc_not_free = (((X_spatial.block(0,0,1,n_ribosomes).array())-0).abs() < R).cast<int>().sum() ;
                            if (loc_not_free == 0){  // kin
                                wn(2) = parameters[2];
                            }

                            if (X_full(0,max_length-1) == 1){  // kout
                                forward_rate_matrix(0,max_length-1) = 0;
                                wn(3) = parameters[3];
                            }
                        '''

        self.eigen_paths = self.find_eigen_path()
        print('eigen instillation found...')

        self.reserved_model_names = ['build', 'model_maker_cpp', 'models',
                                     'rsnapsim_model_maker']
        self.find_models()




    def find_models(self):
        '''
        function to find all models (compiled or failed to compile)

        This function walks the directory and finds all names
        used in ./models/*/

        Returns
        -------
        models : list of str
            strings of all model names found (folders in ./models).

        '''
        fpath = self.file_path
        build_dirs = [y for y in [x[0] for x in os.walk(fpath)] if 'build' in y]
        split_paths = []
        for pathstr in build_dirs:
            path = os.path.normpath(pathstr)
            split_paths = split_paths + [path.split(os.sep), ]
        models = []
        for pathstr in split_paths:
            try:
                if pathstr[-2] == 'build':
                    if pathstr[-3] not in models:
                        models.append(pathstr[-3])
            except:
                pass

        for model in models:
            if model in self.reserved_model_names:
                models.remove(model)
        self.available_models = models

        return models


    def find_eigen_path(self):
        '''
        find an instillation and fpath of eigen in the current python
        enviroment to use to compile models.

        Raises
        ------
        EigenMissingError
            Raised when Eigen is not installed, we cannot compile c++ models.

        Returns
        -------
        eigen_paths : list of str
            list of all file paths to eigen instillations.

        '''

        paths = sys.path
        potential_paths = []  #check the <env>/lib/ folder
        for path in paths:

            if path[-3:] == 'lib':
                potential_paths.append(path)
        eigen_paths = []
        for path in potential_paths:  #in each of these try to find an eigen instillation
            base, _ = os.path.split(path)
            if os.path.exists(os.path.join(base, 'Library',
                                           'include', 'eigen3')):

                eigen_paths.append(os.path.join(base, 'Library',
                                                'include', ''))
        if len(eigen_paths) == 0:
            raise EigenMissingError('Eigen is missing, please provide a path'\
                                    ' or if using a conda instillation, use'\
                                        ' conda install eigen')
        return eigen_paths


    def generate_metadata(self, prop_str, model_name):
        '''
        Generate some meta data to write to the __init__.py of the model

        Parameters
        ----------
        prop_str : str
            propensity string used to generate the model.
        model_name : str
            model name being generated.

        Returns
        -------
        metadata : dict
            dictionary of the metadata.
        meta_str : str
            string version of the metadata
            (this is written to the top of the __init.py__).

        '''
        metadata = GenericMetaData().get() #get all the metadata
        # convert the prop string to a hash ID
        metadata['id'] = hashlib.sha256(prop_str.encode('ascii')).hexdigest()
        # convert to a metadata string
        meta_str = ''
        meta_str += 'model name : ' + model_name + '\n'
        meta_str += 'model ID : ' + metadata['id'] + '\n'
        meta_str += 'files created at : ' + metadata['created_at']+ '\n'
        meta_str += 'platform : ' + metadata['platform']+ '\n'
        meta_str += 'python version : ' + metadata['python_version']+ '\n'

        return metadata, meta_str


    def edit_cpp_files(self, model_name, cpp_file_path, prop_str):
        '''
        edit the c++ files to add the parsed propensity strings

        Parameters
        ----------
        model_name : str
            model.
        cpp_file_path : str
            path to the c++ file to edit.
        prop_str : str
            parsed propensity function string.

        Returns
        -------
        None.

        '''
        strfrnt = '// Start of autogenerated propensity function \n'
        strback = '\n// End of autogenerated propensity function \n'
        prop_str = strfrnt + prop_str + strback
        fstr = ''
        with open(cpp_file_path, 'r') as fname:
            for line in fname:
                if "//INSERT_GENERATED_PROPENSITY_HERE" in line:
                    print('replacing_prop...')
                    line = line.replace("//INSERT_GENERATED_PROPENSITY_HERE",
                                        prop_str)

                fstr += line
        with open(cpp_file_path, 'w') as fname:
            fname.write(fstr)

    def edit_pyx_files(self, model_name, min_length, rules_str,
                       original_rules_str, pyx_file_path):
        '''
        edit the pyx files for the model being made

        Parameters
        ----------
        model_name : str
            name of the model being made.
        min_length : int
            minimum index length in the rules strings provided.
        rules_str : str
            parsed rules string (in c++).
        original_rules_str : str
            original rules string that was passed (python-like).
        pyx_file_path : str
            path to the pyx file to edit (model_name.pyx).

        Returns
        -------
        None.

        '''

        #edit each string to edit
        cdef_str = "cdef extern from 'model_%s.h':"%model_name
        rules_str_r = "rules_str = '''%s"%rules_str + "'''"
        original_rules_str_r = "original_rules_str = ''' %s"%original_rules_str + "'''"
        min_length_str = "min_length = %i"%min_length

        #edit all the lines in the .pyx file
        fstr = ''
        with open(pyx_file_path, 'r') as fname:
            for line in fname:

                if "#cdef goes here" in line:
                    line = line.replace("#cdef goes here", cdef_str)
                if "#min length goes here" in line:
                    line = line.replace("#min length goes here",
                                        min_length_str)
                if "#parsed rules go here" in line:
                    line = line.replace("#parsed rules go here",
                                        rules_str_r)
                if "#original rules go here" in line:
                    line = line.replace("#original rules go here",
                                        original_rules_str_r)
                fstr += line

        with open(pyx_file_path, 'w') as fname:
            fname.write(fstr)


    def edit_setup_files(self, model_name, setup_file_path, eigen_path=''):
        '''
        Edit the setup.py files for the model being made

        Parameters
        ----------
        model_name : str
            name of the model being generated.
        setup_file_path : str
            path of the new setup file to rewrite.
        eigen_path : str, optional
            path to eigen instillation, this has to be added
            to the setup include list. The default is ''.

        Returns
        -------
        None.

        '''
        if len(eigen_path) == 0:
            eigen_path = self.eigen_paths[0]

        #eigen_path_str = r'{}'.format(eigen_path) #convert to raw string

    #generate all the lines that need to be edited in the setup_model_name.py
        #we need to replace a couple things:
            #the source names (.cpp and .pyx)
            #the cythonize command
            #the setup = "model_name"
            #eigen include path
            #model name string
        src_str = "sources = ['%s.pyx', 'model_%s.cpp']"%(model_name,
                                                          model_name)
        cythonize_str = "cythonize('%s.pyx')"%model_name
        setup_str = "setup(name='%s',"%model_name
        eigen_str = r"include_list = include_list + [%s ,]"%repr(eigen_path)
        model_name_str = "model_name = '%s'"%model_name

        #go through and edit each line to edit
        fstr = ''
        with open(setup_file_path, 'r') as fname:

            for line in fname:

                if "#sources_go_here" in line:
                    line = line.replace("#sources_go_here", src_str)
                if "#cythonize_goes_here" in line:
                    line = line.replace("#cythonize_goes_here", cythonize_str)
                if "#setup_goes_here" in line:
                    line = line.replace("#setup_goes_here", setup_str)
                if "#eigen_path_goes_here" in line:
                    line = line.replace("#eigen_path_goes_here", eigen_str)
                if "#model_name_goes_here" in line:
                    line = line.replace('#model_name_goes_here',
                                        model_name_str)
                fstr += line
        with open(setup_file_path, 'w') as fname:
            fname.write(fstr)


    def make_init_file(self, folder_path, model_name, rules):
        '''
        make the init_file for the model maker

        Parameters
        ----------
        folder_path : str
            The location of where these files are being made
            (where the model maker is installed).
        model_name : str
            name of the model being made.
        rules : str
            additional rule string, this is converted to a hash
            for the model id metadata

        Returns
        -------
        None.

        '''
        #make metadata
        _, meta_str = self.generate_metadata(rules, model_name)

        init_str = ''
        init_str += "'''\n"
        init_str += meta_str
        init_str += "'''\n\n"

        #add the model name to the string to write
        init_str += 'from . import %s'%model_name

        #make the new init file
        initpath = os.path.join(folder_path, 'models',
                                model_name, '__init__.py')
        with open(initpath, "w") as fname:
            fname.write(init_str)



    def generate_model_files(self, model_name, overwrite=False,
                             eigen_path='', propensity_str=None,
                             original_rules_str=None, min_length=None):
        '''
        Generate all model files required to compile a new model and put
        them in a new model folder

        Parameters
        ----------
        model_name : str
            name of the model to generate.
        overwrite : bool, optional
            overwrite any existing files in ./models/<model_name>. The default is False.
        eigen_path : str, optional
            location of an eigen instillation. The default is ''.
        propensity_str : str, optional
            parsed additional rules to paste into the c++ files. The default is None.
        original_rules_str : str, optional
            original unparsed additional rules to paste into the c++ files.
            The default is None.
        min_length : int, optional
            the maximum index for X_lattice detected (minimum length for a tasep).
            The default is None.

        Raises
        ------
        ExistenceError
            Model folder already exists and overwrite == False.

        Returns
        -------
        None.

        '''
        setup_name = 'setup_' + model_name
        model_c_file_name = 'model_' + model_name
        model_h_file_name = 'model_' + model_name

       # model_base_found = False
        #setup_file_found = False
        #model_c_file_found = False
        #model_h_file_found = False

        for root, dirs, files in os.walk('.'):
            for fname in files:
                old_name = os.path.join(os.path.abspath(root), fname)
                base, extension = os.path.splitext(fname)

                if extension == '.pyx':
                    if base == 'blank_model':
                        folder_path = os.path.abspath(root)
                        model_file_path = old_name
                        #model_base_found = True

                        new_model_file = os.path.join(os.path.abspath(root),
                                                      model_name + extension)

                if extension == '.py':
                    if base == 'setup_blank':
                        setup_file_path = old_name
                        #setup_file_found = True

                        new_setup_file = os.path.join(os.path.abspath(root),
                                                      setup_name + extension)


                if extension == '.cpp':
                    if base == 'blank_model':
                        model_c_file_path = old_name
                        #model_c_file_found = True
                        new_model_c_file = os.path.join(os.path.abspath(root),
                                                        model_c_file_name + extension)

                if extension == '.h':
                    if base == 'blank_model':
                        model_h_file_path = old_name
                        #model_h_file_found = True
                        new_model_h_file = os.path.join(os.path.abspath(root),
                                                        model_h_file_name + extension)

        files_exist = np.sum(np.array([os.path.exists(new_model_file),
                                       os.path.exists(new_setup_file),
                                       os.path.exists(new_model_c_file),
                                       os.path.exists(new_model_h_file)]))

        folder_exist = os.path.exists(os.path.join(folder_path, 'models',
                                                   model_name))
        if folder_exist == 0:
            os.makedirs(os.path.join(folder_path, 'models', model_name))
        else:
            if overwrite:
                shutil.rmtree(os.path.join(folder_path, 'models', model_name))
                os.makedirs(os.path.join(folder_path, 'models', model_name))
            else:
                message = 'model folder by the name "%s" exists at %s already'\
                    ' exists, if you wish to overwrite '\
                        'use overwrite == True'% (model_name, folder_path)
                raise ExistenceError(message)

        if files_exist == 0:
            shutil.copy(model_file_path, new_model_file)
            shutil.copy(setup_file_path, new_setup_file)
            shutil.copy(model_c_file_path, new_model_c_file)
            shutil.copy(model_h_file_path, new_model_h_file)


        else:
            if overwrite:
                shutil.copy(model_file_path, new_model_file)
                shutil.copy(setup_file_path, new_setup_file)
                shutil.copy(model_c_file_path, new_model_c_file)
                shutil.copy(model_h_file_path, new_model_h_file)


            else:
                message = 'model files: by the name "%s" exists at %s already'\
                    ' exists, if you wish to overwrite '\
                        'use overwrite == True'% (model_name, new_model_file)
                raise ExistenceError(message)
        time.sleep(.1)
        self.edit_setup_files(model_name, new_setup_file,
                              eigen_path=eigen_path)
        self.edit_pyx_files(model_name, min_length, propensity_str,
                            original_rules_str, new_model_file)


        self.edit_cpp_files(model_name, new_model_c_file, propensity_str)
        shutil.move(new_model_file, os.path.join(folder_path,
                                                 'models', model_name))
        shutil.move(new_setup_file, os.path.join(folder_path,
                                                 'models', model_name))
        shutil.move(new_model_c_file, os.path.join(folder_path,
                                                   'models', model_name))
        shutil.move(new_model_h_file, os.path.join(folder_path,
                                                   'models', model_name))
        self.make_init_file(folder_path, model_name, propensity_str)

        #attempt to compile


    def compile_model(self, model_name, rules=None, overwrite=False,
                      eigen_path='', verbose=True):
        '''
        Compile a custom Nascent Chain tracking (NCT) TASEP model.

        This command takes a set of rules and will replace the propensity function
        inside a set of "blank" c++/pyx files, shuttle these files to a new model folder
        and attempt to compile this model.

        The model takes the following:
            ./blank_model.cpp
            ./blank_model.pyx
            ./setup_blank.py
            ./blank_model.h

        and will make a new folder containing:
            ./models/<model_name>/__init__.py
            ./models/<model_name>/model_<model_name>.cpp
            ./models/<model_name>/model_<model_name>.h
            ./models/<model_name>/<model_name>.pyx
            ./models/<model_name>/setup_<model_name>.py

        and will then run "python setup_<model_name>.py build_ext --inplace"
        and attempt to compile the model. If successful the model can then be imported
        via:

            from models import <model_name>

        Parameters
        ----------
        model_name : str
            name of the model to try to compile.
        rules : str, optional
            The additional rules containing the propensity function
            to replace the one in the c++ source. The default is None.
        overwrite : bool, optional
            option to overwrite an existing model, otherwise if the model
            name is already in use, it will not delete or recompile the model
            folder. The default is False.
        eigen_path : str, optional
            location to the path of the eigen instillation
            https://eigen.tuxfamily.org/. The default is ''.
        verbose : bool, optional
            print out the steps of compiling and errors if it goes
            wrong. The default is True.

        Raises
        ------
        ModelNameError
            Raised if a model name is already in use and overwrite == False.

        Returns
        -------
        None.

        '''
        #if the user did not provide a rule string, use a dummy one
        if rules == None:
            rules = self.test_prop
        if model_name in self.reserved_model_names:
            msg = 'The model name requested is a reserved keyword and'\
                ' cannot be used to build the model files, please rename'\
                    ' the model.'
            raise ModelNameError(msg)

        # convert the rules passed
        if verbose:
            print('parsing rules...')
            parsed_rules = RuleConverter().parse_rules(rules)
            min_length = RuleConverter().get_min_length(rules)

        if verbose:
            print('generating model files...')

        #generate all the appropriate files
        cwd = os.getcwd()
        os.chdir(self.file_path) #change to the cwd where this file is stored
        self.generate_model_files(model_name, overwrite=overwrite,
                                  eigen_path=eigen_path,
                                  propensity_str=parsed_rules,
                                  min_length=min_length,
                                  original_rules_str=rules)
        #change dir to where the new files are
        os.chdir(os.path.join('models', model_name))
        #attempt to compile
        if verbose:
            print('compiling model...')

        if verbose:
            compile_run = subprocess.run(["python", "setup_%s.py"%model_name,
                                          "build_ext", "--inplace"],
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)

        else:
            compile_run = subprocess.run(["python", "setup_%s.py"%model_name,
                                          "build_ext", "--inplace"])

        # check if the model compiled
        if compile_run.returncode != 0:
            if verbose:
                print('compilation failed...')
                print(compile_run.stdout.decode('utf-8'))
                print(compile_run.stderr.decode('utf-8'))
        else:
            if verbose:
                print('model compiled!')
        os.chdir(cwd)




class RuleConverter():
    '''
    Class that converts the python-like rules string into c++ to paste into
    the blank model files



    Full example:

        additional_rules = ```

        #propensity function for hairpin model
        int hairpin_location = cast_to_int(parameters[4]) # location of the hairpin
        if state[0] == 1:
            #if the hairpin formed, state = off
            wn[0] = parameters[0]
            step[ hairpin_location-1  ] = 0

        if state[1] == 1:
            #if the hairpin formed, state = on
            if sum(X[hairpin_location:hairpin_location+10  ]) < 1:
                wn[1] = parameters[1]

        if free[0]:
            #if the front location is free, allow kin
            wn[2] = parameters[2]

        if X[149] == 1:
            # if its in the final location, add kout
            step[149] = 0
            wn[3] = parameters[3]



        ```

        >>>parsed_rules = RuleConverter().parse_rules(additional_rules)


        parsed_rules = ```

            //propensity function for hairpin model
            int hairpin_location = (static_cast<int>(parameters[4] < 0 ? parameters[4] - 0.5 : parameters[4] + 0.5)); // location of the hairpin
            if( X_states(0) == 1){
            //if the hairpin formed, state = off
                wn[0] = parameters[0];
                forward_rate_matrix( hairpin_location-1  ) = 0;
                }

            if( X_states(1) == 1){
            //if the hairpin formed, state = on
                if( X_full.block(0,hairpin_location,1,hairpin_location+10  ).sum() < 1){
                    wn[1] = parameters[1];
                    }
                }

            if(   ((((X_spatial.block(0,0,1,n_ribosomes).array())-0).abs() < R).cast<int>().sum() == 0) ){
            //if the front location is free, allow kin
                wn[2] = parameters[2];
                }

            if( X_full(149) == 1){
            // if its in the final location, add kout
                forward_rate_matrix(149) = 0;
                wn[3] = parameters[3];
                }



        ```




    '''
    def __init__(self):

        #functions denoted by () to find and convert
        self.keyfuns = {'sum': '.sum()',
                        'cast_to_int':'static_cast<int>(%s < 0 ? %s - 0.5 :'\
                                      ' %s + 0.5)'}

        #keywords with [] to find and convert
        self.keywords = {'state': 'X_states',
                         'wn':'wn',
                         'step':'forward_rate_matrix',
                         'X':'X_full',
                         'occupied':'  ((((X_spatial.block(0,0,1,n_ribosomes)'\
                                    '.array())-%s).abs() < R).cast<int>().'\
                                    'sum() != 0) ',
                         'free':'  ((((X_spatial.block(0,0,1,n_ribosomes).array()'\
                                ')-%s).abs() < R).cast<int>().sum() == 0) ',}

        #if operators to parse
        self.if_operators = ['<', '>', '==', '!=', '>=', '<=']


    '''
    these are not ready yet
    def make_in(self,in_loc, reaction_num, par_num):

        in_str = 'int loc_not_free = (((X_spatial.block(0,0,1,n_ribosomes).array())-0).abs() < R).cast<int>().sum() ;\nif (loc_not_free == 0){  // kin\n'


        in_str = in_str + 'wn(%i) = parameters[%i];'%(reaction_num, par_num)
        in_str = in_str + '}'

        return in_str

    def make_out(self, out_loc, reaction_num, par_num, run_through):

        out_str = 'if (X_full(0,max_length-1) == 1){  // kout\n '
        if run_through:
            out_str = 'forward_rate_matrix(0,max_length-1) = 0;\n'

        out_str = out_str + 'wn(%i) = parameters[%i];'%(reaction_num, par_num)
        out_str = out_str + '}'

        return out_str

    '''

    @staticmethod
    def check_brackets(string):

        '''
        Check if there are mismatched parenthesis in a string,
        returns none if it doesnt error out.

        Parameters
        ----------
        string : str
            string to check.

        Raises
        ------
        MisMatchedParenthesis
            Error if mismatched (.
        MisMatchedBrackets
            Error if mismatched [.

        Returns
        -------
        None.

        '''
        for bracket in ['[', '(']:
            restr = {'[': '\[|\]', '(':'\(|\)'}[bracket]
            par_locations = [(m.start(0), m.end(0)) for m in re.finditer(restr, string)]
            #exception here if unclosed
            total_pars = int(len(par_locations))
            if total_pars%2 != 0:
                if bracket == '(':
                    raise MisMatchedParenthesis('MisMatched Parenthes'/
                                                'is () in the stri'/
                                                'ng: %s'% string)
                if bracket == '[':
                    raise  MisMatchedBrackets('MisMatched Brackets [] '/
                                              'in the string: %s'% string)


    @staticmethod
    def get_parenthesis_dict(string):
        '''
        Get all matching pairs of parenthesis from a string, error if mismatched

        Parameters
        ----------
        string : str
        string to find the parenthesis pairs of

        Raises
        ------
        MisMatchedParenthesis
        Error if there is an unbalanced ( or )

        Returns
        -------
        list of tuples
        list of tuples of matches from a string. Example:
        'hello(3) + hello(5)' would return
        [(5,7), (16,18)]

        '''
        istart = []
        bracket_dict = {}

        for i, char in enumerate(string):
            if char == '(':
                istart.append(i+1)
            if char == ')':
                try:
                    bracket_dict[istart.pop()] = i
                except:
                    raise MisMatchedParenthesis('Missing "("  with closed '\
                                                ' parenthesis at %i'%i)

        if len(istart) != 0:  # still a ( that cant be matched
            raise MisMatchedParenthesis('Mismatched ")" ')
        return list(bracket_dict.items())

    @staticmethod
    def get_square_bracket_dict(string):
        '''
        Get all matching pairs of square brackets from a string, error if mismatched

        Parameters
        ----------
        string : str
        string to find the square bracket pairs of

        Raises
        ------
        MisMatchedBrackets
        Error if there is an unbalanced [ or ]

        Returns
        -------
        list of tuples
        list of tuples of matches from a string. Example:
        'hello[3] + hello[5]' would return
        [(5,7), (16,18)]

        '''
        istart = []
        bracket_dict = {}

        for i, char in enumerate(string):
            if char == '[':
                istart.append(i+1)
            if char == ']':
                try:
                    bracket_dict[istart.pop()] = i
                except:
                    raise MisMatchedBrackets('Missing "["  with closed '\
                                             'parenthesis at %i'%i)

        if len(istart) != 0:  # still a ( that cant be matched
            raise MisMatchedBrackets('Mismatched "]" ')
        return list(bracket_dict.items())

    @staticmethod
    def invert_pair_locations(pair_locations, length):
        '''
        invert a set of pair locations of a string so we have all indexes of all substrings

        example:

            s = 'abc[defg]hijk'
            pair_locations = [(3,8), ] #pair strings of [*]
            length = 12

            result: [(0,3),(8,11) ] #inverted pair strings!

        Parameters
        ----------
        pair_locations : list of tuples of ints
            substring pair locations of brackets in a string in a format [(x,y), (z,w)...] .
        length : int
            total length of the string the pair locations come from.

        Returns
        -------
        inverse_pair_locations : list of tuples of ints
            inverted substring pair locations of brackets (indexes of substrings outside brackets).

        '''
        order = ''
        inverse_pair_locations = []

        if pair_locations[0][0] == 0:
            order += 'm'
        else:
            order += 'i'
            inverse_pair_locations = inverse_pair_locations + [(0, pair_locations[0][0]),]
        for i in range(len(pair_locations)-1):
            inverse_pair_locations = inverse_pair_locations + [(pair_locations[i][1], pair_locations[i+1][0]),]

        if pair_locations[-1][1] != length:
            inverse_pair_locations = inverse_pair_locations + [(pair_locations[-1][1], length),]
            order += 'i'
        else:
            order += 'm'
        return inverse_pair_locations

    @staticmethod
    def get_orders_of_substrings(pair_locations, inverse_pair_locations,):
        '''
        given all pair locations inside and outside of brackets, get the order
        they should be reconstructed in.

        example:

         s = 'abc[defg]hijk'
         pair_locations = [(3,8), ] #pair strings of [*]
         length = 12
         inverse_locations = [(0,3),(8,11)]

         result: [(0,3),(3,8),(8,11)],  'imi'

         imi for inverse, match, inverse

        Parameters
        ----------
        pair_locations : list of tuples of ints
            substring pair locations of brackets in a string in a format [(x,y), (z,w)...].
        inverse_pair_locations : list of tuples of ints
            inverted substring pair locations of brackets (indexes of substrings outside brackets).

        Returns
        -------
        combined_strings : list of tuples of ints
            all pair locations.
        orders : str
            string denoting which pairs are outside ('i') and inside brackets ('m').

        '''
        n_substrings = len(pair_locations) + len(inverse_pair_locations)
        k = 0
        m = 0
        n = 0
        orders = ''
        combined_strings = []
        for i in range(n_substrings):
            if m < len(pair_locations):
                if pair_locations[m][0] == k:
                    k = pair_locations[m][1]
                    combined_strings = combined_strings + [pair_locations[m],]
                    orders += 'm'
                    m += 1
            if n < len(inverse_pair_locations):
                if inverse_pair_locations[n][0] == k:
                    k = inverse_pair_locations[n][1]
                    combined_strings = combined_strings + [inverse_pair_locations[n],]
                    orders += 'i'
                    n += 1
        return combined_strings, orders


    def get_min_length(self, rule_string):
        '''
        Function that gets the minimum length construct needed for a model.

        Parameters
        ----------
        rule_string : str
            The full additional rules string

        Returns
        -------
        int
            Maximum value of an index for the stepping matrix or
            X_lattice matrix. This corresponds to the minimum length
            construct needed by the model when passing values to it.

        '''

        #remove all comments, trailing_whitespace / leading whitespace
        rules_split, _, _, _, _, _, _ = self.get_valid_rules(rule_string)
        all_indexes = []
        for i in range(len(rules_split)):
            a = [(m.start(0), m.end(0)) for m in re.finditer('step\[(.*?)\]',
                                                             rules_split[i])]
            b = [(m.start(0), m.end(0)) for m in re.finditer('X\[(.*?)\]',
                                                             rules_split[i])]

            a = [re.findall(r'\d+', rules_split[i][x[0]:x[1]]) for x in a]
            b = [re.findall(r'\d+', rules_split[i][x[0]:x[1]]) for x in b]
            all_indexes = all_indexes + a + b
        all_indexes = [item for sublist in all_indexes for item in sublist]
        return max([int(x) for x in all_indexes])


    def get_valid_rules(self, rule_string):
        '''
        master function that parses out stand alone comments and
        trailing comment along with pairs of line orders.

        Parameters
        ----------
        rule_string : str
            the string to parse.

        Returns
        -------
        rules_split : list of str
            list of split up lines without trailing comments.
        rules_line_ids : list of ints
            ids of the order of the lines that were rules.
        comments : list of str
            list of stand alone comments (denoted by #).
        comment_line_ids : list of ints
            ids of the order of comments.
        trailing_comments : list of strs
            list of trailing comment strings.
        trailing_comment_ids : list of ints
            ids of the lines containing trailing comments.
        total_lines : int
            N of total lines.

        '''
        rules_split_all = [x for x in rule_string.split('\n') if len(x.lstrip()) != 0]
        total_lines = len(rules_split_all)
        rules_line_ids = [x[0] for x in enumerate(rules_split_all)]

        comments = [x for x in rules_split_all if x.lstrip()[0] == '#']
        comment_line_ids = [x[0] for x in enumerate(rules_split_all) if x[1].lstrip()[0] == '#']

        trailing_comments = []
        trailing_comment_ids = []
        for i in range(len(rules_split_all)):
            if i not in comment_line_ids:
                if '#' in rules_split_all[i]:
                    trailing_comments.append(rules_split_all[i][rules_split_all[i].find('#'):])
                    trailing_comment_ids = trailing_comment_ids + [i,]

        for idx in comment_line_ids:
            rules_line_ids.remove(idx)
        rules_split = []
        for idx in rules_line_ids:
            if idx in trailing_comment_ids:
                rules_split = rules_split + [rules_split_all[idx][:rules_split_all[idx].find('#')],]
            else:
                rules_split = rules_split + [rules_split_all[idx],]

        return rules_split, rules_line_ids, comments, comment_line_ids, trailing_comments, trailing_comment_ids, total_lines

    def parse_rules(self, rule_string):
        '''
        Master function that attempts to convert the given rule set to a
        parsed ruleset

        Parameters
        ----------
        rule_string : str
            the rule string to convert to c++.

        Returns
        -------
        parsed_rule_full_str : str
            full rule string parsed to c++.

        '''

        #split up the rules by lines, remove stand alone line comments
        #(lines starting with #) and trailing commments, (anything after a #)
        rules_split, rules_line_ids, comments, comment_line_ids, trailing_comments, trailing_comment_ids, total_lines = self.get_valid_rules(rule_string)

        indents = [len(x) - len(x.lstrip()) for x in rules_split] #get the indent levels

        #outer levels (where indents = 0)
        #inner levels (where indents != 0)
        outer_levels = np.where(np.array(indents) == 0)[0]
        inner_levels = np.where(np.array(indents) != 0)[0]

        #convert all of these to consistent levels and parse the tabs
        max_outer = np.max(outer_levels)
        max_inner = np.max(inner_levels)
        max_total = np.max([max_outer, max_inner])
        outer_levels = np.append(outer_levels, max_total + 1)

        levels = []
        for i in range(0, len(outer_levels)-1):
            tosort = [indents[x] for x in range(outer_levels[i]+1, outer_levels[i+1])]
            orderkey = list(set(tosort))
            orderkey.sort()
            levels.append([0,])
            levels.append([orderkey.index(tosort[x]) + 1 for x in range(len(tosort))])

        levels.append([0,])

        #strip out all whitespace now that we have the indentation levels
        levels = [item for sublist in levels for item in sublist]
        rules_stripped = [x.lstrip().rstrip() for x in rules_split]
        #parse all the comments (just convert the front # to a //)
        comments_stripped = [x.lstrip().rstrip() for x in comments]
        trailing_comments_stripped = [x.lstrip().rstrip() for x in trailing_comments]
        comments_parsed = ['//' + x[1:] for x in comments_stripped]
        trailing_comments_parsed = ['//' + x[1:] for x in trailing_comments_stripped]

        #now for each valid rule line, convert those as well
        rules_parsed = []
        for i in range(len(rules_stripped)):
            #check if its an if statement to parse
            if rules_stripped[i][:2] == 'if':
                rules_parsed = rules_parsed + [self.convert_funs(
                    self.convert_if(rules_stripped[i])),]

            elif rules_stripped[i][:4] == 'for ':
                rules_parsed = rules_parsed + [self.convert_funs(
                    self.convert_for(rules_stripped[i])),]
            else:
                rules_parsed = rules_parsed + [self.convert_funs(
                    self.convert_assignment(rules_stripped[i])),]

        #step through based on the indentation levels and add the appropriate
        #{} for the c++
        parsed_rule_full_str = ''
        current_level = 0
        for i in range(total_lines):

            if i in rules_line_ids:
                j = rules_line_ids.index(i)

                front_str = '    '*levels[j]
                current_level = levels[j]
                if i in trailing_comment_ids:
                    m = trailing_comment_ids.index(i)
                    parsed_rule_full_str += front_str +  rules_parsed[j] + ' ' + trailing_comments_parsed[m] +  '\n'
                else:
                    parsed_rule_full_str += front_str +  rules_parsed[j] + '\n'
                add_newline = False
                n = levels[j] -levels[j+1] + 1
                for k in range(0, levels[j] -levels[j+1]):
                    add_newline = True
                    n -= 1
                    parsed_rule_full_str += '    '*n + '}\n'

                if add_newline:
                    parsed_rule_full_str += '\n'

            if i in comment_line_ids:
                j = comment_line_ids.index(i)
                front_str = '    '*current_level
                parsed_rule_full_str += front_str +  comments_parsed[j] + '\n'

        #return the final parsed c++ rules
        return parsed_rule_full_str



    def convert_var(self, var, key):
        '''
        convert a keyword variable to its c++ counterpart

        for example:

        X[0:34] will be converted to X_full.block(0,0,1,34)

        Parameters
        ----------
        var : str
        full string to convert.
        key : str
        the python-like key to convert.

        Returns
        -------
        newvar : str
        the variable string converted to c++.


        '''
        newvar = var[:]

        index_str = newvar[len(key)+1:-1] #remove key and brackets

        if key in ['state', 'step', 'X']:
            if ':' not in index_str:
                new_indexes = '(' + index_str +')'

            else:
                j, k = str(index_str.split(':')[0]), str(index_str.split(':')[1])
                new_indexes = '.block(0,%s,1,%s)'%(j, k)

            newvar = self.keywords[key] + new_indexes

        if key in ['free', 'occupied']:
            newvar = self.keywords[key]%index_str

        return newvar


    def convert_keyfuns(self, string, bracket='('):
        '''
        convert key funs (sets bracket to '(' for keystring func )
        '''
        key_list = list(self.keyfuns.keys())
        return self.convert_keystrings(key_list, string, bracket, fun=True)

    def convert_keywords(self, string, bracket='['):
        '''
        convert key words (sets bracket to '[' for keystring func )
        '''
        key_list = list(self.keywords.keys())
        return self.convert_keystrings(key_list, string, bracket)


    def convert_range(self, string):
        '''
        gets the substrings x, y, z from "range(x,y,z)"

        Parameters
        ----------
        string : str
            string to convert to a for loop in c++.

        Returns
        -------
        parsed substrings.

        '''
        strcopy = string[:]
        iters = strcopy[6:-1].split(',')
        if len(iters) == 2:
            range_iter = '1'
            range_start, range_stop = iters
        else:
            range_start, range_stop, range_iter = iters

        return range_start, range_stop, range_iter

    def convert_for(self, string):
        '''
        convert the statement in the format

            "for i in range(x,y,z):"

        to

            for(int i = x; i < y; i+=z){

        Parameters
        ----------
        string : str
            for loop start to convert.

        Returns
        -------
        forstr : str
            for loop statement converted to c++.

        '''
        key = 'range'
        self.check_brackets(string)
        #get pair locations
        pair_locations = self.get_parenthesis_dict(string)
        statements_to_edit = [string[x[0]-len(key)-1:x[1]+1] for x in pair_locations]

        statements_to_edit2 = []
        bracket_to_edit = []
        for i in range(len(statements_to_edit)):

            if key+'(' == statements_to_edit[i][:len(key)+1]:

                statements_to_edit2 = statements_to_edit2 + [statements_to_edit[i],]
                bracket_to_edit = bracket_to_edit + [(pair_locations[i][0] -len(key)-1, pair_locations[i][1]+1),]

            inverse_brackets = self.invert_pair_locations(bracket_to_edit, len(string))
            #get the order of these substring pairs
            substring_pairs, _ = self.get_orders_of_substrings(bracket_to_edit, inverse_brackets)

            #loop through these pairs and rebuild the string, converting the key + bracket function as appropriate
            all_substrings = [string[x[0]:x[1]] for x in  substring_pairs]
            for i in range(len(all_substrings)):

                if key + '(' == all_substrings[i][:len(key)+1]:
                    substring1, substring2, substring3 = self.convert_range(all_substrings[i])


        itername = string[4:string.index(' in ')]
        forstr = 'for(int %s = %s; %s < %s; %s+= %s){'%(itername,
                                                        substring1,
                                                        itername,
                                                        substring2,
                                                        itername,
                                                        substring3)

        return forstr

    def convert_keystrings(self, key_list, string, bracket, fun=False):
        '''
        take any string that has a key word (X, free, occupied) or
        a key function (sum, cast_to_int) + a desired bracket, and then
        convert them to their c++ equivelents and return the string.

        Parameters
        ----------
        key_list : list
        list of str keys to check for.
        string : str
        original string to convert.
        bracket : str
        which bracket to find and convert, ( or [.
        fun : bool, optional
        are we parsing functions? if so we will use edit_statement.
        The default is False.

        Returns
        -------
        string : str
        The converted string.

        '''


        for key in key_list:

            if key + bracket in string:

                self.check_brackets(string)  #check if the brackets are matching
                #if they match get the matching pairs to split by
                if bracket == '[':
                    pair_locations = self.get_square_bracket_dict(string)
                if bracket == '(':
                    pair_locations = self.get_parenthesis_dict(string)

                #split out the statements split by the desired bracket
                statements_to_edit = [string[x[0]-len(key)-1:x[1]+1] for x in pair_locations]

                # find which of these statements contain key + bracket as the start
                statements_to_edit2 = []
                bracket_to_edit = []
                for i in range(len(statements_to_edit)):

                    if key+bracket == statements_to_edit[i][:len(key)+1]:

                        statements_to_edit2 = statements_to_edit2 + [statements_to_edit[i],]
                        bracket_to_edit = bracket_to_edit + [(pair_locations[i][0] -len(key)-1, pair_locations[i][1]+1),]

                #for each of those matching key + bracket, invert the pair locations
                inverse_brackets = self.invert_pair_locations(bracket_to_edit, len(string))
                #get the order of these substring pairs
                substring_pairs, _ = self.get_orders_of_substrings(bracket_to_edit, inverse_brackets)

                #loop through these pairs and rebuild the string, converting the key + bracket function as appropriate
                all_substrings = [string[x[0]:x[1]] for x in  substring_pairs]
                string_to_rebuild = []
                for i in range(len(all_substrings)):
                    #split by the statement
                    if key + bracket == all_substrings[i][:len(key)+1]:
                        if fun:
                            string_to_rebuild = string_to_rebuild + [self.edit_statement(all_substrings[i], key),]
                        else:

                            string_to_rebuild = string_to_rebuild + [self.convert_var(all_substrings[i], key),]
                    else:
                        string_to_rebuild = string_to_rebuild + [all_substrings[i], ]

                string = ''.join(string_to_rebuild) #rebuild the full string
        return string

    def convert_assignment(self, statement_to_parse):
        '''
        convert an assignment such as

        state[0] = 1 to X_state(0) = 1;


        Parameters
        ----------
        statement_to_parse : str
            string to convert.

        Returns
        -------
        statement_copy : TYPE
            the converted assignment statement.

        '''

        statement_copy = statement_to_parse[:] #copy input string
        statement_copy = self.convert_keywords(statement_copy, bracket='[')

        if statement_copy[-1] != ';':
            statement_copy = statement_copy + ';'

        return statement_copy



    def convert_if(self, statement_to_parse):
        '''
        convert an if statement to c++

        Parameters
        ----------
        statement_to_parse : str
            if statment to convert.

        Returns
        -------
        statement_copy : str
            if statement in c++ (before function conversions).

        '''

        statement_copy = statement_to_parse[:] #copy input string
        statement_copy = statement_copy.replace(' and ', ' && ')
        statement_copy = statement_copy.replace(' or ', ' || ')



        statement_copy = self.convert_keywords(statement_copy)
        found_operator = False
        for operator in ['||', '&&']:
            found_operator =  True
            statement_copy = statement_copy.replace(operator,
                                                    ')' + operator + '(')

        statement_copy = 'if(' + statement_copy[2:]
        if statement_copy[-1] == ':':
            statement_copy = statement_copy[:-1]
        statement_copy = statement_copy + '){'
        #we found an operator so now we have to add an extra set of ()
        if found_operator:
            statement_copy = statement_copy.replace('if(', 'if((')
            statement_copy = statement_copy.replace('){', ')){')

        return statement_copy


    def convert_funs(self, statement_to_parse):
        '''
        convert all functions within a string

        Parameters
        ----------
        statement_to_parse : str
            statement string to convert.

        Returns
        -------
        statement_copy : str
            function converted string.

        '''
        statement_copy = statement_to_parse[:]
        statement_copy = self.convert_keyfuns(statement_copy)
        return statement_copy


    def edit_statement(self, string, key):
        '''
        convert the key functions to their c++ equivelant

        for example:

            sum(X_full.block(0,0,1,34)) to X_full.block(0,0,1,34).sum()

        Parameters
        ----------
        string : str
            the string to convert.
        key : str
            key to convert.

        Returns
        -------
        TYPE
            converted statement string.

        '''
        if key == 'sum':

            if string != 'sum()':
                string = string.replace('sum(', '')[:-1] + '.sum()'

        if key == 'cast_to_int':
            string = string.replace('cast_to_int(', '')[:-1]
            string = '(' + self.keyfuns['cast_to_int']%(string, string, string) + ')'
        return string

