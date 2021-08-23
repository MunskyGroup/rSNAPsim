# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:57:33 2020

@author: willi
"""


import json, codecs
from json import encoder
import numpy as np
from . import GenericMetaData
GenericMetaData = GenericMetaData.GenericMetaData


class SSA_Soln():
    '''
    SSA container class

    holds intensity / ribosome data as well as the propensities used

    __.n_traj = number of trajectories
    __.k = propensities used for the simulation
    __.rib_density = ribosome density per mRNA strand
    __.ribosome_means


    '''
    def __init__(self):
        self.n_traj = 0   #number trajectories
        self.k = []       #propensities
        self.no_rib_per_mrna = 0    #number of ribosomes per mrna strand
        self.rib_density = 0      #ribosome density
        self.ribosome_means = 0  #mean ribosomes
        self.rib_vec = 0          #actual vector of ribosome locations
        self.intensity_vec = []   #intensity vectors per SSA trajectory
        self.time_vec_fixed = []   #time array
        self.__meta = GenericMetaData().get()



    def save(self, filename, precision='.4f'):
        ext = filename.split('.')[-1]

        if ext == 'txt':
            self.__save_txt(filename)
        if ext == 'json':
            self.__save_json(filename, precision=precision)


    def load(self, filename):
        ext = filename.split('.')[-1]

        if ext == 'txt':
            self.__load_from_txt(filename)
        if ext == 'json':
            self.__load_from_json(filename)

    def __save_txt(self, filename):

        if '.txt' in filename:
            filehandle = open(filename, 'a')
            for key in self.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    filehandle.write((key + '\r\n'))
                    np.savetxt(filehandle, np.atleast_2d(self.__dict__[key]),
                               delimiter=',', fmt='%s')
                    filehandle.write(('\r\n'))

        else:
            filename = filename + '.txt'
            filehandle = open(filename, 'a')
            for key in self.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    filehandle.write((key + '\r\n'))
                    np.savetxt(filehandle, np.atleast_2d(self.__dict__[key]),
                               delimiter=',', fmt='%s')
                    filehandle.write(('\r\n'))
        filehandle.close()


    def __load_from_txt(self, filename):
        if '.txt' in filename:
            ssa_obj = np.loadtxt(filename, dtype=str, delimiter='\n')
            solutions = []
            for i in range(0, len(ssa_obj)-1):
                label = ssa_obj[i]


                if label in ['rib_means',
                             'rib_vec',
                             'n_traj',
                             'start_time',
                             'k',
                             'time_vec_fixed',
                             'dwelltime',
                             'mean_autocorr',
                             'no_rib_per_mrna',
                             'ke_sim',
                             'autocorr_vec',
                             'ribosome_means',
                             'error_autocorr',
                             'rib_density',
                             'time',
                             'ke_true',
                             'evaluating_inhibitor',
                             'time_inhibit',
                             'evaluating_frap']:

                    if label in ['start_time', 'no_rib_per_mrna',
                                 'ke_sim', 'dwelltime', 'ke_true',
                                 'time_inhibit']:

                        array = np.fromstring(ssa_obj[i+1], dtype=float, sep=',')[0]
                        exec(('self.'+label+ '=array'))
                    elif label in ['n_traj']:
                        array = int(np.fromstring(ssa_obj[i+1], dtype=float, sep=',')[0])
                        exec(('self.'+label+ '=array'))
                    else:
                        array = np.fromstring(ssa_obj[i+1], dtype=float, sep=',')
                        exec(('self.'+label+ '=array'))

                if label in ['evaluating_inhibitor', 'evaluating_frap']:
                    if 'False' in ssa_obj[i+1]:
                        exec(('self.'+label+ '=False'))
                    if 'True' in ssa_obj[i+1]:
                        exec(('self.'+label+ '=True'))


            for i in range(0, len(ssa_obj)-1):
                label = ssa_obj[i]

                if label == 'intensity_vec':

                    #tvec = self.time_vec_fixed[np.where(self.time_vec_fixed >= self.start_time)]
                    i_vec = np.zeros((self.n_traj, len(self.time)))

                    for j in range(self.n_traj):
                        array = np.fromstring(ssa_obj[i+j+1],
                                              dtype=float, sep=',')
                        i_vec[j] = array

                    exec(('self.'+label+ '=i_vec'))

                if label == 'solutions':
                    for j in range(self.n_traj):
                        array = np.fromstring(ssa_obj[i+j+1],
                                              dtype=float, sep=',')
                        solutions.append(array)

                    exec(('self.'+label+ '=solutions'))




    def make_dict(self):
        ssadict = {}
        for key in self.__dict__.keys():
            print(key)
            if key != 'rib_vec' and key != 'ribosome_means':
                try:
                    ssadict[key] = self.__dict__[key].tolist()
                except:
                    ssadict[key] = self.__dict__[key]

            if key == 'col_points':
                col_pt = [x.tolist() for x in self.__dict__[key]]
                ssadict[key] = col_pt


        return ssadict


    def __set_float_pres(self, precision='.4f'):
        print(precision)
        encoder.FLOAT_REPR = lambda o: format(o, precision)


    def __make_float(self, float_like, precision):
        return float(('%' + precision) % float_like)

    def __format_floats(self, arraylike, precision='.4f'):

        n_decimals = int(precision.split('.')[1][:-1])
        if isinstance(arraylike, np.ndarray):
            tmp_arr = np.around(arraylike, decimals=n_decimals) .tolist()
            return tmp_arr
        else:
            tmp_arr = arraylike
            if isinstance(tmp_arr, float):
                return self.__make_float(tmp_arr, precision)
            if isinstance(tmp_arr, list):
                return [self.__make_float(x, precision) for x in tmp_arr]


    def __save_json(self, filename, precision='.4f'):

        if '.json' in filename:

            ssadict = {}
            for key in self.__dict__.keys():

                if key != 'col_points':

                    if type(self.__dict__[key]) in [float, list, np.ndarray]:

                        ssadict[key] = self.__format_floats(self.__dict__[key],
                                                            precision=precision)
                    else:
                        ssadict[key] = self.__dict__[key]

                if key == 'col_points':
                    col_pt = [x.tolist() for x in self.__dict__[key]]
                    ssadict[key] = col_pt

            json.dump(
                ssadict, codecs.open(
                    filename, 'w', encoding='utf-8'),
                separators=(',', ':'), sort_keys=True, indent=4)

        else:
            filename = filename + '.json'

            ssadict = {}
            for key in self.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa_harr.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa_harr.__dict__[key]

            json.dump(
                ssadict, codecs.open(
                    filename, 'w', encoding='utf-8'),
                separators=(',', ':'), sort_keys=True, indent=4)


    def __load_from_json(self, filename):
        if '.json' in filename:

            obj_text = codecs.open(filename, 'r', encoding='utf-8').read()
            ssadict = json.loads(obj_text)

            for key in ssadict.keys():
                if key in ['solutions', 'all_trna_results',
                           'rib_means', 'time_vec_fixed',
                           'mean_autocorr', 'autocorr_vec',
                           'error_autocorr', 'rib_density',
                           'intensity_vec', 'I']:

                    self.__dict__[key] = np.array(ssadict[key])

                elif key in ['colpoints']:

                    cpts = [np.array(x) for x in ssadict[key]]
                    self.__dict__[key] = cpts
                else:
                    self.__dict__[key] = ssadict[key]
