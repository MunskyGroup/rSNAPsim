# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:56:11 2020

@author: willi
"""


from . import GenericMetaData
GenericMetaData = GenericMetaData.GenericMetaData
import numpy as np
import json, codecs

class ODE_Soln():
    def __init__(self):
        
        self.__arrays = [    'probe_loc',
                             'x0',
                             'k',
                             'N',
                             'time',
                             'mu_It',
                             'intensity_acc',
                             'intensity_acc_norm',]
        
        self.__arrays_2d = ['var_state_ss']
        
        self.__vals = ['fi','kb','N','solve_time','mu_I_ss','var_I_ss']
        self.__meta = GenericMetaData().get()
        
        
        pass
    
    def save(self,filename):
        ext = filename.split('.')[-1]
        
        if 'txt' == ext:
            self.__save_txt(filename)
        if 'json' == ext:
            self.__save_json(filename)
            

    def load(self,filename):
        ext = filename.split('.')[-1]
        
        if 'txt' == ext:
            self.__load_from_txt(filename)
        if 'json' == ext:
            self.__load_from_json(filename)
            
    def __save_txt(self,filename):

        if '.txt' not in filename:    
            filename = filename + '.txt'


      
        f = open(filename,'w')
        for key in self.__dict__.keys():
            
            if key in self.__arrays + self.__arrays_2d + self.__vals:

                f.write((key + '\r\n'))
                np.savetxt(f, np.atleast_2d(self.__dict__[key]), delimiter=',', fmt='%s')
                f.write(('\r\n'))
                
        f.close()

    def __load_from_txt(self, filename):
        if '.txt' in filename:
            ode_obj = np.loadtxt(filename, dtype=str,delimiter='\n')
          
            solutions = []
            for i in range(0,len(ode_obj)-1):
                label = ode_obj[i]                
                
                if label in self.__arrays:
                    
                    array = np.fromstring(ode_obj[i+1], dtype=float, sep=',')
                    exec(('self.'+label+ '=array'))  
                
                if label in self.__vals:
                    value = np.fromstring(ode_obj[i+1], dtype=float, sep=',')
                    exec(('self.'+label+'=value'))
                    
                if label in self.__arrays_2d:
                    j = 1
                    array2d = np.array([[]])
                    while ode_obj[i+j] not in self.__arrays + self.__vals:
                        if j == 1:
                            array2d = np.array([np.fromstring(ode_obj[i+j], dtype=float, sep=',')   ])
                        else:
                            str_2d = np.fromstring(ode_obj[i+j], dtype=float, sep=',')    
                            array2d = np.vstack((array2d,str_2d))
                        
                        j+=1
                        
                   
                    exec(('self.'+label+'=array2d'))
                        
                


    def __save_from_json(self, filename):

        if '.json' not in filename:
            filename =  filename + '.json'

        odedict = {}
        for key in self.__dict__.keys():
            if key in self.__arrays:
                
                odedict[key] = self.ssa_harr.__dict__[key].tolist()
            else:
                odedict[key] = self.ssa_harr.__dict__[key]

        json.dump(odedict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)


    def __load_json(self,filename):
        if '.json' in filename:

            obj_text = codecs.open(filename, 'r', encoding='utf-8').read()
            ode_dict = json.loads(obj_text)

            for key in ode_dict.keys():
                if key in self.__arrays:

                    self.__dict__[key] = np.array(ode_dict[key])
                else:
                    self.__dict__[key] = ode_dict[key]