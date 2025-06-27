# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 12:00:58 2024

@author: willi
"""

from . import custom_errors as custom_err
import numpy as np
from itertools import product
import inspect
import re
import os
import shutil
import time
import sys
import hashlib
import platform
import subprocess

class RuleConverterLambda():
    def __init__(self):
        
        self.keyfuns = {'sum': '.sum()',
                        'cast_to_int':'static_cast<int>(%s < 0 ? %s - 0.5 :'\
                                      ' %s + 0.5)'}

        #keywords with [] to find and convert

        self.key_words = dict(zip('k,t,p,ke,o,l,pr,s,r,nr'.split(','),
                          ['parameters', 'tc', 'rib_arr', 'elong_mat','occupied',
                           'lattice_arr', 'probe_mat', 'state_arr', 'resource_arr', 'NR']))

        self.special_func = {'np.any':'.any()', 'sum': '.sum()',
                             }
        
        self.special_func_just_name = {
                             'np.max':'std::max',
                             'np.min':'std::min',
                             'np.sum':'std::accumulate',
                             'min':'std::min',
                             'max':'std::max',
                             '.max':'.max',
                             '.min':'.min',
                             '.mean':'.mean',
                             '.sum':'.sum',
                             '.any':'.any'}
        
        self.brackets = ['[',']','(',')']

        #if operators to parse
        self.operators = ['+', '*', '-', '**', '^','%', '&', '/']
        self.if_operators = ['<', '>', '==', '!=', '>=', '<=', '~','!']
        self.delimiters = self.operators + self.if_operators + self.brackets
        self.re_splitw ='|'.join('(?<={})'.format(re.escape(delim)) for delim in self.delimiters)
        self.re_splitw_all ='|'.join('(?<={})'.format(re.escape(delim)) for delim in self.delimiters + [',',':'])
        self.re_splitwo_all ='|'.join(map(re.escape, self.delimiters + [':', ',']))
        self.re_splitwo = '|'.join(map(re.escape, self.delimiters))
        self.re_splitwo_comma = '|'.join(map(re.escape, self.delimiters + [',']))
        self.re_splitw_comma = '|'.join('(?<={})'.format(re.escape(delim)) for delim in self.delimiters + [','])
        
        pass
    
    def split_out_comments(self, rule_string):
        if '#' in rule_string:
            return rule_string.split('#')[0],rule_string.split('#')[1]
        else:
            return rule_string, ''
        
    def convert_comment(self, comment_str):
        if len(comment_str) != 0:
            return '// ' + comment_str
        else:
            return ''
        
        
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
                    raise custom_err.MisMatchedParenthesis('Missing "("  with closed '\
                                                ' parenthesis at %i'%i)

        if len(istart) != 0:  # still a ( that cant be matched
            raise custom_err.MisMatchedParenthesis('Mismatched ")" ')
        return list(bracket_dict.items())

    @staticmethod
    def get_curly_bracket_dict(string):
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
            if char == '{':
                istart.append(i+1)
            if char == '}':
                try:
                    bracket_dict[istart.pop()] = i
                except:
                    raise custom_err.MisMatchedBrackets('Missing "{"  with closed '\
                                             'parenthesis at %i'%i)

        if len(istart) != 0:  # still a ( that cant be matched
            raise custom_err.MisMatchedBrackets('Mismatched "}" ')
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
                    raise custom_err.MisMatchedBrackets('Missing "["  with closed '\
                                             'parenthesis at %i'%i)

        if len(istart) != 0:  # still a ( that cant be matched
            raise custom_err.MisMatchedBrackets('Mismatched "]" ')
        return list(bracket_dict.items())


    def split(self, string, maxsplit=0):
        delimiters = self.operators + self.if_operators + self.brackets + [',',':']
        elements = re.split(self.re_splitw_all, string, maxsplit)
        
        elem = [] 
        for i in range(len(elements)):
            element = elements[i]
            if element != '':
                if re.search(self.re_splitwo_all,element) != None:
                    if len(re.search(self.re_splitwo_all,element).group(0)) > 0:
                        e1 =  re.split(self.re_splitwo_all, element, maxsplit)[0] 
                        e2 = re.split(self.re_splitw_all, element, maxsplit)[0].replace(e1, '')
                        elem = elem + [e1,] + [e2,]
                    else:
                        elem = elem + element
                else:
                    elem = elem + [element]            
        return [x for x in elem if x != '']

    def get_elements(self, propensity_string, defined_dict = {'footprint':'9'}, operator_dict = {'~':'!'}):
        '''
        Deconstruct a propensity string into its consituent elements
        '''
        propensity_string = propensity_string.replace(' ','')
        elements = self.split(propensity_string)
        keyword_elements = [int(x in self.key_words.keys()) for x in elements]
        function_elements = [int(x in self.special_func.keys()) for x in elements]
        functionname_elements = [int(x in self.special_func_just_name.keys()) for x in elements]
        index_elements = [0] + [int(elements[x-1] == '[') for x in range(1, len(elements))  ]
        operator_elements = [elements[i] in self.operators + self.if_operators for i in range(len(elements)) ]
        bracket_elements = [elements[i] in self.brackets for i in range(len(elements)) ]
        number_elements = [  (~index_elements[i] + 2)  * elements[i].isnumeric() for i in range(len(elements)) ]
        commacolon_elements = [int(x in [',',':']) for x in elements]
        predefined_elements = [int(x in defined_dict.keys()) for x in elements]
        
        
        print(number_elements)
        print(operator_elements)
        print(elements)
        unknown_elements = []
        
        esum = np.sum(np.array([predefined_elements, commacolon_elements, keyword_elements, functionname_elements, function_elements, index_elements, operator_elements , bracket_elements, number_elements ]), axis=0)
        unknown_elements = [elements[int(x)] for x in np.where(esum==0)[0]]
        print('UNKNOWN ELEMENTS: ')
        print(unknown_elements)
        print(esum)
        if not np.all(esum):
            print(esum)
            raise custom_err.UnknownElementError('The following elements are undefined or unimplemented yet: ' + str(unknown_elements) + ' in the following propensity function: ' + propensity_string)
            
        
        # convert keywords here
        new_elem = []
        for i in range(len(elements)):
            if keyword_elements[i] > 0:
                new_elem = new_elem + [self.key_words[elements[i]]]
            else:
                new_elem = new_elem + [elements[i]]
                
        # convert special defined variables here like footprint
        new_elem2 = []
        for i in range(len(new_elem)):
            if predefined_elements[i] > 0:
                new_elem2 = new_elem2 + [defined_dict[new_elem[i]]]
            else:
                new_elem2 = new_elem2 + [new_elem[i]]
                
        # SPECIAL CONVERSIONS FOR OPERATORS
        new_elem3 = []
        for i in range(len(new_elem)):
            if operator_elements[i] > 0:
                if new_elem2[i] in operator_dict.keys():
                    new_elem3 = new_elem3 + [operator_dict[new_elem2[i]]]
                else:
                    new_elem3 = new_elem3 + [new_elem2[i]]      
            else:
                new_elem3 = new_elem3 + [new_elem2[i]]        
               
        # convert indexes
        new_str = ''.join(new_elem3)
        square_brackets = self.get_square_bracket_dict(new_str) # square bracket pairs
        fun_brackets = self.get_parenthesis_dict(new_str) # normal parenthesis pairs (functions)
        inverted_brackets = self.invert_pair_locations(square_brackets, len(new_str))
        pair_orders, pair_orders_str = self.get_orders_of_substrings(square_brackets, inverted_brackets)
        
        new_str = new_str.replace('[','{')
        new_str = new_str.replace(']','}')
        '''
        substrs = []
        for pair in pair_orders:
            substr = new_str[pair[0]:pair[1]]
            print(substr)
            # insane python index parsing conversion to eigen c++ 
            # TODO: ANYTHING DEPENDENT ON KNOWING MAX COL OR ROW SIZE IS NOT IMPLEMENTED
            if '[' not in substr and ']' not in substr:        
                if ':' not in substr: # no : in index automatically convert
                    new_indexes = '{' + substr + '}'
                elif ',' in substr:       # , and : in index
                    new_indexes = '〈' + substr + '〉' # MARK THESE AS SPECIAL, MULTIDIMENSIONALS

            else:
                substr = substr.replace('[','{').replace(']','}')
            substrs.append(substr)
            
        substrs2 = []
        for substr in substrs:
            if '~' in substr:
                substrs2.append(substr.replace('~','!'))
            else:
                substrs2.append(substr)
        '''
        
        # convert functions finally here
        '''
        new_elem2 = []
        for i in range(len(elements)):
            if function_elements[i]:
                new_elem2 = new_elem2 + [self.special_func[new_elem[i]]]
            else:
                new_elem2 = new_elem2 + [new_elem[i]]
        '''
        
        return new_str


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


    def convert_2d_indexes(self, substr, max_col = 'max_col', max_row = 'max_row'):
        defined_inds = [0,0,0,0]
        
        # HAVE TO PARSE CASES LIKE [max(3,4), 1:max(a[:3])]
        #have to find the center comma: 
        curly_brackets = self.get_curly_bracket_dict(substr)
        square_brackets = self.get_square_bracket_dict(substr)
        parenthesis = self.get_parenthesis_dict(substr)
        ops = curly_brackets + parenthesis + square_brackets
        true_comma=-1
        print(substr)
        if len(ops) !=0:
            inverted_ops = self.invert_pair_locations(ops, len(substr))
            true_comma = ''
            for i in range(len(substr)):
                if substr[i] == ',':
                    print(i, [print(i) for pair in inverted_ops])
                    if any([i in range(*pair) for pair in inverted_ops]):
                        true_comma = i
            ind1 = substr[:true_comma]
            ind2 = substr[true_comma+1:]
        else:
            if ',' in substr:
                true_comma = substr.index(',')
                ind1 = substr[:true_comma]
                ind2 = substr[true_comma+1:]
            else:
                ind1 = substr
                ind2 = ''
        
        # get the true : for ind1
        true_colon1 = -1
        if len(ops) !=0:
            curly_brackets = self.get_curly_bracket_dict(ind1)
            square_brackets = self.get_square_bracket_dict(ind1)
            parenthesis = self.get_parenthesis_dict(ind1)
            ops = curly_brackets + parenthesis + square_brackets
            inverted_ops = self.invert_pair_locations(ops, len(ind1))
            true_colon1=-1
            for i in range(len(ind1)):
                if ind1[i] == ':':
                    if any([i in range(*pair) for pair in inverted_ops]):
                        true_colon1 = i
        else:
            if ':' in substr[:true_comma]:
                true_colon1 = substr.index(':')
            else:
                true_colon1 = -1
                
        # get the true : for ind2
        true_colon2 = -1
        if len(ops) != 0:
            curly_brackets = self.get_curly_bracket_dict(ind2)
            square_brackets = self.get_square_bracket_dict(ind2)
            parenthesis = self.get_parenthesis_dict(ind2)
            ops = curly_brackets + parenthesis + square_brackets
            inverted_ops = self.invert_pair_locations(ops, len(ind2))
            true_colon2=-1
            for i in range(len(ind2)):
                if ind1[i] == ':':
                    if any([i in range(*pair) for pair in inverted_ops]):
                        true_colon2 = i
        else:
            if ':' in substr[true_comma:] and true_comma !=-1:
                true_colon2 = substr[true_comma:].index(':') + len(substr[:true_comma])
            
                
            
        
        # SPECIAL CASE: if its  just a single value dont make a .block(), just return it
        if true_colon1==-1 and true_colon2==-1 and true_comma !=-1:
            istr = '(' + ind1 + ',' + ind2 + ')'
            return istr
        
        
        if true_comma != -1: # TWO DEFINED INDICES [xx:xx,xx:xx]
            front, back = ind1,ind2
            if len(back) == 0:  #ONLY ROW DEFINED  [xx:xx,]
                b1,b2 = '0',max_col
            else:
                if true_colon2 != -1: #[xx:xx,yy:yy]
                    b1,b2 = substr[true_comma:true_colon2+1], substr[true_colon2:]  
                    if len(b1) == 0: #[xx:xx,:yy]
                        b1 = '0'
                    if len(b2) == 0: #[xx:xx,yy:]
                        b2 = max_col                   
                else:  #[xx:xx, N]
                    b1,b2 = back, back + '+1'
                    
        else: #ROWS ONLY [xx:xx] or [xx]
            front = substr
            back = '0:' + max_col
            b1,b2 = '0', max_col
        
        if len(front) == 0:
            f1,f2 = '0',max_row
        else:
            if true_colon1 != -1: # two indeces in front [xx:xx, yy:yy]
                f1,f2  = substr[:true_colon1], substr[true_colon1+1:true_comma] #[xx:xx, yy:yy]
                if len(f1) == 0: #[:xx, yy:yy]
                    f1 = '0'
                if len(f2) == 0: #[xx:, yy:yy]
                    f2 = max_row                   
            else:
                f1,f2 = front, front + '+1' #[xx, yy:yy]
                
        istr = '.block('  + f1 + ',' + b1 + ',' +'std::min(' + f2 + '-' + f1 + ',' + max_row + '-' + f1 +')' + ',' +'std::min(' + b2 + '-' + b1 + ',' + max_col + '-' + b1 +'))'
        if istr == '.block(0,0,'+max_row+'-0,' + max_col + '-0)': # SPECIAL CASE [:,:] DELETE THIS ITS NOT NEEDED
            istr = ''
            
        # SPECIAL CASE FLAG TO KEEP ARRAY SIZE CONSTANT FOR VARIABLE LENGTH INDICES    
        
        return istr


    def convert_2D_ind_kepr(self, function_str, vector_str):
        
        square_bracket_dict = self.get_curly_bracket_dict(function_str) #contents of {} brackets
        arr_inds = [(m.start(), m.start()+len(vector_str)) for m in re.finditer(vector_str, function_str)]

        la_bstart = [[x for x in square_bracket_dict if x[0]-1==y[1]] for y in arr_inds]
        la_bstart = [x for xs in la_bstart for x in xs]
        
        arr_substrings = [(arr_inds[i][0],la_bstart[i][1]+1) for i in range(len(arr_inds))]
        inverted_arr_substrings = self.invert_pair_locations(arr_substrings, len(function_str))
        ordered_substrings, orderimi = self.get_orders_of_substrings(arr_substrings, inverted_arr_substrings)
        substrings = [function_str[x[0]:x[1]] for x in ordered_substrings]
        
        tochange = []
        for i in range(len(orderimi)):
            if orderimi[i] == 'm':
                tochange.append(i)
        

        k = 0
        for ind in tochange:
            print(substrings)
            print('---')
            substr = substrings[ind]
            print(la_bstart)
            ind_str = function_str[la_bstart[k][0]:la_bstart[k][1]] 
            
            ind_str = self.convert_2d_indexes(ind_str, max_row = '3', max_col = 'L')
            

            substrings[ind] = vector_str + ind_str
            k+=1
        
        return ''.join(substrings)
        
    def convert_2D_ind_P(self, function_str):
        
        vector_str = 'rib_arr'
        square_bracket_dict = self.get_curly_bracket_dict(function_str) #contents of {} brackets
        arr_inds = [(m.start(), m.start()+len(vector_str)) for m in re.finditer(vector_str, function_str)]

        la_bstart = [[x for x in square_bracket_dict if x[0]-1==y[1]] for y in arr_inds]
        la_bstart = [x for xs in la_bstart for x in xs]
        
        arr_substrings = [(arr_inds[i][0],la_bstart[i][1]+1) for i in range(len(arr_inds))]
        inverted_arr_substrings = self.invert_pair_locations(arr_substrings, len(function_str))
        ordered_substrings, orderimi = self.get_orders_of_substrings(arr_substrings, inverted_arr_substrings)
        substrings = [function_str[x[0]:x[1]] for x in ordered_substrings]
        
        tochange = []
        for i in range(len(orderimi)):
            if orderimi[i] == 'm':
                tochange.append(i)
        
        k = 0
        for ind in tochange:
            print(substrings)
            print('---')
            substr = substrings[ind]
            print(la_bstart)
            ind_str = function_str[la_bstart[k][0]:la_bstart[k][1]] 
            
            ind_str = self.convert_2d_indexes(ind_str, max_row = 'max_rib', max_col = 'rib_arr_col_size')
            

            substrings[ind] = vector_str + ind_str
            k+=1
        
        return ''.join(substrings)
    
        return function_str 
        
        
   # test_inds = ['0','1:',':1','0,:','1:2', '0,1:2','0,:1','0,1:',':,:','0:2,0:2',':2,:2','1:3,1:3']

    def convert_colon_indexes_1D(self,substr, max_col, max_row):
        new_indexes = ''
        notated_inds = [int(len(x) > 0) for x in [item for sublist in [x.split(':') for x in substr.split(',')] for item in sublist]]
        print(notated_inds)
        if substr.count(':') == 1:
            notated_inds = [int(len(x) > 0) for x in [item for sublist in [x.split(':') for x in substr.split(',')] for item in sublist]]
            # [1,:] columnwise
            if notated_inds == [1, 0, 0]:
                new_indexes = '.col{%s}'%substr[:-2]
            # [:,1] rowwise
            if notated_inds == [0, 0, 1]:
                new_indexes = '.row{%s}'%substr[2:]    
            
            # [1, x:y] or  # [x:y, 1]
            if notated_inds == [1,1,1]:
                # [x:y, 1]
                if ':' in substr.split(',')[0]:
                    i, j = int(substr.split(':')[0]), int(substr.split(':')[1].split(',')[1])
                    k, l = int(substr.split(':')[1].split(',')[0]), j+1
                    p = k - i
                    q = l - j
                    new_indexes = '.block{%s,%s,%s,%s}'%(i, j, p, q)                                
                    
                # [1, x:y]
                else:
                    i, j = int(substr.split(':')[0].split(',')[0]),   int(substr.split(':')[0].split(',')[1])
                    k, l = i+1, int(substr.split(':')[1])
                    p = k - i
                    q = l - j
                    new_indexes = '.block{%s,%s,%s,%s}'%(i, j, p, q)                                
                                                

            # [1, :y] or # [y:, 1]
            if notated_inds == [1,0,1]:
                # [y:, 1]
                if ':' in substr.split(',')[0]:
                    i, j = int(substr.split(':')[0]), int(substr.split(':')[1].split(',')[1])
                    k, l = max_row - i, j+1
                    p = k - i
                    q = l - j
                    new_indexes = '.block{%s,%s,%s,%s}'%(i, j, p, q)
                    
                # [1, :y] 
                else:
                    i, j = int(substr.split(':')[0].split(',')[0]), 0
                    k, l = i+1, int(substr.split(':')[1])
                    p = k - i
                    q = l - j
                    new_indexes = '.block{%s,%s,%s,%s}'%(i, j, p, q)                                
                
            # [1, y:]
            if notated_inds == [1,1,0]:
                i, j = int(substr.split(':')[0].split(',')[0]), int(substr.split(':')[0].split(',')[1])
                k, l = i+1, max_col - j
                p = k - i
                q = l - j
                new_indexes = '.block{%i,%i,%i,%i}'%(i, j, p, q)
                
            # [:y, 1]
            if notated_inds == [0,1,1]:
                i, j = 0, int(substr.split(':')[1].split(',')[1])
                k, l = int(substr.split(':')[1].split(',')[0]), j+1
                p = k - i
                q = l - j
                new_indexes = '.block{%i,%i,%i,%i}'%(i, j, p, q)
                
                
        else: # DOUBLE : INDEXING
            notated_inds = [int(len(x) > 0) for x in [item for sublist in [x.split(':') for x in substr.split(',')] for item in sublist]]
            # [u:v, x:y]
            
            if notated_inds == [1]:
                i = int(substr)
                new_indexes = '.row{%i}'%(i)
                
            if sum(notated_inds) == 4:
                i, j = int(substr.split(':')[0]), int(substr.split(':')[1].split(',')[1])
                k, l = int(substr.split(':')[2].split(',')[0]), int(substr.split(':')[2])
                p = k - i
                q = l - j
                new_indexes = '.block{%i,%i,%i,%i}'%(i, j, p, q)
                
            # [:, x:y]
            if notated_inds == [0,0,1,1]:
                new_indexes = '.middleCols{%s, %s}'% (substr.split(',')[1].split(':')[0], substr.split(',')[1].split(':')[1]  )
                
            # [x:y, :]
            if notated_inds == [1,1,0,0]:
                new_indexes = '.middleRows{%s, %s}'% (substr.split(':')[0], substr.split(':')[1].split(',')[0])
            
            if notated_inds == [1,0,0,0]:
            # [2:, :]
                i = str(max_row - int(substr.split(':')[0]))
                new_indexes= '.bottomRows{%s}'% (i)
            
            if notated_inds == [0,0,1,0]:
            # [:, 2:]
                j = str(max_col - int(substr.split(':')[1].split(',')[1]))
                new_indexes= '.rightCols{%s}'% (j)
            
            if notated_inds == [0,1,0,0]:
            # [:2, :]
                i = substr.split(':')[1].split(',')[0]
                new_indexes= '.topRows{%s}'% (i)
                
            if notated_inds == [0,0,0,1]:
            # [:, :2]
                j = substr.split(':')[2]
                new_indexes= '.topRows{%s}'% (j)

            if notated_inds == [0,1,0,1]:
            # [:2, :2]
                new_indexes= '.bottomLeftCorner{%s, %s}'% (substr.split(':')[1].split(',')[0], substr.split(':')[2] )
                
            if notated_inds == [1,0,0,1]:
                # [2:, :2]
                i = str(max_row - int(substr.split(':')[0]))
                new_indexes= '.bottomLeftCorner{%s, %s}'% (i, substr.split(':')[2])

            if notated_inds == [0,1,1,0]:
            # [:2, 2:]
                j = str(max_col - int(substr.split(':')[1].split(',')[1]))
                new_indexes= '.topRightCorner{%s, %s}'% (substr.split(':')[1].split(',')[0], j)
                
            if notated_inds == [1,0,1,0]:
            # [2:, 2:]
                i = str(max_row - int(substr.split(':')[0]))
                j = str(max_col - int(substr.split(':')[1].split(',')[1]))
                new_indexes= '.bottomRightCorner{%s, %s}'% (i, j)


            # [:, :]
            if notated_inds == [0,0,0,0]:
                new_indexes = ''
        return new_indexes
    
    def convert_1d_array_inds(self, function_str, vector_str, max_row):
        square_bracket_dict = self.get_curly_bracket_dict(function_str) #contents of {} brackets
        arr_inds = [(m.start(), m.start()+len(vector_str)) for m in re.finditer(vector_str, function_str)]

        la_bstart = [[x for x in square_bracket_dict if x[0]-1==y[1]] for y in arr_inds]
        la_bstart = [x for xs in la_bstart for x in xs]
        
        print(function_str)
        print(vector_str)
        print(arr_inds)
        print(la_bstart)
        if len(la_bstart) == 0: # dont convert if theres no brackets to convert ex "lattice_arr*3" vs "lattice_arr{1,2}*3"
            return function_str
        arr_substrings = [(arr_inds[i][0],la_bstart[i][1]+1) for i in range(len(arr_inds))]
        inverted_arr_substrings = self.invert_pair_locations(arr_substrings, len(function_str))
        ordered_substrings, orderimi = self.get_orders_of_substrings(arr_substrings, inverted_arr_substrings)
        substrings = [function_str[x[0]:x[1]] for x in ordered_substrings]
        
        tochange = []
        for i in range(len(orderimi)):
            if orderimi[i] == 'm':
                tochange.append(i)
        
        k = 0
        for ind in tochange:
            print(substrings)
            print('---')
            substr = substrings[ind]
            print(la_bstart)
            ind_str = function_str[la_bstart[k][0]:la_bstart[k][1]] 

            
            if ':' in ind_str:
                first_index, second_index = ind_str.split(':')
                if '-' in second_index:
                    #TODO:
                    x=1 # PARSE NEGATIVE INDS TODO
                if '-' in first_index:
                    x=1 # PARSE NEGATIVE INDS
                newstr = vector_str + '(Eigen::seq(%s,std::min(%s, %s)))'%(first_index,second_index, max_row)
                substrings[ind] = newstr
            else:
                substrings[ind] = vector_str + '[' + ind_str +  ']'
            k+=1
        
        return ''.join(substrings)


    def convert_special_function_name_only(self, string):
        elements = re.split(self.re_splitw, string, 0)
        print(elements)
        for special_func in self.special_func_just_name.keys():
            if special_func + '(' in elements:
                for i in range(len(elements)):
                    if elements[i] == special_func + '(':
                        elements[i] = self.special_func_just_name[special_func] + '('
            
        return ''.join(elements)
    

    
    def convert_special_function(self, string):
        elements = re.split(self.re_splitwo, string, 0)
        print(elements)
        for special_func in self.special_func.keys():
            total_to_change = elements.count(special_func)
            print(total_to_change)
            for k in range(total_to_change):
                parts = string.split(special_func)
                parenthesis = self.get_parenthesis_dict(string)
                match_inds = [(m.start(),m.end()) for m in re.finditer(special_func,string)]
                pairs_to_change = []

                
                for i in range(len(match_inds)):
                    for j in range(len(parenthesis)):
                        if parenthesis[j][0] == match_inds[i][1]+1:
                            pairs_to_change = [(match_inds[i], parenthesis[j])]
                            
                string = string[:pairs_to_change[0][0][0]] + string[pairs_to_change[0][0][1]:pairs_to_change[0][1][1]+1] + self.special_func[special_func] + string[pairs_to_change[0][1][1]+1:]
        
        return string
        

    def convert_list_comprehension(self,function_str):
        ind_letter = function_str.split(' ')[-3]
        function_substr = ''.join(function_str.split(' ')[:-4])[1:]
        elements = re.split(self.re_splitwo, function_substr, 0)
        elements_w = re.split(self.re_splitwo, function_substr, 0)
        elements_w_comma = re.split(self.re_splitw_comma, function_substr,0)

        new_elements = []
        if function_substr !='i':
            for element in elements_w_comma:
                if ind_letter in element:
                    element.replace(ind_letter, 'i')
                else:
                    pass
                new_elements.append(element)
            return ''.join(new_elements)
            
        else:
            return function_substr
            

    def convert_lambda_function_to_c(self, rule_string):
        
        
        print(rule_string)
        # first parse out comment if applicable and convert it
        function_str, comment_str = self.split_out_comments(rule_string)
        c_comment = self.convert_comment(comment_str)
        print('Detected comment:')
        print(c_comment)
        
        print('fstring:')
        print(function_str)
        # DETECT IF ITS A LIST COMPREHENSION
        if 'in range(nr)]' in function_str:
            function_str = self.convert_list_comprehension(function_str)
        print(function_str)
        
        # convert keywords and ~operator
        elements = self.get_elements(function_str)
        print('*** GETTING ELEMENTS***')
        print(elements)
        
        function_str = ''.join(elements)
        print(function_str)
        square_bracket_dict = self.get_curly_bracket_dict(function_str) #contents of {} brackets
        parenthesis_dict = self.get_parenthesis_dict(function_str) #contents of parenthesis
        
        
       # k,t,p,ke,o,l,pr,s,r,nr'
       
        # SPECIAL values: t, nr
        ## CONVERT SPECIAL VALUES, doesnt need to be done
       
       # vectors: k, l, o, s, r
        print('*** CONVERTING 1D ARRAYS***')
        ## CONVERT 1D INDEXES 
        # convert lattice_arr
        if 'lattice_arr' in function_str:
            function_str = self.convert_1d_array_inds(function_str,'lattice_arr', 'L-1')

        # convert state_arr
        if 'state_arr' in function_str:
            function_str = self.convert_1d_array_inds(function_str,'state_arr', 'n_states')
    
        # convert resource_arr
        if 'resource_arr' in function_str:
            function_str = self.convert_1d_array_inds(function_str,'resource_arr', 'n_resources')

        # convert parameters
        if 'parameters' in function_str:
            function_str = self.convert_1d_array_inds(function_str,'parameters', 'n_parameters')

        # convert occupied
        if 'occupied' in function_str:
            function_str = self.convert_1d_array_inds(function_str,'occupied', 'max_rib')    

        print(function_str)

        
        ## CONVERT 2D INDEXED MATRICES
        # matrices: p, ke, pr   
        
        # P is [Max_ribosomes by Nrxns + Ncolors + 4]
        # KE is 3 x L
        # PR is 3 x L
        if 'elong_mat' in function_str:
            function_str = self.convert_2D_ind_kepr(function_str,'elong_mat' )
        if 'probe_mat' in function_str:
            function_str = self.convert_2D_ind_kepr(function_str,'probe_mat' )
        if 'rib_arr' in function_str:
            function_str = self.convert_2D_ind_P(function_str)
        
        ## CONVERT SPECIAL FUNCTIONS
        print('*** CONVERTING SPECIAL FUNCTIONS***')
        function_str = self.convert_special_function_name_only(function_str)
        print(function_str)

        print('*** CONVERTING SPECIAL FUNCTIONS***')
        function_str = self.convert_special_function(function_str)
        print(function_str)


        #self.key_words = dict(zip('k,t,p,ke,o,l,pr,s,r,nr'.split(','),
        #                  ['parameters', 'tc', 'rib_arr', 'elong_mat','occupied',
        #                   'lattice_arr', 'probe_mat', 'state_arr', 'resource_arr', 'NR']))

        # convert previous 1 dimensional square indexes based on the designated array:
        
        ## ADD COMMENT BACK IN
        function_str = function_str + '; ' + c_comment
        print(function_str)

        return function_str
        
    
    
    def make_c_propensities(self, propensity_lambdas, ribosome=0):

        propensities_str_list = [y[0].replace('\n','') for y in [inspect.getsourcelines(x)[0] for x in propensity_lambdas]]
        propensity_names_list = [y[0].replace('\n','').split('=')[0].replace(' ','') for y in [inspect.getsourcelines(x)[0] for x in propensity_lambdas]]
        
        # =.join is to handle stuff like s == 1 in the lambda string since we have to split on the first =
        propensity_function_list = ['='.join(y[0].replace('\n','').split('=')[1:]).replace('lambda k,t,p,ke,o,l,pr,s,r,nr: ','') for y in [inspect.getsourcelines(x)[0] for x in propensity_lambdas]]
        
        rules = []
        for i in range(len(propensities_str_list)):
            propensities_str_list[i]
            print(propensity_function_list[i])
            if not ribosome:
                rules.append( 'wn[%s] = '%str(i) + self.convert_lambda_function_to_c(propensity_function_list[i]) + '\n')
            else:
                rules.append( 'wn[k+%s] = '%str(i) + self.convert_lambda_function_to_c(propensity_function_list[i]) + '\n')
        return ''.join(rules)
        



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
                    
        try:
            self.eigen_path = self.find_eigen_path()
            print('eigen instillation found...')
            print(self.eigen_path)
        except:
            
            self.eigen_path=''
        self.reserved_model_names = ['build', 'model_maker_cpp', 'models',
                                     'rsnapsim_model_maker']
        print(self.eigen_path)
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
        potential_paths = []  #check the <env>/lib/ folder and check the <env>/Library/folder
                
        for path in paths:

            if path[-4:] == 'DLLs':
                potential_paths.append(path)
                potential_paths.append(path[:-3] + 'Library')
            if path[-3:] == 'Lib':
                potential_paths.append(path)
                potential_paths.append(path[:-3] + 'Library')
            if path[-7:] == 'include':
                potential_paths.append(path)

        
        eigen_paths = []
        for path in potential_paths:  #in each of these try to find an eigen instillation
            base, _ = os.path.split(path)
            if os.path.exists(os.path.join(base, 'Library',
                                           'include', 'eigen3')):

                eigen_paths.append(os.path.join(base, 'Library',
                                                'include', ''))
            if os.path.exists(os.path.join(base, 'Eigen')):
                eigen_paths.append(os.path.join(base,''))
            if os.path.exists(os.path.join(base, 'eigen3')):
                eigen_paths.append(os.path.join(base,''))
        
        
        if len(eigen_paths) == 0:
            raise custom_err.EigenMissingError('Eigen is missing, please provide a path'\
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
        metadata = {} #get all the metadata
        
        metadata['created_at'] = time.strftime('%Y-%m-%d %H:%M:%S',
                                        time.localtime(time.time()))
        metadata['user'] = os.path.expanduser("~")
        metadata['platform'] = platform.platform()
        metadata['python_version'] = sys.version

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


    def edit_cpp_files(self, model_name, cpp_file_path, prop_str_constant, prop_str_ribosome):
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
        #prop_str = strfrnt + prop_str + strback
        fstr = ''
        with open(cpp_file_path, 'r') as fname:
            for line in fname:

                if "//INSERT_GENERATED_CONSTANT_PROPENSITY_HERE" in line:
                    print('replacing_constant_prop...')
                    line = line.replace("//INSERT_GENERATED_CONSTANT_PROPENSITY_HERE",
                                        prop_str_constant)
                if "//INSERT_GENERATED_RIBOSOME_PROPENSITY_HERE" in line:
                    print('replacing_ribosome_prop...')
                    line = line.replace("//INSERT_GENERATED_RIBOSOME_PROPENSITY_HERE",
                                        prop_str_ribosome)

                fstr += line

        with open(cpp_file_path, 'w') as fname:
            fname.write(fstr)

    def edit_pyx_files(self, model_name, prop_str_constant, prop_str_ribosome,
                       original_rules_str, model_id, pyx_file_path):
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
        rules_str_r = "rules_str = '''%s"%prop_str_constant + '/n' + prop_str_ribosome + "'''"
        original_rules_str_r = "original_rules_str = ''' %s"%original_rules_str + "'''"
        #min_length_str = "min_length = %i"%min_length

        #edit all the lines in the .pyx file
        fstr = ''
        with open(pyx_file_path, 'r') as fname:
            for line in fname:
            
                if "#model id goes here" in line:
                    line = line.replace("#model id goes here" , 'model_id = ' + "'" + model_id + "'")

                if "#cdef goes here" in line:
                    line = line.replace("#cdef goes here", cdef_str)
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
        
        if len(eigen_path) > 0:
            eigen_path = self.eigen_path[0]

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


    def make_init_file(self, folder_path, model_name, constant_propensity_str, ribosome_propensity_str):
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
        _, meta_str = self.generate_metadata(constant_propensity_str + '/n' + ribosome_propensity_str, model_name)

        init_str = ''
        init_str += "'''\n"
        init_str += meta_str
        init_str += "'''\n\n"

        #add the model name to the string to write
        init_str += 'from . import %s'%model_name#'from . import %s'%model_name

        #make the new init file
        initpath = os.path.join(folder_path, 'models',
                                model_name, '__init__.py')
        with open(initpath, "w") as fname:
            fname.write(init_str)



    def generate_model_files(self, model_name, model_id, overwrite=False,
                             eigen_path='', constant_propensity_str=None, ribosome_propensity_str=None,
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
                    if base == 'blank_model_src':
                        model_c_file_path = old_name
                        #model_c_file_found = True
                        new_model_c_file = os.path.join(os.path.abspath(root),
                                                        model_c_file_name + extension)

                if extension == '.h':
                    if base == 'blank_model_src':
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
        print(folder_exist)
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
                raise custom_err.ExistenceError(message)
        
        print(model_file_path)
        print(new_model_file)
        print(model_c_file_path)
        print(new_model_c_file)
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
                raise custom_err.ExistenceError(message)
        time.sleep(.1)
        self.edit_setup_files(model_name, new_setup_file,
                              eigen_path=eigen_path)
        self.edit_pyx_files(model_name, constant_propensity_str, ribosome_propensity_str,
                            original_rules_str, model_id, new_model_file)


        self.edit_cpp_files(model_name, new_model_c_file, constant_propensity_str, ribosome_propensity_str)
        shutil.move(new_model_file, os.path.join(folder_path,
                                                 'models', model_name))
        shutil.move(new_setup_file, os.path.join(folder_path,
                                                 'models', model_name))
        shutil.move(new_model_c_file, os.path.join(folder_path,
                                                   'models', model_name))
        shutil.move(new_model_h_file, os.path.join(folder_path,
                                                   'models', model_name))
        self.make_init_file(folder_path, model_name, constant_propensity_str, ribosome_propensity_str)

        #attempt to compile


    def compile_model(self, model_name, model_id, constant_prop, ribosome_prop, overwrite=False,
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

        if model_name in self.reserved_model_names:
            msg = 'The model name requested is a reserved keyword and'\
                ' cannot be used to build the model files, please rename'\
                    ' the model.'
            raise custom_err.ModelNameError(msg)
            
        if eigen_path == '':
            eigen_path = self.eigen_path

        # convert the rules passed
        if verbose:
            print('converting propensities to c...')
            print(constant_prop)
            
        constant_propensities_c = RuleConverterLambda().make_c_propensities(constant_prop)
        ribosome_propensities_c = RuleConverterLambda().make_c_propensities(ribosome_prop, ribosome=1)
        propensities_str_list = [y[0].replace('\n','') for y in [inspect.getsourcelines(x)[0] for x in constant_prop + ribosome_prop]]
        original_rules = '/n'.join(propensities_str_list)
        min_length = 0

        if verbose:
            print('generating model files...')

        #generate all the appropriate files
        cwd = os.getcwd()
        os.chdir(self.file_path) #change to the cwd where this file is stored
        
        self.generate_model_files( model_name, model_id, overwrite=overwrite,
                             eigen_path=eigen_path, constant_propensity_str=constant_propensities_c, ribosome_propensity_str=ribosome_propensities_c,
                             original_rules_str=original_rules, min_length=1)

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

            print('compilation failed...')
            print(compile_run.stdout.decode('utf-8'))
            print(compile_run.stderr.decode('utf-8'))
        else:
            if verbose:
                print('model compiled!')
        os.chdir(cwd)
        
        return compile_run.returncode

