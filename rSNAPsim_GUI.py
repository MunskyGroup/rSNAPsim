# -*- coding: utf-8 -*-
"""
Created on Fri Jul 06 13:40:30 2018

@author: William
"""

import os
import sys
os.chdir('ssa_cpp')
print(os.getcwd())
print(os.listdir(os.getcwd()))



import ssa_translation

try:


    import ssa_translation

except:
    pass
os.chdir('..')

import re                         #import regex
import matplotlib                 #matplotlib
#matplotlib.use("TkAgg")          #switch backends for matplotlib for tkinter




'''
try:
   # sys.path.append('C:\\Users\\wsraymon\\Github\\ssa_cpp\\translation_ssa')

    import ssa_translation
except:



    try:
        sys.path.append('C:\\Users\\wsraymon\\Github\\ssa_cpp\\translation_ssa')

        import ssa_translation
    except:
        pass

    pass

'''
import rSNAPsim





import platform
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
from matplotlib import cm
import matplotlib.patches as mpatches
import matplotlib.animation as animation
from matplotlib.collections import PatchCollection
import json, codecs
try:
    from tkColorChooser import askcolor
except:
    from tkinter.colorchooser import *
from scipy.stats import norm

try:
    from Bio import SeqIO
    from Bio import Entrez
except:
    pass



import ast
import operator as op



import PIL.Image
import PIL.ImageTk
import PIL
try:
    import Tkinter as tk
    import ttk
    from Tkinter import OptionMenu as OM


    from ttk import Frame, Label, Entry, Progressbar, Notebook, Separator, Widget, Combobox #import necessary tkk functions
    import numpy as np
    from Tkinter import Entry as tkEntry
    import tkMessageBox as tmb  #import tkitner file dialogs for opening and errors
    import tkFileDialog as tfd

except ImportError:  # Python 3
    import tkinter as tk
    from tkinter import ttk

    from tkinter import OptionMenu as OM


    from tkinter.ttk import Frame, Label, Entry, Progressbar, Notebook, Separator, Widget, Combobox #import necessary tkk functions
    import numpy as np
    from tkinter import Entry as tkEntry

    import tkinter.messagebox as tmb
    import tkinter.filedialog as tfd





import matplotlib.backends.tkagg as tkagg
try:
    from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg, NavigationToolbar2TkAgg)
except:
    from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg)    


import copy





class AutocompleteComboBox(Combobox):
        '''
        Autocompleting Combo Box adpated
        '''
        def set_completion_list(self, completion_list):
            self.completion_list = completion_list
            self.matches = []
            self.match_locs = []
            self.posistion = 0
            self.select_index = 0
            self.bind('<KeyRelease>', self.handle_keys)
            self['values'] = self.completion_list  # Setup our popup menu
            self.typed = ''
            self.prev_len = 0
            


        def fill(self):
    
            self.matches = []
            current_str = self.typed

                
            init = len(current_str)

            
            for i in range(len(self.completion_list)):
                if self.completion_list[i].lower().startswith(current_str):
                   
                    self.matches.append(self.completion_list[i])
                    self.match_locs.append(i)
                

            if len(self.matches) > 0 :
                self.delete(0,tk.END)
                self.insert(0,self.matches[0])
                self.select_range(init,tk.END)
                self.posistion = init
                
            self.prev_len = init
                

    
        def handle_keys(self, event):

            if self.get().lower() == '':
                self.typed = ''
            
            if event.keysym == 'BackSpace':
                self.posistion = self.posistion-1 
                self.typed = self.get().lower()
                
            if event.keysym == "Left":
                self.posistion = self.posistion-1                               
            if event.keysym == "Right":
                self.posistion = self.index(tk.END) 
            
            if len(event.keysym) == 1:
                self.typed = (self.typed + event.keysym)
                self.fill()
                    
            if event.keysym == "Return":
                self.select_range(tk.END,tk.END)
                self.typed = self.get().lower()
              


class ScrollingFrame(Frame):
    '''
    Y vertically scrolling frame
    '''
    def __init__(self, parent, *args, **kw):
        Frame.__init__(self, parent)

        scrollbar = tk.Scrollbar(self, orient="vertical")
        inside_canvas = tk.Canvas(self, yscrollcommand=scrollbar.set, bd=1.0, highlightthickness=0, selectborderwidth=0, relief=tk.GROOVE, *args, **kw)
        inside_canvas.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
        scrollbar.pack(fill=tk.Y,expand=tk.FALSE,side=tk.RIGHT,anchor=tk.E)
        scrollbar.configure(command=inside_canvas.yview)
        self.inside_frame = Frame(inside_canvas)


        canvas_window = inside_canvas.create_window(0, 0, window=self.inside_frame, anchor=tk.NW)

        def bind_mouse_wheel(event):
            inside_canvas.bind_all("<MouseWheel>", handle_scroll_event)

        def unbind_mouse_wheel(event):
            inside_canvas.unbind_all("<MouseWheel>")

        self.bind('<Enter>', bind_mouse_wheel)
        self.bind('<Leave>', unbind_mouse_wheel)

        def handle_scroll_event(event):
            inside_canvas.yview_scroll(int(-1*(event.delta/(parent.winfo_height()/4))), "units")

        def update_inside_frame(event):
            size = (self.inside_frame.winfo_reqwidth(), self.inside_frame.winfo_reqheight())
            inside_canvas.config(scrollregion="0 0 %s %s" % size)
            if self.inside_frame.winfo_reqwidth() != inside_canvas.winfo_width():
                inside_canvas.config(width=self.inside_frame.winfo_reqwidth())

        def update_inside_canvas(event):
            if self.inside_frame.winfo_reqwidth() != inside_canvas.winfo_width():
                inside_canvas.itemconfigure(canvas_window, width=inside_canvas.winfo_width())

        self.inside_frame.bind('<Configure>', update_inside_frame)
        inside_canvas.bind('<Configure>', update_inside_canvas)




class tooltip():
    '''
    general tooltip class for tkinter widgets

    *args

        tkwidget = the tkwidget the tool tip is attached too
        text = text to display
        wraplength = how many characters to wrap
        timer = how long does the mouse have to hover til the popup appears
        bg = background color

    '''

    def __init__(self, tkwidget, text='info', wraplength=200, timer=1000, bg='#fcf2de'):
        self.text = text
        self.tkwidget = tkwidget
        self.timer = timer #ms
        self.tkwidget.bind("<Enter>", self.enter_widget)
        self.tkwidget.bind("<Leave>", self.leave_widget)
        self.event_id = None
        self.tool_window = None
        self.justify = 'left'
        self.bg = bg
        self.wrap_len = wraplength

    def enter_widget(self,event):
        '''
        start the timer when the mouse enters the widget
        '''
        self.start_timer()


    def leave_widget(self, event):
        '''
        stop the timer if the user removes the mouse before the timer ends
        '''
        self.stop_timer()
        self.destroy_tip()

    def start_timer(self):
        self.stop_timer()
        id = self.event_id

        self.event_id = None
        if id == None:
            self.event_id = self.tkwidget.after(self.timer, self.popup_tooltip)


    def stop_timer(self):
        id = self.event_id
        self.event_id = None
        if id:
            self.event_id = self.tkwidget.after_cancel(id)

    def popup_tooltip(self):
        '''
        popup a new small window on which to place text
        '''
        xy = [0,0]
        xy[0] += self.tkwidget.winfo_rootx()+ self.tkwidget.winfo_width() -9
        xy[1] += self.tkwidget.winfo_rooty()+ self.tkwidget.winfo_height() -9


        self.tool_window = tk.Toplevel(self.tkwidget)
        self.tool_window.overrideredirect(True)
        self.tool_window.config(bg='#FFFFFF')
        self.tool_window.geometry("+%d+%d"%(xy[0], xy[1]))

        tt_label = tk.Label(self.tool_window, text=self.text, justify=self.justify, bg = self.bg, bd=1, wraplength=self.wrap_len)
        tt_label.pack(padx=1, pady=1, fill=tk.BOTH)

    def destroy_tip(self):
        '''
        destroy the tooltip window
        '''
        try:
            self.tool_window.destroy()
        except:
            pass







class GUI(Frame):

    def __init__(self, parent): #define the intial setup of the GUI


        Frame.__init__(self, parent)    #hold overs for now from old gui
        self.parent = parent
        self.initUI()
        self.sms = rSNAPsim.rSNAPsim()
        self.gb_rec = None
        self.seqtext = None
        self.seq_map_open = False
        self.seq_text_open = False
        self.ssas=[]
        self.ssasaves = []


        self.icons = 1

        self.icondata =tk.PhotoImage(data='''
                        R0lGODlhIwAhAHAAACH5BAEAAPwALAAAAAAjACEAhwAAAAAAMwAAZgAAmQAAzAAA/
                        wArAAArMwArZgArmQArzAAr/wBVAABVMwBVZgBVmQBVzABV/wCAAACAMwCAZgCAmQ
                        CAzACA/wCqAACqMwCqZgCqmQCqzACq/wDVAADVMwDVZgDVmQDVzADV/wD/AAD/MwD
                        /ZgD/mQD/zAD//zMAADMAMzMAZjMAmTMAzDMA/zMrADMrMzMrZjMrmTMrzDMr/zNV
                        ADNVMzNVZjNVmTNVzDNV/zOAADOAMzOAZjOAmTOAzDOA/zOqADOqMzOqZjOqmTOqz
                        DOq/zPVADPVMzPVZjPVmTPVzDPV/zP/ADP/MzP/ZjP/mTP/zDP//2YAAGYAM2YAZm
                        YAmWYAzGYA/2YrAGYrM2YrZmYrmWYrzGYr/2ZVAGZVM2ZVZmZVmWZVzGZV/2aAAGa
                        AM2aAZmaAmWaAzGaA/2aqAGaqM2aqZmaqmWaqzGaq/2bVAGbVM2bVZmbVmWbVzGbV
                        /2b/AGb/M2b/Zmb/mWb/zGb//5kAAJkAM5kAZpkAmZkAzJkA/5krAJkrM5krZpkrm
                        ZkrzJkr/5lVAJlVM5lVZplVmZlVzJlV/5mAAJmAM5mAZpmAmZmAzJmA/5mqAJmqM5
                        mqZpmqmZmqzJmq/5nVAJnVM5nVZpnVmZnVzJnV/5n/AJn/M5n/Zpn/mZn/zJn//8wA
                        AMwAM8wAZswAmcwAzMwA/8wrAMwrM8wrZswrmcwrzMwr/8xVAMxVM8xVZsxVmcxVzMx
                        V/8yAAMyAM8yAZsyAmcyAzMyA/8yqAMyqM8yqZsyqmcyqzMyq/8zVAMzVM8zVZszVmcz
                        VzMzV/8z/AMz/M8z/Zsz/mcz/zMz///8AAP8AM/8AZv8Amf8AzP8A//8rAP8rM/8rZv8
                        rmf8rzP8r//9VAP9VM/9VZv9Vmf9VzP9V//+AAP+AM/+AZv+Amf+AzP+A//+qAP+qM/+
                        qZv+qmf+qzP+q///VAP/VM//VZv/Vmf/VzP/V////AP//M///Zv//mf//zP///wAAAAA
                        AAAAAAAAAAAj/ANUIHEiwoMGDCBMqXMiwocOHEAXu2xex4cSJAComvLgPgMeMGg1i/Ph
                        xI0eKJkl6DBVqI8mJCEeuDAWK5UGZH1EWlMmyJs2WIlV6vEhQZs2joUyF0imx40qVJ53
                        OPIUUFMymJGsKJQkPHgCfYFkCmJgGZ89QW8vBUxevK7yfoKiiHYsTwE+xJMmVU3tR3du
                        wacuRewsXr9py6jh2PQtKqF548fz+vbu3nL6TfgvP9bh3rWKWcvWSo3cSsrpQcu+GQtz
                        1M9LK+fIlu7y2K1LGaiNfnni6sF6/ktctPvtTaee18/T5vb26MjzhXU/fDvu7reSkd0E
                        JPi6cJfaaSn02ble7GKxczs4Ze6+Kdvu39YxnFg57V+lH+OxXTj9FM3XVmeD1B5Z+crE
                        XoE/3sVfcXEj5R9xd8imlmVZ2zQcXVbfpZyFSHjnYn2pxUUhcTQVuaKJ4KJqI2m0Shjd
                        idqm5aOGMIcJooI31ZadjKAEBADs=''')






    def initUI(self):

        ############################
        # This command sets up the buttons and user interface
        ############################
        self.parent.title("rSNAPsim")      #Set the title
        self.pack(fill=tk.BOTH, expand=True)     #expand the gui to full window

        self.icondict = {'popout':tk.PhotoImage(data='''R0lGODlhCwALAHAAACwAAAAACwALAIcAAAAAADMAAGYAAJkAAMwAAP8AKwAAKzMAK2YAK5kAK8wAK
                             /8AVQAAVTMAVWYAVZkAVcwAVf8AgAAAgDMAgGYAgJkAgMwAgP8AqgAAqjMAqmYAqpkAqswAqv8A1Q
                             AA1TMA1WYA1ZkA1cwA1f8A/wAA/zMA/2YA/5kA/8wA//8zAAAzADMzAGYzAJkzAMwzAP8zKwAzKzM
                             zK2YzK5kzK8wzK/8zVQAzVTMzVWYzVZkzVcwzVf8zgAAzgDMzgGYzgJkzgMwzgP8zqgAzqjMzqmYz
                             qpkzqswzqv8z1QAz1TMz1WYz1Zkz1cwz1f8z/wAz/zMz/2Yz/5kz/8wz//9mAABmADNmAGZmAJlmA
                             MxmAP9mKwBmKzNmK2ZmK5lmK8xmK/9mVQBmVTNmVWZmVZlmVcxmVf9mgABmgDNmgGZmgJlmgMxmgP
                             9mqgBmqjNmqmZmqplmqsxmqv9m1QBm1TNm1WZm1Zlm1cxm1f9m/wBm/zNm/2Zm/5lm/8xm//+ZAAC
                             ZADOZAGaZAJmZAMyZAP+ZKwCZKzOZK2aZK5mZK8yZK/+ZVQCZVTOZVWaZVZmZVcyZVf+ZgACZgDOZ
                             gGaZgJmZgMyZgP+ZqgCZqjOZqmaZqpmZqsyZqv+Z1QCZ1TOZ1WaZ1ZmZ1cyZ1f+Z/wCZ/zOZ/2aZ/
                             5mZ/8yZ///MAADMADPMAGbMAJnMAMzMAP/MKwDMKzPMK2bMK5nMK8zMK//MVQDMVTPMVWbMVZnMVc
                             zMVf/MgADMgDPMgGbMgJnMgMzMgP/MqgDMqjPMqmbMqpnMqszMqv/M1QDM1TPM1WbM1ZnM1czM1f/
                             M/wDM/zPM/2bM/5nM/8zM////AAD/ADP/AGb/AJn/AMz/AP//KwD/KzP/K2b/K5n/K8z/K///VQD/
                             VTP/VWb/VZn/Vcz/Vf//gAD/gDP/gGb/gJn/gMz/gP//qgD/qjP/qmb/qpn/qsz/qv//1QD/1TP/1
                             Wb/1Zn/1cz/1f///wD//zP//2b//5n//8z///8AAAAAAAAAAAAAAAAILAD3CRxIUCCAgwD2IUyoEK
                             HBgw8HLozYMCFDiA4JYrzIsSDEigs/hpxYsGRAADs='''),
                        'pause':tk.PhotoImage(data='''R0lGODlhDwAPAPcAAAAAAAAAMwAAZgAAmQAAzAAA/wArAAArMwArZgArmQArzAAr/wBVAABV
                        MwBVZgBVmQBVzABV/wCAAACAMwCAZgCAmQCAzACA/wCqAACqMwCqZgCqmQCqzACq/wDVAADVMwDVZgDVmQD
                        VzADV/wD/AAD/MwD/ZgD/mQD/zAD//zMAADMAMzMAZjMAmTMAzDMA/zMrADMrMzMrZjMrmTMrzDMr/zNVAD
                        NVMzNVZjNVmTNVzDNV/zOAADOAMzOAZjOAmTOAzDOA/zOqADOqMzOqZjOqmTOqzDOq/zPVADPVMzPVZjPVm
                        TPVzDPV/zP/ADP/MzP/ZjP/mTP/zDP//2YAAGYAM2YAZmYAmWYAzGYA/2YrAGYrM2YrZmYrmWYrzGYr/2ZV
                        AGZVM2ZVZmZVmWZVzGZV/2aAAGaAM2aAZmaAmWaAzGaA/2aqAGaqM2aqZmaqmWaqzGaq/2bVAGbVM2bVZmb
                        VmWbVzGbV/2b/AGb/M2b/Zmb/mWb/zGb//5kAAJkAM5kAZpkAmZkAzJkA/5krAJkrM5krZpkrmZkrzJkr/5
                        lVAJlVM5lVZplVmZlVzJlV/5mAAJmAM5mAZpmAmZmAzJmA/5mqAJmqM5mqZpmqmZmqzJmq/5nVAJnVM5nVZ
                        pnVmZnVzJnV/5n/AJn/M5n/Zpn/mZn/zJn//8wAAMwAM8wAZswAmcwAzMwA/8wrAMwrM8wrZswrmcwrzMwr
                        /8xVAMxVM8xVZsxVmcxVzMxV/8yAAMyAM8yAZsyAmcyAzMyA/8yqAMyqM8yqZsyqmcyqzMyq/8zVAMzVM8z
                        VZszVmczVzMzV/8z/AMz/M8z/Zsz/mcz/zMz///8AAP8AM/8AZv8Amf8AzP8A//8rAP8rM/8rZv8rmf8rzP
                        8r//9VAP9VM/9VZv9Vmf9VzP9V//+AAP+AM/+AZv+Amf+AzP+A//+qAP+qM/+qZv+qmf+qzP+q///VAP/VM
                        //VZv/Vmf/VzP/V////AP//M///Zv//mf//zP///wAAAAAAAAAAAAAAACH5BAEAAPwALAAAAAAPAA8AAAg6
                        APcJHEiwoMGDBQEoHKgQgMGGDBcmlLgP4kSHAi0S1KgxIsaKFD2KvDhyI8WOGU+GTPkRJcKXMPcFBAA7'''),
                        'play':tk.PhotoImage(data='''R0lGODlhDwAPAHAAACwAAAAADwAPAIcAAAAAADMAAGYAAJkAAMwAAP8AKwAAKzMAK2YAK5kAK
                        8wAK/8AVQAAVTMAVWYAVZkAVcwAVf8AgAAAgDMAgGYAgJkAgMwAgP8AqgAAqjMAqmYAqpkAqswAqv8A1QAA
                        1TMA1WYA1ZkA1cwA1f8A/wAA/zMA/2YA/5kA/8wA//8zAAAzADMzAGYzAJkzAMwzAP8zKwAzKzMzK2YzK5k
                        zK8wzK/8zVQAzVTMzVWYzVZkzVcwzVf8zgAAzgDMzgGYzgJkzgMwzgP8zqgAzqjMzqmYzqpkzqswzqv8z1Q
                        Az1TMz1WYz1Zkz1cwz1f8z/wAz/zMz/2Yz/5kz/8wz//9mAABmADNmAGZmAJlmAMxmAP9mKwBmKzNmK2ZmK
                        5lmK8xmK/9mVQBmVTNmVWZmVZlmVcxmVf9mgABmgDNmgGZmgJlmgMxmgP9mqgBmqjNmqmZmqplmqsxmqv9m
                        1QBm1TNm1WZm1Zlm1cxm1f9m/wBm/zNm/2Zm/5lm/8xm//+ZAACZADOZAGaZAJmZAMyZAP+ZKwCZKzOZK2a
                        ZK5mZK8yZK/+ZVQCZVTOZVWaZVZmZVcyZVf+ZgACZgDOZgGaZgJmZgMyZgP+ZqgCZqjOZqmaZqpmZqsyZqv
                        +Z1QCZ1TOZ1WaZ1ZmZ1cyZ1f+Z/wCZ/zOZ/2aZ/5mZ/8yZ///MAADMADPMAGbMAJnMAMzMAP/MKwDMKzPMK
                        2bMK5nMK8zMK//MVQDMVTPMVWbMVZnMVczMVf/MgADMgDPMgGbMgJnMgMzMgP/MqgDMqjPMqmbMqpnMqszMq
                        v/M1QDM1TPM1WbM1ZnM1czM1f/M/wDM/zPM/2bM/5nM/8zM////AAD/ADP/AGb/AJn/AMz/AP//KwD/KzP/K
                        2b/K5n/K8z/K///VQD/VTP/VWb/VZn/Vcz/Vf//gAD/gDP/gGb/gJn/gMz/gP//qgD/qjP/qmb/qpn/qsz/q
                        v//1QD/1TP/1Wb/1Zn/1cz/1f///wD//zP//2b//5n//8z///8AAAAAAAAAAAAAAAAIMQD3CRxIsKDBgwQBK
                        ESYUOFChg4jIoxI0SDFiw0vShSocSPHjhYxHqwI8SHDkygRBgQAOw==''')

                        }


        #small safer evaluation of expressions for entry boxes

        #allow the following operators:
        ops = {ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor, ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
                     ast.USub: op.neg, ast.RShift: op.rshift, ast.LShift: op.lshift, ast.Eq: op.eq, ast.NotEq: op.ne, ast.Lt: op.lt, ast.LtE: op.le, ast.Gt: op.gt, ast.GtE: op.ge}

        def eval_expression(expression):
            return safer_eval(ast.parse(expression, mode='eval').body)

        def safer_eval(node):
            if isinstance(node, ast.Num):
                return node.n

            #evaluate expressions like 1 num = op 1 num
            elif isinstance(node, ast.BinOp):
                return ops[type(node.op)](safer_eval(node.left), safer_eval(node.right))
            #evaluate negatives
            elif isinstance(node, ast.UnaryOp):
                return ops[type(node.op)](safer_eval(node.operand))

            #evaluate booleans
            elif isinstance(node, ast.BoolOp):
                return ops[type(node.op)](safer_eval(node.left), safer_eval(node.right))

            #evaluate comparisons? why are you multiplying values by a true false
            elif isinstance(node, ast.Compare):

                return ops[type(node.ops[0])](safer_eval(node.left), safer_eval(node.comparators[0]))

            else:
                raise TypeError(node)





        self.tooltips = {'Select File Sequence':['Select a text or genbank (.gb) file to open'],
                         'Select POI (multiple detected POI)':['If there are multiple possible proteins detected, select between the proteins of interest by sequence location with the drop down menu'],
                                     'Colorize':['Color the gene sequence, this will take some time'],
                                     'popout':['popout the plot to a new figure'],
                                     'Gene Sequence':['Displays the gene sequence of nucleotides'],
                                     'Sequence Viewer':['Displays the selected tag in white and the protein of interest in green. The vertical orange bars are locations of tag epitopes'],
                                     'Simulation Time':['Time to run the simulation in seconds'],
                                     'Repetitions':['Number of times to run the simulation'],
                                     'Exp. Elongation Rate':['Experimental or Expected elongation rate to normalize the rates based on tRNA copy numbers'],
                                     'Exp. Initiation Rate':['Experimental or Expected initiation rate of ribosomes initiating translation'],
                                     'Time Courses':['Plots of all the run simulation trajectories (intensities of fluorescence)'],
                                     'Simulated Cell':['Gif display of a simulated cell with current parameters'],
                                     'Simulate':['Open dialog to generate a new simulation file'],
                                     'simfile':['Name of the simulated cell gif file to play'],
                                     'Autocorrelation':['Autocorrelation of fluorescence'],
                                     'Ribosome Density':['Ribosome density by codon position, displayed as probability of finding a ribosome at that codon position'],
                                     'Codon Usage':['Codon usage by position on the nucleotide sequence'],
                                     'DC (um^2/s)':['Diffusion Constant for brownian motion'],
                                     'kRNA':['Rate of RNA creation'],
                                     'kOn_RNA':['Rate for RNA to begin translation'],
                                     'kOff_RNA':['Rate for RNA to turn off (not used currently)'],
                                     'kDecay_RNA':['Rate for RNA to be degraded'],
                                     'Cell size':['Size of the cell in micrometers'],
                                     'Perturbation Time':['The time at which the perturbation is applied to the simulation'],
                                     'Harringtonine':['A molecule treatment that stops initiation of ribosomes into translation but does not inhibit ribosomes currently translating'],
                                     'FRAP':['Fluorescence Recovery After Photobleaching, manual photobleaching of all current fluorophores, but new fluorophores that fold or are translated are not affected'],
                                     'Dynamic Property Values':['Calculated parameters based on the perturbation simulation'],
                                     'Sensitivity':['A measurement of how much changing the codons copy number would affect fluorescence'],
                                     'Dynamic Properties':['Calculated parameters from the simulated wild type and modified sequences autocorrelations'],
                                     'Autocorrelations_codon':['The autocorrelations of the Wild Type sequence and the Modified sequence'],


}




        #plots have different background colors by default, so setting the face color of the plot to the default will make it look
        #like it blends in to the application, try to set to windows default, if not make a button and grab its default background
        try:

            color = self.parent.winfo_rgb('systembuttonface')  #get the default face color for blank windows
            dc=[0,0,0]  #convert that face color to rgb hex
            dc[0] = hex(int(color[0]/256))[-2:]
            dc[1] = hex(int(color[1]/256))[-2:]
            dc[2] = hex(int(color[2]/256))[-2:]
            self.default_color = '#' + dc[0] + dc[1] + dc[2]
        except:
            self.default_color = '#ececec'

            colorbutton = tk.Button(self.parent, text='color')
            color = colorbutton.cget('bg')
            color = self.parent.winfo_rgb(color)
            colorbutton.pack()
            colorbutton.destroy()

            dc=[0,0,0]  #convert that face color to rgb hex
            dc[0] = hex(int(color[0]/256))[-2:]
            dc[1] = hex(int(color[1]/256))[-2:]
            dc[2] = hex(int(color[2]/256))[-2:]
            self.default_color = '#' + dc[0] + dc[1] + dc[2]

        #self.tk_setPalette(background = self.default_color,disabledforeground='#aaaaaa',inactivebackground='#aaaaaa',disabledbackground='#aaaaaa')
        self.main_color = '#fcaa0a'
        self.global_dpi = 70 #global dpi setting for plots
        self.loaded_gif = ''

        menubar = tk.Menu(self.parent) #the top file, info dropdowns

        self.parent.config(menu=menubar)
        filemenu = tk.Menu(menubar, tearoff=0)

        filemenu.add_separator()
        filemenu.add_command(label="Exit",command=self.parent.quit)
        menubar.add_cascade(label='File',menu=filemenu)

        self.steady_state = tk.BooleanVar()
        self.steady_state.set(True)
        settingsmenu = tk.Menu(menubar, tearoff=0)
        settingsmenu.add_command(label="Change tRNA Copynumbers", command=self.change_trna_copys)
        #settingsmenu.add_checkbutton(label="Steady State", onvalue=True, offvalue=False, variable=self.steady_state)

        menubar.add_cascade(label='Settings', menu=settingsmenu)


        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Guide", command=self.open_guide)

        filemenu.add_command(label="Info", command=self.open_info)
        menubar.add_cascade(label='Help', menu=filemenu)

        style  = self.initstyle()                #get the custom style
        style.theme_use("new")                  #set the style as the custom
        self.Main_Nb = ttk.Notebook(self)       #set up the notebook tabs
        self.Main_Nb.pack(fill='both', expand=False)

        seq_frame = tk.Frame(self.Main_Nb, name="seq")  #first notebook Frame
        seq_frame.pack(expand=True, side='top', fill='x')
#
#        ss_frame = tk.Frame(self.Main_Nb,name="ss")   #SSA notebook Frame
#        ss_frame.pack(expand=True,side='top',fill='both')

        #ss_frame.bind("<Visibility>", self.__update_ss_frame)

        newseq_frame = tk.Frame(self.Main_Nb,name="nseq")  #first notebook Frame
        newseq_frame.pack(expand=True, side='top', fill='x')


        newcodon_frame = tk.Frame(self.Main_Nb,name="ncodon")
        newcodon_frame.pack(expand=True, side='top', fill='x')

        runsims_frame = tk.Frame(self.Main_Nb,name="psims")
        runsims_frame.pack(expand=True, side='top', fill='x')

        #elon_frame.bind("<Visibility>", self.__update_elon_frame)

        #codon_frame = tk.Frame(self.Main_Nb,name="codon")        #codon notebook Frame
        #codon_frame.pack(expand=True,side='top',fill='both')
        #codon_frame.bind("<Visibility>", self.__update_codon_frame)


        self.Main_Nb.add(newseq_frame, text="   Sequence Selector   ")
        self.Main_Nb.add(seq_frame, text="   Sequence Properties   ")
        self.Main_Nb.add(newcodon_frame, text= "  Codon Adjustments  ")
        self.Main_Nb.add(runsims_frame, text= "   Previous Simulations  " )
       # self.Main_Nb.add(ss_frame,text = "   Sequence Properties  ")
        #self.Main_Nb.add(elon_frame,text= "   Elongation Assays   ")
        #self.Main_Nb.add(codon_frame,text="   Codon Optimization  ")

        self.stoc_Nb = ttk.Notebook(self)
        self.stoc_Nb.pack(expand=True,side='top', fill='both')


        ss_frame = tk.Frame(self.stoc_Nb,name="ss")   #SSA notebook Frame
        ss_frame.pack(expand=True,side='top', fill='both')

        #elon_frame = tk.Frame(self.stoc_Nb,name="elon") #elongation notebook Frame
       # elon_frame.pack(expand=True,side='top',fill='both')

        fcs_frame = tk.Frame(self.stoc_Nb,name="fcs") #elongation notebook Frame
        fcs_frame.pack(expand=True,side='top', fill='both')

        sim_frame = tk.Frame(self.stoc_Nb,name="sim") #elongation notebook Frame
        sim_frame.pack(expand=True,side='top', fill='both')
        
        kym_frame = tk.Frame(self.stoc_Nb,name='kym')
        kym_frame.pack(expand=True,side='top', fill='both')


        #ccodon_frame = tk.Frame(self.stoc_Nb,name="codon")        #codon notebook Frame
        #ccodon_frame.pack(expand=True,side='top',fill='both')

        self.stoc_Nb.add(ss_frame, text="  Stochastic Simulations  ")
        self.stoc_Nb.add(fcs_frame, text="   FCS  ")
        #self.stoc_Nb.add(elon_frame,text= "   Elongation Assays   ")
        #self.stoc_Nb.add(ccodon_frame,text="   Codon Optimization  ")
        
        self.stoc_Nb.add(kym_frame, text="   Kymographs   ")
        self.stoc_Nb.add(sim_frame, text="   Simulated Cell  ")
        
        global_font_size = 9

        '''
        Previous simulation frame
        '''
        self.psimscroll = ScrollingFrame(runsims_frame,height=100)
        self.psimscroll.pack(anchor=tk.NW,expand=True,fill='both',side=tk.LEFT)

        loadbutton = tk.Button(runsims_frame,text='Load',command= self.load_ss_data,font=('SystemButtonText',global_font_size))
        loadbutton.pack(side=tk.LEFT)



        


        '''
        Sequence Frame

        This frame contains the information from the loaded sequence file / text
        '''

        leftcodonframe = tk.Frame(newcodon_frame)
        leftcodonframe.pack(side=tk.LEFT)

        rightcodonframe=tk.Frame(newcodon_frame)
        rightcodonframe.pack(side=tk.LEFT)
        outside_codon_frame_top = tk.Frame(leftcodonframe)
        outside_codon_frame_top.pack(anchor=tk.NW)
        codonscroll = ScrollingFrame(leftcodonframe,height=100)
        codonscroll.pack(anchor=tk.NW,expand=True,fill='both')

        codon_frame = codonscroll.inside_frame

        outside_codon_frame_bottom = tk.Frame(leftcodonframe)
        outside_codon_frame_bottom.pack(anchor=tk.NW)
        gpy = 2
        gpx = 2

        #gs_input_l = tk.Label(codon_frame,text='Codon Optimization Design',font=('SystemButtonText',13,'bold'),bg='#888888',fg='#FFFFFF')
        #gs_input_l.grid(row=4,column=0,columnspan=10,sticky=tk.W+tk.E+tk.N+tk.S)

        aa_key = ['Ala (A)','Arg (R)','Asn (N)','Asp (D)','Cys (C)','Gln (Q)','Glu (E)','Gly (G)','His (H)','Ile (I)','Leu (L)','Lys (K)','Met (M)','Phe (F)','Pro (P)','Ser (S)','Thr (T)','Trp (W)','Tyr (Y)','Val (V)']
        self.aa_entries = []
        self.aa_vars = []
        self.sensitivity_entries = []
        aa_options = [('Natural',0),('Common',1),('Rare',2)]


        rbframe = tk.Frame(codon_frame)
        rbframe.grid(row=5,column=0,columnspan=10)
        '''
        titles = ['AA','No.','WT','Common','Rare','Sensitivity']
        for i in range(len(titles)):
            l = tk.Label(outside_codon_frame_top,text=titles[i],font=('SystemButtonText',8,'bold'))
            l.grid(row=1,column=i,padx=2,pady=2)
            if i == 5:
                ss_tt = tooltip(l,self.tooltips['Sensitivity'][0])

        l = tk.Label(outside_codon_frame_top,text='            ')
        l.grid(row=1,column=6,pady=2)
        '''

        '''
        for i in range(len(titles)):

            l = tk.Label(outside_codon_frame_top,text=titles[i],font=('SystemButtonText',10,'bold'))
            l.grid(row=1,column=i+7,padx=2,pady=2)

            if i == 5:
                ss2_tt = tooltip(l,self.tooltips['Sensitivity'][0])
        '''

        k = 0
        for aa in aa_key:
            k+=1
            if k < 11:
                if k%2 == 0:
                    cc = '#DDDDDD'
                else:
                    cc = self.default_color

                cl = tk.Label(rbframe,text=' ',bg=cc)
                cl.grid(row=5+k,column=0,columnspan=8,sticky=tk.N+tk.E+tk.W+tk.S)

                l = tk.Label(rbframe,text=aa,font=('SystemButtonText',8,'bold'),bg=cc)
                l.grid(row=5+k,column=0)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+k,column=1,pady=0)
                self.aa_entries.append(e)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+k,column=5,pady=0)
                self.sensitivity_entries.append(e)


                v = tk.IntVar(value=0)
                j=0
                m = 0
                mo = ['WT','Common','Rare']
                for text,mode in aa_options:

                    b = tk.Radiobutton(rbframe,text=mo[m],variable=v,value=mode,indicatoron=1,bg=cc,font=('SystemButtonText',global_font_size))
                    m+=1
                    b.grid(row=5+k,column=2+j,pady=0)
                    j+=1
                self.aa_vars.append(v)

                #spacer = tk.Label(rbframe,text='            ')
                #spacer.grid(row=5+k,column=6,pady=0)

            if k >=11:
                if k%2 == 0:
                    cc = '#DDDDDD'
                else:
                    cc = self.default_color

                cl = tk.Label(rbframe,text=' ',bg=cc)
                cl.grid(row=5+k,column=0,columnspan=10,sticky=tk.N+tk.E+tk.W+tk.S)

                l = tk.Label(rbframe,text=aa,font=('SystemButtonText',8,'bold'),bg=cc)
                l.grid(row=5+(k),column=0)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+(k),column=1,pady=0)
                self.aa_entries.append(e)

                v = tk.IntVar(value=0)
                j = 0
                m = 0
                mo = ['WT','Common','Rare']
                for text,mode in aa_options:
                    b = tk.Radiobutton(rbframe,text=mo[m],variable=v,value=mode,indicatoron=1,bg=cc,font=('SystemButtonText',global_font_size))
                    m+=1
                    b.grid(row=5+(k),column=2+j,pady=0)
                    j+=1
                self.aa_vars.append(v)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+(k),column=5,pady=0)
                self.sensitivity_entries.append(e)

        #l = tk.Label(outside_codon_frame_bottom,text='* Sensitivity of the codon manipulation can vary between 0 and 1, where 1 is the most sensitive amino acid for manipulation',font=('SystemButtonText',8,))
        #l.grid(row=0,column=0,columnspan=8)

        allnatural = tk.Button(rightcodonframe,text='All Wild Type Codons',command=self.set_natural_codons,font=('SystemButtonText',global_font_size))
        allnatural.grid(row=0,column=1,columnspan=3,pady=gpy,padx=gpx)

        allcommon = tk.Button(rightcodonframe,text='All Common Codons',command=self.set_common_codons,font=('SystemButtonText',global_font_size))
        allcommon.grid(row=1,column=1,columnspan=4,pady=gpy,padx=gpx)

        allrare = tk.Button(rightcodonframe,text='All Rare Codons',command=self.set_rare_codons,font=('SystemButtonText',global_font_size))
        allrare.grid(row=2,column=1,columnspan=4,pady=gpy,padx=gpx)


        fullsizecodon = tk.Button(rightcodonframe,text= 'Full size Codon editor',command=self.fullsize_codon_popout,font=('SystemButtonText',global_font_size))
        fullsizecodon.grid(row=1,column=0,pady=gpy,padx=gpx)

        edittrna = tk.Button(rightcodonframe,text='Edit tRNA Abundance',command=self.change_trna_copys,font=('SystemButtonText',global_font_size))
        edittrna.grid(row=0,column=0,pady=gpy,padx=gpx)



        gpx = 2
        gpy = 2
        #gs_input_l = tk.Label(newseq_frame,text='Gene Sequence Input',font=('SystemButtonText',13,'bold'),bg='#888888',fg='#FFFFFF')
        #gs_input_l.grid(row=0,column=0,columnspan=12,sticky=tk.W+tk.E+tk.N+tk.S)
        gbid_label = tk.Label(newseq_frame,text='GenBank ID',font=('SystemLabelText',global_font_size))
        gbid_label.grid(row=1,column=0,sticky=tk.E)

        gbid_tt = tooltip(gbid_label,text='GenBank Accession number, If you would like to poll genbank please input a valid accession number for your gene sequence',wraplength=200)



        #seq_frame.grid_columnconfigure(9,weight=1)

        self.gbid_entry = tk.Entry(newseq_frame,width=30,justify=tk.RIGHT)   #GeneBank ID entry
        self.gbid_entry.grid(row=1,column=1,sticky=tk.W,padx=gpx,pady=gpy)
        self.gbid_entry.bind('<FocusIn>',self.__bind_genbank)
        self.gbid_entry.bind('<FocusOut>',self.__unbind_genbank)


        ginput_tt = tooltip(self.gbid_entry,text='GenBank Accession number, If you would like to poll genbank please input a valid accession number for your gene sequence')


        pullasc = tk.Button(newseq_frame,text='Pull',command=self.poll_genebank,font=('SystemButtonText',global_font_size))
        pullasc.grid(row=1,column=2,sticky=tk.W,padx=gpx,pady=gpy,rowspan=2)

        self.H2B = tk.IntVar(value=0)  #checkbutton for H2B variable
        h2bcheck = tk.Checkbutton(newseq_frame,text='Use Default (H2B)?',variable=self.H2B)
        #h2bcheck.grid(row=1,column=2,sticky=tk.W,padx=gpx,pady=gpy)

        orl = tk.Label(newseq_frame,text='     OR     ',font=('SystemButtonText',14,'bold'))
        orl.grid(row=1,column=4,padx=gpx,pady=gpy,rowspan=2)

        seqbutton = tk.Button(newseq_frame,text='Select File Sequence',command=self.open_seq_file,font=('SystemButtonText',global_font_size))
        seqbutton.grid(row=1,column=5,padx=gpx,pady=gpy,rowspan=2)

        seqbutton_tt = tooltip(seqbutton,self.tooltips['Select File Sequence'][0])

        self.seqdisp = tk.Entry(newseq_frame,width=30)    #sequence display for just loading a text file
        self.seqdisp.grid(row=1,column=6,padx=gpx,pady=gpy,rowspan=2)
        self.seqdisp.config(state='disabled')

        #set up the gene sequence properties labels and entries





        topframe = tk.Frame(newseq_frame)
        topframe.grid(row=3,column=0,columnspan=13,sticky=tk.W+tk.E+tk.N+tk.S)



        gbid = tk.Label(topframe,text='Simulation Time',font=('SystemLabelText',global_font_size))
        gbid.grid(row=1,column=7,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid5_tt = tooltip(gbid,self.tooltips['Simulation Time'][0])

        gbid = tk.Label(topframe,text='Repetitions',font=('SystemLabelText',global_font_size))
        gbid.grid(row=3,column=7,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid4_tt = tooltip(gbid,self.tooltips['Repetitions'][0])

        gbid = tk.Label(topframe,text='Perturbation Time',font=('SystemLabelText',global_font_size))
        gbid.grid(row=2,column=7,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid3_tt = tooltip(gbid,self.tooltips['Perturbation Time'][0])


        pbte = tk.Entry(topframe,justify=tk.RIGHT)
        pbte.grid(row=2,column=9,sticky=tk.W,padx=gpx,pady=gpy,columnspan=1)

        sste = tk.Entry(topframe,justify=tk.RIGHT)
        sste.grid(row=1,column=9,sticky=tk.W,padx=gpx,pady=gpy,columnspan=1)

        repe = tk.Entry(topframe,justify=tk.RIGHT)
        repe.grid(row=3,column=9,sticky=tk.W,padx=gpx,pady=gpy,columnspan=2)

        #self.ss_gs_inputs = [sste,pbte,repe]

        sec = tk.Label(topframe,text='sec',font=('SystemLabelText',global_font_size))
        sec.grid(row=1,column=10,sticky=tk.W)
        sec = tk.Label(topframe,text='sec',font=('SystemLabelText',global_font_size))
        sec.grid(row=2,column=10,sticky=tk.W)

        run_ss = tk.Button(topframe,text='Run',command=self.run_ssa,font=('SystemButtonText',14))
        run_ss.grid(row=1,column=13,rowspan=3,sticky=tk.W,padx=gpx,pady=gpy)



        gl = tk.Label(topframe,text='Burn In',font=('SystemLabelText',global_font_size))
        gl.grid(row=1,column=3,sticky=tk.E,padx=gpx,pady=gpy)

        gl = tk.Label(topframe,text='sec',font=('SystemLabelText',global_font_size))
        gl.grid(row=1,column=5,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(topframe,text='Exp. Elongation Rate',font=('SystemLabelText',global_font_size))
        gl.grid(row=2,column=3,sticky=tk.E,padx=gpx,pady=gpy)

        gl3_tt = tooltip(gl,self.tooltips['Exp. Elongation Rate'][0])

        gl = tk.Label(topframe,text='aa/sec',font=('SystemLabelText',global_font_size))
        gl.grid(row=2,column=5,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(topframe,text='Exp. Initiation Rate',font=('SystemLabelText',global_font_size))
        gl.grid(row=3,column=3,sticky=tk.E,padx=gpx,pady=gpy)

        gl4_tt = tooltip(gl,self.tooltips['Exp. Initiation Rate'][0])

        gl = tk.Label(topframe,text='1/sec',font=('SystemLabelText',global_font_size))
        gl.grid(row=3,column=5,sticky=tk.W,padx=gpx,pady=gpy)

        gl_e = tk.Entry(newseq_frame,state=tk.DISABLED,justify=tk.RIGHT)
        gl_e.grid(row=1,column=7,padx=3,pady=8)
        
        backend_entry = tk.Entry(newseq_frame,justify=tk.RIGHT,width=8)
        backend_entry.grid(row=1,column=8,padx=3,pady=8)
        
        try:
            ssa_translation.__file__
            backend_entry.insert(0,'C++')
            
        except:
            backend_entry.insert(0,'Python')
        backend_entry.config(state='disabled')
        

        burn_e = tk.Entry(topframe,justify=tk.RIGHT)
        burn_e.grid(row=1,column=4,sticky=tk.W,padx=gpx,pady=gpy)


        eer_e = tk.Entry(topframe,justify=tk.RIGHT)
        eer_e.grid(row=2,column=4,sticky=tk.W,padx=gpx,pady=gpy)

        eir_e = tk.Entry(topframe,justify=tk.RIGHT)
        eir_e.grid(row=3,column=4,sticky=tk.W,padx=gpx,pady=gpy)

        self.ss_gs_inputs = [pbte,sste,repe,burn_e]

        self.ssa_inputs = [pbte,sste,repe,gl_e,eer_e,eir_e,burn_e]    #perturb time, simtime, repetitions, experiemental rate, initiation rate

        #topframe.grid_columnconfigure(9,weight=2)


        save_seq= tk.Button(newseq_frame,text='Save file', command = self.save_seq,font=('SystemButtonText',global_font_size))
        save_seq.grid(row=1,column=3,columnspan=1)

        select_protein = tk.Label(topframe,text='Select POI',font=('SystemLabelText',global_font_size))
        select_protein.grid(row=1,column=0,columnspan=1)

        time_res = tk.Label(topframe,text='time steps',font=('SystemLabelText',global_font_size))
        time_res.grid(row=3,column=0)

        self.time_res_e = tk.Entry(topframe,justify=tk.RIGHT)
        self.time_res_e.grid(row=3,column=1,sticky=tk.EW)
        #select_protein_tt = tooltip(select_protein,self.tooltips['Select POI (multiple detected POI)'][0])

        #self.ss_gs_inputs = [gbid,sste,repe]
        ginput_tt = tooltip(select_protein,text='If the genbank file given has more than one protein sequence, the dropdown menu provides them by site location on the gene for user selection')


        tag = tk.Label(topframe,text='Tag ',font=('SystemLabelText',global_font_size))
        tag.grid(row=2,column=0,sticky=tk.E,padx=2,pady=2)

        tagtypes = ['T_Flag','T_Flag','T_Hemagglutinin','T_SunTag']           #current filetypes for dropdown

        self.tagtype = tk.StringVar(self,value='T_Flag')
        giflabel = ttk.OptionMenu(topframe,self.tagtype,*tagtypes)    #make the dropdown
        giflabel.grid(row=2,column=1,sticky=tk.EW)
        giflabel.config(width=20)

        self.poi_list = ['---','---']
        self.poi_var = tk.StringVar(self,value='---')
        self.poi_menu = ttk.OptionMenu(topframe,self.poi_var,*self.poi_list)
        self.poi_var.trace('w',self.change_poi)
        self.poi_menu.grid(row=1,column=1,sticky=tk.EW)
        self.poi_menu.config(width=20)


        self.ss_rates = [gl_e,eer_e,eir_e]



        perturbs = [("None","none"),("FRAP","frap"),("Harringtonine","harr")]
        self.perturb = tk.StringVar(value="none")
        nonebutton = tk.Radiobutton(topframe,text='No Perturbation',value="none",variable=self.perturb,font=('SystemButtonText',global_font_size))
        frapbutton = tk.Radiobutton(topframe,text='FRAP',value="frap",variable=self.perturb,font=('SystemButtonText',global_font_size))
        harrbutton = tk.Radiobutton(topframe,text='Harringtonine',value="harr",variable=self.perturb,font=('SystemButtonText',global_font_size))

        nonebutton.grid(row=1,column=12,sticky=tk.W,padx=4)
        frapbutton.grid(row=2,column=12,sticky=tk.W,padx=4)
        harrbutton.grid(row=3,column=12,sticky=tk.W,padx=4)



        #self.ss_rates_elon = [gl_e,eer_e,eir_e]



        '''
        gpx = 4
        gpy = 6
        gs_input_l = tk.Label(seq_frame,text='Gene Sequence Input',font=('SystemButtonText',13,'bold'),bg='#888888',fg='#FFFFFF')
        gs_input_l.grid(row=0,column=0,columnspan=12,sticky=tk.W+tk.E+tk.N+tk.S)
        gbid_label = tk.Label(seq_frame,text='GenBank ID')
        gbid_label.grid(row=1,column=0,sticky=tk.E)

        gbid_tt = tooltip(gbid_label,text='GenBank Accession number, If you would like to poll genbank please input a valid accession number for your gene sequence',wraplength=200)



        seq_frame.grid_columnconfigure(9,weight=1)

        self.gbid_entry = tk.Entry(seq_frame,width=30,justify=tk.RIGHT)   #GeneBank ID entry
        self.gbid_entry.grid(row=1,column=1,sticky=tk.W,padx=gpx,pady=gpy)


        ginput_tt = tooltip(self.gbid_entry,text='GenBank Accession number, If you would like to poll genbank please input a valid accession number for your gene sequence')


        pullasc = tk.Button(seq_frame,text='Pull',command=self.poll_genebank)
        pullasc.grid(row=1,column=2,sticky=tk.W,padx=gpx,pady=gpy,rowspan=2)



        self.H2B = tk.IntVar(value=0)  #checkbutton for H2B variable
        h2bcheck = tk.Checkbutton(seq_frame,text='Use Default (H2B)?',variable=self.H2B)
        #h2bcheck.grid(row=1,column=2,sticky=tk.W,padx=gpx,pady=gpy)

        orl = tk.Label(seq_frame,text='        OR        ',font=('SystemButtonText',14,'bold'))
        orl.grid(row=1,column=3,padx=gpx,pady=gpy,rowspan=2)

        seqbutton = tk.Button(seq_frame,text='Select File Sequence',command=self.open_seq_file)
        seqbutton.grid(row=1,column=4,padx=gpx,pady=gpy,rowspan=2)

        seqbutton_tt = tooltip(seqbutton,self.tooltips['Select File Sequence'][0])

        self.seqdisp = tk.Entry(seq_frame,width=37)    #sequence display for just loading a text file
        self.seqdisp.grid(row=1,column=5,padx=gpx,pady=gpy,rowspan=2)
        self.seqdisp.config(state='disabled')

        #set up the gene sequence properties labels and entries
        seq_info= tk.Label(seq_frame,text='Sequence Information',font=('SystemButtonText',13,'bold'),bg='#888888',fg='#FFFFFF')
        seq_info.grid(row=3,column=0,columnspan=10,sticky=tk.W+tk.E+tk.N+tk.S)

        seq_prop_l = tk.Label(seq_frame,text='Sequence Properties',font=('SystemButtonText',12,'bold'))
        seq_prop_l.grid(row=5,column=0,columnspan=3,padx=gpx,pady=gpy)

        gene_seq = tk.Label(seq_frame,text='Gene Sequence',font=('SystemButtonText',12,'bold'))
        gene_seq.grid(row=5,column=3,columnspan=2,padx=gpy,pady=gpy,sticky=tk.E)

        gene_seq_tt = tooltip(gene_seq,self.tooltips['Gene Sequence'][0])

        '''



        #select_protein_tt = tooltip(select_protein,self.tooltips['Select POI (multiple detected POI)'][0])


        ginput_tt = tooltip(select_protein,text='If the genbank file given has more than one protein sequence, the dropdown menu provides them by site location on the gene for user selection')



        prop_frame = tk.Frame(seq_frame,bg='#FFFFFF')
        prop_frame.grid(row=6,column=0,columnspan=3,padx=10,pady=6,sticky=tk.W + tk.E+tk.N+tk.S,rowspan=2)





        gn = tk.Label(prop_frame,text='Gene Name',font=('SystemButtonText',9,'bold'),bg='#FFFFFF')
        tt = tk.Label(prop_frame,text='Tag Type',font=('SystemButtonText',9,'bold'),bg='#FFFFFF')
        ne = tk.Label(prop_frame,text='Number of Epitopes',font=('SystemButtonText',9,'bold'),bg='#FFFFFF')
        ep = tk.Label(prop_frame,text='Epitopes Position',font=('SystemButtonText',9,'bold'),bg='#FFFFFF')

        ne_t = tooltip(ne,text='Number of epitopes for binding on the Tag sequence')
        ep_t = tooltip(ep,text='Epitope codon locations along the Tag sequence')



        gn.grid(row=0,column=0,padx=2,pady=2,sticky=tk.E)
        tt.grid(row=1,column=0,padx=2,pady=2,sticky=tk.E)
        ne.grid(row=2,column=0,padx=2,pady=2,sticky=tk.E)
        ep.grid(row=3,column=0,padx=2,pady=2,sticky=tk.E)

        self.Seq_props = []

        gne = tk.Entry(prop_frame,width=25,bg='#FFFFFF')
        tte = tk.Entry(prop_frame,width=25,bg='#FFFFFF')
        nee = tk.Entry(prop_frame,width=25,bg='#FFFFFF')
        epe = tk.Entry(prop_frame,width=25,bg='#FFFFFF')

        gne.grid(row=0,column=1,padx=2,pady=2,sticky=tk.W)
        tte.grid(row=1,column=1,padx=2,pady=2,sticky=tk.W)
        nee.grid(row=2,column=1,padx=2,pady=2,sticky=tk.W)
        epe.grid(row=3,column=1,padx=2,pady=2,sticky=tk.W)





        prop_frame2 = tk.Frame(seq_frame,bg='#FFFFFF')
        prop_frame2.grid(row=6,column=4,columnspan=3,padx=10,pady=6,sticky=tk.W + tk.E+tk.N+tk.S,rowspan=2)

        gn = tk.Label(prop_frame2,text='            Gene Length',font=('SystemButtonText',9,'bold'),bg='#FFFFFF')
        tt = tk.Label(prop_frame2,text='         Tag Length',font=('SystemButtonText',9,'bold'),bg='#FFFFFF')
        ne = tk.Label(prop_frame2,text='         Total Length',font=('SystemButtonText',9,'bold'),bg='#FFFFFF')
        gn.grid(row=0,column=0,padx=2,pady=2,sticky=tk.E)
        tt.grid(row=1,column=0,padx=2,pady=2,sticky=tk.E)
        ne.grid(row=2,column=0,padx=2,pady=2,sticky=tk.E)

        gle = tk.Entry(prop_frame2,width=25,bg='#FFFFFF')
        tle = tk.Entry(prop_frame2,width=25,bg='#FFFFFF')
        ttle = tk.Entry(prop_frame2,width=25,bg='#FFFFFF')
        gle.grid(row=0,column=1,padx=2,pady=2,sticky=tk.E)
        tle.grid(row=1,column=1,padx=2,pady=2,sticky=tk.E)
        ttle.grid(row=2,column=1,padx=2,pady=2,sticky=tk.E)

        cd = tk.Label(prop_frame2,text='codons',font=('SystemLabelText',global_font_size),bg='#FFFFFF')
        cd.grid(row=0,column=2,padx=2,pady=2,sticky=tk.W)

        cd = tk.Label(prop_frame2,text='codons',font=('SystemLabelText',global_font_size),bg='#FFFFFF')
        cd.grid(row=1,column=2,padx=2,pady=2,sticky=tk.W)

        cd = tk.Label(prop_frame2,text='codons',font=('SystemLabelText',global_font_size),bg='#FFFFFF')
        cd.grid(row=2,column=2,padx=2,pady=2,sticky=tk.W)


        #Sequence property entries,
        #gne = Gene name, sq[0]
        #tte = Tag Type, sq[1]
        #nee = Number of epitopes, sq[2]
        #epe = epitopes position, seq[3]
        #gle = gene length, seq[4]
        #tle = tag length, seq[5]
        #ttle = total length, seq[6]

        self.Seq_props = [gne,tte,nee,epe,gle,tle,ttle]


        veiwframeshow = tk.Button(seq_frame,text='View Sequence',command=self.view_sequence_newwin,font=('SystemButtonText',global_font_size))
        veiwframeshow.grid(row=6, column=9, padx=2)

        plotseq = tk.Button(seq_frame,text='Sequence Map',command=self.seqmap,font=('SystemButtonText',global_font_size))
        plotseq.grid(row=7,column=9,padx=2)

        #Text entry for the ATA... nucleotide sequence
        self.textframe = tk.Frame()
        self.viewframe = tk.Frame()

        colorize = tk.Button(self.textframe,text='Colorize',command=self.color_text,font=('SystemButtonText',global_font_size))
        colorize.grid(row=6,column=8)

        colorize_tt = tooltip(colorize,text='Color the sequence display by G, C, T, and A')


        self.gs_text = tk.Text(self.textframe,height=20,width=71,bg='#FFFFFF',font=('SystemButtonText',8))
        self.gs_text.grid(row=0,column=0,columnspan=4,rowspan=2,sticky=tk.W+tk.E+tk.N+tk.S,padx=gpx,pady=gpy)

        gs_info = tk.Label(self.viewframe,text='Gene sequence displays the Tag region + POI',font = ('SystemButtonText',8))
        gs_info.grid(row=8,column=3,columnspan=4,sticky=tk.W+tk.E+tk.N+tk.S,padx=gpx,pady=2)

        seq_viewer = tk.Label(self.viewframe,text='Sequence Viewer',font=('SystemButtonText',14,'bold'))
        seq_viewer.grid(row=9,column=0,columnspan=3,sticky=tk.E)

        seq_viewer_tt = tooltip(seq_viewer,self.tooltips['Sequence Viewer'][0])
        #photo1 = tk.PhotoImage(file=("icons" + backslash +"popout.gif"))
        #photo1 = tk.PhotoImage(file=(os.path.join("icons","popout.gif")))

        photo1 = self.icondict['popout']
        seqpop = tk.Button(self.viewframe,image=self.icondict['popout'],command=self.popout_fig_seq) #button for expanding out the sequence plot
        seqpop.grid(row=9,column=3,columnspan=1,sticky=tk.W)
        seqpop.image = self.icondict['popout']

        seqpop_tt = tooltip(seqpop,self.tooltips['popout'][0])

        seq_frame.grid_rowconfigure(10,weight=1)
        figframe = tk.Frame(self.viewframe)
        figframe.grid(row=10,column=0,columnspan=10,padx=gpx,pady=gpy,sticky=tk.W+tk.N+tk.E+tk.S)


        # Create the figure we desire to add to an existing canvas
        self.protein_fig = mpl.figure.Figure(figsize=(1,1))#figsize=(12,5),dpi=self.global_dpi)
        self.protein_fig.set_tight_layout(True)
        self.ax = self.protein_fig.add_subplot(111)


        self.vis_canvas = FigureCanvasTkAgg(self.protein_fig,master=figframe)  #tkinter specific matplotlib canvas
        self.vis_canvas.draw()
        self.vis_canvas.get_tk_widget().pack(fill='both',side='top',expand=True) #place the figure widget
        self.protein_fig.patch.set_facecolor(self.default_color)  #make the edges of the figure the same as the background so it looks like it blends in

        # the blank protien example using plasmid_H2B
        sample_protein = [146,318,464,[2,11,20,196,206,218,228,300,309,318]]  #plasmid h2b example
        self.plot_sequence(self.ax,sample_protein[0],sample_protein[1],sample_protein[2],sample_protein[3],'Flag','Plasmid_H2B')




        poi = tk.Label(self.viewframe,text='POI = Protein of intrest',font=('SystemLabelText',8))
        poi.grid(row=12,column=0,columnspan=2)



        '''
        Stochastic Simulation Frame
        '''
        '''
        topframe = tk.Frame(ss_frame)
        topframe.grid(row=0,column=0,columnspan=11,sticky=tk.W+tk.E+tk.N+tk.S)
        ss_frame.grid_rowconfigure(0,weight=2)
        ss_frame.grid_columnconfigure(0,weight=2)

        #set up all the top entries and labels
        gs_input_l = tk.Label(topframe,text='Gene Sequence Input',font=('SystemButtonText',11,'bold'),bg='#888888',fg='#FFFFFF')
        gs_input_l.grid(row=0,column=0,columnspan=11,sticky=tk.W+tk.E+tk.N+tk.S)

        gbid = tk.Label(topframe,text='GenBank ID')
        gbid.grid(row=1,column=0,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid = tk.Label(topframe,text='Selected File')
        gbid.grid(row=2,column=0,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)


        gbid_e = tk.Entry(topframe,justify=tk.RIGHT)
        gbid_e.grid(row=1,column=2,sticky=tk.EW,padx=gpx,pady=gpy,columnspan=2)

        sfe = tk.Entry(topframe,justify=tk.RIGHT)
        sfe.grid(row=2,column=2,sticky=tk.EW,padx=gpx,pady=gpy,columnspan=2)


        gbid = tk.Label(topframe,text='Simulation Time')
        gbid.grid(row=1,column=4,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid_tt = tooltip(gbid,self.tooltips['Simulation Time'][0])

        gbid = tk.Label(topframe,text='Repetitions')
        gbid.grid(row=2,column=4,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid2_tt = tooltip(gbid,self.tooltips['Repetitions'][0])


        sste = tk.Entry(topframe,justify=tk.RIGHT)
        sste.grid(row=1,column=6,sticky=tk.EW,padx=gpx,pady=gpy,columnspan=1)

        repe = tk.Entry(topframe,justify=tk.RIGHT)
        repe.grid(row=2,column=6,sticky=tk.EW,padx=gpx,pady=gpy,columnspan=1)

        self.ss_gs_inputs = [gbid_e,sfe,sste,repe]  #gene id, selected file, simulation time, repeittions
        sec = tk.Label(topframe,text='sec')
        sec.grid(row=1,column=7,sticky=tk.W)


        #run the stochastic simulation
        run_ss = tk.Button(topframe,text='Run',command=self.run_ssa,font=('SystemButtonText',14),fg='#FFFFFF',bg='#33BB33')
        run_ss.grid(row=1,column=9,rowspan=3,sticky=tk.W,padx=gpx,pady=gpy)


        gl = tk.Label(topframe,text='Gene Length')  #parameters for ssa run
        gl.grid(row=3,column=0,sticky=tk.E,padx=gpx,pady=gpy)

        gl = tk.Label(topframe,text='Codons')
        gl.grid(row=3,column=2,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(topframe,text='Exp. Elongation Rate')
        gl.grid(row=3,column=3,sticky=tk.E,padx=gpx,pady=gpy)

        gl_tt = tooltip(gl,self.tooltips['Exp. Elongation Rate'][0])


        gl = tk.Label(topframe,text='aa/sec')
        gl.grid(row=3,column=5,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(topframe,text='Exp. Initiation Rate')
        gl.grid(row=3,column=6,sticky=tk.E,padx=gpx,pady=gpy)

        gl2_tt = tooltip(gl,self.tooltips['Exp. Initiation Rate'][0])

        gl = tk.Label(topframe,text='1/sec')
        gl.grid(row=3,column=8,sticky=tk.W,padx=gpx,pady=gpy)

        gl_e = tk.Entry(topframe,state=tk.DISABLED,justify=tk.RIGHT)
        gl_e.grid(row=3,column=1,sticky=tk.EW,padx=gpx,pady=gpy)

        eer_e = tk.Entry(topframe,justify=tk.RIGHT)
        eer_e.grid(row=3,column=4,sticky=tk.EW,padx=gpx,pady=gpy)

        eir_e = tk.Entry(topframe,justify=tk.RIGHT)
        eir_e.grid(row=3,column=7,sticky=tk.EW,padx=gpx,pady=gpy)

        self.ss_rates = [gl_e,eer_e,eir_e]  #gene length, exp elongation, exp init

        topframe.grid_columnconfigure(9,weight=2)
        #topframe.grid_columnconfigure(2,weight=2)

        '''
        
        
        

        kyminfoframe = tk.Frame(kym_frame)
        kyminfoframe.grid(row=0,column=0,sticky=tk.W+tk.E,padx=gpx,pady=gpy)
        kymspacerframe = tk.Frame(kym_frame)
        kymspacerframe.grid(row=0,column=1,sticky=tk.W+tk.E,padx=gpx,pady=gpy)
        kymspacer = tk.Label(kymspacerframe,text = '                                                               ')
        kymspacer.pack()

        kymtopframeL = tk.Frame(kym_frame)
        kymtopframeL.grid(row=1,column=0,columnspan=1,sticky=tk.W+tk.E+tk.N+tk.S,padx=gpx,pady=gpy)     
        
        kymtopframeR = tk.Frame(kym_frame)
        kymtopframeR.grid(row=1,column=1,sticky=tk.W+tk.E+tk.N+tk.S,padx=gpx,pady=gpy)     
        
        kymbottomframe = tk.Frame(kym_frame)
        kymbottomframe.grid(row=2,column=0,columnspan=1,sticky=tk.W+tk.E+tk.N+tk.S,padx=gpx,pady=gpy)

        kym_cmap_list=np.array(cm.datad.keys()).astype(str).tolist()
        
        self.kym_menu = AutocompleteComboBox(kyminfoframe)
        self.kym_menu.set_completion_list(kym_cmap_list)
        
        
        
        self.kym_menu.grid(row=0,column=1,sticky=tk.W,pady=2,padx=2)
        self.kym_menu.config(width=20)        
        
        kymcmaplabel = tk.Label(kyminfoframe, text = 'Intensity color map:',font=('SystemLabelText',global_font_size) )
        kymcmaplabel.grid(row=0,column=0,sticky = tk.W)
        
        
        kym_traj_color = tk.Label(kyminfoframe,text='Line color:',font=('SystemLabelText',global_font_size) )
        kym_traj_color.grid(row=0,column=3,sticky=tk.W)
        
        self.kymlinecolor = tk.Button(kyminfoframe,text='      ',command=self.kymline_color_picker,bg='#550c7c',relief=tk.FLAT)
        self.kymbgcolor = tk.Button(kyminfoframe,text='      ',command=self.kymbg_color_picker,bg='#FFFFFF',relief=tk.FLAT)
        self.kymlinecolor.grid(row=0,column=4,pady=2,padx=2,sticky=tk.W)
        self.kymbgcolor.grid(row=0,column=6,pady=2,padx=2,sticky=tk.W)
        
        kymcmaplabel2 = tk.Label(kyminfoframe,text= 'Background color:',font=('SystemLabelText',global_font_size))
        kymcmaplabel2.grid(row=0,column=5,pady=2,padx=2,sticky=tk.W)
        
        self.kymintense = tk.BooleanVar(value=True)
        
        intensity_on = tk.Checkbutton(kyminfoframe,text='Show Intensity', var =self.kymintense,font=('SystemButtonText',global_font_size) )
        intensity_on.grid(row=0,column=7,padx=2,pady=2,sticky=tk.W)

        self.traj_select = tk.Entry(kyminfoframe,width=10)
        self.traj_select.grid(row=0,column=9,padx=2,pady=2)
        
        self.insert_entry_readable(self.traj_select,'0')
        
        kymlabel3 = tk.Label(kyminfoframe,text='Trajectory: ',font=('SystemLabelText',global_font_size))
        kymlabel3.grid(row=0,column=8,padx=2,pady=2)
        
        plot_kym = tk.Button(kyminfoframe,text='Plot',command=self.kymograph,font=('SystemButtonText',global_font_size))
        plot_kym.grid(row=0,column=10,padx=2,pady=2)
                
        

        
        self.kym_fig = mpl.figure.Figure(figsize=(1,1))#figsize=(2,5),dpi=60)
        self.kym_fig.set_tight_layout(True)
        self.kymax = self.kym_fig.add_subplot(111)
        self.kym_fig.patch.set_facecolor(self.default_color)
        self.kym_fig.tight_layout(h_pad=1.0)

        #self.kymax.set_xlabel('time (sec)')
        #self.kymax.set_ylabel('Intensity (a.u.)')
        self.kymax.set_title("")



        self.kym_canvas = FigureCanvasTkAgg(self.kym_fig,master=kymtopframeL)
        self.kym_canvas.draw()
        self.kym_canvas.get_tk_widget().pack(expand=True,fill='both',side='left') #stickying this way it makes it fill all avaliable space

        
        self.kym_fig2 = mpl.figure.Figure(figsize=(1,1))#figsize=(2,5),dpi=60)
        self.kym_fig2.set_tight_layout(True)
        self.kymax2 = self.kym_fig2.add_subplot(111)
        self.kym_fig2.patch.set_facecolor(self.default_color)
        self.kym_fig2.tight_layout(h_pad=1.0)

        #self.kymax.set_xlabel('time (sec)')
        #self.kymax.set_ylabel('Intensity (a.u.)')
        self.kymax2.set_title("")
        
       # figframe = tk.Frame(kymtopframe,bg='#00ff00')
        #figframe.pack(expand=True,fill='both',side='left')


        self.kym_canvas2 = FigureCanvasTkAgg(self.kym_fig2,master=kymtopframeR)
        self.kym_canvas2.draw()
        self.kym_canvas2.get_tk_widget().pack(expand=True,fill='both',side='left') 
        

        self.kym_fig3 = mpl.figure.Figure(figsize=(5,2))#figsize=(2,5),dpi=60)
        self.kym_fig3.set_tight_layout(True)
        self.kymax3 = self.kym_fig3.add_subplot(111)
        self.kym_fig3.patch.set_facecolor(self.default_color)
        self.kym_fig3.tight_layout(h_pad=1.0)

        #self.kymax.set_xlabel('time (sec)')
        #self.kymax.set_ylabel('Intensity (a.u.)')
        self.kymax3.set_title("")
        
       # figframe = tk.Frame(kymtopframe,bg='#00ff00')
        #figframe.pack(expand=True,fill='both',side='left')


        self.kym_canvas3 = FigureCanvasTkAgg(self.kym_fig3,master=kymbottomframe)
        self.kym_canvas3.draw()
        self.kym_canvas3.get_tk_widget().pack(expand=True,fill='both',side='left') 

        kymsmallframe = tk.Frame(kym_frame)
        kymsmallframe.grid(row=2,column=0,sticky=tk.NW)

        kym_frame.rowconfigure(1,weight=3)
        kym_frame.columnconfigure(0,weight=3)
        kym_frame.columnconfigure(1,weight=4)
        
        
    
        
        
        gs_input_l = tk.Label(ss_frame,text='Stochastic Simulations',font=('SystemButtonText',11,'bold'),bg='#888888',fg='#FFFFFF')
        gs_input_l.grid(row=4,column=0,columnspan=11,sticky=tk.W+tk.E+tk.N+tk.S)


        tccframe = tk.Frame(ss_frame)
        tccframe.grid(row=5,column=0,columnspan=1,sticky=tk.EW,rowspan=1)
        timecourse = tk.Label(tccframe,text='Time Courses',font=('SystemButtonText',11,'bold'))
        timecourse.grid(row=0,column=0,columnspan=1,sticky=tk.NW,rowspan=1)

        timecourse_tt = tooltip(timecourse,self.tooltips['Time Courses'][0])

        seqpop = tk.Button(tccframe,image=photo1,command=self.popout_fig_ssa,font=('SystemButtonText',global_font_size)) #button for expanding out the sequence plot
        seqpop.grid(row=0,column=1,columnspan=1,sticky=tk.W)
        seqpop.image = photo1 #keep photo alive

        self.tcplottype = tk.StringVar(value='All')
        tt_dropdown = tk.OptionMenu(ss_frame,self.tcplottype,"All","Average Time Course","PDF")
        tt_dropdown.config(font=('SystemButtonText',global_font_size))
        tt_dropdown['menu'].config(font=('SystemButtonText',global_font_size))

        tt_dropdown.grid(row=5,column=4,sticky=tk.EW,padx=3)

        tt_dropdown.config(width=20)
        
        pdf_time_label = tk.Label(ss_frame,text='PDF time slice',font=('SystemLabelText',global_font_size))
        pdf_time_label.grid(row=5,column=2,sticky=tk.EW,padx=3)
        self.pdf_time_entry = tk.Entry(ss_frame,width=5)
        self.pdf_time_entry.grid(row=5,column=3,sticky=tk.EW,padx=3)

        self.tcplottype.trace('w',lambda name, index, mode,tcplottype= self.tcplottype: self.update_timecourse_plot(self.tcplottype))




        simbuttonframe = tk.Frame(sim_frame)
        simbuttonframe.grid(row=7,column=0,columnspan=4,sticky=tk.EW,pady=2,padx=gpx)

        # set up the updating animation canvas and buttons / dialog
        simulate = tk.Button(simbuttonframe,text = 'Simulate',font=('SystemButtonText',9,'bold'),command=self.simulate_cell)
        simulate.grid(row=1,column=0,sticky=tk.NW,rowspan =1)

        simulate_tt = tooltip(simulate,self.tooltips['Simulate'][0])

        self.simulate_fname = tk.Entry(simbuttonframe,width=20,justify=tk.RIGHT)       #file name entry
        self.simulate_fname.grid(row=1,column=3,sticky=tk.NW,pady=2,padx=gpx)


        files = os.listdir(os.getcwd())
        firstfile=False
        for f in files:
            if firstfile == False:
                if '.gif' in f:
                    im = PIL.Image.open(f)
                    from itertools import count

                    try:
                        for i in range(1):
                            testframe = PIL.ImageTk.PhotoImage(im.copy())
                            if testframe.width() == 450 and testframe.height() == 450:
                                self.simulate_fname.insert(0, f[:-4])
                                firstfile = True
                            im.seek(i)




                    except EOFError:
                        pass


        simulate_fname_tt = tooltip(self.simulate_fname,self.tooltips['simfile'][0])


        self.simscale = tk.Scale(simbuttonframe,orient=tk.HORIZONTAL,showvalue=0)
        self.simscale.grid(row=0,column=0,columnspan=7,sticky=tk.EW)
        filetypes = ['.gif','.gif','.html']           #current filetypes for dropdown
        self.filetype = tk.StringVar(self,value='.gif')
        giflabel = ttk.OptionMenu(simbuttonframe,self.filetype,*filetypes)    #make the dropdown
        giflabel.grid(row=1,column=4)
        
        photo2 = self.icondict['play']
        #photo2 = tk.PhotoImage(file=os.path.join("icons","play.gif"))    #load the play button icon

        #delaylabel = tk.Label(simbuttonframe,text='delay speed (ms)' )
        #delaylabel.grid(row=1,column=2)
        self.delayentry = tk.Entry(simbuttonframe,width=7,justify=tk.RIGHT)
        self.delayentry.grid(row=1,column=2,sticky=tk.NW,pady=2,padx=gpx)




        play = tk.Button(simbuttonframe,image=photo2,command=self.play_sim,font=('SystemButtonText',global_font_size))  #play button
        play.grid(row=1,column=5,sticky=tk.NW,pady=2,padx=gpx)
        play.image=photo2  #be sure to keep the object handle alive so its not garbage collected
        
        photo3 = self.icondict['pause']
        #photo2 = tk.PhotoImage(file=os.path.join("icons","pause.gif"))
        pause = tk.Button(simbuttonframe,image=photo3,command=self.pause_sim,font=('SystemButtonText',global_font_size))
        pause.grid(row=1,column=6,sticky=tk.NW,pady=2,padx=gpx)
        pause.image=photo3  #be sure to keep the object handle alive so its not garbage collected
        try:
            fr = tk.PhotoImage(file=os.path.join("icons",'blank_frame.gif'))   #load a test gif for the sim file
        except:
            pass

        im2 = PIL.Image.open(os.path.join("icons",'blank_frame.gif'))
        animationframe = tk.Frame(sim_frame)
        animationframe.grid(row=5,column=0,columnspan=4,rowspan=2,sticky=tk.W+tk.E+tk.N+tk.S)

        #animation canvas
        w,h = im2.size
        self.animation = tk.Canvas(animationframe,width=w-40,height=h-40,bg='#FFFFFF')
        self.animation.pack(fill=tk.BOTH)
        try:
            self.animation.create_image(-40,-20,anchor=tk.NW,image=fr) #-40 to remove left side deadspace
            self.animation.update()

            self.animation.image = fr
        except:
            pass

        simulate_cell_tt = tooltip(self.animation,self.tooltips['Simulated Cell'][0])

        figframe = tk.Frame(ss_frame)
        figframe.grid(row=6,column=0,columnspan=7,rowspan=2,sticky=tk.W+tk.E+tk.N+tk.S,padx=gpx,pady=gpy)



        dialogframe = tk.Frame(sim_frame)
        dialogframe.grid(row=5,column=4)




        im = PIL.Image.open(os.path.join("icons",("simulation_flowchart" + ".png")))
        photo = PIL.ImageTk.PhotoImage(im)
        piclabel = tk.Label(dialogframe,image=photo)
        piclabel.grid(row=0,column=0,columnspan=6)
        piclabel.image = photo


        dlabel = tk.Label(dialogframe,text='DC (um^2/s)')
        krnal = tk.Label(dialogframe,text='kRNA')
        konl = tk.Label(dialogframe,text='kOn_RNA')
        koffl = tk.Label(dialogframe,text='kOff_RNA')
        kdecayl = tk.Label(dialogframe,text='kDecay_RNA')

        dlabel.grid(row=1,column=0)
        krnal.grid(row=2,column=0)
        konl.grid(row=3,column=0)
        koffl.grid(row=4,column=0)
        kdecayl.grid(row=5,column=0)

        self.krna = tk.StringVar(value='')
        self.krna.trace('w',lambda name, index, mode,krna= self.krna: self.update_rna_pred(self.krna))

        de = tk.Entry(dialogframe,width=15,justify=tk.RIGHT)
        krnae = tk.Entry(dialogframe,width=15,justify=tk.RIGHT,textvariable=self.krna)
        kone= tk.Entry(dialogframe,width=15,justify=tk.RIGHT)
        koffe= tk.Entry(dialogframe,width=15,justify=tk.RIGHT)
        kdecaye= tk.Entry(dialogframe,width=15,justify=tk.RIGHT)

        de.grid(row=1,column=1,pady=2,padx=2)
        krnae.grid(row=2,column=1,pady=2,padx=2)
        kone.grid(row=3,column=1,pady=2,padx=2)
        koffe.grid(row=4,column=1,pady=2,padx=2)
        kdecaye.grid(row=5,column=1,pady=2,padx=2)

        de.insert(tk.END,'1')
        krnae.insert(tk.END,'0.08')
        kone.insert(tk.END,'3')
        koffe.insert(tk.END,'1')
        kdecaye.insert(tk.END,'.005')




        sep = ttk.Separator(dialogframe,orient='vertical')
        sep.grid(row=1,column=2,sticky=tk.NS,rowspan=5)

        frame2 = tk.Frame(dialogframe,bg='#000000')
        frame2.grid(row=1,column=3,rowspan=3,sticky=tk.E + tk.W + tk.S + tk.N,columnspan=2)
        dialogframe.grid_columnconfigure(3,weight=2)
        rnacl = tk.Label(frame2,text='RNA color',fg='#FFFFFF',bg='#000000')
        fluorcl = tk.Label(frame2,text='Tag color',fg='#FFFFFF',bg='#000000')

        self.rnacb = tk.Button(frame2,text='      ',command=self.rna_color_picker,bg='#ff0000',relief=tk.FLAT,font=('SystemButtonText',global_font_size))
        self.fluorcb = tk.Button(frame2,text='      ',command=self.fluor_color_picker,bg='#00ff00',relief=tk.FLAT,font=('SystemButtonText',global_font_size))
        self.rnacb.grid(row=1,column=4,pady=2,padx=2,sticky=tk.E)
        self.fluorcb.grid(row=2,column=4,pady=2,padx=2,sticky=tk.E)
        rnacl.grid(row=1,column=3,sticky=tk.W,pady=2,padx=2)
        fluorcl.grid(row=2,column=3,sticky=tk.W,pady=2,padx=2)

        csl = tk.Label(dialogframe,text='Cell size')
        csl.grid(column=3,row=5,pady=2,padx=2)
        cse = tk.Entry(dialogframe,width=8,justify=tk.RIGHT)
        cse.grid(column=4,row=5,pady=2,padx=2)

        cse.insert(tk.END,'50')

        sep = ttk.Separator(dialogframe,orient='vertical')
        sep.grid(row=6,column=0,sticky=tk.EW,columnspan=5)

        fnamel = tk.Label(dialogframe,text='Filename')
        fnamel.grid(row=7,column=0,columnspan=1,padx=2)
        fnamee = tk.Entry(dialogframe,justify=tk.RIGHT)
        fnamee.grid(row=7,column=1,columnspan=4,sticky=tk.E+tk.N+tk.W+tk.S,pady=2,padx=2)

        fnamee.insert(tk.END,'simulated_cell')

        self.sim_entries = [de,krnae,kone,koffe,kdecaye,cse,fnamee]


        dc_tt = tooltip(dlabel,self.tooltips['DC (um^2/s)'][0])
        krna_tt = tooltip(krnal,self.tooltips['kRNA'][0])
        kon_tt = tooltip(konl,self.tooltips['kOn_RNA'][0])
        koff_tt = tooltip(koffl,self.tooltips['kOff_RNA'][0])
        kd_tt = tooltip(kdecayl,self.tooltips['kDecay_RNA'][0])
        cs_tt = tooltip(csl,self.tooltips['Cell size'][0])


        dce_tt = tooltip(de,self.tooltips['DC (um^2/s)'][0])
        krnae_tt = tooltip(krnae,self.tooltips['kRNA'][0])
        kone_tt = tooltip(kone,self.tooltips['kOn_RNA'][0])
        koffe_tt = tooltip(koffe,self.tooltips['kOff_RNA'][0])
        kde_tt = tooltip(kdecaye,self.tooltips['kDecay_RNA'][0])
        cse_tt = tooltip(cse,self.tooltips['Cell size'][0])

        self.expected_rna = tk.Label(dialogframe,text='Expected RNA: ')
        self.expected_rna.grid(row=8,column=0,columnspan=3,sticky=tk.E)


        runsim = tk.Button(dialogframe,text='RUN',command=self.gui_sim_cell,font=('SystemButtonText',global_font_size))
        runsim.grid(row=8,column=3,columnspan=2,sticky=tk.E)

        try:


            t = int(self.ss_gs_inputs[2].get())
            e = int(float(self.krna.get()) * t)
            std = 2*int(np.sqrt(e))


            if e > 9999 or std > 9999:

                self.expected_rna.config(text=('Expected RNA: ' + "{:.2e}".format(e) + ' +- '+"{:.2e}".format(std)  ))
            else:
                self.expected_rna.config(text=('Expected RNA: ' + str(e) + ' +- '+str(std)   ))

            self.update()
        except:
            pass

        try:
            self.ssa
        except:
            #tmb.showerror('Error','No SSA trajectory run, please run a simulation')
            #return
            pass











        # The ssa trajectory figure
        self.time_course_fig = mpl.figure.Figure(figsize=(1,1))#figsize=(2,5),dpi=60)
        self.time_course_fig.set_tight_layout(True)
        self.tcax = self.time_course_fig.add_subplot(111)
        self.time_course_fig.patch.set_facecolor(self.default_color)
        self.time_course_fig.tight_layout(h_pad=1.0)

        self.tcax.set_xlabel('time (sec)')
        self.tcax.set_ylabel('Intensity (a.u.)')
        self.tcax.set_title("")



        self.tc_canvas = FigureCanvasTkAgg(self.time_course_fig,master=figframe)
        self.tc_canvas.draw()
        self.tc_canvas.get_tk_widget().pack(expand=True,fill='both',side='top') #stickying this way it makes it fill all avaliable space

        ss_frame.columnconfigure(0,weight=3)
        ss_frame.rowconfigure(6,weight=3)




        gs_input_l = tk.Label(fcs_frame,text='Simulated Fluorescence Correlation Spectroscopy',font=('SystemButtonText',11,'bold'),bg='#888888',fg='#FFFFFF')
        gs_input_l.pack(expand=False,fill='x',side='top',anchor=tk.N)   ##.grid(row=8,column=0,columnspan=11,sticky=tk.W+tk.E+tk.N+tk.S)


        #ss_frame.grid_rowconfigure(6,weight=2)





        acframe = tk.Frame(fcs_frame)

        aclabel = tk.Label(acframe,text='Autocorrelation',font=('SystemButtonText',11,'bold'))
        acframe.pack(expand=False,fill='x',side='top',anchor=tk.N) #(row=9,column=0,columnspan=3,sticky=tk.NW,padx=10)
        aclabel.grid(row=0,column=0)

        '''
        self.CAIlabel = tk.Entry(ss_frame,width=5)
        self.CAIlabel.grid(row=9,column=4,sticky=tk.W,pady=1)
        self.CAIlabel.config(state='disabled')

        cailab = tk.Label(ss_frame,text='CAI')
        cailab.grid(row=9,column=3,sticky=tk.E,pady=1)
        '''

        aclabel_tt = tooltip(aclabel,self.tooltips['Autocorrelation'][0])

        photo1 = self.icondict['popout']

        seqpop = tk.Button(acframe,image=photo1,command=self.popout_fig_acc) #button for expanding out the acc plot
        seqpop.grid(row=0,column=1,columnspan=1,sticky=tk.W)
        seqpop.image = photo1

        figframe = tk.Frame(fcs_frame)
        figframe.pack(expand=True,fill='both',side='top',anchor=tk.N)   #.grid(row=10,column=0,columnspan=5,padx=gpx,pady=gpy,sticky=tk.W+tk.N+tk.E+tk.S)
        #ss_frame.grid_rowconfigure(10,weight=2)
        self.autocor_fig = mpl.figure.Figure(figsize=(4,3))  #autocorrelation figure
        self.acax = self.autocor_fig.add_subplot(111)
        self.autocor_fig.patch.set_facecolor(self.default_color)
        self.autocor_fig.set_tight_layout(True)


        self.ac_canvas = FigureCanvasTkAgg(self.autocor_fig,master=figframe)
        self.ac_canvas.draw()
        self.ac_canvas.get_tk_widget().pack(expand=True,side="top",fill="both")


        #datatable.grid_columnconfigure(1,weight=2)





        #gs_input_l = tk.Label(ss_frame,text='Gene Sequence Analysis',font=('SystemButtonText',11,'bold'),bg='#888888',fg='#FFFFFF')
        #gs_input_l.grid(row=11,column=0,columnspan=11,sticky=tk.W+tk.E+tk.N+tk.S)
        rb_frame = tk.Frame(fcs_frame)
        rb =  tk.Label(rb_frame,text='Ribosome Density',font=('SystemButtonText',9,'bold'))
        rb_frame.pack(expand=False,fill='x',side='top',anchor=tk.NW) #.grid(row=11,column=0,columnspan=2,sticky=tk.EW,padx=10)
        rb.grid(row=0,column=0,columnspan=1,sticky=tk.E)

        seqpop = tk.Button(rb_frame,image=photo1,command=self.popout_fig_caicodon) #button for expanding out the plot
        seqpop.grid(row=0,column=1,sticky=tk.W)
        seqpop.image = photo1

        rb_tt = tooltip(rb,self.tooltips['Ribosome Density'][0])

        figframe = tk.Frame(fcs_frame)
        figframe.pack(expand=True,fill='both',side='top',anchor=tk.NW)   #.grid(row=12,column=0,columnspan=5,padx=0,pady=0,sticky=tk.W+tk.N+tk.E+tk.S,rowspan=4)

        self.rib_fig = mpl.figure.Figure(figsize=(4,2))#figsize=(3,2),dpi = 60)    #ribosome figures
        self.rbax = self.rib_fig.add_subplot('111')

        self.rb_canvas = FigureCanvasTkAgg(self.rib_fig,master=figframe)
        self.rb_canvas.draw()
        self.rb_canvas.get_tk_widget().pack(side='top',fill='both',expand=True)
        self.rib_fig.patch.set_facecolor(self.default_color)
        self.rib_fig.set_tight_layout(True)


        rb_frame2 = tk.Frame(fcs_frame)
        rb =  tk.Label(rb_frame2,text='Codon Usage',font=('SystemButtonText',9,'bold'))
        rb_frame2.pack(expand=False,fill='x',side='top',anchor=tk.NW)   #(row=11,column=6,columnspan=1,sticky=tk.EW)
        rb.grid(row=0,column=0)


        seqpop = tk.Button(rb_frame2,image=photo1,command=self.popout_fig_ribdense) #button for expanding out the plot
        seqpop.grid(row=0,column=1,columnspan=1,sticky=tk.W,rowspan=1)
        seqpop.image = photo1



        rb2_tt = tooltip(rb,self.tooltips['Codon Usage'][0])



        figframe = tk.Frame(fcs_frame)
        figframe.pack(expand=True,fill='both',side='top')   #.grid(row=12,column=5,columnspan=5,padx=0,pady=0,sticky=tk.W+tk.N+tk.E+tk.S,rowspan=4)

        #CAI plot
        self.codonuse_fig = mpl.figure.Figure(figsize=(4,2))#figsize=(2,2),dpi=60)
        self.cuax = self.codonuse_fig.add_subplot('111',xmargin=.05,ymargin=.05)

        self.cuax.title.set_visible(False)


        self.codonuse_fig.set_facecolor(self.default_color)
        self.codonuse_fig.set_tight_layout(True)



        self.cu_canvas = FigureCanvasTkAgg(self.codonuse_fig,master=figframe)
        self.cu_canvas.draw()
        self.cu_canvas.get_tk_widget().pack(side='top',fill='both',expand=True)


        upperdata_frame = tk.Frame(fcs_frame)
        upperdata_frame.pack(expand=False,fill='x',side='left',padx=2,pady=2)

        leftdata_frame = tk.Frame(fcs_frame)
        leftdata_frame.pack(expand=False,fill='x',side='left',padx=2,pady=2)

        datatable = tk.Frame(upperdata_frame)  #generate the parameter table (tk.entries)
        datatable.pack(expand=False,fill='both',side='left')    #.grid(row=9,column=5,columnspan=8,sticky=tk.E,padx=gpx,pady=gpy,rowspan=2)

        datatable2 = tk.Frame(leftdata_frame)  #generate the parameter table (tk.entries)
        datatable2.pack(expand=False,fill='both',side='left')    #.grid(row=9,column=5,columnspan=8,sticky=tk.E,padx=gpx,pady=gpy,rowspan=2)

        #ss_frame.grid_columnconfigure(5,weight=2)


        height = 5
        width = 3
        lsa = np.array([['Codon Adaptation Index','','--'],
                        #['Gene Length','','codons'],
                        ['Elongation rate','','aa/sec'],
                        ['Dist. btw Rib.','','codons'],
                        ['Elongation ','','sec'],
                        ['# Rib. per mRNA','','Rib/mRNA']],dtype=object)

        self.calculated_parameters = []
        for i in range(height): #Rows
            for j in range(width): #Columns
                if i%2 ==0:

                    b = tk.Entry(datatable, text="",bd=1,bg='#DDDDDD',relief=tk.RIDGE,font=('SystemButtonText',10))
                else:
                    b = tk.Entry(datatable, text="",bd=1,relief= tk.RIDGE,font=('SystemButtonText',10),readonlybackground='#FFFFFF')
                if i == 0:
                    b.grid(row=i, column=j,ipady=1,sticky=tk.EW)

                else:
                    b.grid(row=i,column=j,ipady=1,sticky=tk.EW)

                if j == 1:
                    b.config(justify=tk.RIGHT)

                b.insert(tk.END,lsa[i,j])



                b.config(state='readonly')
                self.calculated_parameters.append(b)





        height = 5
        width = 3
        lsa = np.array([['True Elongation rate','','aa/sec'],
                        #['Gene Length','','codons'],
                        ['Run off time (Harr)','','sec'],
                        ['Recovery Time (FRAP)','','sec'],
                        ['Perturbed Elong. Rate','','aa/sec'],
                        ['---','','---']],dtype=object)


        for i in range(height): #Rows
            for j in range(width): #Columns
                if i%2 ==0:

                    b = tk.Entry(datatable2, text="",bd=1,bg='#DDDDDD',relief=tk.RIDGE,font=('SystemButtonText',10))
                else:
                    b = tk.Entry(datatable2, text="",bd=1,relief= tk.RIDGE,font=('SystemButtonText',10),readonlybackground='#FFFFFF')
                if i == 0:
                    b.grid(row=i, column=j,ipady=1,sticky=tk.EW)

                else:
                    b.grid(row=i,column=j,ipady=1,sticky=tk.EW)

                if j == 1:
                    b.config(justify=tk.RIGHT)

                b.insert(tk.END,lsa[i,j])



                b.config(state='readonly')
                self.calculated_parameters.append(b)



        dataframe = tk.Frame(leftdata_frame)
        dataframe.pack(expand=False,fill='both',side='left')

        #save_data = tk.Button(dataframe,text='Save Data',fg='#FFFFFF',bg='#33BB33',command=self.save_ss_data)
        #save_data.pack(side='top',fill='both',pady=2,padx=2)

        #load_data = tk.Button(dataframe,text='Load Data',fg='#FFFFFF',bg='#33BB33',command=self.load_ss_data)
        #load_data.pack(side='top',fill='both',pady=2,padx=2)






        '''
        elongation frame
        '''
        '''
        topframe = tk.Frame(elon_frame)
        topframe.grid(row=0,column=0,columnspan=10,sticky=tk.W+tk.E+tk.N+tk.S)
        gs_input_l = tk.Label(topframe,text='Gene Sequence Input',font=('SystemButtonText',13,'bold'),bg='#888888',fg='#FFFFFF')
        gs_input_l.grid(row=0,column=0,columnspan=10,sticky=tk.W+tk.E+tk.N+tk.S)

        gbid = tk.Label(topframe,text='GenBank ID')
        gbid.grid(row=1,column=0,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid = tk.Label(topframe,text='Selected File')
        gbid.grid(row=2,column=0,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)


        gbid_e = tk.Entry(topframe,justify=tk.RIGHT)
        gbid_e.grid(row=1,column=2,sticky=tk.W,padx=gpx,pady=gpy,columnspan=2)

        sfe = tk.Entry(topframe,justify=tk.RIGHT)
        sfe.grid(row=2,column=2,sticky=tk.W,padx=gpx,pady=gpy,columnspan=2)


        gbid = tk.Label(topframe,text='Simulation Time')
        gbid.grid(row=1,column=4,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid5_tt = tooltip(gbid,self.tooltips['Simulation Time'][0])

        gbid = tk.Label(topframe,text='Repetitions')
        gbid.grid(row=3,column=4,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid4_tt = tooltip(gbid,self.tooltips['Repetitions'][0])

        gbid = tk.Label(topframe,text='Perturbation Time')
        gbid.grid(row=2,column=4,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid3_tt = tooltip(gbid,self.tooltips['Perturbation Time'][0])


        pbte = tk.Entry(topframe,justify=tk.RIGHT)
        pbte.grid(row=2,column=6,sticky=tk.W,padx=gpx,pady=gpy,columnspan=1)

        sste = tk.Entry(topframe,justify=tk.RIGHT)
        sste.grid(row=1,column=6,sticky=tk.W,padx=gpx,pady=gpy,columnspan=1)

        repe = tk.Entry(topframe,justify=tk.RIGHT)
        repe.grid(row=3,column=6,sticky=tk.W,padx=gpx,pady=gpy,columnspan=2)

        self.elon_gs_inputs = [gbid_e,sfe,sste,pbte,repe]
        sec = tk.Label(topframe,text='sec')
        sec.grid(row=1,column=7,sticky=tk.W)
        sec = tk.Label(topframe,text='sec')
        sec.grid(row=2,column=7,sticky=tk.W)

        run_ss = tk.Button(topframe,text='Run',command=self.run_ssa_frap,font=('SystemButtonText',14),fg='#FFFFFF',bg='#33BB33')
        run_ss.grid(row=1,column=9,rowspan=3,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(topframe,text='Gene Length')
        gl.grid(row=4,column=0,sticky=tk.E,padx=gpx,pady=gpy)

        gl = tk.Label(topframe,text='Codons')
        gl.grid(row=4,column=2,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(topframe,text='Exp. Elongation Rate')
        gl.grid(row=4,column=3,sticky=tk.E,padx=gpx,pady=gpy)

        gl3_tt = tooltip(gl,self.tooltips['Exp. Elongation Rate'][0])

        gl = tk.Label(topframe,text='aa/sec')
        gl.grid(row=4,column=5,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(topframe,text='Exp. Initiation Rate')
        gl.grid(row=4,column=6,sticky=tk.E,padx=gpx,pady=gpy)

        gl4_tt = tooltip(gl,self.tooltips['Exp. Initiation Rate'][0])

        gl = tk.Label(topframe,text='1/sec')
        gl.grid(row=4,column=8,sticky=tk.W,padx=gpx,pady=gpy)

        gl_e = tk.Entry(topframe,state=tk.DISABLED,justify=tk.RIGHT)
        gl_e.grid(row=4,column=1,sticky=tk.W,padx=gpx,pady=gpy)

        eer_e = tk.Entry(topframe,justify=tk.RIGHT)
        eer_e.grid(row=4,column=4,sticky=tk.W,padx=gpx,pady=gpy)

        eir_e = tk.Entry(topframe,justify=tk.RIGHT)
        eir_e.grid(row=4,column=7,sticky=tk.W,padx=gpx,pady=gpy)


        topframe.grid_columnconfigure(9,weight=2)


        self.ss_rates_elon = [gl_e,eer_e,eir_e]
        '''

        '''
        gs_input_l = tk.Label(elon_frame,text='Simulations with Perturbations',font=('SystemButtonText',13,'bold'),bg='#888888',fg='#FFFFFF')
        gs_input_l.grid(row=5,column=0,columnspan=10,sticky=tk.W+tk.E+tk.N+tk.S)

        perframe = tk.Frame(elon_frame)
        perframe.grid(row=6,column=0,columnspan=10,sticky=tk.W+tk.E+tk.N+tk.S)

        elon_frame.grid_rowconfigure(6,weight=1)
        elon_frame.grid_columnconfigure(0,weight=1)
        perframe.grid_columnconfigure(0,weight=1)
        perframe.grid_columnconfigure(2,weight=1)


        perframe.grid_rowconfigure(3,weight=2)

        self.harring_fig = mpl.figure.Figure(figsize=(1,1))#figsize=(1,7))
        self.hrax = self.harring_fig.add_subplot('111')
        self.harring_fig.set_facecolor(self.default_color)
        self.harring_fig.set_tight_layout(True)

        # Keep this handle alive, or else figure will disappear

        self.hr_canvas = FigureCanvasTkAgg(self.harring_fig,master=perframe)
        self.hr_canvas.draw()
        self.hr_canvas.get_tk_widget().grid(row=3,column=0,columnspan=2,sticky=tk.W+tk.E+tk.N+tk.S,padx=gpx,pady=gpy)


        self.frap_fig = mpl.figure.Figure(figsize=(1,1))
        self.frax = self.frap_fig.add_subplot('111')
        self.frap_fig.set_facecolor(self.default_color)
        self.frap_fig.set_tight_layout(True)


        # Keep this handle alive, or else figure will disappear

        self.fr_canvas = FigureCanvasTkAgg(self.frap_fig,master=perframe)
        self.fr_canvas.draw()
        self.fr_canvas.get_tk_widget().grid(row=3,column=2,columnspan=2,sticky=tk.W+tk.E+tk.N+tk.S,padx=gpx,pady=gpy)



        harframe = tk.Frame(perframe)
        harframe.grid(row=0,column=0,columnspan=2,sticky=tk.W+tk.E+tk.N+tk.S)

        frframe = tk.Frame(perframe)
        frframe.grid(row=0,column=2,columnspan=2,sticky=tk.W+tk.E+tk.N+tk.S)

        frlabel = tk.Label(harframe,text='                                        Harringtonine',font=('SystemButtonText',12,'bold'))
        hrlabel = tk.Label(frframe,text='                                                      FRAP',font=('SystemButtonText',12,'bold'))



        hrlabel.pack(side=tk.LEFT,anchor=tk.E)
        frlabel.pack(side=tk.LEFT,anchor=tk.E)


        fr_tt = tooltip(hrlabel,self.tooltips['FRAP'][0])
        hr_tt = tooltip(frlabel,self.tooltips['Harringtonine'][0])


        seqpop = tk.Button(harframe,image=photo1,command=self.popout_fig_harr) #button for expanding out the sequence plot
        seqpop.pack(side=tk.LEFT,anchor=tk.W)
        seqpop.image = photo1


        seqpop = tk.Button(frframe,image=photo1,command=self.popout_fig_frap) #button for expanding out the sequence plot
        seqpop.pack(side=tk.LEFT,anchor=tk.W)
        seqpop.image = photo1


        calcpframe = tk.Frame(perframe,bg='#FFFFFF')
        calcpframe.grid(row=4,column=0,columnspan=1,padx=gpx,pady=gpy)

        clabel = tk.Label(calcpframe,text='Dynamic Property Values',font=('SystemButtonText',9,'bold'),bg='#FFFFFF')
        clabel.grid(row=0,column=0,columnspan=3,pady=8)

        clabel2_tt = tooltip(clabel,self.tooltips['Dynamic Property Values'][0])

        er = tk.Label(calcpframe,text='Elongation Rate',bg='#FFFFFF')
        er.grid(row=1,column=0,pady=8)
        rt = tk.Label(calcpframe,text='Run-off time',bg='#FFFFFF')
        rt.grid(row=2,column=0,pady=8)

        ere = tk.Entry(calcpframe,width=20)
        ere.grid(row=1,column=1,pady=8)
        rte = tk.Entry(calcpframe,width=20)
        rte.grid(row=2,column=1,pady=8)

        er = tk.Label(calcpframe,text='aa/sec',bg='#FFFFFF')
        er.grid(row=1,column=2,pady=8,sticky=tk.W)
        rt = tk.Label(calcpframe,text='sec',bg='#FFFFFF')
        rt.grid(row=2,column=2,pady=8,sticky=tk.W)


        saveframe = tk.Frame(perframe)
        saveframe.grid(row=4,column=1,padx=gpx,pady=gpy)

        save_data = tk.Button(saveframe,text='Save H Data',fg='#FFFFFF',bg='#33BB33',command=self.save_harr_data)
        save_data.grid(row=0,column=0,columnspan=1,pady=gpy,padx=gpx,sticky=tk.W+tk.S+tk.N)

        save_data = tk.Button(saveframe,text='Load H Data',fg='#FFFFFF',bg='#33BB33',command=self.load_harr_data)
        save_data.grid(row=1,column=0,columnspan=1,pady=gpy,padx=gpx,sticky=tk.W+tk.S+tk.N)




        calcpframe = tk.Frame(perframe,bg='#FFFFFF')
        calcpframe.grid(row=4,column=2,columnspan=1,padx=gpx,pady=gpy)

        clabel = tk.Label(calcpframe,text='Dynamic Property Values',font=('SystemButtonText',9,'bold'),bg='#FFFFFF')
        clabel.grid(row=0,column=0,columnspan=3,pady=8)

        clabel_tt = tooltip(clabel,self.tooltips['Dynamic Property Values'][0])

        er = tk.Label(calcpframe,text='Elongation Rate',bg='#FFFFFF')
        er.grid(row=1,column=0,pady=8)
        rt = tk.Label(calcpframe,text='Recovery time',bg='#FFFFFF')
        rt.grid(row=2,column=0,pady=8)

        ere2 = tk.Entry(calcpframe,width=20)
        ere2.grid(row=1,column=1,pady=8)
        rte2 = tk.Entry(calcpframe,width=20)
        rte2.grid(row=2,column=1,pady=8)

        er = tk.Label(calcpframe,text='aa/sec',bg='#FFFFFF')
        er.grid(row=1,column=2,pady=8,sticky=tk.W)
        rt = tk.Label(calcpframe,text='sec',bg='#FFFFFF')
        rt.grid(row=2,column=2,pady=8,sticky=tk.W)


        saveframe = tk.Frame(perframe)
        saveframe.grid(row=4,column=3,padx=gpx,pady=gpy)

        save_data = tk.Button(saveframe,text='Save F Data',fg='#FFFFFF',bg='#33BB33',command=self.save_frap_data)
        save_data.grid(row=0,column=0,columnspan=1,pady=gpy,padx=gpx,sticky=tk.W+tk.S+tk.N)


        save_data = tk.Button(saveframe,text='Load F Data',fg='#FFFFFF',bg='#33BB33',command=self.load_frap_data)
        save_data.grid(row=1,column=0,columnspan=1,pady=gpy,padx=gpx,sticky=tk.W+tk.S+tk.N)


        self.cpv_entries = [ere,rte,ere2,rte2]
        '''

        '''
        Codon Frame
        '''
        '''
        gs_input_l = tk.Label(codon_frame,text='Gene Sequence Input',font=('SystemButtonText',13,'bold'),bg='#888888',fg='#FFFFFF')
        gs_input_l.grid(row=0,column=0,columnspan=10,sticky=tk.W+tk.E+tk.N+tk.S)

        gbid = tk.Label(codon_frame,text='GenBank ID')
        gbid.grid(row=1,column=0,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid = tk.Label(codon_frame,text='Selected File')
        gbid.grid(row=2,column=0,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)


        gbid_e = tk.Entry(codon_frame,width=20,justify=tk.RIGHT)
        gbid_e.grid(row=1,column=2,sticky=tk.W,padx=gpx,pady=gpy,columnspan=2)

        sfe = tk.Entry(codon_frame,width=20,justify=tk.RIGHT)
        sfe.grid(row=2,column=2,sticky=tk.W,padx=gpx,pady=gpy,columnspan=2)


        gbid = tk.Label(codon_frame,text='Simulation Time')
        gbid.grid(row=1,column=4,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)
        gbidcodon_tt = tooltip(gbid,self.tooltips['Simulation Time'][0])

        gbid = tk.Label(codon_frame,text='Repetitions')
        gbid.grid(row=2,column=4,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbidcodon2_tt = tooltip(gbid,self.tooltips['Repetitions'][0])

        sste = tk.Entry(codon_frame,width=20,justify=tk.RIGHT)
        sste.grid(row=1,column=6,sticky=tk.W,padx=gpx,pady=gpy,columnspan=1)

        repe = tk.Entry(codon_frame,width=20,justify=tk.RIGHT)
        repe.grid(row=2,column=6,sticky=tk.W,padx=gpx,pady=gpy,columnspan=2)

        self.co_gs_inputs = [gbid_e,sfe,sste,repe]
        sec = tk.Label(codon_frame,text='sec')
        sec.grid(row=1,column=7,sticky=tk.W)

        run_ss = tk.Button(codon_frame,text='Run',command=self.run_ssa_codons,font=('SystemButtonText',14),fg='#FFFFFF',bg='#33BB33')
        run_ss.grid(row=1,column=9,rowspan=3,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(codon_frame,text='Gene Length')
        gl.grid(row=3,column=0,sticky=tk.E,padx=gpx,pady=gpy)

        gl = tk.Label(codon_frame,text='Codons')
        gl.grid(row=3,column=2,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(codon_frame,text='Exp. Elongation Rate')
        gl.grid(row=3,column=3,sticky=tk.E,padx=gpx,pady=gpy)

        glcodon_tt = tooltip(gl,self.tooltips['Exp. Elongation Rate'][0])

        gl = tk.Label(codon_frame,text='aa/sec')
        gl.grid(row=3,column=5,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(codon_frame,text='Exp. Initiation Rate')
        gl.grid(row=3,column=6,sticky=tk.E,padx=gpx,pady=gpy)

        glcodon2_tt = tooltip(gl,self.tooltips['Exp. Initiation Rate'][0])

        gl = tk.Label(codon_frame,text='1/sec')
        gl.grid(row=3,column=8,sticky=tk.W,padx=gpx,pady=gpy)

        gl_e = tk.Entry(codon_frame,state=tk.DISABLED,width=15,justify=tk.RIGHT)
        gl_e.grid(row=3,column=1,sticky=tk.W,padx=gpx,pady=gpy)

        eer_e = tk.Entry(codon_frame,width=16,justify=tk.RIGHT)
        eer_e.grid(row=3,column=4,sticky=tk.W,padx=gpx,pady=gpy)

        eir_e = tk.Entry(codon_frame,width=16,justify=tk.RIGHT)
        eir_e.grid(row=3,column=7,sticky=tk.W,padx=gpx,pady=gpy)

        self.ss_rates_codon = [gl_e,eer_e,eir_e]

        '''

        '''
        gs_input_l = tk.Label(codon_frame,text='Codon Optimization Design',font=('SystemButtonText',13,'bold'),bg='#888888',fg='#FFFFFF')
        gs_input_l.grid(row=4,column=0,columnspan=10,sticky=tk.W+tk.E+tk.N+tk.S)

        aa_key = ['Ala (A)','Arg (R)','Asn (N)','Asp (D)','Cys (C)','Gln (Q)','Glu (E)','Gly (G)','His (H)','Ile (I)','Leu (L)','Lys (K)','Met (M)','Phe (F)','Pro (P)','Ser (S)','Thr (T)','Trp (W)','Tyr (Y)','Val (V)']
        self.aa_entries = []
        self.aa_vars = []
        self.sensitivity_entries = []
        aa_options = [('Natural',0),('Common',1),('Rare',2)]

        rbframe = tk.Frame(codon_frame)
        rbframe.grid(row=5,column=0,columnspan=10)
        titles = ['Amino Acid','No.','Wild Type','Common','Rare','Sensitivity']
        for i in range(len(titles)):
            l = tk.Label(rbframe,text=titles[i],font=('SystemButtonText',10,'bold'))
            l.grid(row=1,column=i,padx=2,pady=2)
            if i == 5:
                ss_tt = tooltip(l,self.tooltips['Sensitivity'][0])
        l = tk.Label(rbframe,text='            ')
        l.grid(row=1,column=6,pady=2)


        for i in range(len(titles)):

            l = tk.Label(rbframe,text=titles[i],font=('SystemButtonText',10,'bold'))
            l.grid(row=1,column=i+7,padx=2,pady=2)

            if i == 5:
                ss2_tt = tooltip(l,self.tooltips['Sensitivity'][0])

        k = 0
        for aa in aa_key:
            k+=1
            if k < 11:
                l = tk.Label(rbframe,text=aa,font=('SystemButtonText',10,'bold'))
                l.grid(row=5+k,column=0)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+k,column=1,pady=2)
                self.aa_entries.append(e)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+k,column=5,pady=2)
                self.sensitivity_entries.append(e)


                v = tk.IntVar(value=0)
                j=0
                for text,mode in aa_options:
                    b = tk.Radiobutton(rbframe,text=' ',variable=v,value=mode)
                    b.grid(row=5+k,column=2+j,pady=2)
                    j+=1
                self.aa_vars.append(v)

                spacer = tk.Label(rbframe,text='            ')
                spacer.grid(row=5+k,column=6,pady=2)

            if k >=11:

                l = tk.Label(rbframe,text=aa,font=('SystemButtonText',10,'bold'))
                l.grid(row=5+(k-10),column=7)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+(k-10),column=8,pady=2)
                self.aa_entries.append(e)

                v = tk.IntVar(value=0)
                j = 0
                for text,mode in aa_options:
                    b = tk.Radiobutton(rbframe,text='',variable=v,value=mode)
                    b.grid(row=5+(k-10),column=9+j,pady=2)
                    j+=1
                self.aa_vars.append(v)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+(k-10),column=12,pady=2)
                self.sensitivity_entries.append(e)

        l = tk.Label(rbframe,text='* Sensitivity of the codon manipulation can vary between 0 and 1, where 1 is the most sensitive amino acid for manipulation',font=('SystemButtonText',8,))
        l.grid(row=16,column=0,columnspan=15)

        allnatural = tk.Button(codon_frame,text='All Wild Type Codons',command=self.set_natural_codons)
        allnatural.grid(row=6,column=0,columnspan=3,pady=gpy,padx=gpx)

        allcommon = tk.Button(codon_frame,text='All Common Codons',command=self.set_common_codons)
        allcommon.grid(row=6,column=3,columnspan=4,pady=gpy,padx=gpx)

        allrare = tk.Button(codon_frame,text='All Rare Codons',command=self.set_rare_codons)
        allrare.grid(row=6,column=7,columnspan=4,pady=gpy,padx=gpx)

        '''

        '''
        calcpframe = tk.Frame(ccodon_frame,bg='#FFFFFF')
        calcpframe.grid(row=7,column=6,columnspan=6,padx=gpx,pady=gpy,rowspan=2)

        clabel = tk.Label(calcpframe,text='Dynamic Properties',font=('SystemButtonText',9,'bold'),bg='#FFFFFF')
        clabel.grid(row=0,column=0,columnspan=3,pady=6)

        clabel_tt = tooltip(clabel,self.tooltips['Dynamic Properties'][0])

        er = tk.Label(calcpframe,text='Modified Elongation Rate',bg='#FFFFFF')
        er.grid(row=1,column=0,pady=6)
        rt = tk.Label(calcpframe,text='Wild Type Elongation Rate',bg='#FFFFFF')
        rt.grid(row=2,column=0,pady=6)

        rt = tk.Label(calcpframe,text='Modified Ribosome Distance',bg='#FFFFFF')
        rt.grid(row=3,column=0,pady=6)

        rt = tk.Label(calcpframe,text='Wild Type Ribosome Distance',bg='#FFFFFF')
        rt.grid(row=4,column=0,pady=6)

        ere2 = tk.Entry(calcpframe,width=20)
        ere2.grid(row=1,column=1,pady=6)
        rte2 = tk.Entry(calcpframe,width=20)
        rte2.grid(row=2,column=1,pady=6)

        drd = tk.Entry(calcpframe,width=20)
        drd.grid(row=3,column=1,pady=6)
        nrd = tk.Entry(calcpframe,width=20)
        nrd.grid(row=4,column=1,pady=6)

        er = tk.Label(calcpframe,text='aa/sec',bg='#FFFFFF')
        er.grid(row=1,column=2,pady=6,sticky=tk.W)
        rt = tk.Label(calcpframe,text='aa/sec',bg='#FFFFFF')
        rt.grid(row=2,column=2,pady=6,sticky=tk.W)

        er = tk.Label(calcpframe,text='codons',bg='#FFFFFF')
        er.grid(row=3,column=2,pady=6,sticky=tk.W)
        rt = tk.Label(calcpframe,text='codons',bg='#FFFFFF')
        rt.grid(row=4,column=2,pady=6,sticky=tk.W)

        self.codon_parameter_entries = [ere2,rte2,drd,nrd]

        acc_label = tk.Label(ccodon_frame,text='Autocorrelations',font=('SystemButtonText',12,'bold'))
        acc_label.grid(row=7,column=1,columnspan=2,sticky=tk.SE,padx=gpx,pady=gpy)

        acc_l_tt = tooltip(acc_label,self.tooltips['Autocorrelations_codon'][0])

        seqpop = tk.Button(ccodon_frame,image=photo1,command=self.popout_fig_codoncorr) #button for expanding out the sequence plot
        seqpop.grid(row=7,column=3,columnspan=1,sticky=tk.W)
        seqpop.image = photo1

        figframe = tk.Frame(ccodon_frame)
        figframe.grid(row=8,column=0,columnspan=6,sticky=tk.W+tk.E+tk.N+tk.S,padx=gpx,pady=gpy)

        codon_frame.grid_rowconfigure(8,weight=1)

        self.accodon_fig = mpl.figure.Figure(figsize=(1,1))
        self.accax = self.accodon_fig.add_subplot('111')
        self.accodon_fig.set_facecolor(self.default_color)
        self.accodon_fig.set_tight_layout(True)

        # Keep this handle alive, or else figure will disappear

        self.acc_canvas = FigureCanvasTkAgg(self.accodon_fig,master=figframe)
        self.acc_canvas.draw()
        self.acc_canvas.get_tk_widget().pack(expand=True,side='top',fill='both')

        save_data = tk.Button(ccodon_frame,text='Save Modified SSA Data',fg='#FFFFFF',bg='#33BB33',command=self.save_design_codon_data)
        save_data.grid(row=9,column=0,columnspan=5,pady=gpy,padx=gpx)

        save_data = tk.Button(ccodon_frame,text='Save Wild type SSA Data',fg='#FFFFFF',bg='#33BB33',command=self.save_natural_codon_data)
        save_data.grid(row=9,column=5,columnspan=5,pady=gpy,padx=gpx)


        '''

    def __bind_genbank(self,event):
     
        self.gbid_entry.bind('<Return>',self.poll_genebank)
    def __unbind_genbank(self,event):
        self.gbid_entry.unbind('<Return>')
        

    def initstyle(self):
        ''' Initalize the custom tkinter style (rounded tabs and spacing)
        '''
        try:
            style = ttk.Style()

            style.theme_create( "new", parent="alt", settings={
                    "TProgressbar":{"configure":{"background":['#fcaa0a'],"foreground":['#fcaa0a']}},
                    "TNotebook": {"configure": {"tabmargins": [2, 5, 2, 0], "tabposition":['nw'],"background":['#777777']} } ,
                    "TNotebook.Tab": {
                        "configure": {"padding": [5, 5], "background": '#D3D3D3' },
                        "map":       {"background": [("selected", '#EEEEEE')] } } ,

                    "TCombobox":{"padding":[5,5]},
                    "TCombobox.Entry":{"padding":[5,5]} ,
                    "TLabel":{'font':'helvetica 50',"bg":self.default_color },
                    "TButton.Label":{"fontsize":[45]}})

            return style
        except:
            return


    def fullsize_codon_popout(self):
        try:
            cmnorm = self.sms.POI.codon_sensitivity/np.max(self.sms.POI.codon_sensitivity).flatten().tolist()
        except:
            return

        self.fullsizecodon = tk.Toplevel(self.parent)

        self.fullsizecodon.title('Full Size Codon Editor')
        codon_frame = tk.Frame(self.fullsizecodon)
        codon_frame.pack(fill=tk.BOTH)

        self.fullsizecodon.lift()


        #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))


        self.fullsizecodon.tk.call('wm','iconphoto',self.fullsizecodon._w,self.icondata)

        gpx=2
        gpy=2
        '''
        gs_input_l = tk.Label(codon_frame,text='Gene Sequence Input',font=('SystemButtonText',13,'bold'),bg='#888888',fg='#FFFFFF')
        gs_input_l.grid(row=0,column=0,columnspan=10,sticky=tk.W+tk.E+tk.N+tk.S)

        gbid = tk.Label(codon_frame,text='GenBank ID')
        gbid.grid(row=1,column=0,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        gbid = tk.Label(codon_frame,text='Selected File')
        gbid.grid(row=2,column=0,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)


        gbid_e = tk.Entry(codon_frame,width=20,justify=tk.RIGHT)
        gbid_e.grid(row=1,column=2,sticky=tk.W,padx=gpx,pady=gpy,columnspan=2)

        sfe = tk.Entry(codon_frame,width=20,justify=tk.RIGHT)
        sfe.grid(row=2,column=2,sticky=tk.W,padx=gpx,pady=gpy,columnspan=2)


        gbid = tk.Label(codon_frame,text='Simulation Time')
        gbid.grid(row=1,column=4,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)
        #gbidcodon_tt = tooltip(gbid,self.tooltips['Simulation Time'][0])

        gbid = tk.Label(codon_frame,text='Repetitions')
        gbid.grid(row=2,column=4,sticky=tk.E,padx=gpx,pady=gpy,columnspan=2)

        #gbidcodon2_tt = tooltip(gbid,self.tooltips['Repetitions'][0])

        sste = tk.Entry(codon_frame,width=20,justify=tk.RIGHT)
        sste.grid(row=1,column=6,sticky=tk.W,padx=gpx,pady=gpy,columnspan=1)

        repe = tk.Entry(codon_frame,width=20,justify=tk.RIGHT)
        repe.grid(row=2,column=6,sticky=tk.W,padx=gpx,pady=gpy,columnspan=2)

        #self.co_gs_inputs = [gbid_e,sfe,sste,repe]
        sec = tk.Label(codon_frame,text='sec')
        sec.grid(row=1,column=7,sticky=tk.W)

        run_ss = tk.Button(codon_frame,text='Run',command=self.run_ssa_codons,font=('SystemButtonText',14),fg='#FFFFFF',bg='#33BB33')
        run_ss.grid(row=1,column=9,rowspan=3,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(codon_frame,text='Gene Length')
        gl.grid(row=3,column=0,sticky=tk.E,padx=gpx,pady=gpy)

        gl = tk.Label(codon_frame,text='Codons')
        gl.grid(row=3,column=2,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(codon_frame,text='Exp. Elongation Rate')
        gl.grid(row=3,column=3,sticky=tk.E,padx=gpx,pady=gpy)

        glcodon_tt = tooltip(gl,self.tooltips['Exp. Elongation Rate'][0])

        gl = tk.Label(codon_frame,text='aa/sec')
        gl.grid(row=3,column=5,sticky=tk.W,padx=gpx,pady=gpy)

        gl = tk.Label(codon_frame,text='Exp. Initiation Rate')
        gl.grid(row=3,column=6,sticky=tk.E,padx=gpx,pady=gpy)

        #glcodon2_tt = tooltip(gl,self.tooltips['Exp. Initiation Rate'][0])

        gl = tk.Label(codon_frame,text='1/sec')
        gl.grid(row=3,column=8,sticky=tk.W,padx=gpx,pady=gpy)

        gl_e = tk.Entry(codon_frame,state=tk.DISABLED,width=15,justify=tk.RIGHT)
        gl_e.grid(row=3,column=1,sticky=tk.W,padx=gpx,pady=gpy)

        eer_e = tk.Entry(codon_frame,width=16,justify=tk.RIGHT)
        eer_e.grid(row=3,column=4,sticky=tk.W,padx=gpx,pady=gpy)

        eir_e = tk.Entry(codon_frame,width=16,justify=tk.RIGHT)
        eir_e.grid(row=3,column=7,sticky=tk.W,padx=gpx,pady=gpy)
        '''
        #self.ss_rates_codon = [gl_e,eer_e,eir_e]

        try:
            cmnorm = self.sms.POI.codon_sensitivity/np.max(self.sms.POI.codon_sensitivity).flatten().tolist()
        except:
            return
        #colors = self.__cmap_map(lambda x: .9*x +.2,cm.Wistia)(cmnorm[0])
        colors = cm.coolwarm(cmnorm[0])


        for i in range(len(self.aa_entries)):

            self.aa_entries[i].config(state='normal')
            self.aa_entries[i].delete(0,tk.END)
            self.aa_entries[i].insert(tk.END,str(self.sms.POI.aa_seq.count(self.sms.aa_keys[i])))
            self.aa_entries[i].config(state='readonly')
            self.aa_entries[i].config(readonlybackground='#FFFFFF')


            self.sensitivity_entries[i].config(state='normal')
            self.sensitivity_entries[i].delete(0,tk.END)
            self.sensitivity_entries[i].insert(tk.END,str(self.sms.POI.codon_sensitivity[0,i]))
            self.sensitivity_entries[i].config(state='readonly')
            color = (colors[i][0:3])*255


            dc=[0,0,0]

            dc[0] = hex(int(color[0]))[2:]
            dc[1] = hex(int(color[1]))[2:]
            dc[2] = hex(int(color[2]))[2:]
            if len(dc[0]) < 2:
                dc[0] = '0'+dc[0]
            if len(dc[1]) < 2:
                dc[1] = '0'+dc[1]
            if len(dc[2]) < 2:
                dc[2] = '0'+dc[2]
            color = '#' + dc[0] + dc[1] + dc[2]

            self.sensitivity_entries[i].config(readonlybackground=color)



       # gs_input_l = tk.Label(codon_frame,text='Codon Optimization Design',font=('SystemButtonText',13,'bold'),bg='#888888',fg='#FFFFFF')
        #gs_input_l.grid(row=4,column=0,columnspan=10,sticky=tk.W+tk.E+tk.N+tk.S)

        aa_key = ['Ala (A)','Arg (R)','Asn (N)','Asp (D)','Cys (C)','Gln (Q)','Glu (E)','Gly (G)','His (H)','Ile (I)','Leu (L)','Lys (K)','Met (M)','Phe (F)','Pro (P)','Ser (S)','Thr (T)','Trp (W)','Tyr (Y)','Val (V)']
        #self.aa_entries = []
        #self.aa_vars = []
        #self.sensitivity_entries = []
        aa_options = [('Natural',0),('Common',1),('Rare',2)]

        rbframe = tk.Frame(codon_frame)
        rbframe.grid(row=5,column=0,columnspan=10)
        titles = ['Amino Acid','No.','Wild Type','Common','Rare','Sensitivity']
        for i in range(len(titles)):
            l = tk.Label(rbframe,text=titles[i],font=('SystemButtonText',10,'bold'))
            l.grid(row=1,column=i,padx=2,pady=2)
            if i == 5:
                ss_tt = tooltip(l,self.tooltips['Sensitivity'][0])
        l = tk.Label(rbframe,text='            ')
        l.grid(row=1,column=6,pady=2)


        for i in range(len(titles)):

            l = tk.Label(rbframe,text=titles[i],font=('SystemButtonText',10,'bold'))
            l.grid(row=1,column=i+7,padx=2,pady=2)

            #if i == 5:
                #ss2_tt = tooltip(l,self.tooltips['Sensitivity'][0])

        k = 0
        for aa in aa_key:
            k+=1
            if k < 11:
                l = tk.Label(rbframe,text=aa,font=('SystemButtonText',10,'bold'))
                l.grid(row=5+k,column=0)



                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+k,column=1,pady=2)

                self.insert_entry(e,self.aa_entries[k-1].get())
                e.config(readonlybackground='#FFFFFF')
                #e.insert(tk.END,self.aa_entries[k-1].get())

                #self.aa_entries.append(e)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+k,column=5,pady=2)

                self.insert_entry(e,self.sensitivity_entries[k-1].get())
                color = (colors[k-1][0:3])*255


                dc=[0,0,0]

                dc[0] = hex(int(color[0]))[2:]
                dc[1] = hex(int(color[1]))[2:]
                dc[2] = hex(int(color[2]))[2:]
                if len(dc[0]) < 2:
                    dc[0] = '0'+dc[0]
                if len(dc[1]) < 2:
                    dc[1] = '0'+dc[1]
                if len(dc[2]) < 2:
                    dc[2] = '0'+dc[2]
                color = '#' + dc[0] + dc[1] + dc[2]

                e.config(readonlybackground=color)


                #self.sensitivity_entries.append(e)


                v = self.aa_vars[k-1]
                j=0
                for text,mode in aa_options:
                    b = tk.Radiobutton(rbframe,text=' ',variable=v,value=mode)
                    b.grid(row=5+k,column=2+j,pady=2)
                    j+=1
                #self.aa_vars.append(v)

                spacer = tk.Label(rbframe,text='            ')
                spacer.grid(row=5+k,column=6,pady=2)

            if k >=11:

                l = tk.Label(rbframe,text=aa,font=('SystemButtonText',10,'bold'))
                l.grid(row=5+(k-10),column=7)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+(k-10),column=8,pady=2)
                self.insert_entry(e,self.aa_entries[k-1].get())
                e.config(readonlybackground='#FFFFFF')
                #e.insert(tk.END,self.aa_entries[k-1].get())
                #self.aa_entries.append(e)

                v = self.aa_vars[k-1]
                j = 0
                for text,mode in aa_options:
                    b = tk.Radiobutton(rbframe,text='',variable=v,value=mode)
                    b.grid(row=5+(k-10),column=9+j,pady=2)
                    j+=1
                #self.aa_vars.append(v)

                e = tk.Entry(rbframe,width=8)
                e.grid(row=5+(k-10),column=12,pady=2)
                self.insert_entry(e,self.sensitivity_entries[k-1].get())

                color = (colors[k-1][0:3])*255


                dc=[0,0,0]

                dc[0] = hex(int(color[0]))[2:]
                dc[1] = hex(int(color[1]))[2:]
                dc[2] = hex(int(color[2]))[2:]
                if len(dc[0]) < 2:
                    dc[0] = '0'+dc[0]
                if len(dc[1]) < 2:
                    dc[1] = '0'+dc[1]
                if len(dc[2]) < 2:
                    dc[2] = '0'+dc[2]
                color = '#' + dc[0] + dc[1] + dc[2]

                e.config(readonlybackground=color)
                #self.sensitivity_entries.append(e)
        global_font_size= 8
        l = tk.Label(rbframe,text='* Sensitivity of the codon manipulation can vary between 0 and 1, where 1 is the most sensitive amino acid for manipulation',font=('SystemButtonText',8,))
        l.grid(row=16,column=0,columnspan=15)

        allnatural = tk.Button(codon_frame,text='All Wild Type Codons',command=self.set_natural_codons,font=('SystemButtonText',global_font_size))
        allnatural.grid(row=6,column=0,columnspan=3,pady=gpy,padx=gpx)

        allcommon = tk.Button(codon_frame,text='All Common Codons',command=self.set_common_codons,font=('SystemButtonText',global_font_size))
        allcommon.grid(row=6,column=3,columnspan=4,pady=gpy,padx=gpx)

        allrare = tk.Button(codon_frame,text='All Rare Codons',command=self.set_rare_codons,font=('SystemButtonText',global_font_size))
        allrare.grid(row=6,column=7,columnspan=4,pady=gpy,padx=gpx)





    def view_sequence_newwin(self):

        self.seq_text_open = True
        info = tk.Toplevel(self.parent)

        info.title('Raw Sequence Text')


        info.lift()
        #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
        info.tk.call('wm','iconphoto',info._w,self.icondata)



        self.textframe = tk.Frame(info)
        self.textframe.pack(expand=True,side='top',fill='both')
        self.gs_text = tk.Text(self.textframe,height=20,width=71,bg='#FFFFFF',font=('SystemButtonText',8))
        self.gs_text.pack(expand=True,side='top',fill='both')


        try:
            self.gs_text.delete(0.0,tk.END)
            self.gs_text.insert(tk.END,self.seqtext)
        except:
            pass


        colorize = tk.Button(self.textframe,text='Colorize',command=self.color_text,font=('SystemButtonText',global_font_size))
        colorize.pack()

        colorize_tt = tooltip(colorize,text='Color the sequence display by G, C, T, and A')

    def seqmap(self):

        try:
            self.seqmapwindow.lift()
            self.seq_map_open = True
        except:
            self.seq_map_open = False

        if self.seq_map_open == False:
            self.seq_map_open = True
            info = tk.Toplevel(self.parent)
            self.seqmapwindow = info
            info.title('Sequence Map')

            info.lift()
            #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
            info.tk.call('wm','iconphoto',info._w,self.icondata)
            info.minsize(400, 150)
            self.viewframe = tk.Frame(info)
            self.viewframe.pack(expand=True,side='top',fill='both')


            figframe = tk.Frame(self.viewframe)
            figframe.pack(expand=True,side='top',fill='both')





        # Create the figure we desire to add to an existing canvas
            self.protein_fig = mpl.figure.Figure(figsize=(1,1))#figsize=(12,5),dpi=self.global_dpi)
            self.protein_fig.set_tight_layout(True)
            self.ax = self.protein_fig.add_subplot(111)


            self.vis_canvas = FigureCanvasTkAgg(self.protein_fig,master=figframe)  #tkinter specific matplotlib canvas
            self.vis_canvas.draw()
            self.vis_canvas.get_tk_widget().pack(fill='both',side='top',expand=True) #place the figure widget
            self.protein_fig.patch.set_facecolor(self.default_color)  #make the edges of the figure the same as the background so it looks like it blends in

            # the blank protien example using plasmid_H2B
            sample_protein = [146,318,464,[2,11,20,196,206,218,228,300,309,318]]  #plasmid h2b example
            self.plot_sequence(self.ax,sample_protein[0],sample_protein[1],sample_protein[2],sample_protein[3],'Flag','Plasmid_H2B')




            poi = tk.Label(self.viewframe,text='POI = Protein of intrest',font=('SystemButtonText',8))
            poi.pack()




    def change_trna_copys(self):

        global_font_size=8
        try:
            self.trna_Nb.winfo_children()
            return
        except:
            pass


        info = tk.Toplevel(self.parent)

        info.title('Change tRNA Abundance')
        infoframe = tk.Frame(info)
        infoframe.pack(fill=tk.BOTH)

        info.lift()
        #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
        info.tk.call('wm','iconphoto',info._w,self.icondata)



        self.trna_Nb = ttk.Notebook(info)
        self.trna_Nb.pack(fill='both',expand=True)



        wt_frame = tk.Frame(self.trna_Nb,name="wildtype")  #first notebook Frame
        wt_frame.pack()

        self.trna_Nb.add(wt_frame,text="   WT   ")

        keys = ['TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATA','ATG','GTT','GTC','GTA','GTG',
                'TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG',
                'TAT','TAC','TAA','TAG','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG',
                'TGT','TGC','TGA','TGG','CGT','CGC','CGA','CGG','AGT','AGC','AGA','AGG','GGT','GGC','GGA','GGG']


        aa_keys = ['Phe-F','Leu-L','Ile-I','Met-M','Val-V','Ser-S','Pro-P','Thr-T','Ala-A','Tyr-Y','STOP','His-H','Gln-Q','Asn-N','Lys-K','Asp-D','Glu-E','Cys-C','STOP','Trp-W','Arg-R','Ser-S','Arg-R','Gly-G']

        sense_dict = {'Phe-F':[],'Leu-L':[],'Ile-I':[],'Met-M':[],'Val-V':[],
                          'Ser-S':[],'Pro-P':[],'Thr-T':[],'Ala-A':[],'Tyr-Y':[],'STOP':[],'His-H':[],'Gln-Q':[],'Asn-N':[],'Lys-K':[],
                          'Asp-D':[],'Glu-E':[],'Cys-C':[],'Trp-W':[],'Arg-R':[],'Gly-G':[]}

        try:
            self.sms.POI.nt_seq
            for i in range(len(self.sms.aa_keys)):
                key = self.sms.aa_keys[i]


                if key != '*':
                    skeys = sense_dict.keys()
                    for skey in skeys:
                        if skey[-1] == key:
                            sense_dict[skey] = [self.sms.POI.codon_sensitivity[0][i]]
        except:
            pass







        aan =    np.array([2,6,3,1,4,4,4,4,4,2,2,2,2,2,2,2,2,2,1,1,4,2,2,4])
        aan2 =    np.array([0,2,2,4,3,1,4,4,4,4,4,2,2,2,2,2,2,2,2,2,1,1,4,2,2,4])
        aansum = np.cumsum(aan)
        aansum2 = np.cumsum(aan2)
        colorpalette = ['#e9efb1','#c2d4dd','#ddeedd','#f9d3ac','#b0aac0','#d8a9b8']
        ukeys = []
        for key in keys:
            ukeys.append(key.replace('T','U'))

        k=0
        m = 0
        span =2
        niter = 1

        self.wt_entries = []
        self.fast_entries = []
        self.slow_entries = []

        for i in range(4):
            for j in range(16):


                if k in aansum:


                    m +=1
                    span = aan[m]
                    c = tk.Label(wt_frame,text =aa_keys[m],bg = colorpalette[m%len(colorpalette)],font=('SystemLabelText',global_font_size))
                    c.grid(row=j,column=5*i+2,rowspan = span,sticky=tk.N+tk.S+tk.W+tk.E)



                    sense_label = tk.Label(wt_frame,text='---',relief=tk.GROOVE,bg = colorpalette[m%len(colorpalette)],font=('SystemLabelText',global_font_size))
                    sense_label.grid(row=j,column=5*i+4,rowspan = span,sticky=tk.N+tk.S+tk.W+tk.E)


                    try:
                        self.sms.POI.nt_seq

                        sense_label.config(text=sense_dict[aa_keys[m]][0])
                    except:
                        pass



                b = tk.Label(wt_frame,text=keys[k],relief=tk.GROOVE,bg = colorpalette[m%len(colorpalette)],font=('SystemLabelText',global_font_size))
                b.grid(row=j,column=5*i+1,sticky=tk.N+tk.S+tk.W+tk.E)


                if m ==0:
                    if niter ==1:
                        span = aan[m]
                        c = tk.Label(wt_frame,text =aa_keys[m],bg = colorpalette[m%len(colorpalette)],font=('SystemLabelText',global_font_size))
                        c.grid(row=j,column=5*i+2,rowspan = span,sticky=tk.N+tk.S+tk.W+tk.E)
                        niter+=1

                        sense_label = tk.Label(wt_frame,text='---',relief=tk.GROOVE,bg = colorpalette[m%len(colorpalette)],font=('SystemLabelText',global_font_size))
                        sense_label.grid(row=j,column=5*i+4,rowspan = span,sticky=tk.N+tk.S+tk.W+tk.E)

                        try:
                            self.sms.POI.nt_seq

                            sense_label.config(text=sense_dict[aa_keys[m]][0])
                        except:
                            pass


                a = tk.Entry(wt_frame,width=10,bg = colorpalette[m%len(colorpalette)],font=('SystemLabelText',global_font_size))
                a.grid(row=j,column=5*i+3,sticky=tk.N+tk.S+tk.W+tk.E)

                self.wt_entries.append(a)


                a.insert(tk.END,self.sms.strGeneCopy[keys[k]])







                k+=1

        for i in range(5):
            spacer = tk.Label(wt_frame,text='  ')
            spacer.grid(row=0,column=5*i,rowspan = 16,sticky=tk.N+tk.S+tk.W+tk.E)

        self.wt_defaults = copy.deepcopy(self.sms.strGeneCopy)
        self.ft_defaults = copy.deepcopy(self.sms.strGeneCopy_fast)
        self.sl_defaults = copy.deepcopy(self.sms.strGeneCopy_slow)

        wtbuttons = tk.Frame(wt_frame,relief=tk.GROOVE)
        wtbuttons.grid(row=16,column=0,columnspan=19,sticky=tk.E)
        wt_default = tk.Button(wtbuttons,text='Default',command=self.set_default_wt,font=('SystemButtonText',global_font_size) )
        wt_default.grid(row=0,column=0,padx=3)

        self.wt_mult_entry = tk.Entry(wtbuttons,width=8)
        self.wt_mult_entry.grid(row=0,column=2)

        wt_label = tk.Label(wtbuttons,text='Multiply by: ',font=('SystemButtonText',global_font_size) )
        wt_label.grid(row=0,column=1)

        wt_mult = tk.Button(wtbuttons,text='Enter',command=self.wt_mult_go,font=('SystemButtonText',global_font_size))
        wt_mult.grid(row=0,column=4,padx=1)

        wt_mult = tk.Button(wtbuttons,text='Update',command=self.update_tRNA_dicts,font=('SystemButtonText',global_font_size))
        wt_mult.grid(row=0,column=5,padx=1)





        '''
        fast_frame = tk.Frame(self.trna_Nb,name="fast")  #first notebook Frame
        fast_frame.pack()

        self.trna_Nb.add(fast_frame,text="   Fast   ")

        k=0
        m = 0
        span =2
        niter = 1

        for i in range(4):
            for j in range(16):


                if k in aansum:


                    m +=1
                    span = aan[m]
                    c = tk.Label(fast_frame,text =aa_keys[m],bg = colorpalette[m%len(colorpalette)])
                    c.grid(row=j,column=4*i+2,rowspan = span,sticky=tk.N+tk.S+tk.W+tk.E)

                b = tk.Label(fast_frame,text=keys[k],relief=tk.GROOVE,bg = colorpalette[m%len(colorpalette)])
                b.grid(row=j,column=4*i+1,sticky=tk.N+tk.S+tk.W+tk.E)


                if m ==0:
                    if niter ==1:
                        span = aan[m]
                        c = tk.Label(fast_frame,text =aa_keys[m],bg = colorpalette[m%len(colorpalette)])
                        c.grid(row=j,column=4*i+2,rowspan = span,sticky=tk.N+tk.S+tk.W+tk.E)
                        niter+=1


                a = tk.Entry(fast_frame,width=5,bg = colorpalette[m%len(colorpalette)])
                a.grid(row=j,column=4*i+3,sticky=tk.N+tk.S+tk.W+tk.E)

                self.fast_entries.append(a)

                a.insert(tk.END,self.sms.strGeneCopy_fast[keys[k]])


                k+=1

        for i in range(5):
            spacer = tk.Label(fast_frame,text='  ')
            spacer.grid(row=0,column=4*i,rowspan = 16,sticky=tk.N+tk.S+tk.W+tk.E)


        ftbuttons = tk.Frame(fast_frame,relief=tk.GROOVE)
        ftbuttons.grid(row=16,column=0,columnspan=19,sticky=tk.E)
        ft_default = tk.Button(ftbuttons,text='Default',command=self.set_default_ft)
        ft_default.grid(row=0,column=0,padx=3)

        self.ft_mult_entry = tk.Entry(ftbuttons,width=8)
        self.ft_mult_entry.grid(row=0,column=2)

        ft_label = tk.Label(ftbuttons,text='Multiply by: ' )
        ft_label.grid(row=0,column=1)

        ft_mult = tk.Button(ftbuttons,text='Enter',command=self.ft_mult_go)
        ft_mult.grid(row=0,column=4,padx=1)



        slo_frame = tk.Frame(self.trna_Nb,name="slow")  #first notebook Frame
        #slo_frame.pack()

        self.trna_Nb.add(slo_frame,text="   Slow   ")

        k=0
        m = 0
        span =2
        niter = 1

        for i in range(4):
            for j in range(16):


                if k in aansum:


                    m +=1
                    span = aan[m]
                    c = tk.Label(slo_frame,text =aa_keys[m],bg = colorpalette[m%len(colorpalette)])
                    c.grid(row=j,column=4*i+2,rowspan = span,sticky=tk.N+tk.S+tk.W+tk.E)

                b = tk.Label(slo_frame,text=keys[k],relief=tk.GROOVE,bg = colorpalette[m%len(colorpalette)])
                b.grid(row=j,column=4*i+1,sticky=tk.N+tk.S+tk.W+tk.E)


                if m ==0:
                    if niter ==1:
                        span = aan[m]
                        c = tk.Label(slo_frame,text =aa_keys[m],bg = colorpalette[m%len(colorpalette)])
                        c.grid(row=j,column=4*i+2,rowspan = span,sticky=tk.N+tk.S+tk.W+tk.E)
                        niter+=1


                a = tk.Entry(slo_frame,width=5,bg = colorpalette[m%len(colorpalette)])
                a.grid(row=j,column=4*i+3,sticky=tk.N+tk.S+tk.W+tk.E)

                self.slow_entries.append(a)

                a.insert(tk.END,self.sms.strGeneCopy_slow[keys[k]])


                k+=1

        for i in range(5):
            spacer = tk.Label(slo_frame,text='  ')
            spacer.grid(row=0,column=4*i,rowspan = 16,sticky=tk.N+tk.S+tk.W+tk.E)
        '''

    def set_default_wt(self):

        keys = ['TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATA','ATG','GTT','GTC','GTA','GTG',
                'TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG',
                'TAT','TAC','TAA','TAG','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG',
                'TGT','TGC','TGA','TGG','CGT','CGC','CGA','CGG','AGT','AGC','AGA','AGG','GGT','GGC','GGA','GGG']
        k=0
        for entry in self.wt_entries:
            entry.delete(0,tk.END)
            entry.insert(tk.END,self.wt_defaults[keys[k]])
            k+=1

    def wt_mult_go(self):

        try:
            mult = float(self.wt_mult_entry.get())
        except:
            try:

                mult = float(eval_expression(self.wt_mult_entry.get()))
            except:
                return

        for entry in self.wt_entries:
            try:
                val = np.round(float(entry.get())*mult,5)

                entry.delete(0,tk.END)
                entry.insert(tk.END, val  )
            except:
                pass
        self.update_tRNA_dicts()


    '''
    def set_default_ft(self):

        keys = ['TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATA','ATG','GTT','GTC','GTA','GTG',
                'TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG',
                'TAT','TAC','TAA','TAG','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG',
                'TGT','TGC','TGA','TGG','CGT','CGC','CGA','CGG','AGT','AGC','AGA','AGG','GGT','GGC','GGA','GGG']
        k=0
        for entry in self.ft_entries:
            entry.delete(0,tk.END)
            entry.insert(tk.END,self.ft_defaults[keys[k]])
            k+=1

    def ft_mult_go(self):
        try:
            mult = float(self.ft_mult_entry.get())
        except:
            try:
                mult = (self.ft_mult_entry.get())
            except:
                return




        for entry in self.ft_entries:
            try:
                val = np.round(float(entry.get())*mult,5)

                entry.delete(0,tk.END)
                entry.insert(tk.END, val  )
            except:
                pass
        self.update_tRNA_dicts()
    '''

    def update_tRNA_dicts(self):
        keys = ['TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATA','ATG','GTT','GTC','GTA','GTG',
                'TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG',
                'TAT','TAC','TAA','TAG','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG',
                'TGT','TGC','TGA','TGG','CGT','CGC','CGA','CGG','AGT','AGC','AGA','AGG','GGT','GGC','GGA','GGG']
        k = 0
        for entry in self.wt_entries:

            self.sms.strGeneCopy[keys[k]] = float(entry.get())
            k+=1

        self.sms._rSNAPsim__update_sensitivity()






    def open_guide(self):
        info = tk.Toplevel(self.parent)

        info.title('Guide')
        infoframe = tk.Frame(info)
        infoframe.pack(fill=tk.BOTH)

        info.lift()
        #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
        info.tk.call('wm','iconphoto',info._w,self.icondata)

        self.guide_Nb = ttk.Notebook(info)       #set up the notebook tabs
        self.guide_Nb.pack(fill='both',expand=True)

        seq_frame = tk.Frame(self.guide_Nb,name="first")  #first notebook Frame
        seq_frame.pack()

        self.guide_Nb.add(seq_frame,text="   Running an SSA   ")


        tc_seq = tk.Frame(seq_frame)
        tc_seq.pack(side=tk.LEFT)

        tc_1 = tk.Label(tc_seq,text='Open a file')
        tc_1.grid(row=0,column=0,padx=2,pady=2)

        self.tc_labels = tk.Frame(seq_frame)
        self.tc_labels.pack(side=tk.LEFT)

        im = PIL.Image.open(os.path.join("icons","seq0.png"))
        photo = PIL.ImageTk.PhotoImage(im)

        self.guide_index = 0

        self.label1seq = tk.Label(self.tc_labels,image=photo)
        self.label1seq.image = photo
        self.label1seq.grid(row=0,column=0,padx=2,pady=2,columnspan=5)

        backbutton = tk.Button(self.tc_labels,text='back',command=self.back_guide)
        nextbutton = tk.Button(self.tc_labels,text='next',command=self.next_guide)
        firstbutton = tk.Button(self.tc_labels,text='first',command=self.first_guide)

        backbutton.grid(row=1,column=0)
        nextbutton.grid(row=1,column=1)
        firstbutton.grid(row=1,column=2)



    def back_guide(self):
        try:
            self.guide_index -=1
            im = PIL.Image.open(os.path.join("icons",("seq" + str(self.guide_index) + ".png")))
            photo = PIL.ImageTk.PhotoImage(im)
            self.label1seq.config(image=photo)
            self.label1seq.image = photo
        except:
            pass


    def next_guide(self):
        try:
            self.guide_index +=1
            im = PIL.Image.open(os.path.join("icons",("seq" + str(self.guide_index) + ".png")))
            photo = PIL.ImageTk.PhotoImage(im)
            self.label1seq.config(image=photo)
            self.label1seq.image = photo
        except:
            pass


    def first_guide(self):
        try:
            self.guide_index = 0
            im = PIL.Image.open(os.path.join("icons",("seq" + str(self.guide_index) + ".png")))
            photo = PIL.ImageTk.PhotoImage(im)
            self.label1seq.config(image=photo)
            self.label1seq.image = photo
        except:
            pass

    def save_seq(self):
        try:
            self.gb_rec.name
            gb = True
        except:
            gb = False

        if gb:
            try:
                filename = tfd.asksaveasfilename(title = "Save as",initialfile  =self.seq_name, defaultextension='.txt',filetypes = (("text files","*.txt"),("genbank file","*.gb")))
            except:
                pass
        else:

            try:
                filename = tfd.asksaveasfilename(title = "Save as",initialfile  =self.seq_name, defaultextension='.txt',filetypes = (("text files","*.txt")))
            except:
                pass

        f = open(filename,'w')
        if '.txt' in filename:


            np.savetxt(f,[self.seq_name, self.seq_str],delimiter=',',fmt='%s')
            f.write(('\r\n'))



        if '.gb' in filename:



            f.write(self.gb_rec.format('gb'))




        f.close()


    def poll_genebank(self,event=None):
        asc = self.gbid_entry.get()

        pollwindow = tk.Toplevel(self.parent)
        pollwindow.title('Accessing Entrez...')
        pollwindow.lift()
        #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
        pollwindow.tk.call('wm','iconphoto',pollwindow._w,self.icondata)
        pollwindow.focus_force()
        pollframe = tk.Frame(pollwindow)

        pollframe.pack(expand=True,fill='both',side='top')
        plabel = tk.Label(pollframe,text='Searching NCBI        ',font=('SystemButtonText',16))
        plabel.grid(row=0,column=1,sticky=tk.W)


        x = 50
        y = 50

        pollcanvas = tk.Canvas(pollframe,width=90,height=90)
        pollcanvas.grid(row=0,column=0)

        r1 = 30
        r2=15
        r3=5

        pollcanvas.create_rectangle([x-22,y-22,x+22,y+22],outline='white',width=1,fill='white')
        rect1 = pollcanvas.create_rectangle([x-r1,y-r1,x+r1,y+r1],outline='#fcaa0a',width=12,tag='rect')

        rect2 = pollcanvas.create_rectangle([x-r2,y-r2,x+r2,y+r2],outline='white',width=0,tag='rect')#,fill='white')

        rect3 = pollcanvas.create_rectangle([x-r3,y-r3,x+r3,y+r3],outline='#fcaa0a',width=12,tag='rect')#,fill='#fcaa0a')



        r1 = 30
        r2=2
        r3=15

        for i in range(1,220):
            time.sleep(.02)
            r1 +=1
            r2 +=1
            r3 +=1

            if r1 >= 30:
                r1 = 1
            if r2 >= 30:
                r2 = 1
            if r3>= 30:
                r3 = 1

            radii = [r1,r2,r3]
            k = 0
            for obj in [rect1,rect2,rect3]:
                pollcanvas.coords(obj,(x-radii[k],y-radii[k],x+radii[k],y+radii[k]))
                k+=1

            if i == 40:
                plabel.configure(text='Searching NCBI        ')
                pollwindow.update()
            if i == 80:
                plabel.configure(text='Searching NCBI        ')
                pollwindow.update()


            if i == 120:
                Entrez.email = "wsraymon@rams.colostate.edu"
                Entrez.tool = 'SingleMoleculeSimulator'

                er = False
                try:
                    handle =  Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=asc)
                    gb_record = SeqIO.read(handle, "genbank") #using "gb" as an alias for "genbank"
                    handle.close()
                except:
                    er = True
                    return
                plabel.configure(text='Searching NCBI        ')
                pollwindow.update()


            if i == 160:
                plabel.configure(text='Searching NCBI        ')
                pollwindow.update()

            if i == 200:
                plabel.configure(text='Searching NCBI        ')
                pollwindow.update()

            pollcanvas.update_idletasks()



        if er == True:
            tmb.showinfo(title='HTTP Error', message='Could not find specified ascession ID')

            pollwindow.destroy()
            return
        pollwindow.destroy()

        self.gb_rec = gb_record

        self.seqdisp.config(state='normal')
        self.seqdisp.delete(0,tk.END)

        self.seq_name = gb_record.name

        self.seqdisp.insert(tk.END,self.seq_name)

        self.seqdisp.config(state='readonly')
        self.seqtext = str(gb_record.seq)

        self.gb_obj = gb_record

        self.gs_text.config(state='normal')
        self.gs_text.delete(0.0,tk.END)

        self.gs_text.insert(tk.END,self.seqtext.upper())

        self.gs_text.config(state='disabled')


        self.sms.sequence_str = str(gb_record.seq)
        self.sms.sequence_name = gb_record.name
        self.sms.gb_obj = gb_record

        self.sms.get_orfs(self.sms.sequence_str, min_codons = 80)
        self.sms.get_temporal_proteins()

        if len(self.sms.pois) !=0:
            self.multiple_poi = False
            self.sms.analyze_poi(self.sms.pois[0],self.sms.pois_seq[0])

        else:
            POIs = []
            for i in range(len(self.gb_obj.features)):

                try:

                    self.gb_obj.features[i].qualifiers['translation']
                    POIs.append(self.gb_obj.features[i])
                except:
                    pass

            if len(POIs) !=0:
                self.POI_gb_obj = POIs
            else:
                self.POI_gb_obj = [self.gb_obj.features[0]]
                POIs = [self.gb_obj.features[0]]

                for i in range(0,8):
                    self.sms.get_orfs(self.sms.sequence_str, min_codons = 80-i*10)
                    self.sms.get_temporal_proteins()

                    if len(self.sms.proteins['1']) != 0:
                        break
                    if len(self.sms.proteins['2']) != 0:
                        break
                    if len(self.sms.proteins['3']) != 0:
                        break



            self.poi_list = ['---']
            for i in range(len(POIs)):

                self.poi_list.append(str(POIs[i].location))


            self.poi_menu['menu'].delete(0,'end')
            for poi in self.poi_list:
                self.poi_menu['menu'].add_command(label=poi,command=tk._setit(self.poi_var,poi))




        if len(self.poi_list) >=2:

            self.poi_var.set(self.poi_list[1])



        if len(self.sms.pois) != 0:
            #self.Seq_props = [gne,tte,nee,epe,gle,tle,ttle]

            self.insert_entry(self.Seq_props[0],self.sms.POI.name)


            tags = ''
            for text in self.sms.POI.tag_types:
                tags = tags + text + ' '
            self.insert_entry(self.Seq_props[1],tags)

            eplen = ''
            for key in self.sms.POI.tag_types:
                eplen = eplen + str(len(self.sms.POI.tag_epitopes[key])) + ' '
            self.insert_entry(self.Seq_props[2],eplen)

            ep_pos = ''
            ep_pos2 = []
            for key in self.sms.POI.tag_types:
                ep_pos = ep_pos + str(self.sms.POI.tag_epitopes[key]) + ' '
                ep_pos2 = ep_pos2 + self.sms.POI.tag_epitopes[key]

            self.insert_entry(self.Seq_props[3],ep_pos)

            self.insert_entry(self.Seq_props[4],str(self.sms.POI.gene_length))

            self.insert_entry(self.Seq_props[5],str(self.sms.POI.total_length - self.sms.POI.gene_length))
            self.insert_entry(self.Seq_props[6],str(self.sms.POI.total_length))
                
            self.fill_defaults()

    
            self.plot_sequence(self.ax,self.sms.POI.gene_length,(self.sms.POI.total_length - self.sms.POI.gene_length)
                , self.sms.POI.total_length ,ep_pos2,tags,self.sms.POI.name)

            self.vis_canvas.draw()
            self.update_idletasks()
            
            
            
            
        

    def flatten_alpha_for_mac(self,img):
        '''
        Strip alpha values for rgba image formats (happens when opening images in mac OS)
        '''
        alpha_layers = img.split()[-1]
        allvalues = alpha_layers.tobytes()
        newalphas = []

        for i in range(0, len(allvalues)):
            if allvalues[i] < max(allvalues) + 30:
                newalphas.append(0)
            else:
                newalphas.append(255)
        mask = PIL.Image.frombytes('L', img.size, bytes(newalphas))
        img.putalpha(mask)

        return img

    def open_info(self):
        info = tk.Toplevel(self.parent)

        info.title('Information')
        infoframe = tk.Frame(info)
        infoframe.pack(fill=tk.BOTH)

        info.lift()
        #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
        info.tk.call('wm','iconphoto',info._w,self.icondata)

        infolabel = tk.Label(infoframe,text='rSNAPsim info',font=('SystemButtonText',12,'bold') )
        infolabel.grid(row=0,column=0,padx=2,pady=2,sticky=tk.W)

        i2label = tk.Label(infoframe,bg='#FFFFFF',wraplength=920,text=' The rSNAPsim GUI allows the user to load a gene sequence / tagged sequence and simulate ribosomal translation trajectorys via a Stochastic Simulation Algorithim (SSA)')
        i2label.grid(row=1,column=0,padx=2,pady=2,sticky=tk.W)

        im = PIL.Image.open(os.path.join("icons","flowchart.png"))
        if platform.system() != 'Windows':
            im = self.flatten_alpha_for_mac(im)


        photo = PIL.ImageTk.PhotoImage(im)

        imageframe = tk.Frame(infoframe,bg='#FFFFFF')
        imageframe.grid(row=2,column=0,sticky=tk.N+tk.E+tk.W+tk.S)

        imagecanvas = tk.Label(imageframe,image=photo)

        imagecanvas.image = photo
        imagecanvas.pack()



    def color_text(self):



        colors = {'A':'#e85f7a','G':'#69e588','T':'#6998db','C':'#ebed87'}




        self.gs_text.tag_config('A',background=colors['A'])
        self.gs_text.tag_config('G',background=colors['G'])
        self.gs_text.tag_config('T',background=colors['T'])
        self.gs_text.tag_config('C',background=colors['C'])

        text = self.seqtext.upper()

        self.gs_text.config(state='normal')
        self.gs_text.delete(0.0,tk.END)
        for i in range(len(text)):

            self.gs_text.insert(tk.END,text[i],(text[i]))

        self.gs_text.config(state='disabled')


    def open_seq_file(self):
        seq = tfd.askopenfilename(defaultextension='.txt')

        self.seq_file = seq

        if '.txt' in seq:
            with open(seq) as f:
                raw = f.readlines()
            
            
            self.sequence_name = os.path.basename(seq)
            print(self.sequence_name)
            
            if len(raw[0].splitlines()[0])< 50:
                if raw[0].splitlines()[0] != '':
                    name = raw[0].splitlines()[0]
                    self.sequence_name = name

            raw = ''.join(raw)
           
            print(self.sequence_name)
            
            self.sequence_str= ''
            
            
            onlychar = re.split(r'[^A-Za-z]', raw)
            validt = ['A', 'G', 'T', 'C']
            validu = ['A', 'G', 'U', 'C']
            namelen = 0
            
            
            
       
            for i in range(len(onlychar)):
                section = onlychar[i]
              
                if set(section.upper()) == set(validt):
                    self.sequence_str += section.upper()
                    
                
                elif set(section.upper()) == set(validu):
                    self.sequence_str += section.upper()                    
                else:
                    pass
                    

    
            self.seq_name = self.sequence_name
            self.seqtext = self.sequence_str
            
            
            self.seqdisp.config(state='normal')
            self.seqdisp.delete(0,tk.END)
            self.seqdisp.insert(0,self.seq_name)

            self.seqdisp.config(state='readonly')




            self.gs_text.config(state='normal')
            self.gs_text.delete(0.0,tk.END)

            self.gs_text.insert(0.0,self.sequence_str)

            self.gs_text.config(state='disabled')

        if '.gb' in seq:
            gb_record = SeqIO.read(open(seq,'r'),'genbank')
            self.gb_rec = gb_record
            self.seqdisp.config(state='normal')
            self.seqdisp.delete(0,tk.END)

            self.seq_name = gb_record.name

            self.seqdisp.insert(tk.END,self.seq_name)

            self.seqdisp.config(state='readonly')
            self.seqtext = str(gb_record.seq)

            self.gb_obj = gb_record

            self.gs_text.config(state='normal')
            self.gs_text.delete(0.0,tk.END)

            self.gs_text.insert(tk.END,self.seqtext.upper())

            self.gs_text.config(state='disabled')


            self.insert_entry(self.gbid_entry,self.gb_rec.id)
            self.gbid_entry.config(state='normal')




        self.sms.open_seq_file(seq)
        self.sms.get_orfs(self.sms.sequence_str, min_codons = 80)
        self.sms.get_temporal_proteins()

        if len(self.sms.pois) != 0:
            self.multiple_poi = False
            self.sms.analyze_poi(self.sms.pois[0],self.sms.pois_seq[0])

        else:
            POIs = []
            for i in range(len(self.gb_obj.features)):

                try:

                    self.gb_obj.features[i].qualifiers['translation']
                    POIs.append(self.gb_obj.features[i])
                except:
                    pass

            if len(POIs) !=0:
                self.POI_gb_obj = POIs
            else:
                self.POI_gb_obj = self.gb_obj.features[0]


            self.poi_list = ['---']
            for i in range(len(POIs)):

                self.poi_list.append(str(POIs[i].location))


            self.poi_menu['menu'].delete(0,'end')
            for poi in self.poi_list:
                self.poi_menu['menu'].add_command(label=poi,command=tk._setit(self.poi_var,poi))




        if len(self.poi_list) ==2:

            self.poi_var.set(self.poi_list[1])



        if len(self.sms.pois) != 0:
            #self.Seq_props = [gne,tte,nee,epe,gle,tle,ttle]

            self.insert_entry(self.Seq_props[0],self.sms.POI.name)


            tags = ''
            for text in self.sms.POI.tag_types:
                tags = tags + text + ' '
            self.insert_entry(self.Seq_props[1],tags)

            eplen = ''
            for key in self.sms.POI.tag_types:
                eplen = eplen + str(len(self.sms.POI.tag_epitopes[key])) + ' '
            self.insert_entry(self.Seq_props[2],eplen)

            ep_pos = ''
            ep_pos2 = []
            for key in self.sms.POI.tag_types:
                ep_pos = ep_pos + str(self.sms.POI.tag_epitopes[key]) + ' '
                ep_pos2 = ep_pos2 + self.sms.POI.tag_epitopes[key]

            self.insert_entry(self.Seq_props[3],ep_pos)

            self.insert_entry(self.Seq_props[4],str(self.sms.POI.gene_length))

            self.insert_entry(self.Seq_props[5],str(self.sms.POI.total_length - self.sms.POI.gene_length))
            self.insert_entry(self.Seq_props[6],str(self.sms.POI.total_length))

            try:
                self.fullsizecodon.destroy()
            except:
                pass
            
        self.fill_defaults()


    def fill_defaults(self):
        #self.insert_entry(self.ss_gs_inputs[1],self.seq_name)
        if self.ss_gs_inputs[0].get() == '':
            self.ss_gs_inputs[0].insert(tk.END,'500')

        if self.ss_gs_inputs[2].get() == '':
            self.ss_gs_inputs[2].insert(tk.END,'50')

        if self.ss_gs_inputs[1].get() == '':
            self.ss_gs_inputs[1].insert(tk.END,'1000')

        if self.ss_gs_inputs[3].get() == '':
            self.ss_gs_inputs[3].insert(tk.END,'1000')

        #self.ss_rates = [gl_e,eer_e,eir_e]

        self.insert_entry(self.ss_rates[0],(str(self.sms.POI.total_length)+ ' codons'))

        if self.ss_rates[1].get() == '':
            self.ss_rates[1].insert(tk.END,'10')

        if self.ss_rates[2].get() == '':
            self.ss_rates[2].insert(tk.END,'0.03')

        if self.time_res_e.get() == '':
            self.time_res_e.insert(tk.END,'1000')

        cmnorm = self.sms.POI.codon_sensitivity/np.max(self.sms.POI.codon_sensitivity).flatten().tolist()
        #colors = self.__cmap_map(lambda x: .9*x +.2,cm.Wistia)(cmnorm[0])
        colors = cm.coolwarm(cmnorm[0])



        for i in range(len(self.aa_entries)):

            self.aa_entries[i].config(state='normal')
            self.aa_entries[i].delete(0,tk.END)
            self.aa_entries[i].insert(tk.END,str(self.sms.POI.aa_seq.count(self.sms.aa_keys[i])))
            self.aa_entries[i].config(state='readonly')
            self.aa_entries[i].config(readonlybackground='#FFFFFF')


            self.sensitivity_entries[i].config(state='normal')
            self.sensitivity_entries[i].delete(0,tk.END)
            self.sensitivity_entries[i].insert(tk.END,str(self.sms.POI.codon_sensitivity[0,i]))
            self.sensitivity_entries[i].config(state='readonly')


            color = (colors[i][0:3])*255


            dc=[0,0,0]

            dc[0] = hex(int(color[0]))[2:]
            dc[1] = hex(int(color[1]))[2:]
            dc[2] = hex(int(color[2]))[2:]
            if len(dc[0]) < 2:
                dc[0] = '0'+dc[0]
            if len(dc[1]) < 2:
                dc[1] = '0'+dc[1]
            if len(dc[2]) < 2:
                dc[2] = '0'+dc[2]
            color = '#' + dc[0] + dc[1] + dc[2]

            self.sensitivity_entries[i].config(readonlybackground=color)


            ep_pos = ''
            ep_pos2 = []
            for key in self.sms.POI.tag_types:
                ep_pos = ep_pos + str(self.sms.POI.tag_epitopes[key]) + ' '
                ep_pos2 = ep_pos2 + self.sms.POI.tag_epitopes[key]
                
            tags = ''
            for text in self.sms.POI.tag_types:
                tags = tags + text + ' '
            self.insert_entry(self.Seq_props[1],tags)

            self.plot_sequence(self.ax,self.sms.POI.gene_length,(self.sms.POI.total_length - self.sms.POI.gene_length)
                , self.sms.POI.total_length ,ep_pos2,tags,self.sms.POI.name)

            self.vis_canvas.draw()
            self.update_idletasks()






    def insert_entry(self,entry,text):
        '''simpler command for readonly entries
        '''
        entry.config(state='normal')
        entry.delete(0,tk.END)
        entry.insert(tk.END,text)
        entry.config(state='readonly')

    def insert_entry_readable(self,entry,text):
        '''simpler command for readonly entries
        '''
        entry.config(state='normal')
        entry.delete(0,tk.END)
        entry.insert(tk.END,text)
       

    def load_data(self):
        x=1


    def save_data(self):
        x=1


    def change_poi(self,*args):
        '''
        POI selection from the multiple POI detection dropdown menu
        '''

        poi = self.poi_var.get()

        try:

            a = self.POI_gb_obj
        except:
            return

        for i in range(len(self.POI_gb_obj)):
            if poi == str(self.POI_gb_obj[i].location):
                poi_obj = self.POI_gb_obj[i]


        inds = str(poi_obj.location)
        inds = inds[:-3]
        inds = inds.split(':')
        inds[0] = int(inds[0].replace('[',''))
        inds[1] = int(inds[1].replace(']',''))

        nt_seq = self.seqtext[inds[0]:inds[1]]


        nttag = self.sms.tag_full[self.tagtype.get()]

        self.sms.analyze_poi(self.sms.nt2aa(nttag) + poi_obj.qualifiers['translation'][0],nttag+ nt_seq )

        #self.insert_entry(self.Seq_props[0],self.sms.POI.name)


        tags = ''
        for text in self.sms.POI.tag_types:
            tags = tags + text + ' '
        self.insert_entry(self.Seq_props[1],tags)

        eplen = ''
        for key in self.sms.POI.tag_types:
            eplen = eplen + str(len(self.sms.POI.tag_epitopes[key])) + ' '
        self.insert_entry(self.Seq_props[2],eplen)

        ep_pos = ''
        ep_pos2 = []
        for key in self.sms.POI.tag_types:
            ep_pos = ep_pos + str(self.sms.POI.tag_epitopes[key]) + ' '
            ep_pos2 = ep_pos2 + self.sms.POI.tag_epitopes[key]

        self.insert_entry(self.Seq_props[3],ep_pos)

        self.insert_entry(self.Seq_props[4],str(self.sms.POI.gene_length))

        self.insert_entry(self.Seq_props[5],str(self.sms.POI.total_length - self.sms.POI.gene_length))
        self.insert_entry(self.Seq_props[6],str(self.sms.POI.total_length))




        self.plot_sequence(self.ax,self.sms.POI.gene_length,(self.sms.POI.total_length - self.sms.POI.gene_length)
            , self.sms.POI.total_length ,ep_pos2,tags,(self.sms.POI.name + ' ' + poi))

        self.vis_canvas.draw()
        self.update_idletasks()

    def play_sim(self):
        '''
        Plays the simulated cell animation if the file in the entry box is a valid gif file

        The way this works is by loading all the frames into memory then drawing the images every 50 ms on the tkinter canvas "self.animation"
        '''
        fname = str(self.simulate_fname.get() +  self.filetype.get())
        self.pause_sim= False
        if '.gif' == self.filetype.get():
            if fname != self.loaded_gif:
                fr = []
                frame = 1

                nframes= 0
                i = 0

                im = PIL.Image.open(fname)
                from itertools import count

                try:
                    for i in count(1):
                        fr.append(PIL.ImageTk.PhotoImage(im.copy()))
                        im.seek(i)
                except EOFError:
                    pass


                for obj in self.animation.find_all():  #find all objections on the animation canvas
                    #if obj != 1:   #if its not the base image delet
                    self.animation.delete(obj)

                self.animation.create_image(-40,-20,anchor=tk.NW,image=fr[0])
                self.animation.image = fr[0]
                self.update()



                try:
                    delay = float(self.delayentry.get())
                except:
                    delay = 30

                self.delay = float(min([delay,100.0]))
                self.loaded_gif = fname
                self.frames = fr

                for i in range(len(fr)):
                    if i > 0:
                        self.animation.delete(last_img) #be sure to clear the canvas of the past frame so as not to memory overflow
                    last_img = self.animation.create_image(-40,-20,anchor=tk.NW,image=fr[i])
                    self.animation.image = fr[i]

                    self.update()

                    time.sleep(float(delay/1000))


            else:
                for obj in self.animation.find_all():  #find all objections on the animation canvas
                    if obj != 1:   #if its not the base image delet
                        self.animation.delete(obj)


                try:
                    delay = float(self.delayentry.get())
                except:
                    delay = 30

                self.delay = float(min([delay,100.0]))


                for i in range(len(self.frames)):
                    if i > 0:
                        self.animation.delete(last_img) #be sure to clear the canvas of the past frame so as not to memory overflow
                    last_img = self.animation.create_image(-40,-20,anchor=tk.NW,image=self.frames[i])
                    self.animation.image = self.frames[i]
                    self.animation.update()
                    self.update()
                    if self.pause_sim == True:
                        return
                    time.sleep(float(self.delay/1000))

            self.simscale.config(command=self.scale_gif)
            self.simscale.config(to=len(fr)-1)

    def pause_sim(self):
        self.pause_sim = True

    def scale_gif(self,event):
        for obj in self.animation.find_all():  #find all objections on the animation canvas
            #if obj != 1:   #if its not the base image delet
            self.animation.delete(obj)

        last_img = self.animation.create_image(-40,-20,anchor=tk.NW,image=self.frames[int(event)])
        self.animation.image = self.frames[int(event)]

        self.update()


    def simulate_cell(self):
        '''
        Simulate cell function, makes sure the file / file path given is valid
        then calls sms.get_simulated_mov() to make the gif / html file
        '''
        fname = self.simulate_fname.get()



        try:
            backslash = '\\'
            im = PIL.image.open('icons' + backslash + 'flowchart.png')
        except:
            backslash = '/'

        if '\\' in self.simulate_fname.get():
            fname = self.simulate_fname.get()
            fname = fname.replace('\\',backslash)

            fname.split(backslash)[:-1]
            fname_dir = backslash.join(fname.split(backslash)[:-1])
            if os.path.isdir(fname_dir) == False:
                tmb.showerror('Directory Error','Invalid Directory')
                return

        if '/' in self.simulate_fname.get():
            fname = self.simulate_fname.get()
            fname = fname.replace('/',backslash)
            fname.split(backslash)[:-1]
            fname_dir = backslash.join(fname.split(backslash)[:-1])
            if os.path.isdir(fname_dir) == False:
                tmb.showerror('Directory Error','Invalid Directory')
                return

        #if theres already a simulation dialog open, do not open a new one
        try:
            self.sim_entries[0].get()
            return
        except:
            pass


        info = tk.Toplevel(self.parent)

        info.title('Simulate')
        dialogframe = tk.Frame(info)
        dialogframe.pack(fill=tk.BOTH)

        info.lift()
        #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
        info.tk.call('wm','iconphoto',info._w,self.icondata)


        im = PIL.Image.open(os.path.join("icons",("simulation_flowchart" + ".png")))
        photo = PIL.ImageTk.PhotoImage(im)
        piclabel = tk.Label(dialogframe,image=photo)
        piclabel.grid(row=0,column=0,columnspan=6)
        piclabel.image = photo


        dlabel = tk.Label(dialogframe,text='DC (um^2/s)')
        krnal = tk.Label(dialogframe,text='kRNA')
        konl = tk.Label(dialogframe,text='kOn_RNA')
        koffl = tk.Label(dialogframe,text='kOff_RNA')
        kdecayl = tk.Label(dialogframe,text='kDecay_RNA')

        dlabel.grid(row=1,column=0)
        krnal.grid(row=2,column=0)
        konl.grid(row=3,column=0)
        koffl.grid(row=4,column=0)
        kdecayl.grid(row=5,column=0)

        self.krna = tk.StringVar(value='')
        self.krna.trace('w',lambda name, index, mode,krna= self.krna: self.update_rna_pred(self.krna))

        de = tk.Entry(dialogframe,width=15,justify=tk.RIGHT)
        krnae = tk.Entry(dialogframe,width=15,justify=tk.RIGHT,textvariable=self.krna)
        kone= tk.Entry(dialogframe,width=15,justify=tk.RIGHT)
        koffe= tk.Entry(dialogframe,width=15,justify=tk.RIGHT)
        kdecaye= tk.Entry(dialogframe,width=15,justify=tk.RIGHT)

        de.grid(row=1,column=1,pady=2,padx=2)
        krnae.grid(row=2,column=1,pady=2,padx=2)
        kone.grid(row=3,column=1,pady=2,padx=2)
        koffe.grid(row=4,column=1,pady=2,padx=2)
        kdecaye.grid(row=5,column=1,pady=2,padx=2)

        de.insert(tk.END,'1')
        krnae.insert(tk.END,'0.08')
        kone.insert(tk.END,'3')
        koffe.insert(tk.END,'1')
        kdecaye.insert(tk.END,'.005')




        sep = ttk.Separator(dialogframe,orient='vertical')
        sep.grid(row=1,column=2,sticky=tk.NS,rowspan=5)

        frame2 = tk.Frame(dialogframe,bg='#000000')
        frame2.grid(row=1,column=3,rowspan=3,sticky=tk.E + tk.W + tk.S + tk.N,columnspan=2)
        dialogframe.grid_columnconfigure(3,weight=2)
        rnacl = tk.Label(frame2,text='RNA color',fg='#FFFFFF',bg='#000000')
        fluorcl = tk.Label(frame2,text='Tag color',fg='#FFFFFF',bg='#000000')

        self.rnacb = tk.Button(frame2,text='      ',command=self.rna_color_picker,bg='#ff0000',relief=tk.FLAT)
        self.fluorcb = tk.Button(frame2,text='      ',command=self.fluor_color_picker,bg='#00ff00',relief=tk.FLAT)
        self.rnacb.grid(row=1,column=4,pady=2,padx=2,sticky=tk.E)
        self.fluorcb.grid(row=2,column=4,pady=2,padx=2,sticky=tk.E)
        rnacl.grid(row=1,column=3,sticky=tk.W,pady=2,padx=2)
        fluorcl.grid(row=2,column=3,sticky=tk.W,pady=2,padx=2)

        csl = tk.Label(dialogframe,text='Cell size')
        csl.grid(column=3,row=5,pady=2,padx=2)
        cse = tk.Entry(dialogframe,width=8,justify=tk.RIGHT)
        cse.grid(column=4,row=5,pady=2,padx=2)

        cse.insert(tk.END,'50')

        sep = ttk.Separator(dialogframe,orient='vertical')
        sep.grid(row=6,column=0,sticky=tk.EW,columnspan=5)

        fnamel = tk.Label(dialogframe,text='Filename')
        fnamel.grid(row=7,column=0,columnspan=1,padx=2)
        fnamee = tk.Entry(dialogframe,justify=tk.RIGHT)
        fnamee.grid(row=7,column=1,columnspan=4,sticky=tk.E+tk.N+tk.W+tk.S,pady=2,padx=2)

        fnamee.insert(tk.END,'simulated_cell')

        self.sim_entries = [de,krnae,kone,koffe,kdecaye,cse,fnamee]


        dc_tt = tooltip(dlabel,self.tooltips['DC (um^2/s)'][0])
        krna_tt = tooltip(krnal,self.tooltips['kRNA'][0])
        kon_tt = tooltip(konl,self.tooltips['kOn_RNA'][0])
        koff_tt = tooltip(koffl,self.tooltips['kOff_RNA'][0])
        kd_tt = tooltip(kdecayl,self.tooltips['kDecay_RNA'][0])
        cs_tt = tooltip(csl,self.tooltips['Cell size'][0])


        dce_tt = tooltip(de,self.tooltips['DC (um^2/s)'][0])
        krnae_tt = tooltip(krnae,self.tooltips['kRNA'][0])
        kone_tt = tooltip(kone,self.tooltips['kOn_RNA'][0])
        koffe_tt = tooltip(koffe,self.tooltips['kOff_RNA'][0])
        kde_tt = tooltip(kdecaye,self.tooltips['kDecay_RNA'][0])
        cse_tt = tooltip(cse,self.tooltips['Cell size'][0])

        self.expected_rna = tk.Label(dialogframe,text='Expected RNA: ')
        self.expected_rna.grid(row=8,column=0,columnspan=3,sticky=tk.E)


        runsim = tk.Button(dialogframe,text='RUN',command=self.gui_sim_cell)
        runsim.grid(row=8,column=3,columnspan=2,sticky=tk.E)

        try:


            t = int(self.ss_gs_inputs[2].get())
            e = int(float(self.krna.get()) * t)
            std = 2*int(np.sqrt(e))


            if e > 9999 or std > 9999:

                self.expected_rna.config(text=('Expected RNA: ' + "{:.2e}".format(e) + ' +- '+"{:.2e}".format(std)  ))
            else:
                self.expected_rna.config(text=('Expected RNA: ' + str(e) + ' +- '+str(std)   ))

            self.update()
        except:
            pass

        try:
            self.ssa
        except:
            #tmb.showerror('Error','No SSA trajectory run, please run a simulation')
            #return
            pass



    def update_timecourse_plot(self,event):
        try:
            self.ssa.intensity_vec
        except:
            return

        ptype = self.tcplottype.get()

        self.tcax.clear()
        #self.plot_ssa_intensity(self.tcax,self.ssa.intensity_vec)




        if ptype == 'All':
            self.plot_ssa_intensity(self.tcax,self.ssa.intensity_vec)

        if ptype == 'Average Time Course':
            self.plot_ssa_average(self.tcax,self.ssa.intensity_vec)

        if ptype == 'PDF':
            self.plot_ssa_pdf(self.tcax,self.ssa.intensity_vec)


        self.tc_canvas.draw()



    def update_rna_pred(self,event):
        try:
            self.expected_rna.config(text=('Expected RNA: '))
        except:
            return


        try:
            t = int(self.ss_gs_inputs[2].get())
            e = int(float(self.krna.get()) * t)
            std = 2*int(np.sqrt(e))
        except:
            return

        if e > 9999 or std > 9999:

            self.expected_rna.config(text=('Expected RNA: ' + "{:.2e}".format(e) + ' +- '+"{:.2e}".format(std)  ))
        else:
            self.expected_rna.config(text=('Expected RNA: ' + str(e) + ' +- '+str(std)   ))

        self.update()


    def kymline_color_picker(self):
        try:
            color = askcolor()

        except:
            pass

        self.kymlinecolor.config(bg=color[1])

    def kymbg_color_picker(self):
        try:
            color = askcolor()

        except:
            pass
        self.kymbgcolor.config(bg=color[1])



    def rna_color_picker(self):
        try:
            color = askcolor()

        except:
            pass

        self.rnacb.config(bg=color[1])

    def fluor_color_picker(self):
        try:
            color = askcolor()

        except:
            pass
        self.fluorcb.config(bg=color[1])


    def gui_sim_cell(self):

        [de,krnae,kone,koffe,kdecaye,cse,fnamee] = self.sim_entries

        rnacolor = self.rnacb.cget('bg')
        fcolor = self.fluorcb.cget('bg')
        filename = fnamee.get()


        try:
            backslash = '\\'
            im = PIL.image.open('icons' + backslash + 'flowchart.png')
        except:
            backslash = '/'

        if '\\' in filename:
            fname = filename.get()
            fname = fname.replace('\\',backslash)

            fname.split(backslash)[:-1]
            fname_dir = backslash.join(fname.split(backslash)[:-1])
            if os.path.isdir(fname_dir) == False:
                tmb.showerror('Directory Error','Invalid Directory')
                return

        if '/' in filename:
            fname = filename.get()
            fname = fname.replace('/',backslash)
            fname.split(backslash)[:-1]
            fname_dir = backslash.join(fname.split(backslash)[:-1])
            if os.path.isdir(fname_dir) == False:
                tmb.showerror('Directory Error','Invalid Directory')
                return

        cell_radius = float(cse.get())
        diffusion_constant = float(de.get())
        kRNA = float(krnae.get())
        kon = float(kone.get())
        koff = float(koffe.get())
        kdecay = float(kdecaye.get())




        tstep = int(self.ss_gs_inputs[2].get())

        t = np.linspace(0,tstep,tstep)


        dna_s = np.array([[ 0,  0],
                          [ 1, -1]])

        dna_w1 = np.array([[kRNA,0],
                           [0, 0]],dtype=float)


        dna_w0 = np.array([[0],[0]])


        dna_si = GenericSSA(type='linear' )
        dna_si.W1 = dna_w1
        dna_si.W0 = dna_w0
        dna_si.S = dna_s

        dna_si.ti = t[0]
        dna_si.tf = t[-1]
        dna_si.n = 1
        xi = np.zeros((2,1))
        xi[0] = 1
        dna_si.xi = xi
        dna_si.ptimes = len(t)

        dna_si.time_variant = False
        dna_si._solve(1)
        rna_creation_data = dna_si.data




        stoich = np.array([[  0,    0,  1],
                           [  -1,  1, -1],
                           [  1, -1, 0]])

        propensity = np.array([
                              [0, kon, 0],
                              [0, 0,koff],
                              [0,kdecay, 0]], dtype=float)

        w0 = np.array([[0],[0],[0]])

        solver_instance = GenericSSA(type='linear' )
        solver_instance.W1 = propensity
        solver_instance.W0 = w0
        solver_instance.S = stoich

        solver_instance.ti = t[0]
        solver_instance.tf = t[-1]
        solver_instance.n = 1
        xi = np.zeros((3,1))
        xi[1] = 1
        solver_instance.xi = xi
        solver_instance.ptimes = len(t)

        solver_instance.time_variant = False




        print('simulating RNA activation....')



        R = cell_radius
        squarelen = float(R/np.sqrt(2))

        n_RNA_t = np.zeros((len(t),int(np.max(rna_creation_data[1]))))

        nRNA =0
        nparticles = (int(np.max(rna_creation_data[1])))
        for i in range(len(t)):
            while nRNA != rna_creation_data[1][i]:
                data = solver_instance._solve(1)

                rnaonoff = data[1] + 1 - data[0]



                n_RNA_t[i:,nRNA   ] = rnaonoff[i:].flatten()
                nRNA +=1


        rna_particles = n_RNA_t.T
        rna_exist = np.where(rna_particles >0,1,0)



        print('simulating RNA motion....')
        rna_locations = np.empty((nparticles,len(t),2  ))

        dt = t[-1]/len(t)

        delta = diffusion_constant
        x = np.empty((2,len(t)))
        centers = np.zeros(x.shape)

        def linecirc(m,b,xc,yc,r):

            if np.isinf(m) == False:
                a = 1+m**2
                e = 2*(m*(b-yc)-xc)
                c = yc**2+xc**2 + b**2-2*yc*b-r**2
                x = np.roots([a,e,c])

                if np.isreal(x).all() == False:
                    x = [np.nan,np.nan]
                    y = [np.nan,np.nan]
                else:
                    y = [b + m*x[0], b+m*x[1]]

            elif abs(xc-b) > r:
                x = [np.nan,np.nan]
            else:
                x = [b,b]
                step = np.sqrt(r**2-(b-xc)**2)
                y = [yc + step,yc-step]

            return [x[0],y[0]],[x[1],y[1]]


        def dist(x1,y1,x2,y2):
            return np.sqrt((x1-x2)**2+(y1-y2)**2)


        for i in range(nparticles):
            x = np.empty((2,len(t)))
            x[:,0] = np.random.random()*squarelen
            x0 = [  ((R+squarelen/4) - (R-squarelen/4))*np.random.random() + (R-squarelen/4)   ,((R+squarelen/4) - (R-squarelen/4))*np.random.random() + (R-squarelen/4) ]


            x0 = x0 - np.array([R,R])
            x[:,0] =x0
            r = norm.rvs(size=np.array(x0).shape + (len(t),), scale=delta*np.sqrt(dt))



            out = np.empty(r.shape)

            np.cumsum(r, axis=-1, out=out)
            out+= np.expand_dims(np.array(x0),axis=-1)

            #out = np.array([[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36],
                            #[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36]])

            centers = np.zeros(out.shape)
            dists = np.zeros((x.shape[1],1)).flatten()

            incirc = np.hypot(out.T[:, 0]-centers.T[:, 0], out.T[:, 1]-centers.T[:, 1])
            dists[np.where(out[0] !=0)] = incirc[np.where(out[0] !=0)]

            while len(np.where(dists>R)[0])!=0:   #trajectory left the cell
                out = out.T
                left_cell = np.where(dists>R)[0][0]



                pts = [[out[left_cell][0],out[left_cell][1] ],[out[left_cell-1][0],out[left_cell-1][1]]   ]

                p = np.polyfit([out[left_cell][0], out[left_cell-1][0]],[out[left_cell][1], out[left_cell-1][1]],1)
                m = p[0]
                b = p[1]

                intercepts = linecirc(m,b,0,0,R)
                if dist(*tuple(intercepts[0])+tuple(pts[0])) > dist(*tuple(intercepts[1])+tuple(pts[0])):
                    inter = np.array(intercepts[1])
                else:
                    inter = np.array(intercepts[0])

                a = out[left_cell] - inter


                out[left_cell-1:] = out[left_cell-1:] - 2*(np.dot(inter,a)/np.linalg.norm(inter)**2)*inter




                dists = np.zeros((x.shape[1],1)).flatten()
                out = out.T
                incirc = np.hypot(out.T[:, 0]-centers.T[:, 0], out.T[:, 1]-centers.T[:, 1])
                dists[np.where(out[0] !=0)] = incirc[np.where(out[0] !=0)]

            data = ((out.T).T*rna_exist[i].T).T
            data[np.where(rna_exist[i]==0)] = -R
            rna_locations[i,:,:] =  data


        rna_loc_compressed = rna_locations[np.where(np.sum(np.sum(rna_locations+R,axis=1),axis=1) >0)]

        try:
            ssa_obj =  self.ssa
        except:
            ssa_obj = None

        if ssa_obj == None:
            print('no ssa data given')
            print('simulating translation....')

            ssa_obj = self.sms.ssa_solver(nRepetitions = int(rna_loc_compressed.shape[0]),tf = tstep,tstep = tstep)

            ivec = ssa_obj.intensity_vec/np.max(ssa_obj.intensity_vec)
            ivec = ivec.T  #get the intensity vec for the "fluorescence"


        else:
            print('Translation data given')
            print('Given ' + str(ssa_obj.n_traj) + ' Needed '+str(int(rna_loc_compressed.shape[0])) )
            if ssa_obj.n_traj  < int(rna_loc_compressed.shape[0]):
                print('simulating ' + str(int(rna_loc_compressed.shape[0]) - ssa_obj.n_traj) + ' additional trajectories....')
                ssa_obj = self.sms.ssa_solver_append(ssa_obj, n=int(rna_loc_compressed.shape[0]) - ssa_obj.n_traj)
                ivec = ssa_obj.intensity_vec/np.max(ssa_obj.intensity_vec)
                ivec = ivec.T  #get the intensity vec for the "fluorescence"

            else:
                ivec = ssa_obj.intensity_vec[0:int(rna_loc_compressed.shape[0])]/np.max(ssa_obj.intensity_vec[0:int(rna_loc_compressed.shape[0])])
                ivec = ivec[0:int(rna_loc_compressed.shape[0])].T  #get the intensity vec for the "fluorescence"







        print('making movie...')
        #simulate brownian motion
        def update_line(num, xpos,ypos, line):  #function for the FuncAnimation
            if num !=0:

                for child in ax.get_children():  #remove the previous patch collection (green spots)

                    if isinstance(child, PatchCollection):
                        child.remove()
                    if isinstance(child, mpatches.Ellipse):
                        child.remove()

            patches = []
            radi = 3*ivec[inds[num]]   #create a max radius of 3 for intensity vecs


            for x1, y1, r in zip(xpos[inds[num]],ypos[inds[num]], radi):   #make circle objects of radius based on ivec
                circle = mpatches.Circle((x1, y1), r,color=fcolor)
                patches.append(circle)
                #fig.gca().add_artist(circle)

            line.set_data(xpos[inds[num]],ypos[inds[num]])
            line.set_linewidth(0)
            line.set_marker('o')
            line.set_markersize(1)
            line.set_color(rnacolor)
            line.set
            p = PatchCollection(patches,zorder=3,facecolors=(fcolor,))  #create a patch collection to add to axis
            ax.add_collection(p)  #adds the circles to axis


            p = mpatches.Circle((0,0),radius=R,color='black')  #add the black circle
            ax.add_patch(p)


            whitep = mpatches.Ellipse((-R,-R),width=7,height=7,color='white',zorder=5)  #add the black circle
            ax.add_patch(whitep)

            plt.xlabel(str(inds[num]))  #update time label


            return line,
        xpos = rna_loc_compressed.T[0]
        ypos = rna_loc_compressed.T[1]


        filetype='.gif'

        if filetype == '.gif':

            Writer = animation.writers['pillow']
        if filetype == '.html':
            Writer = animation.writers['html']

        #print(Writer.__dict__)
        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
        imagesize = 5
        dpi=90
        fig1 = plt.figure(figsize=(imagesize,imagesize),dpi=dpi)  #make figure
        fig1.tight_layout()

        ax= fig1.add_subplot('111')
        plt.yticks([])
        plt.xticks([])
        p = mpatches.Circle((0,0),radius=R,color='black')  #add the black circle
        ax.add_patch(p)
        plt.gca().set_aspect('equal', adjustable='box')

        l, = plt.plot([], [], 'r-')
        plt.xlim(-R-10, R+10)
        plt.ylim(-R-10, R+10)
        plt.xlabel('0')
        plt.title('Simulated Cell')

        inds = np.linspace(0,len(t)-1,len(t)).astype(int)
        #creates the animation
        line_ani = animation.FuncAnimation(fig1, update_line, tstep, fargs=(xpos,ypos, l),
                                           interval=50, blit=True)

        line_ani.save((filename + filetype), writer=writer)  #save the animation



    def ssa_solver_gui2(self, nt_seq=None, all_k=None, k_elong_mean=10, k_initiation=.03, probePosition=[], n_traj=100, tf=1000, start_time=0, tstep=1000, time_inhibit=0, evaluating_frap=False, evaluating_inhibitor=False,force_python = False):
        '''
        Solve stochastic simulation algorithms (SSA) for the translation simulation.

        *keyword args*

            **nt_seq**, nucleotide sequence to simulate

            **all_k**, the propensity rates for each codon location (obtained via get_k)

            **k_elong_mean**, average elongation rate to normalize by

            **k_initiation**, rate of mRNA translation initiation

            **probePosition**, binary vector of probe positions, i.e. where the tag epitopes start by codon position

            **n_traj**, number of trajectories

            **tf**, final time point

            **tstep**, number of time steps to record from 0 to tf

            **time_inhibit**, inhibition time of translation either, harringtonine assay or FRAP

            **evaluating_frap**, true or false for evaluating frap assay at time_inhibit

            **evaluating_inhibitor**, true or false for evaluating harringtonine at time_inhibit

        *returns*

            **ssa_obj**, a ssa() class containing the raw ribosome posistions simulated and statistics such as intensity vectors from the SSA trajectory group

        '''

        if len(probePosition) == 0:
            '''
            try:
                probePosition = []
                for key in self.POI.tag_epitopes.keys():
                    probePosition = probePosition + self.POI.tag_epitopes[key]
                probePosition = np.unique(probePosition).tolist()
            except:
                print('No POI found')
                #nt_seq = self.tag_full['T_flag'] + nt_seq
            '''
            
            pv,probePosition = self.sms.get_probvec()
        

        if nt_seq == None:
          nt_seq = self.POI.nt_seq
        genelength = int(len(nt_seq)/3)

        if all_k == None:


            codons = nt_seq
            genelength = int(len(codons)/3)
            seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
            k_elongation = np.zeros((1, genelength))
            tRNA_copynumber = np.zeros((1, genelength))

            for i in range(len(seperated_codons)):
                tRNA_copynumber[0, i] = self.strGeneCopy[seperated_codons[i]]

            mean_tRNA_copynumber = np.mean(list(self.strGeneCopy.values()))

            k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean
            all_k = [k_initiation] + k_elongation.flatten().tolist()[:-1] + [10]

        
        if isinstance(probePosition,list):
            probePosition = np.array([probePosition]).astype(int)
            
     

        non_consider_time = start_time
      
        '''
        if probePosition.shape[0] <= 1:
            pv = np.zeros((1, genelength+1)).astype(int).flatten()
            
            for i in range(len(probePosition[0])):
                pv[probePosition[0][i]:] = i+1
        else:
            pv = np.zeros((probePosition.shape[0], genelength+1)).astype(int)
            for j in range(probePosition.shape[0]):
                for i in range(len(probePosition)):
                    pv[j][probePosition[j][i]:] = i+1      
        '''

        npoints = tstep #non_consider_time + tstep
        
        time_vec_fixed = np.linspace(0, npoints-1, npoints, dtype=np.float64)
        truetime = np.linspace(0, tf, tstep, dtype=np.float64)

        rib_vec = []

        solutions = []
        


        evf = int(evaluating_frap)
        evi = int(evaluating_inhibitor)
        try:
            intime = float(time_inhibit)
        except:
            intime = 0

#        if evaluating_frap == True or evaluating_inhibitor == True:
#            for i in range(nRepetitions):
#
#                soln = self.SSA(all_k,time_vec_fixed,inhibit_time=time_inhibit+non_consider_time,FRAP=evaluating_frap,Inhibitor=evaluating_inhibitor)
#                solutions.append(soln)
#        else:

        solutionssave = []
        
        st = time.time() 
        
        try:
            if force_python == True:
                st[0]
                
            rib_vec = []
    
            solutions = []            
            solutionssave = []
            N_rib = 200
            all_results = np.zeros((n_traj, N_rib*len(time_vec_fixed)), dtype=np.int32)
            all_ribtimes = np.zeros((n_traj,int(1.3*all_k[0]*truetime[-1])),dtype=np.float64)
            result = np.zeros((len(time_vec_fixed)*N_rib), dtype=np.int32)
            nribs = np.array([0],dtype=np.int32)
            k = np.array(all_k)
            seeds = np.random.randint(0, 0x7FFFFFF, n_traj)
            all_frapresults = np.zeros((n_traj,N_rib*len(time_vec_fixed)),dtype=np.int32)
            all_collisions = np.zeros((n_traj,int(1.3*all_k[0]*truetime[-1])),dtype=np.int32)
            all_nribs = np.zeros((n_traj,1))
            
            for i in range(n_traj):
                result = np.zeros((len(time_vec_fixed)*N_rib), dtype=np.int32)
                ribtimes = np.zeros((int(1.3*k[0]*truetime[-1])),dtype=np.float64)
                frapresult = np.zeros((len(time_vec_fixed)*N_rib),dtype=np.int32)
                coltimes = np.zeros((int(1.3*k[0]*truetime[-1])),dtype=np.int32)
                nribs = np.array([0],dtype=np.int32)
                
                ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
                #ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
                all_results[i, :] = result
                all_frapresults[i,:] = frapresult
                all_ribtimes[i,:] = ribtimes
                all_collisions[i,:] = coltimes
                all_nribs[i,:] = nribs
    
            for i in range(n_traj):
                soln = all_results[i, :].reshape((N_rib, len(time_vec_fixed)))
                validind = np.where(np.sum(soln,axis=1)!=0)[0]
                if np.max(validind) != N_rib-1:
                    validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)
            
                so = soln[(validind,)]
                
                solutionssave.append(so)
                solutions.append(soln)
            
            collisions = np.array([[]])
            for i in range(n_traj):
                totalrib = all_nribs[0][0]
            
                if totalrib > all_collisions.shape[1]:
                    collisions = np.append(collisions, all_collisions[0][:totalrib])
            
                else:
                   
                    collisions = np.append(collisions, all_collisions[0][:])
            
            sttime = time.time() - st
        except:
            
            print('C++ library failed, Using Python Implementation')
            rib_vec = []
    
            solutions = []            
            solutionssave = []
            N_rib = 200
            collisions = np.array([[]])
            all_results = np.zeros((n_traj, N_rib*len(time_vec_fixed)), dtype=np.int32)
            all_col_points = []
            for i in range(n_traj):
                
                soln,all_ribtimes,Ncol,col_points = self.SSA(all_k, truetime, inhibit_time=time_inhibit+non_consider_time, FRAP=evaluating_frap, Inhibitor=evaluating_inhibitor)
                #soln = soln.reshape((1, (len(time_vec_fixed)*N_rib)))
                
                collisions = np.append(collisions,Ncol)
                validind = np.where(np.sum(soln,axis=1)!=0)[0]
                all_col_points.append(np.array(col_points))
                if np.max(validind) != N_rib-1:
                    validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)
                
                so = soln[(validind,)]
               
                solutionssave.append(so)

                solutions.append(soln)
            
                result = soln.reshape((1, (len(time_vec_fixed)*N_rib)))
                all_results[i, :] = result
            
            sttime = time.time() - st


                #rb = sparse.lil_matrix((len(time_vec_fixed),genelength),dtype=int)
                #for j in range(soln.shape[1]):

                    #if len(np.where(soln[:,j]!=0)[0]) !=0:
                    #print(np.where(soln[:,j]!=0)[0])


                    #rb[j,np.where(soln[:,j]!=0)[0]] = 1


                        #for value in soln[:,j][np.where(soln[:,j]!=0)[0]].astype(int):

                            #rb[j, value-1] = 1

                #rib_vec.append(rb)

        



        no_ribosomes = np.zeros((n_traj, (genelength+1)))
        
        startindex = np.where(truetime >= non_consider_time)[0][0]
        
        #all_results = all_results[:,startindex*N_rib:]

        for i in range(len(solutions)):
            for j in range(len(solutions[0][0][startindex:])):
                rib_pos = solutions[i][startindex:, j][np.nonzero(solutions[i][startindex:, j])]
            
                no_ribosomes[i, rib_pos.astype(int)] += 1
        no_ribosomes = no_ribosomes[:, 1:]

        ribosome_means = np.mean(no_ribosomes, axis=0)
        ribosome_density = ribosome_means/npoints

        no_ribosomes_per_mrna = np.mean(no_ribosomes)
        
 

        if probePosition.shape[0] <=1:
            I = np.zeros((n_traj, len(time_vec_fixed[startindex:])))
         
            
        else:
            I = np.zeros((int(probePosition.shape[0]),n_traj, len(time_vec_fixed[startindex:])))
         

        #I = np.zeros((1,tstep+1))
        
        if evaluating_frap == False:
            if probePosition.shape[0] <=1:
                for i in range(n_traj):
        
                    traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
        
                    I[i, :] = np.sum(pv[traj], axis=1)[startindex:].T
            else:
                for j in range(probePosition.shape[0]):
                    for i in range(n_traj):
            
                        traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
            
                        I[j,i, :] = np.sum(pv[j][traj], axis=1)[startindex:].T                
    
    
            intensity_vec = I
        
        else:
            fraptime = time_inhibit
            
            inds = np.where(truetime > fraptime)

            inds2 = np.where(truetime  < fraptime+20)
            inds = np.intersect1d(inds,inds2)
            endfrap = inds[-1]-1
         
            
            for i in range(n_traj):
    
                traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
                
                nribs = np.sum(solutionssave[i][:,endfrap]!=0)
             
                #ribloc = solutionssave[i][:,endfrap]
                
                #adj_pv = pv[solutionssave[i][:,inds[-1]][:nribs]]
                revI = self.__get_negative_intensity(traj,genelength,pv,truetime,fraptime,fraptime+20)
                print(revI)

                I[i, :inds[0]-startindex] = np.sum(pv[traj], axis=1)[startindex:inds[0]].T
                I[i,inds[0]-startindex:] = 0
                I[i,endfrap-startindex:] = np.sum(pv[traj],axis=1)[endfrap-startindex:].T
                I[i,endfrap-startindex:len(revI)+endfrap-startindex] = I[i,endfrap-startindex:len(revI)+endfrap-startindex] + revI
      
                
                
    
    
            intensity_vec = I




        ssa_obj = ssa()
        ssa_obj.no_ribosomes = no_ribosomes
        ssa_obj.n_traj = n_traj
        ssa_obj.k = all_k
        ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        ssa_obj.rib_density = ribosome_density
        ssa_obj.rib_means = ribosome_means
        ssa_obj.rib_vec = rib_vec
        ssa_obj.intensity_vec = intensity_vec
        ssa_obj.time_vec_fixed = time_vec_fixed
        ssa_obj.time = truetime
        ssa_obj.time_rec = truetime[startindex:]
        ssa_obj.start_time = non_consider_time
        try:
            ssa_obj.col_points = all_col_points
        except:
            pass


        ssa_obj.evaluating_inhibitor = evaluating_inhibitor
        ssa_obj.evaluating_frap = evaluating_frap
        ssa_obj.time_inhibit = time_inhibit
        ssa_obj.solutions = solutionssave
        ssa_obj.solvetime = sttime
        ssa_obj.collisions = collisions
        
        
        try:
            ssa_obj.ribtimes = all_ribtimes[np.where(all_ribtimes > 0)]
        except:
            pass


        #solt = solutions.T

        fragmented_trajectories = []
        fragtimes = []
        maxlen = 0
    
        fragmentspertraj= []
        for k in range(n_traj):
            ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
            changes = ind[1:] - ind[:-1]
            addindexes = np.where(changes > 0)[0]
            subindexes = np.where(changes < 0)[0]
            
            sub = solutions[k][:,1:] - solutions[k][:,:-1]
            neutralindexes = np.unique(np.where(sub < 0)[1])
            neutralindexes = np.setxor1d(neutralindexes, subindexes)
            
            for index in neutralindexes:
                pre = solutions[k][:,index]
                post = solutions[k][:,index+1]
                changecount = 0
                while len(np.where(post - pre < 0)[0]) > 0:
    
                    post = np.append([genelength],post)
                    pre = np.append(pre,0)
                    
                    changecount+=1
                
                for i in range(changecount):
                    addindexes = np.sort(np.append(addindexes,index))
                    subindexes = np.sort(np.append(subindexes,index))
                    
                changes[index] = -changecount
                ind[index] += changecount
             
                
            for index in np.where(np.abs(changes)>1)[0]:
                if changes[index] < 0:
                    for i in range(np.abs(changes[index])-1):
                        subindexes = np.sort(np.append(subindexes,index))
                else:
                    for i in range(np.abs(changes[index])-1):
                        addindexes = np.sort(np.append(addindexes,index))   
                
            truefrags = len(subindexes)
     
                
        
           
            if len(subindexes) < len(addindexes):
                subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(len(truetime)-1)).astype(int))
                
            
            fragmentspertraj.append(len(subindexes))
            
            for m in range(min(len(subindexes),len(addindexes))):
                traj = solutions[k][:, addindexes[m]:subindexes[m]+1]
                traj_ind = changes[addindexes[m]:subindexes[m]+1]
                
                startind = ind[addindexes[m]]
                minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
                fragment = np.array([])
            
                    
                
                iterind = startind
                
                if subindexes[m]-addindexes[m] > 0:
                    if len(minusloc) > 1:
                        if m <= truefrags:
                            for n in range(len(minusloc)-1):
                                iterind = iterind + min(0,traj_ind[minusloc[n]])
                                fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                                
                                
                                
                  
                
                      
                            
                            fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                            
                        else:
                            for n in range(len(minusloc)-1):

                                iterind = iterind + min(0,traj_ind[minusloc[n]])
                                
                                fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                  
                                
                            fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
          
                        
                    
                    else:

                        fragment = solutions[k][startind][addindexes[m]:subindexes[m]+1].flatten()
                   
                
                    
                    fragtimes.append(addindexes[m]+1)
                       
                    
                    fragmented_trajectories.append(fragment)
                    #if m <= truefrags:
                        #kes.append(genelength/truetime[len(fragment)])
            
                    if len(fragment) > maxlen:
                        maxlen = len(fragment)
                    
    
            fragarray = np.zeros((len(fragmented_trajectories), maxlen))
            for i in range(len(fragmented_trajectories)):
                fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
            
        ssa_obj.fragments = fragarray
        ssa_obj.fragtimes = fragtimes
        ssa_obj.frag_per_traj = fragmentspertraj
        ssa_obj.full_frags = truefrags
        ssa_obj.all_results = all_results
        
        if probePosition.shape[0] > 1:
            for i in range(probePosition.shape[0]):
                if i > 0:
                    autocorr_vec2, mean_autocorr2, error_autocorr2, dwelltime2, ke_sim2  = self.get_autocorr(intensity_vec[i], truetime, 0, genelength)
                    autocorr_vec = np.vstack((autocorr_vec,autocorr_vec2))
                    mean_autocorr = np.vstack((mean_autocorr,mean_autocorr2))
                    error_autocorr = np.vstack((error_autocorr,error_autocorr2))
                    dwelltime.append(dwelltime2)
                    ke_sim.append(ke_sim2)
                else:
                    autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.get_autocorr(intensity_vec[i], truetime, 0, genelength)
                    dwelltime = [dwelltime]
                    ke_sim = [ke_sim]
            
        else:
            autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.get_autocorr(intensity_vec, truetime, 0, genelength)
        
        
        ssa_obj.autocorr_vec = autocorr_vec
        ssa_obj.mean_autocorr = mean_autocorr
        ssa_obj.error_autocorr = error_autocorr
        ssa_obj.dwelltime = dwelltime
        ssa_obj.ke_sim = ke_sim
        ssa_obj.ke_true = float(genelength)/np.mean(ssa_obj.ribtimes)
        ssa_obj.probe = probePosition

        return ssa_obj




    def ssa_solver_gui(self,nt_seq=None,all_k=None,k_elong_mean=10,k_initiation=.03,probePosition=[],start_time = 0, nRepetitions=100,tf = 1000,tstep = 1000,multi=False,time_inhibit = False,evaluating_frap = False,evaluating_inhibitor = False,numssa=1):

        '''
        Solve stochastic simulation algorithms (SSA) for the translation simulation.

        *keyword args*

            **nt_seq**, nucleotide sequence to simulate

            **all_k**, the propensity rates for each codon location (obtained via get_k)

            **k_elong_mean**, average elongation rate to normalize by

            **k_initiation**, rate of mRNA translation initiation

            **probePosition**, binary vector of probe positions, i.e. where the tag epitopes start by codon position

            **n_traj**, number of trajectories

            **tf**, final time point

            **tstep**, number of time steps to record from 0 to tf

            **time_inhibit**, inhibition time of translation either, harringtonine assay or FRAP

            **evaluating_frap**, true or false for evaluating frap assay at time_inhibit

            **evaluating_inhibitor**, true or false for evaluating harringtonine at time_inhibit

        *returns*

            **ssa_obj**, a ssa() class containing the raw ribosome posistions simulated and statistics such as intensity vectors from the SSA trajectory group

        '''
        
        force_python = False
        
        n_traj = nRepetitions
        if len(probePosition) == 0:
            '''
            try:
                probePosition = []
                for key in self.POI.tag_epitopes.keys():
                    probePosition = probePosition + self.POI.tag_epitopes[key]
                probePosition = np.unique(probePosition).tolist()
            except:
                print('No POI found')
                #nt_seq = self.tag_full['T_flag'] + nt_seq
            '''
            
            pv,probePosition = self.sms.get_probvec()
                #nt_seq = self.tag_full['T_flag'] + nt_seq
        pv,probePosition = self.sms.get_probvec()
        if nt_seq == None:
            nt_seq = self.POI.nt_seq
        genelength = int(len(nt_seq)/3)

        if all_k == None:


            codons = nt_seq
            genelength = int(len(codons)/3)
            seperated_codons = [codons[i:i+3] for i in range(0, len(codons), 3)] #split codons by 3
            k_elongation = np.zeros((1, genelength))
            tRNA_copynumber = np.zeros((1, genelength))

            for i in range(len(seperated_codons)):
                tRNA_copynumber[0, i] = self.sms.strGeneCopy[seperated_codons[i]]

            mean_tRNA_copynumber = np.mean(list(self.sms.strGeneCopy.values()))

            k_elongation = (tRNA_copynumber / mean_tRNA_copynumber) * k_elong_mean
            all_k = [k_initiation] + k_elongation.flatten().tolist()[:-1] + [10]

        if isinstance(probePosition,list):
            probePosition = np.array([probePosition]).astype(int)


        non_consider_time = start_time
        #pv = np.zeros((1, genelength+1)).astype(int).flatten()

        #for i in range(len(probePosition)):
            #pv[probePosition[i]:] = i+1

        npoints = tstep #non_consider_time + tstep
        offset=self.prog['value']

        time_vec_fixed = np.linspace(0, npoints-1, npoints, dtype=np.float64)
        
        
        truetime = np.linspace(0,tf,tstep,dtype=np.float64)

        rib_vec = []

        solutions = []


        evf = int(evaluating_frap)
        evi = int(evaluating_inhibitor)
        try:
            intime = float(time_inhibit)
        except:
            intime = 0

#        if evaluating_frap == True or evaluating_inhibitor == True:
#            for i in range(nRepetitions):
#
#                soln = self.SSA(all_k,time_vec_fixed,inhibit_time=time_inhibit+non_consider_time,FRAP=evaluating_frap,Inhibitor=evaluating_inhibitor)
#                solutions.append(soln)
#        else:


        solutionssave = []
        
        st = time.time() 
            
        
        try:
            if force_python:
                st[0]
        
            N_rib = 200
                        
            rib_vec = []
    
            solutions = []            
            solutionssave = []
            
            all_results = np.zeros((n_traj, N_rib*len(time_vec_fixed)), dtype=np.int32)
            result = np.zeros((len(time_vec_fixed)*N_rib), dtype=np.int32)
            all_ribtimes = np.zeros((n_traj,int(1.3*all_k[0]*truetime[-1])),dtype=np.float64)
            nribs = np.array([0])
    
            all_frapresults = np.zeros((n_traj,N_rib*len(time_vec_fixed)),dtype=np.int32)
            all_collisions = np.zeros((n_traj,int(1.3*all_k[0]*truetime[-1])),dtype=np.int32)
            all_nribs = np.zeros((n_traj,1))
    
            k = np.array(all_k)
            seeds = np.random.randint(0, 0x7FFFFFF, n_traj)
    
            updatetime = time.time()
    
            for i in range(n_traj):
    
                stime = time.time()
                sstime = time.time()
                
                frapresult = np.zeros((len(time_vec_fixed)*N_rib),dtype=np.int32)
                coltimes = np.zeros((int(1.3*k[0]*truetime[-1])),dtype=np.int32)
                nribs = np.array([0],dtype=np.int32)
                result = np.zeros((len(time_vec_fixed)*N_rib), dtype=np.int32)
                ribtimes = np.zeros((int(1.3*k[0]*truetime[-1])),dtype=np.float64)
                
                ssa_translation.run_SSA(result, ribtimes, coltimes, k[1:-1],frapresult, truetime, k[0], k[-1], evf, evi, intime, seeds[i],nribs)
                
                ttime = time.time()
                simtime = ttime-stime
                if i ==0:
                    ptime = simtime
                else:
                    ptime = (ptime*i + simtime)/(i+1)
    
                if ttime-updatetime > .02:
                    self.prog['value'] = i+offset
                    self.nssa.config(text=str(i+offset))
                    if i < int(.1*n_traj):
                        self.time.config(text=('Av. time per sim: ' + str(np.round(ptime,3)) + 's ETA: ' + str(np.round(ptime*n_traj,3))+ 's'))
    
                    self.newwin.update()
                    updatetime = ttime
    
    
    
                all_results[i, :] = result
                all_ribtimes[i,:] = ribtimes
    
            for i in range(n_traj):
                soln = all_results[i, :].reshape((N_rib, len(time_vec_fixed)))
    
                validind = np.where(np.sum(soln,axis=1)!=0)[0]
                if np.max(validind) != N_rib-1:
                    validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)
    
                so = soln[validind]
    
                solutions.append(so)
                solutionssave.append(so)
                
            collisions = np.array([[]])
                
            for i in range(n_traj):
                totalrib = all_nribs[0][0]
            
                if totalrib > all_collisions.shape[1]:
                    collisions = np.append(collisions, all_collisions[0][:totalrib])
            
                else:
                   
                    collisions = np.append(collisions, all_collisions[0][:])
    
            self.prog['value'] = i+offset
            self.nssa.config(text=str(i+offset))
            #self.time.config(text=('Av. time per sim: ' + str(np.round(ptime,3)) + 's ETA: ' + str(np.round(ptime*n_traj,6))+ 's'))
    
            self.newwin.update()
            updatetime = ttime
            
            sttime = time.time() - st

        
        except:
            print('C++ library failed, Using Python Implementation')
            rib_vec = []
    
            solutions = []            
            solutionssave = []
            N_rib = 200
            collisions = np.array([[]])
            all_results = np.zeros((n_traj, N_rib*len(time_vec_fixed)), dtype=np.int32)
            sstime = time.time()
            updatetime = time.time()
            for i in range(n_traj):
                stime = time.time()
                soln,all_ribtimes,Ncol,colpoints = self.sms.SSA(all_k, truetime, inhibit_time=time_inhibit+non_consider_time, FRAP=evaluating_frap, Inhibitor=evaluating_inhibitor)

                collisions = np.append(collisions,Ncol)
                validind = np.where(np.sum(soln,axis=1)!=0)[0]
                if np.max(validind) != N_rib-1:
                    validind = np.append(np.where(np.sum(soln,axis=1)!=0)[0],np.max(validind)+1)

                so = soln[validind]

                solutions.append(so)
                solutionssave.append(so)
                result = soln.reshape((1, (len(time_vec_fixed)*N_rib)))
                all_results[i, :] = result

                ttime = time.time()
                elap = ttime-sstime
                simtime = ttime-stime
                if i ==0:
                    ptime = simtime
                else:
                    ptime = (ptime*i + simtime)/(i+1)

                if ttime-updatetime > .02:
                    self.prog['value'] = i+offset
                    self.nssa.config(text=str(i+offset))
                    self.time.config(text=('Av. time per traj: ' + str(np.round(ptime,3)) + 's    ETA: ' + str(np.round(ptime*n_traj,2))+ 's   Elapsed: ' + str(np.round(elap,2)) + 's'))
                    self.newwin.update()
                    updatetime = ttime

                #rb = sparse.lil_matrix((len(time_vec_fixed),genelength),dtype=int)
                #for j in range(soln.shape[1]):

                    #if len(np.where(soln[:,j]!=0)[0]) !=0:
                    #print(np.where(soln[:,j]!=0)[0])


                    #rb[j,np.where(soln[:,j]!=0)[0]] = 1


                        #for value in soln[:,j][np.where(soln[:,j]!=0)[0]].astype(int):

                            #rb[j, value-1] = 1

                #rib_vec.append(rb)
        

        self.prog['value'] = self.prog['value']+1
        self.nssa.config(text=('Complete, Collecting... '+ str(self.prog['value'])))
        self.newwin.update()

        no_ribosomes = np.zeros((n_traj, (genelength+1)))
        startindex = np.where(truetime >= non_consider_time)[0][0]

        
        
        #all_results = all_results[:,startindex*N_rib:]

        for i in range(len(solutions)):
            for j in range(len(solutions[0][0][startindex:])):

                #print( solutions[i][startindex:, j])
                rib_pos = solutions[i][:, startindex+j][np.nonzero(solutions[i][:, startindex+j])]

                no_ribosomes[i, rib_pos.astype(int)] += 1
        no_ribosomes = no_ribosomes[:, 1:]

        ribosome_means = np.mean(no_ribosomes, axis=0)
        ribosome_density = ribosome_means/npoints

        no_ribosomes_per_mrna = np.mean(no_ribosomes)




        if probePosition.shape[0] <=1:
            I = np.zeros((n_traj, len(time_vec_fixed[startindex:])))
         
            
        else:
            I = np.zeros((int(probePosition.shape[0]),n_traj, len(time_vec_fixed[startindex:]) ))
        #I = np.zeros((1,tstep+1))

        if evaluating_frap == False:
            if probePosition.shape[0] <=1:
                for i in range(n_traj):
        
                    traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
        
                    I[i, :] = np.sum(pv[traj], axis=1)[startindex:].T
            else:
                for j in range(probePosition.shape[0]):
                    for i in range(n_traj):
            
                        traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
            
                        I[j,i, :] = np.sum(pv[j][traj], axis=1)[startindex:].T                
    
    
            intensity_vec = I
        
        else:
            fraptime = time_inhibit
            
       
            

            inds = np.where(truetime > fraptime)

            inds2 = np.where(truetime  < fraptime+20)
            inds = np.intersect1d(inds,inds2)
            
          
            endfrap = inds[-1]-1

         
            
            for i in range(n_traj):
    
                traj = all_results[i, :].reshape((N_rib, len(time_vec_fixed))).T
                
                nribs = np.sum(solutionssave[i][:,endfrap]!=0)
             
                #ribloc = solutionssave[i][:,endfrap]
                
                #adj_pv = pv[solutionssave[i][:,inds[-1]][:nribs]]
                frap_app = 20

                revI = self.sms.get_negative_intensity(traj,genelength,pv,truetime,fraptime+start_time,fraptime+start_time+frap_app)
                

                I[i, :] = np.sum(pv[traj], axis=1)[startindex:].T
                
                             
                I[i,inds[0]:inds[0]+20] = 0
                #I[i,endfrap-startindex:] = np.sum(pv[traj],axis=1)[endfrap-startindex:].T

                I[i,inds[0]+frap_app:len(revI)+inds[0]+frap_app] = I[i,inds[0]+frap_app:len(revI)+inds[0]+frap_app] + revI
                
               
                
                
            print(np.min(I))
    
            intensity_vec = I



    



        ssa_obj = ssa()
        ssa_obj.no_ribosomes = no_ribosomes
        ssa_obj.n_traj = n_traj
        ssa_obj.k = all_k
        ssa_obj.no_rib_per_mrna = no_ribosomes_per_mrna
        ssa_obj.rib_density = ribosome_density
        ssa_obj.rib_means = ribosome_means
        ssa_obj.rib_vec = rib_vec
        ssa_obj.intensity_vec = intensity_vec
        ssa_obj.time_vec_fixed = time_vec_fixed
        ssa_obj.time = truetime
        ssa_obj.time_rec = truetime[startindex:]
        ssa_obj.start_time = non_consider_time
        try:
            ssa_obj.col_points = all_col_points
        except:
            pass


        ssa_obj.evaluating_inhibitor = evaluating_inhibitor
        ssa_obj.evaluating_frap = evaluating_frap
        ssa_obj.time_inhibit = time_inhibit
        ssa_obj.solutions = solutionssave
        ssa_obj.solvetime = sttime
        ssa_obj.collisions = collisions
        
        
        try:
            ssa_obj.ribtimes = all_ribtimes[np.where(all_ribtimes > 0)]
        except:
            pass


        #solt = solutions.T
        fragmented_trajectories = []
        fragtimes = []
        maxlen = 0
    
        fragmentspertraj= []
        for k in range(n_traj):
            ind = np.array([next(j for j in range(0,solutions[k].shape[0]) if int(solutions[k][j, i]) == 0 or int(solutions[k][j, i]) == -1) for i in range(0, solutions[k].shape[1])])
            changes = ind[1:] - ind[:-1]
            addindexes = np.where(changes > 0)[0]
            subindexes = np.where(changes < 0)[0]
            
            sub = solutions[k][:,1:] - solutions[k][:,:-1]
            neutralindexes = np.unique(np.where(sub < 0)[1])
            neutralindexes = np.setxor1d(neutralindexes, subindexes)
            
            for index in neutralindexes:
                pre = solutions[k][:,index]
                post = solutions[k][:,index+1]
                changecount = 0
                while len(np.where(post - pre < 0)[0]) > 0:
    
                    post = np.append([genelength],post)
                    pre = np.append(pre,0)
                    
                    changecount+=1
                
                for i in range(changecount):
                    addindexes = np.sort(np.append(addindexes,index))
                    subindexes = np.sort(np.append(subindexes,index))
                    
                changes[index] = -changecount
                ind[index] += changecount
             
                
            for index in np.where(np.abs(changes)>1)[0]:
                if changes[index] < 0:
                    for i in range(np.abs(changes[index])-1):
                        subindexes = np.sort(np.append(subindexes,index))
                else:
                    for i in range(np.abs(changes[index])-1):
                        addindexes = np.sort(np.append(addindexes,index))   
                
            truefrags = len(subindexes)
     
                
        
           
            if len(subindexes) < len(addindexes):
                subindexes = np.append(subindexes, (np.ones((len(addindexes)-len(subindexes)))*(len(truetime)-1)).astype(int))
                
            
            fragmentspertraj.append(len(subindexes))
            
            for m in range(min(len(subindexes),len(addindexes))):
                traj = solutions[k][:, addindexes[m]:subindexes[m]+1]
                traj_ind = changes[addindexes[m]:subindexes[m]+1]
                
                startind = ind[addindexes[m]]
                minusloc = [0] + np.where(traj_ind < 0)[0].astype(int).tolist()
                fragment = np.array([])
            
                    
                
                iterind = startind
                
                if subindexes[m]-addindexes[m] > 0:
                    if len(minusloc) > 1:
                        if m <= truefrags:
                            for n in range(len(minusloc)-1):
                                iterind = iterind + min(0,traj_ind[minusloc[n]])
                                fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                                
                                
                                
                  
                
                      
                            
                            fragment = np.append(fragment, traj[0, minusloc[-1]+1:].flatten())
                            
                        else:
                            for n in range(len(minusloc)-1):

                                iterind = iterind + min(0,traj_ind[minusloc[n]])
                                
                                fragment = np.append(fragment, traj[iterind, minusloc[n]+1:minusloc[n+1]+1].flatten()) 
                  
                                
                            fragment = np.append(fragment, traj[m-truefrags, minusloc[-1]+1:].flatten())
          
                        
                    
                    else:

                        fragment = solutions[k][startind][addindexes[m]:subindexes[m]+1].flatten()
                   
                
                    
                    fragtimes.append(addindexes[m]+1)
                       
                    
                    fragmented_trajectories.append(fragment)
                    #if m <= truefrags:
                        #kes.append(genelength/truetime[len(fragment)])
            
                    if len(fragment) > maxlen:
                        maxlen = len(fragment)
                    
    
            fragarray = np.zeros((len(fragmented_trajectories), maxlen))
            for i in range(len(fragmented_trajectories)):
                fragarray[i][0:len(fragmented_trajectories[i])] = fragmented_trajectories[i]
            
        ssa_obj.fragments = fragarray
        ssa_obj.fragtimes = fragtimes
        ssa_obj.frag_per_traj = fragmentspertraj
        ssa_obj.full_frags = truefrags
        ssa_obj.all_results = all_results
        
        if probePosition.shape[0] > 1:
            for i in range(probePosition.shape[0]):
                if i > 0:
                    autocorr_vec2, mean_autocorr2, error_autocorr2, dwelltime2, ke_sim2  = self.sms.get_autocorr(intensity_vec[i], truetime, 0, genelength)
                    autocorr_vec = np.vstack((autocorr_vec,autocorr_vec2))
                    mean_autocorr = np.vstack((mean_autocorr,mean_autocorr2))
                    error_autocorr = np.vstack((error_autocorr,error_autocorr2))
                    dwelltime.append(dwelltime2)
                    ke_sim.append(ke_sim2)
                else:
                    autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.sms.get_autocorr(intensity_vec[i], truetime, 0, genelength)
                    dwelltime = [dwelltime]
                    ke_sim = [ke_sim]
            
        else:
            autocorr_vec, mean_autocorr, error_autocorr, dwelltime, ke_sim = self.sms.get_autocorr(intensity_vec, truetime, 0, genelength)
        
        
        ssa_obj.autocorr_vec = autocorr_vec
        ssa_obj.mean_autocorr = mean_autocorr
        ssa_obj.error_autocorr = error_autocorr
        ssa_obj.dwelltime = dwelltime
        ssa_obj.ke_sim = ke_sim
        ssa_obj.ke_true = float(genelength)/np.mean(ssa_obj.ribtimes)
        ssa_obj.probe = probePosition
        
        ssa.solvetime = np.round(time.time()-sstime,3)
        return ssa_obj



    def popout_ssabar(self):

        global_font_size=8
        #self.ssa_inputs = [pbte,sste,repe,gl_e,eer_e,eir_e]    #perturb time, simtime, repetitions, experiemental rate, initiation rate
        self.newwin = tk.Toplevel(self.parent)
        self.newwin.lift()
        #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
        self.newwin.tk.call('wm','iconphoto',self.newwin._w,self.icondata)
        self.newwin.title(' ')
        newframe = tk.Frame(self.newwin)
        newframe.grid(row=0,column=0)
        lb = tk.Label(newframe,text='Running SSA...',font=('SystemLabelText',8))
        lb.grid(row=0,column=0,sticky=tk.W)
        self.prog = ttk.Progressbar(newframe,orient="horizontal",length=600, mode="determinate")
        self.prog.grid(row=1,column=0,columnspan=7)
        self.prog["value"] = 0
        self.prog["maximum"] = int(self.ssa_inputs[2].get())


        self.time = tk.Label(newframe,text='Time per traj: 0',font=('SystemLabelText',global_font_size))
        self.time.grid(row=0,column=3,sticky=tk.E)
        self.nssa = tk.Label(newframe,text='Traj: 0',font=('SystemLabelText',global_font_size))
        self.nssa.grid(row=0,column=5,sticky=tk.E)
        nssa2 = tk.Label(newframe,text=( '/' + str(int(self.ssa_inputs[2].get()))  ),font=('SystemLabelText',global_font_size) )
        nssa2.grid(row=0,column=6,sticky=tk.W)

        self.newwin.focus_force()
        self.newwin.grab_set()

        tRNA_norm = np.zeros((1,self.sms.POI.total_length))
        nt_seq = self.sms.POI.nt_seq
        aa_seq = self.sms.POI.aa_seq
        seperated_codons = [nt_seq[i:i+3] for i in range(0, len(nt_seq), 3)] #split codons by 3



        for i in range(self.sms.POI.total_length):
            tRNA_norm[0,i] = self.sms.strGeneCopy[seperated_codons[i]]



        tRNA_design = np.zeros((1,self.sms.POI.total_length))
        for i in range(len(self.aa_vars)):

            fs = self.aa_vars[i].get()
            indexes = [m.start() for m in re.finditer(self.sms.aa_keys[i],aa_seq)]
            for index in indexes:

                if fs ==2:
                    tRNA_design[0,index] = self.sms.slow_codons_value[i]
                if fs ==1:
                    tRNA_design[0,index] = self.sms.fast_codons_value[i]
                if fs == 0:
                    tRNA_design[0,index] = tRNA_norm[0,index]


        k_elongation_mean = float(self.ss_rates[1].get())

        k_init = float(self.ss_rates[2].get())

        mean_tRNA_copynumber = np.mean(list(self.sms.strGeneCopy.values()))


        k_elongation_design = (tRNA_design / mean_tRNA_copynumber) * k_elongation_mean


        all_k_design = [k_init] + k_elongation_design.flatten().tolist() + [10]


        tvec = np.linspace(0,  float(self.ssa_inputs[6].get())+float(self.ssa_inputs[1].get()), int(self.time_res_e.get())  + int(self.ssa_inputs[1].get())+1)

      
        pt = self.perturb.get()
        if pt == 'none':
            self.ssa = self.ssa_solver_gui(all_k = all_k_design,nt_seq = self.sms.POI.nt_seq,nRepetitions=int(self.ssa_inputs[2].get()),
                                       probePosition = self.sms.POI.tag_epitopes['T_Flag'],tf = tvec[-1], tstep = len(tvec),
                                       k_elong_mean=float(self.ssa_inputs[4].get()),k_initiation=float(self.ssa_inputs[5].get()),start_time=float(self.ssa_inputs[6].get()))
        if pt == 'harr':

            self.ssa = self.ssa_solver_gui(all_k = all_k_design,nt_seq = self.sms.POI.nt_seq,nRepetitions=int(self.ssa_inputs[2].get()),
                                       probePosition = self.sms.POI.tag_epitopes['T_Flag'],tf = tvec[-1], tstep = len(tvec),
                                       k_elong_mean=float(self.ssa_inputs[4].get()),k_initiation=float(self.ssa_inputs[5].get()),
                                       time_inhibit=float(self.ssa_inputs[0].get()) + float(self.ssa_inputs[6].get()),evaluating_frap=False,evaluating_inhibitor=True,start_time=float(self.ssa_inputs[6].get()))
        if pt == 'frap':

            self.ssa = self.ssa_solver_gui(all_k = all_k_design,nt_seq = self.sms.POI.nt_seq,nRepetitions=int(self.ssa_inputs[2].get()),
                                       probePosition = self.sms.POI.tag_epitopes['T_Flag'],tf = tvec[-1], tstep = len(tvec),
                                       k_elong_mean=float(self.ssa_inputs[4].get()),k_initiation=float(self.ssa_inputs[5].get()),
                                       time_inhibit=float(self.ssa_inputs[0].get()),evaluating_frap=True,evaluating_inhibitor=False,start_time=float(self.ssa_inputs[6].get()) )

        self.newwin.destroy()
        self.parent.focus()
        self.Main_Nb.select(3)




    def popout_ssabar_elon(self):

        self.newwin = tk.Toplevel(self.parent)
        self.newwin.lift()
        #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
        self.newwin.tk.call('wm','iconphoto',self.newwin._w,self.icondata)
        self.newwin.title(' ')
        newframe = tk.Frame(self.newwin)
        newframe.grid(row=0,column=0)
        lb = tk.Label(newframe,text='Running SSA...',font=('SystemLabelText',8))
        lb.grid(row=0,column=0,sticky=tk.W)
        self.prog = ttk.Progressbar(newframe,orient="horizontal",length=600, mode="determinate")
        self.prog.grid(row=1,column=0,columnspan=7)
        self.prog["value"] = 0
        self.prog["maximum"] = 2*int(self.elon_gs_inputs[4].get())

        self.nssa = tk.Label(newframe,text='Traj: 0',font=('SystemLabelText',8))
        self.nssa.grid(row=0,column=5,sticky=tk.E)
        nssa2 = tk.Label(newframe,text=( '/' + str(2*int(self.elon_gs_inputs[4].get()))  ),font=('SystemLabelText',8) )
        nssa2.grid(row=0,column=6,sticky=tk.W)

        self.newwin.focus_force()
        self.newwin.grab_set()



        self.ssa_frap = self.ssa_solver_gui(self.sms.POI.nt_seq,nRepetitions=int(self.elon_gs_inputs[4].get()),
                                            probePosition = self.sms.POI.tag_epitopes['T_Flag'],tstep = int(self.elon_gs_inputs[2].get()),
                                            k_elongationMean=float(self.ss_rates_elon[1].get()),k_initiation=float(self.ss_rates_elon[2].get()),
                                            time_inhibit=float(self.elon_gs_inputs[3].get()),evaluating_frap=True,evaluating_inhibitor=False,numssa = 2)


        self.ssa_harr = self.ssa_solver_gui(self.sms.POI.nt_seq,nRepetitions=int(self.elon_gs_inputs[4].get()),
                                            probePosition = self.sms.POI.tag_epitopes['T_Flag'],tstep = int(self.elon_gs_inputs[2].get()),
                                            k_elongationMean=float(self.ss_rates_elon[1].get()),k_initiation=float(self.ss_rates_elon[2].get()),
                                            time_inhibit=float(self.elon_gs_inputs[3].get()),evaluating_frap=False,evaluating_inhibitor=True,numssa = 2)


        self.newwin.destroy()
        self.parent.focus()


    def popout_ssabar_codon(self):

        self.newwin = tk.Toplevel(self.parent)
        self.newwin.lift()
        #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
        self.newwin.tk.call('wm','iconphoto',self.newwin._w,self.icondata)
        self.newwin.title(' ')
        newframe = tk.Frame(self.newwin)
        newframe.grid(row=0,column=0)
        lb = tk.Label(newframe,text='Running SSA...')
        lb.grid(row=0,column=0,sticky=tk.W)
        self.prog = ttk.Progressbar(newframe,orient="horizontal",length=600, mode="determinate")
        self.prog.grid(row=1,column=0,columnspan=7)
        self.prog["value"] = 0
        self.prog["maximum"] = 2*int(self.co_gs_inputs[3].get())

        self.nssa = tk.Label(newframe,text='Traj: 0')
        self.nssa.grid(row=0,column=5,sticky=tk.E)
        nssa2 = tk.Label(newframe,text=( '/' + str(2*int(self.co_gs_inputs[3].get()))  ) )
        nssa2.grid(row=0,column=6,sticky=tk.W)

        self.newwin.focus_force()
        self.newwin.grab_set()


        nt_seq = self.sms.POI.nt_seq
        aa_seq = self.sms.POI.aa_seq


        tRNA_design = np.zeros((1,self.sms.POI.total_length))
        tRNA_norm = np.zeros((1,self.sms.POI.total_length))

        seperated_codons = [nt_seq[i:i+3] for i in range(0, len(nt_seq), 3)] #split codons by 3



        for i in range(self.sms.POI.total_length):
            tRNA_norm[0,i] = self.sms.strGeneCopy[seperated_codons[i]]




        for i in range(len(self.aa_vars)):

            fs = self.aa_vars[i].get()
            indexes = [m.start() for m in re.finditer(self.sms.aa_keys[i],aa_seq)]
            for index in indexes:

                if fs ==2:
                    tRNA_design[0,index] = self.sms.slow_codons_value[i]
                if fs ==1:
                    tRNA_design[0,index] = self.sms.fast_codons_value[i]
                if fs == 0:
                    tRNA_design[0,index] = tRNA_norm[0,index]


        k_elongation_mean = float(self.ss_rates_codon[1].get())

        k_init = float(self.ss_rates[2].get())

        mean_tRNA_copynumber = np.mean(list(self.sms.strGeneCopy.values()))

        k_elongation = (tRNA_norm / mean_tRNA_copynumber) * k_elongation_mean
        k_elongation_design = (tRNA_design / mean_tRNA_copynumber) * k_elongation_mean

        all_k_norm = [k_init] + k_elongation.flatten().tolist() + [10]
        all_k_design = [k_init] + k_elongation_design.flatten().tolist() + [10]


        self.ssa_norm =  self.ssa_solver_gui(all_k = all_k_norm,  nRepetitions = int(self.co_gs_inputs[3].get()),evaluating_frap=False,evaluating_inhibitor=False,numssa=2 )
        #self.co_gs_inputs = [gbid_e,sfe,sste,repe]

        self.ssa_design = self.ssa_solver_gui(all_k = all_k_design, nRepetitions = int(self.co_gs_inputs[3].get()),evaluating_frap=False,evaluating_inhibitor=False,numssa=2  )



        self.newwin.destroy()
        self.parent.focus()



    def run_ssa(self):

        try:
            self.sms.POI
        except:
            return
        self.popout_ssabar()


        ptype = self.tcplottype.get()
        self.tcax.clear()

        if ptype == 'All':

            self.plot_ssa_intensity(self.tcax,self.ssa.intensity_vec)
        if ptype == 'Average Time Course':
            self.plot_ssa_average(self.tcax,self.ssa.intensity_vec)

        if ptype == 'PDF':
            self.plot_ssa_pdf(self.tcax,self.ssa.intensity_vec)

        self.tc_canvas.draw()

        self.acax.clear()
        self.acax.cla()
        self.plot_ssa_acc(self.acax,self.ssa.mean_autocorr,self.ssa.error_autocorr)
        self.ac_canvas.draw()

        self.update_ss_frame_data()

        self.rbax.clear()
        self.plot_rib_dense(self.rbax,self.ssa.rib_density)
        self.rb_canvas.draw()

        self.cuax.clear()
        self.plot_cai_codons(self.cuax,self.sms.POI.CAI_codons)
        self.cu_canvas.draw()



        self.update_psim()


    def update_psim(self):

        ssaruns = len(self.psimscroll.inside_frame.winfo_children())
        self.ssas.append(self.ssa)
        if ssaruns %2 == 1:
            bgc = self.default_color
        else:
            bgc = '#DDDDDD'
        newf = tk.Frame(self.psimscroll.inside_frame,bg=bgc)
        newf.pack(fill='x',side='top',expand=True)
        ssalabel = tk.Label(newf,text=('ssa ' + str(ssaruns+1)),bg=bgc,width=4)
        ssalabel.grid(row=0,column=0,padx=2,pady=2)
        numlabel = tk.Label(newf,text=(str(self.ssa.n_traj) + ' traj '),bg=bgc,width=4)
        numlabel.grid(row=0,column=1,padx=2,pady=2)
        nameentry = tk.Entry(newf,width=12)
        nameentry.grid(row=0,column=7,padx=2,pady=2)
        nameentry.insert(tk.END,('ssa' + str(ssaruns+1)))

        timelabel = tk.Label(newf,bg=bgc,text=('time ' + str(self.ssa.time_vec_fixed[0]) + ' --> ' + str(self.ssa.time_vec_fixed[-1]) + '  ' + str(len(self.ssa.time_vec_fixed)   ) + ' steps'),width=25 )
        timelabel.grid(row=0,column=2,padx=2,pady=2)
        harr = 'X'
        if self.ssa.evaluating_inhibitor == True:
            harr = '✓'

        harrlabel = tk.Label(newf, text=('Harr: ' + harr ),bg=bgc,width=6)
        harrlabel.grid(row=0,column=4,padx=2,pady=2)


        harr = 'X'
        if self.ssa.evaluating_frap == True:
            harr = '✓'

        harrlabel = tk.Label(newf, text=('FRAP: ' + harr ),bg=bgc,width=6)
        harrlabel.grid(row=0,column=5,padx=2,pady=2)


        if self.ssa.start_time == self.ssa.time_vec_fixed[0]:
            ss = 'False'
        else:
            ss = 'True'

        starttimelabel = tk.Label(newf,text=('start time: ' + str(self.ssa.start_time) + '     steady state: ' + ss + '  Sim time: ' + str(self.ssa.solvetime)),bg=bgc,width=40)
        starttimelabel.grid(row=0,column=3,padx=2,pady=2)

        savebutton = tk.Button(newf,text='Save',command=lambda j=ssaruns: self.save_ss_data(j),font=('SystemButtonText',8))
        savebutton.grid(row=0,column=8,padx=2,pady=2)

        savebutton = tk.Button(newf,text='Plot',command= lambda j=ssaruns: self.plot_ssa_psim(j),font=('SystemButtonText',8))
        savebutton.grid(row=0,column=9,padx=2,pady=2)
        
        savebutton = tk.Button(newf,text='Del',command= lambda j=ssaruns: self.delete_ss_data(j),font=('SystemButtonText',8))
        savebutton.grid(row=0,column=10,padx=2,pady=2)







    def save_ssa_psim(self,num):
        x=1
        
    def delete_ss_data(self,num):
        
        self.ssa = self.ssas[max(num-1,0)]
        self.ssas.pop(num)
        
        
        for child in self.psimscroll.inside_frame.winfo_children():
            
            child.destroy()
        
        
        for k in range(len(self.ssas)):
            
            
            if k %2 == 1:
                bgc = self.default_color
            else:
                bgc = '#DDDDDD'
            newf = tk.Frame(self.psimscroll.inside_frame,bg=bgc)
            newf.pack(fill='x',side='top',expand=True)
            ssalabel = tk.Label(newf,text=('ssa ' + str(k+1)),bg=bgc,width=4)
            ssalabel.grid(row=0,column=0,padx=2,pady=2)
            numlabel = tk.Label(newf,text=(str(self.ssas[k].n_traj) + ' traj '),bg=bgc,width=4)
            numlabel.grid(row=0,column=1,padx=2,pady=2)
            nameentry = tk.Entry(newf,width=12)
            nameentry.grid(row=0,column=7,padx=2,pady=2)
            nameentry.insert(tk.END,('ssa' + str(k+1)))
    
            timelabel = tk.Label(newf,bg=bgc,text=('time ' + str(self.ssa.time_vec_fixed[0]) + ' --> ' + str(self.ssas[k].time_vec_fixed[-1]) + '  ' + str(len(self.ssas[k].time_vec_fixed)   ) + ' steps'),width=25 )
            timelabel.grid(row=0,column=2,padx=2,pady=2)
            harr = 'X'
            if self.ssas[k].evaluating_inhibitor == True:
                harr = '✓'
    
            harrlabel = tk.Label(newf, text=('Harr: ' + harr ),bg=bgc,width=6)
            harrlabel.grid(row=0,column=4,padx=2,pady=2)
    
    
            harr = 'X'
            if self.ssas[k].evaluating_frap == True:
                harr = '✓'
    
            harrlabel = tk.Label(newf, text=('FRAP: ' + harr ),bg=bgc,width=6)
            harrlabel.grid(row=0,column=5,padx=2,pady=2)
    
    
            if self.ssas[k].start_time == self.ssas[k].time_vec_fixed[0]:
                ss = 'False'
            else:
                ss = 'True'
    
            starttimelabel = tk.Label(newf,text=('start time: ' + str(self.ssas[k].start_time) + '     steady state: ' + ss + '  Sim time: ' + str(self.ssas[k].solvetime)),bg=bgc,width=40)
            starttimelabel.grid(row=0,column=3,padx=2,pady=2)
    
            savebutton = tk.Button(newf,text='Save',command=lambda j=k: self.save_ss_data(j),font=('SystemButtonText',8))
            savebutton.grid(row=0,column=8,padx=2,pady=2)
    
            savebutton = tk.Button(newf,text='Plot',command= lambda j=k: self.plot_ssa_psim(j),font=('SystemButtonText',8))
            savebutton.grid(row=0,column=9,padx=2,pady=2)
            
            savebutton = tk.Button(newf,text='Del',command= lambda j=k: self.delete_ss_data(j),font=('SystemButtonText',8))
            savebutton.grid(row=0,column=10,padx=2,pady=2)

    def plot_ssa_psim(self,num):

        self.ssa = self.ssas[num]


        ptype = self.tcplottype.get()
        self.tcax.clear()

        if ptype == 'All':

            self.plot_ssa_intensity(self.tcax,self.ssa.intensity_vec)
        if ptype == 'Average Time Course':
            self.plot_ssa_average(self.tcax,self.ssa.intensity_vec)

        if ptype == 'PDF':
            self.plot_ssa_pdf(self.tcax,self.ssa.intensity_vec)

        self.tc_canvas.draw()

        self.acax.clear()
        self.acax.cla()
        self.plot_ssa_acc(self.acax,self.ssa.mean_autocorr,self.ssa.error_autocorr)
        self.ac_canvas.draw()

        self.update_ss_frame_data()

        self.rbax.clear()
        self.plot_rib_dense(self.rbax,self.ssa.rib_density)
        self.rb_canvas.draw()

        self.cuax.clear()
        self.plot_cai_codons(self.cuax,self.sms.POI.CAI_codons)
        self.cu_canvas.draw()

        self.update_idletasks()






    def plot_loaded_ssa(self):
        self.tcax.clear()
        self.plot_ssa_intensity(self.tcax,self.ssa.intensity_vec)


        self.tc_canvas.draw()

        self.acax.clear()
        self.acax.cla()
        self.plot_ssa_acc(self.acax,self.ssa.mean_autocorr,self.ssa.error_autocorr)
        self.ac_canvas.draw()

        self.update_ss_frame_data()

        self.rbax.clear()
        self.plot_rib_dense(self.rbax,self.ssa.rib_density)
        self.rb_canvas.draw()

        self.cuax.clear()
        self.plot_cai_codons(self.cuax,self.sms.POI.CAI_codons)
        self.cu_canvas.draw()

    def plot_loaded_frap(self):
        self.frax.clear()

        self.plot_ssa_frap(self.frax,self.ssa_frap.intensity_vec)
        self.fr_canvas.draw()
        time_inds = np.where(self.ssa_frap.time_vec_fixed- 500  > 0)
        tvec = self.ssa_frap.time_vec_fixed[time_inds]-500
        tempind = 50
        recoverythresh = .9
        recoverytime = 1
        application_time = float(self.elon_gs_inputs[3].get())


        tvec_inds = np.where(tvec>=application_time)

        meaniv = np.mean(self.ssa_frap.intensity_vec,axis=0)
        meaniv_norm = meaniv/np.max(meaniv)
        try:
            recoverytime = tvec[np.where(meaniv_norm[tvec_inds] >= recoverythresh )][0]+50

        except:
            recoverytime = tvec[np.where(meaniv_norm[tvec_inds] >= .7 )][0]+50


        elongation_rate = np.round(float(self.sms.POI.total_length) / recoverytime,2)

        self.cpv_entries[2].delete(0,tk.END)
        self.cpv_entries[2].insert(tk.END,str(elongation_rate))

        self.cpv_entries[3].delete(0,tk.END)
        self.cpv_entries[3].insert(tk.END,str(recoverytime))

    def plot_loaded_harr(self):
        self.hrax.clear()
        self.plot_ssa_frap(self.hrax,self.ssa_harr.intensity_vec)
        self.hr_canvas.draw()
        meaniv = np.mean(self.ssa_harr.intensity_vec,axis=0)
        meaniv_norm = meaniv/np.max(meaniv)
        norm_iv = self.ssa_harr.intensity_vec/np.max(self.ssa_harr.intensity_vec)
        sem_iv = np.std(self.ssa_harr.intensity_vec,axis=0) / int(self.elon_gs_inputs[4].get())
        runoffindex = np.where(meaniv_norm < .005)[0][0]
        time_inds = np.where(self.ssa_harr.time_vec_fixed- 500  > 0)
        application_time = float(self.elon_gs_inputs[3].get())
        tvec = self.ssa_harr.time_vec_fixed[time_inds]-500
        runoff_time = (tvec)[runoffindex] - application_time

        elongation_rate = np.round(float(self.sms.POI.total_length) / runoff_time,2)

        self.ssa_har_data = [norm_iv,sem_iv,meaniv,meaniv_norm,runoff_time]
        #self.cpv_entries = [ere,rte,ere2,rte2]
        self.cpv_entries[1].delete(0,tk.END)
        self.cpv_entries[1].insert(tk.END,str(runoff_time))
        self.cpv_entries[0].delete(0,tk.END)
        self.cpv_entries[0].insert(tk.END,str(elongation_rate))


    def run_ssa_frap(self):
        #self.elon_gs_inputs = [gbid_e,sfe,sste,pbte,repe]
        #time_inhibit = False,evaulating_frap = False,evaulating_inhibitor = False

        self.popout_ssabar_elon()




        self.frax.clear()


        self.plot_ssa_frap(self.frax,self.ssa_frap.intensity_vec,self.ssa_frap.time_inhibit)
        self.fr_canvas.draw()

        self.hrax.clear()
        self.plot_ssa_frap(self.hrax,self.ssa_harr.intensity_vec,self.ssa_harr.time_inhibit)
        self.hr_canvas.draw()



        meaniv = np.mean(self.ssa_harr.intensity_vec,axis=0)
        meaniv_norm = meaniv/np.max(meaniv)
        norm_iv = self.ssa_harr.intensity_vec/np.max(self.ssa_harr.intensity_vec)
        sem_iv = np.std(self.ssa_harr.intensity_vec,axis=0) / int(self.elon_gs_inputs[4].get())
        ti = self.ssa_harr.time_inhibit
        time_inds = np.where(self.ssa_harr.time_vec_fixed  > 0)
        tvec = self.ssa_harr.time_vec_fixed[time_inds]

        timeindex = np.where(tvec >= ti)[0][0]
        runoffindex = np.where(meaniv_norm[timeindex:] <.005   )[0][0] + timeindex





        application_time = float(self.elon_gs_inputs[3].get())


        runoff_time = (tvec)[runoffindex] - application_time

        elongation_rate = np.round(float(self.sms.POI.total_length) / runoff_time,2)

        self.ssa_har_data = [norm_iv,sem_iv,meaniv,meaniv_norm,runoff_time]
        #self.cpv_entries = [ere,rte,ere2,rte2]
        self.cpv_entries[1].delete(0,tk.END)
        self.cpv_entries[1].insert(tk.END,str(runoff_time))
        self.cpv_entries[0].delete(0,tk.END)
        self.cpv_entries[0].insert(tk.END,str(elongation_rate))

        tempind = 50
        recoverythresh = .9
        recoverytime = 1
        application_time = float(self.elon_gs_inputs[3].get())


        tvec_inds = np.where(tvec>=application_time)

        meaniv = np.mean(self.ssa_frap.intensity_vec,axis=0)
        meaniv_norm = meaniv/np.max(meaniv)




        try:
            #recoverytime = tvec[np.where(meaniv_norm[tvec_inds] >= recoverythresh )][0]+50
            timeindex = np.where(tvec >= ti+20)[0][0]
            runoffindex = np.where(meaniv_norm[timeindex:] >= recoverythresh  )[0][0] + timeindex


            recoverytime = (tvec)[runoffindex] - application_time



        except:
            #recoverytime = tvec[np.where(meaniv_norm[tvec_inds] >= .7 )][0]+50
            timeindex = np.where(tvec >= ti)[0][0]
            runoffindex = np.where(meaniv_norm[timeindex+1:] >= .7  )[0][0] + timeindex

            recoverytime = (tvec)[runoffindex] - application_time


        elongation_rate = np.round(float(self.sms.POI.total_length) / recoverytime,2)

        self.cpv_entries[2].delete(0,tk.END)
        self.cpv_entries[2].insert(tk.END,str(elongation_rate))

        self.cpv_entries[3].delete(0,tk.END)
        self.cpv_entries[3].insert(tk.END,str(recoverytime))


    def run_ssa_codons(self):
        #self.ss_rates_codon = [gl_e,eer_e,eir_e]

        self.popout_ssabar_codon()

        self.accax.clear()
        self.accax.cla()

        index1 = np.where(self.ssa_norm.mean_autocorr < 0.001)[0][0]
        index2 = np.where(self.ssa_design.mean_autocorr<.001)[0][0]

        #index1 = np.where(self.ssa_norm.mean_autocorr !=0)[0][0]
        #index2 = np.where(self.ssa_design.mean_autocorr !=0)[0][0]

        if index1 >= index2:
            index = index1
        if index2 >= index1:
            index = index2

        self.plot_ssa_acc_codon(self.accax,self.ssa_norm.mean_autocorr[:index],self.ssa_norm.error_autocorr[:index],name='Natural')

        self.plot_ssa_acc_codon(self.accax,self.ssa_design.mean_autocorr[:index],self.ssa_design.error_autocorr[:index],color='#3a7ce8',name='Design')
        self.accax.set_ylabel('Autocorrelation')
        self.acc_canvas.draw()
        #self.codon_parameter_entries = [ere2,rte2,drd,nrd]
        params = [str(self.ssa_design.ke_sim),str(self.ssa_norm.ke_sim),str(int(np.round(self.sms.POI.total_length/self.ssa_design.no_rib_per_mrna,0))),str(int(np.round(self.sms.POI.total_length/self.ssa_norm.no_rib_per_mrna,0))) ]
        i = 0
        for par in params:

            self.codon_parameter_entries[i].config(state='normal')
            self.codon_parameter_entries[i].delete(0,tk.END)
            self.codon_parameter_entries[i].insert(tk.END,par)
            self.codon_parameter_entries[i].config(state='readonly')
            i+=1


    def plot_ssa_frap(self,ax,i_vec,itime):
        ax.cla()
        ax.axis([0,int(i_vec.shape[1]),0,1.1] )
        ticks = np.linspace(0,int(i_vec.shape[1]) ,6).astype(int)
        yticks = [0,.25,.5,.75,1]
        ax.set_xticks(ticks)
        ax.set_yticks(yticks)
        meaniv = np.mean(i_vec,axis=0)
        meaniv = meaniv / np.max(meaniv)

        normiv = i_vec /np.max(i_vec)
        semiv = np.std(normiv,axis=0) #/ i_vec.shape[0]

        ax.plot([itime,itime],[0,1.1],color='r',linestyle='--'   )
        ax.plot(meaniv,color=self.main_color,linewidth=2)
        #ax.plot(meaniv-semiv,color='red',linewidth=1,ls='--',alpha=.5)
        #ax.plot(meaniv+semiv,color='red',linewidth=1,ls='--',alpha=.5)
        ax.fill_between(np.linspace(0,int(i_vec.shape[1])-1,int(i_vec.shape[1])), meaniv-semiv,meaniv+semiv,color=self.main_color,alpha=.2   )
        ax.set_xlabel('time (sec)')
        ax.set_ylabel('Intensity (a.u.)')


    def plot_ssa_frap_non(self,ax,i_vec):
        ax.cla()
        ax.axis([0,int(i_vec.shape[1]),0,np.max(i_vec)+5] )
        ticks = np.linspace(0,int(i_vec.shape[1]) ,6).astype(int)
        yticks = np.linspace(0,np.max(i_vec)+5,4).astype(int)
        ax.set_xticks(ticks)
        ax.set_yticks(yticks)
        meaniv = np.mean(i_vec,axis=0)


        for i in range(0,i_vec.shape[0]):

            ax.plot(i_vec[i],alpha=.5,color='gray')


        ax.plot(i_vec[0],color=self.main_color,linewidth=3)
        ax.plot(meaniv,color='blue',linewidth=2)
        ax.set_xlabel('time (sec)')
        ax.set_ylabel('Intensity (a.u.)')

    def plot_ssa_intensity(self,ax,i_vec):
        ax.cla()
        ax.axis([0,int(i_vec.shape[1]),0,np.max(i_vec)+5] )
        ticks = np.linspace(0,int(i_vec.shape[1]) ,6).astype(int)
        yticks = np.linspace(0,np.max(i_vec)+5,4).astype(int)
        ax.set_xticks(ticks)
        ax.set_yticks(yticks)




        for i in range(0,i_vec.shape[0]):

            ax.plot(i_vec[i],alpha=.6,color=self.main_color)

        #ax.plot(i_vec[0],color=self.main_color,linewidth=3)
        ax.set_xlabel('time (sec)')
        ax.set_ylabel('Intensity (a.u.)')


    def plot_ssa_average(self,ax,i_vec):
        ax.cla()
        ax.axis([0,int(i_vec.shape[1]),0,np.max(i_vec)+5] )
        ticks = np.linspace(0,int(i_vec.shape[1]) ,6).astype(int)
        yticks = np.linspace(0,np.max(i_vec)+5,4).astype(int)
        ax.set_xticks(ticks)
        ax.set_yticks(yticks)




        for i in range(0,i_vec.shape[0]):

            ax.plot(i_vec[i],alpha=.1,color='gray')

        ax.plot(np.mean(i_vec,axis=0),color=self.main_color,linewidth=3)

        ax.plot(np.mean(i_vec,axis=0) - np.std(i_vec,axis=0  ), color=self.main_color,linewidth=1)
        ax.plot(np.mean(i_vec,axis=0) + np.std(i_vec,axis=0  ), color=self.main_color,linewidth=1)
        ax.set_xlabel('time (sec)')
        ax.set_ylabel('Intensity (a.u.)')


    def __is_number(self, string):
        '''
        Is this string a number?
        '''
        try:
            complex(string)
        except ValueError:
            return False
        return True
    
    
    def __get_closest_val(self,value,array):
       
        idx = (np.abs(array - value)).argmin()
        return idx


    def plot_ssa_pdf(self,ax,i_vec):
        ax.cla()
        nbins = int(np.max(i_vec))

        stime = self.ssa.time_rec-self.ssa.start_time
        
        timeslice = self.pdf_time_entry.get()
        
        if self.__is_number(timeslice):
            idx = self.__get_closest_val(int(timeslice),stime)
            ivec_slice = i_vec.T[idx]
        else:
            ivec_slice = i_vec

        argmin = int(np.min(np.floor(ivec_slice/len(self.ssa.probe))))
        argmax = int(np.max(np.ceil(ivec_slice/len(self.ssa.probe))))

        histdata = np.histogram(ivec_slice/len(self.ssa.probe), bins = np.linspace(argmin,argmax, argmax+1).astype(int))

   
     
        '''
        ax.axis([0,np.max(i_vec)+5,0,np.max(hist[0])+5] )
        ticks = np.linspace(0,np.max(i_vec)+5 ,6).astype(int)
        yticks = np.linspace(0,np.max(hist[0])+5,4).astype(int)
        ax.set_xticks(ticks)
        ax.set_yticks(yticks)
        '''


        ax.bar(histdata[1][:-1],histdata[0],color=self.main_color,width=.96)
        
        maxhist = np.sum(histdata[0])
        ticks = ax.get_yticks()
        p_ticks = []
        for tick in ticks:
            val = float(tick)
            p_ticks.append(str(np.round(val/maxhist,2)))
            
        
        ax.set_yticklabels(p_ticks)
        
        '''
        for i in range(0,i_vec.shape[1])
            nbins = int(np.max(ivect))
            ax.hist(ivect[i],bins = nbins,alpha=.5,color=self.main_color,histtype='step',fill=self.main_color)
        '''

        ax.set_xlabel('intensity (ump)')
        ax.set_ylabel('Probability ')


    def plot_ssa_acc_codon(self,ax,mean_acc,error_acc,color=None,name=None):


        if color == None:
            color = self.main_color

        ax.axis([0,int(len(mean_acc)),-.5,1.1] )
        t = np.linspace(0,len(mean_acc),len(mean_acc)+1)
        ticks = np.linspace(0,int(len(mean_acc)),6).astype(int)
        yticks = np.linspace(-.5,1.1,4).astype(int)
        ax.set_xticks(ticks)
        ax.set_yticks(yticks)
        ax.set_xlabel('time (sec)')
        ax.plot(mean_acc,color=color)
        ax.plot(mean_acc-error_acc,'--',color=color,alpha=.3)
        ax.plot(mean_acc+error_acc,'--',color=color,alpha=.3)
        ax.plot([0,int(len(mean_acc))],[0,0],color='r',alpha=.5)


        p = []
        i = 0
        label = ['Wild-Type','Modified']
        k=0
        for line in ax.lines:
            if i%4 == 0:

                p.append(mpatches.Patch(color=line.get_c(), label=label[k]))
                k+=1
            i+=1


        ax.legend(handles=p)
        #ax.text(int(2*len(t)/3),-.4,'time (sec)')
        
    def kymograph(self):
        
        self.kymax.clear()
        self.kymax2.clear()
        n_traj = int(self.traj_select.get())
        self.kymograph_internal(self.kymax,self.kymax2,self.ssa,n_traj,color=self.main_color,bg_intense=self.kymintense.get())
        
        self.kymax3.clear()
        self.plot_rib_dense_kym(self.kymax3,self.ssa.rib_density)

        self.kym_canvas.draw()
        self.kym_canvas2.draw()
        self.kym_canvas3.draw()
        
        
    def kymograph_internal(self,ax,ax1,ssa_obj,n_traj,bg_intense=True,show_intense = True, *args,**kwargs):
        '''
        Constructs a kymograph of ribosome locations
        '''
        
        

        startfrags = 0
        for i in range(n_traj):
            startfrags += ssa_obj.frag_per_traj[i]
            
        endfrags = startfrags + ssa_obj.frag_per_traj[n_traj]
        fragments = ssa_obj.fragments[startfrags:endfrags]
        
        nfrag = fragments.shape[0]
 
        maxlen= fragments.shape[1]
        time  = ssa_obj.time
        
        ivec = ssa_obj.intensity_vec[n_traj]
        
        ftimes = ssa_obj.fragtimes[startfrags:startfrags+endfrags]
        #plt.figure(figsize=(5,10))

        #plt.subplot(gs[0])
        lenplot = np.max(fragments)
        maxin = np.max(ivec)
        
        #ax.set_facecolor('black')
        
        
        cmapguess = self.kym_menu.get()
        if cmapguess not in cm.datad.keys():
            cmap = cm.summer
        else:
            exec(('cmap = cm.' + cmapguess))
            
        bgcolor = self.kymbgcolor.cget('bg')
        linecolor = self.kymlinecolor.cget('bg')
            
        if bg_intense == True:
            for i in range(len(ssa_obj.time_rec)):
                ax.plot([0,lenplot],[ssa_obj.time_rec[i],ssa_obj.time_rec[i]],color = cmap(1.*ivec[i]/maxin),lw=1)
            
        for i in range(nfrag):
            
           
            
            if maxlen <= np.where(fragments[i] > 0 )[0][-1]:       
                timeseg = time[ftimes[i]:ftimes[i]+maxlen]
                ax.plot(fragments[i][0:len(timeseg)] ,timeseg[::-1],color=linecolor )
                
            else:
                timeseg = time[ftimes[i]:]
                stop = np.where(fragments[i] > 0 )[0][-1]+1
                timelen = len(fragments[i][0:stop]) 
                
                ax.plot(fragments[i][0:stop]   ,timeseg[0:timelen],color=linecolor)

        segtime = ssa_obj.time[0:len(ssa_obj.time_rec)]
        

        ax.set_xlabel('Ribosome position (residue)')
        ax.set_ylabel('Time')
        ax.set_ylim(ssa_obj.time_rec[-1], ssa_obj.time_rec[0])
        ax.set_xlim(-10,max(ssa_obj.rib_density.shape)+10 )
                

            
        
        ax.set_facecolor(bgcolor )
    
        ax1.plot(ivec.T,segtime,color=linecolor)
        ax1.set_xlabel('Intensity (AU)')
        ax1.set_facecolor(bgcolor)
        ax1.set_ylim(segtime[-1], segtime[0])
        ax1.set_yticks([])
        ax1.set_xlim(0,maxin+5)
       
            


    def plot_ssa_acc(self,ax,mean_acc,error_acc,color=None,name=None):


        if color == None:
            color = self.main_color

        ax.axis([0,int(len(mean_acc)),-.5,1.1] )
        t = np.linspace(0,len(mean_acc),len(mean_acc)+1)
        ticks = np.linspace(0,int(len(mean_acc)),6).astype(int)
        yticks = np.linspace(-.5,1.1,4).astype(int)
        ax.set_xticks(ticks)
        ax.set_yticks(yticks)
        ax.set_xlabel('time (sec)')
        ax.plot(mean_acc,color=color)
        ax.plot(mean_acc-error_acc,'--',color=color,alpha=.3)
        ax.plot(mean_acc+error_acc,'--',color=color,alpha=.3)
        ax.plot([0,int(len(mean_acc))],[0,0],color='r',alpha=.5)



    def plot_rib_dense_kym(self,ax,rib_density):
        ax.cla()
        ax.axis([0,int(len(rib_density)),0,np.max(rib_density)+.1*np.max(rib_density)])
        ax.set_xlabel('Codon Position')
        ax.set_ylabel('P')
        
        ax.plot(rib_density,color=self.kymlinecolor.cget('bg'))
        bgcolor = self.kymbgcolor.cget('bg')
        ax.set_facecolor(bgcolor )
        
        ax.set_ylim(np.max(rib_density)+.1*np.max(rib_density),np.min(rib_density)+.1*np.min(rib_density) )
        ax.set_xlim(-10,rib_density.shape[0] + 10)

    def plot_rib_dense(self,ax,rib_density):
        ax.cla()
        ax.axis([0,int(len(rib_density)),0,np.max(rib_density)+.1*np.max(rib_density)])
        ax.set_xlabel('Codon Position')
        ax.set_ylabel('P')
        ax.plot(rib_density,self.main_color)

    def plot_cai_codons(self,ax,cai_codon):
        ax.cla()
        ax.axis([0,int(len(cai_codon)),0,np.max(cai_codon)+.1*np.max(cai_codon)])
        ax.set_xlabel('Codon Position')
        ax.set_ylabel('CAI')
        ax.plot(cai_codon,color=self.main_color)




    def save_ss_data(self,num):
        try:
            filename = tfd.asksaveasfilename(title = "Save as",defaultextension='.txt',filetypes = (("text files","*.txt"),("csv","*.csv"),('JSON','*.json'),('XML','*.xml'),("all files","*.*")))
        except:
            pass
        self.ssa = self.ssas[num]
        f = open(filename,'w')
        if '.txt' in filename:
            for key in self.ssa.__dict__.keys():
                print(key)
                #if key != 'rib_vec' and key != 'ribosome_means':

                if key == 'solutions':
                    f.write((key + '\r\n'))
                    for i in range(len(self.ssa.solutions)):

                        small_solution = self.ssa.solutions[0][np.where(np.sum(self.ssa.solutions[0],axis=1)!=0)]
                        np.savetxt(f,np.atleast_2d(small_solution),delimiter=',',fmt='%s')
                        f.write(('\r\n'))

                else:
                    f.write((key + '\r\n'))
                    np.savetxt(f,np.atleast_2d(self.ssa.__dict__[key]),delimiter=',',fmt='%s')
                    f.write(('\r\n'))

        if '.json' in filename:

            ssadict = {}
            for key in self.ssa.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa.__dict__[key]

            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)
        f.close()


    def load_ss_data(self):
        try:
            filename = tfd.askopenfilename(title = "Open",filetypes = (("all files","*.*"),("text files","*.txt"),("csv","*.csv"),('JSON','*.json'),('XML','*.xml')))
        except:
            pass

        if '.txt' in filename:


            self.ssa = ssa()
            self.ssa.load_from_txt(filename)


            self.plot_loaded_ssa()


        if '.json' in filename:
            self.ssa = ssa()
            obj_text = codecs.open(filename, 'r', encoding='utf-8').read()
            ssadict = json.loads(obj_text)

            for key in ssadict.keys():
                if key in ['rib_means','time_vec_fixed','intensity_vec','mean_autocorr','autocorr_vec','error_autocorr','rib_density']:

                    self.ssa.__dict__[key] = np.array(ssadict[key])
                else:
                    self.ssa.__dict__[key] = ssadict[key]

            self.plot_loaded_ssa()

    def save_frap_data(self):
        try:
            filename = tfd.asksaveasfilename(title = "Save as",defaultextension='.txt',filetypes = (("text files","*.txt"),("csv","*.csv"),('JSON','*.json'),('XML','*.xml'),("all files","*.*")))
        except:
            pass

        f = open(filename,'a')
        if '.txt' in filename:
            for key in self.ssa_frap.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    f.write((key + '\r\n'))
                    np.savetxt(f,np.atleast_2d(self.ssa_frap.__dict__[key]),delimiter=',',fmt='%s')
                    f.write(('\r\n'))


        if '.json' in filename:

            ssadict = {}
            for key in self.ssa_frap.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa_frap.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa_frap.__dict__[key]

            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

        f.close()

    def save_harr_data(self):
        try:
            filename = tfd.asksaveasfilename(title = "Save as",defaultextension='.txt',filetypes = (("text files","*.txt"),("csv","*.csv"),('JSON','*.json'),('XML','*.xml'),("all files","*.*")))
        except:
            pass

        f = open(filename,'a')
        if '.txt' in filename:
            for key in self.ssa_harr.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    f.write((key + '\r\n'))
                    np.savetxt(f,np.atleast_2d(self.ssa_harr.__dict__[key]),delimiter=',',fmt='%s')
                    f.write(('\r\n'))

        if '.json' in filename:

            ssadict = {}
            for key in self.ssa_harr.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa_harr.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa_harr.__dict__[key]

            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

        f.close()

    def load_frap_data(self):
        try:
            filename = tfd.askopenfilename(title = "Open",filetypes = (("all files","*.*"),("text files","*.txt"),("csv","*.csv"),('JSON','*.json'),('XML','*.xml')))
        except:
            pass

        if '.txt' in filename:

            self.ssa_frap = ssa()
            self.ssa_frap.load_from_txt(filename)
            self.plot_loaded_frap()

        if '.json' in filename:
            self.ssa_frap = ssa()
            obj_text = codecs.open(filename, 'r', encoding='utf-8').read()
            ssadict = json.loads(obj_text)

            for key in ssadict.keys():
                if key in ['rib_means','time_vec_fixed','intensity_vec','mean_autocorr','autocorr_vec','error_autocorr','rib_density']:

                    self.ssa_frap.__dict__[key] = np.array(ssadict[key])
                else:
                    self.ssa_frap.__dict__[key] = ssadict[key]

            self.plot_loaded_frap()


    def load_harr_data(self):
        try:
            filename = tfd.askopenfilename(title = "Open",filetypes = (("all files","*.*"),("text files","*.txt"),("csv","*.csv"),('JSON','*.json'),('XML','*.xml')))

        except:
            pass
        if '.txt' in filename:

            self.ssa_harr = ssa()
            self.ssa_harr.load_from_txt(filename)
            self.plot_loaded_harr()

        if '.json' in filename:
            self.ssa_harr = ssa()
            obj_text = codecs.open(filename, 'r', encoding='utf-8').read()
            ssadict = json.loads(obj_text)

            for key in ssadict.keys():
                if key in ['rib_means','time_vec_fixed','intensity_vec','mean_autocorr','autocorr_vec','error_autocorr','rib_density']:

                    self.ssa_harr.__dict__[key] = np.array(ssadict[key])
                else:
                    self.ssa_harr.__dict__[key] = ssadict[key]

            self.plot_loaded_harr()


    def save_natural_codon_data(self):
        try:
            filename = tfd.asksaveasfilename(title = "Save as",defaultextension='.txt',filetypes = (("text files","*.txt"),("csv","*.csv"),('JSON','*.json'),('XML','*.xml'),("all files","*.*")))
        except:
            pass

        f = open(filename,'a')
        if '.txt' in filename:
            for key in self.ssa_norm.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    f.write((key + '\r\n'))
                    np.savetxt(f,np.atleast_2d(self.ssa_norm.__dict__[key]),delimiter=',',fmt='%s')
                    f.write(('\r\n'))

        if '.json' in filename:

            ssadict = {}
            for key in self.ssa_norm.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa_norm.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa_norm.__dict__[key]

            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

        f.close()

    def save_design_codon_data(self):
        try:
            filename = tfd.asksaveasfilename(title = "Save as",defaultextension='.txt',filetypes = (("text files","*.txt"),("csv","*.csv"),('JSON','*.json'),('XML','*.xml'),("all files","*.*")))
        except:
            pass

        f = open(filename,'a')
        if '.txt' in filename:
            for key in self.ssa_design.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    f.write((key + '\r\n'))
                    np.savetxt(f,np.atleast_2d(self.ssa_design.__dict__[key]),delimiter=',',fmt='%s')
                    f.write(('\r\n'))
        if '.json' in filename:

            ssadict = {}
            for key in self.ssa_design.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa_design.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa_design.__dict__[key]

            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

        f.close()





    def save_elon_data(self):
        x=1

    def save_codon_data(self):
        x=1

    def set_natural_codons(self):
        for var in self.aa_vars:
            var.set(0)

    def set_common_codons(self):
        for var in self.aa_vars:
            var.set(1)

    def set_rare_codons(self):
        for var in self.aa_vars:
            var.set(2)


    def plot_sequence(self,ax,gene_length,tag_length,total_length,epitopes_pos,tag_name,gene_name):
        ax.cla()
        ax.plot([0,total_length],[0,0],color='k')
        ax.plot([0,total_length],[1,1],color='k')
        ax.plot([0,0],[0,1],color='k')
        ax.plot([total_length,total_length],[0,1],color='k')
        ax.axis([-10,total_length+10,-1,2])
        ax.plot([tag_length,tag_length],[0,1],color='k',linewidth=1)

        for i in range(len(epitopes_pos)):
            ax.plot([epitopes_pos[i],epitopes_pos[i]],[0,1],color='orange',linewidth=2)
        ax.fill_between([tag_length,tag_length,total_length,total_length],[0,1,1,0],color='#00FF00')

        ticks = np.linspace(0,total_length,10).astype(int)
        ax.set_xticks(ticks)
        ax.set_xlabel('Codon Position')
        ax.get_yaxis().set_visible(False)
        ax.text(tag_length/3,1.3,('Tag: ' + tag_name), fontdict = dict(fontsize=13) )
        ax.text(tag_length,1.3,('POI: ' + gene_name),fontdict = dict(fontsize=13))



    def popout_fig_ssa(self):
        '''
        newwin = tk.Toplevel(self.parent)
        newwin.lift()
        img = tk.Image("photo", file=os.path.join("icons","icon.gif"))
        newwin.tk.call('wm','iconphoto',newwin._w,img)
        newwin.title(' ')
        newframe = tk.Frame(newwin)
        newframe.grid(row=0,column=0)
        '''



        fig = mpl.figure.Figure(figsize=(14,5))#figsize=(12,5),dpi=self.global_dpi)
        fig.set_tight_layout(True)
        ax = fig.add_subplot(111)

        fig = plt.figure(figsize=(14,5))
        ax = fig.add_subplot('111')
        self.plot_ssa_intensity(ax,self.ssa.intensity_vec)

        fig.show()

        '''
        vis_canvas = FigureCanvasTkAgg(fig,master=newframe)  #tkinter specific matplotlib canvas
        vis_canvas.draw()
        vis_canvas.get_tk_widget().pack(fill='both',side='top',expand=True) #place the figure widget
        #protein_fig.patch.set_facecolor(self.default_color)  #make the edges of the figure the same as the background so it looks like it blends in

        '''


    def popout_fig_acc(self):
        fig = plt.figure(figsize=(14,5))
        ax = fig.add_subplot('111')


        self.plot_ssa_acc(ax,self.ssa.mean_autocorr,self.ssa.error_autocorr)


    def popout_fig_caicodon(self):
        fig = plt.figure(figsize=(14,5))
        ax = fig.add_subplot('111')


        self.plot_cai_codons(ax,self.sms.POI.CAI_codons)


    def popout_fig_ribdense(self):
        fig = plt.figure(figsize=(14,5))
        ax = fig.add_subplot('111')


        self.plot_rib_dense(ax,self.ssa.rib_density)


    def popout_fig_harr(self):
        fig = plt.figure(figsize=(14,5))
        ax = fig.add_subplot('111')
        self.plot_ssa_frap(ax,self.ssa_harr.intensity_vec)


    def popout_fig_codoncorr(self):
        fig = plt.figure(figsize=(14,5))
        ax = fig.add_subplot('111')

        index1 = np.where(self.ssa_norm.mean_autocorr < 0.001)[0][0]
        index2 = np.where(self.ssa_design.mean_autocorr<.001)[0][0]

        if index1 >= index2:
            index = index1
        if index2 >= index1:
            index = index2

        self.plot_ssa_acc_codon(ax,self.ssa_norm.mean_autocorr[:index],self.ssa_norm.error_autocorr[:index],name='Wild Type')

        self.plot_ssa_acc_codon(ax,self.ssa_design.mean_autocorr[:index],self.ssa_design.error_autocorr[:index],color='#3a7ce8',name='Modified')



    def popout_fig_frap(self):
        fig = plt.figure(figsize=(14,5))
        ax = fig.add_subplot('111')
        self.plot_ssa_frap(ax,self.ssa_frap.intensity_vec)

    def popout_fig_seq(self):
        '''
        Pull out all data from sequence viewer and replot it in a new pop out figure
        '''

        fig = plt.figure(figsize=(14,5))
        ax = fig.add_subplot('111')

        chiles = self.ax.get_children()
        for i in range(len(chiles)):
            if  isinstance(chiles[i],matplotlib.text.Text):


                if 'POI:' in chiles[i]._text:
                    gene_name = chiles[i]._text.replace('POI: ', '')

                if 'Tag:' in chiles[i]._text:
                    tag_name = chiles[i]._text.replace('Tag: ', '')

        ep_pos = []
        for i in range(len(self.ax.lines)):
            if i == 3:
                total_length=self.ax.lines[i].get_xdata()[0]
            if i == 4:
                tag_length=self.ax.lines[i].get_xdata()[0]
            if i > 4:
                ep_pos.append(self.ax.lines[i].get_xdata()[0])

        gene_length = total_length-tag_length
        self.plot_sequence(ax,gene_length,tag_length,total_length,ep_pos,tag_name,gene_name)

    def update_ss_frame_data(self):
        calc = [1,4,7,10,13]

        
        params = [str(np.round(self.sms.POI.CAI,3)),str(self.ssa.ke_sim),str(int(np.round(self.sms.POI.total_length/self.ssa.no_rib_per_mrna,0))),str(np.round(self.sms.POI.total_length/self.ssa.ke_sim,2))  ,str(np.round(self.ssa.no_rib_per_mrna,3))]
        
        i = 0
        for index in calc:
            self.calculated_parameters[index].config(state='normal')
            self.calculated_parameters[index].delete(0,tk.END)
            self.calculated_parameters[index].insert(tk.END,params[i])
            self.calculated_parameters[index].config(state='readonly')
            i+=1
        self.calculated_parameters[16].config(state='normal')
        self.calculated_parameters[16].delete(0,tk.END)
        self.calculated_parameters[16].insert(tk.END,str(np.round(self.ssa.ke_true,2)))
        self.calculated_parameters[16].config(state='readonly')

    


        if self.ssa.evaluating_inhibitor == True:
            meaniv = np.mean(self.ssa.intensity_vec,axis=0)
            meaniv_norm = meaniv/np.max(meaniv)
            norm_iv = self.ssa.intensity_vec/np.max(self.ssa.intensity_vec)
            sem_iv = np.std(self.ssa.intensity_vec,axis=0) / float(self.ssa.n_traj)
            ti = self.ssa.time_inhibit - self.ssa.start_time
            time_inds = np.where(self.ssa.time_vec_fixed  > 0)
            tvec = self.ssa.time_vec_fixed[time_inds]

            timeindex = np.where(tvec >= ti)[0][0]

            runoffindex = np.where(meaniv_norm[timeindex:] <.005   )[0][0] 
            
  

            application_time = float(self.ssa_inputs[6].get())


            runoff_time = (tvec)[runoffindex]

            elongation_rate = np.round(float(self.sms.POI.total_length) / runoff_time,2)

            self.calculated_parameters[22].config(state='normal')
            self.calculated_parameters[22].delete(0,tk.END)
            self.calculated_parameters[22].insert(tk.END,'---')
            self.calculated_parameters[22].config(state='readonly')

            self.calculated_parameters[25].config(state='normal')
            self.calculated_parameters[25].delete(0,tk.END)
            self.calculated_parameters[25].insert(tk.END,elongation_rate)
            self.calculated_parameters[25].config(state='readonly')

            self.calculated_parameters[19].config(state='normal')
            self.calculated_parameters[19].delete(0,tk.END)
            self.calculated_parameters[19].insert(tk.END,runoff_time)
            self.calculated_parameters[19].config(state='readonly')

        if self.ssa.evaluating_frap == True:

            meaniv = np.mean(self.ssa.intensity_vec,axis=0)
            meaniv_norm = meaniv/np.max(meaniv)
            norm_iv = self.ssa.intensity_vec/np.max(self.ssa.intensity_vec)
            sem_iv = np.std(self.ssa.intensity_vec,axis=0) / float(self.ssa.n_traj)
            ti = self.ssa.time_inhibit
            
            tvec = self.ssa.time_rec-self.ssa.start_time

            timeindex = np.where(tvec >= ti)[0][0]
            
            print(tvec)
            print(ti)
            print(np.where(meaniv_norm[timeindex:] <.005   ))
            print(timeindex)
            runoffindex = np.where(meaniv_norm[timeindex:] <.005   )[0][0] + timeindex




            #self.ssadat = [norm_iv,sem_iv,meaniv,meaniv_norm,runoff_time]
            #self.cpv_entries = [ere,rte,ere2,rte2]


            tempind = 50
            recoverythresh = .9
            recoverytime = 1

            application_time = float(self.ssa_inputs[0].get()) 
            print(application_time)


            tvec_inds = np.where(tvec>=application_time)

            meaniv = np.mean(self.ssa.intensity_vec,axis=0)
            meaniv_norm = meaniv/np.max(meaniv)




            try:
                #recoverytime = tvec[np.where(meaniv_norm[tvec_inds] >= recoverythresh )][0]+50
                timeindex = np.where(tvec >= ti+20)[0][0]
                runoffindex = np.where(meaniv_norm[timeindex:] >= recoverythresh  )[0][0] + timeindex


                recoverytime = (tvec)[runoffindex] - application_time



            except:
                #recoverytime = tvec[np.where(meaniv_norm[tvec_inds] >= .7 )][0]+50
                timeindex = np.where(tvec >= ti)[0][0]
                runoffindex = np.where(meaniv_norm[timeindex+1:] >= .7  )[0][0] + timeindex

                recoverytime = (tvec)[runoffindex] - application_time

            elongation_rate = np.round(float(self.sms.POI.total_length) / recoverytime,2)
            self.calculated_parameters[22].config(state='normal')
            self.calculated_parameters[22].delete(0,tk.END)
            self.calculated_parameters[22].insert(tk.END,recoverytime)
            self.calculated_parameters[22].config(state='readonly')

            self.calculated_parameters[25].config(state='normal')
            self.calculated_parameters[25].delete(0,tk.END)
            self.calculated_parameters[25].insert(tk.END,elongation_rate)
            self.calculated_parameters[25].config(state='readonly')

            self.calculated_parameters[19].config(state='normal')
            self.calculated_parameters[19].delete(0,tk.END)
            self.calculated_parameters[19].insert(tk.END,'---')
            self.calculated_parameters[19].config(state='readonly')



    def __update_ss_frame(self,event):
        '''
        If there is a protien of intrest update the dialog boxes on the the stochastic sim frame
        '''
        try:
            self.sms.POI.name
        except:
            return

        #self.ss_gs_inputs = [gbid_e,sfe,sste,repe]

        try:
            self.insert_entry(self.ss_gs_inputs[0],self.gb_rec.id)
        except:
            pass

        self.insert_entry(self.ss_gs_inputs[1],self.seq_name)
        if self.ss_gs_inputs[2].get() == '':
            self.ss_gs_inputs[2].insert(tk.END,'1000')

        if self.ss_gs_inputs[3].get() == '':
            self.ss_gs_inputs[3].insert(tk.END,'50')

        #self.ss_rates = [gl_e,eer_e,eir_e]

        self.insert_entry(self.ss_rates[0],str(self.sms.POI.total_length))

        if self.ss_rates[1].get() == '':
            self.ss_rates[1].insert(tk.END,'10')

        if self.ss_rates[2].get() == '':
            self.ss_rates[2].insert(tk.END,'0.03')


    def __update_elon_frame(self,event):
        '''
        If there is a protien of intrest update the dialog boxes on the the stochastic sim frame
        '''
        try:
            self.sms.POI.name
        except:
            return

        #self.ss_gs_inputs = [gbid_e,sfe,sste,repe]

        try:
            self.insert_entry(self.elon_gs_inputs[0],self.gb_rec.id)
        except:
            pass


        self.insert_entry(self.elon_gs_inputs[1],self.seq_name)
        if self.elon_gs_inputs[2].get() == '':
            self.elon_gs_inputs[2].insert(tk.END,'600')

        if self.elon_gs_inputs[3].get() == '':
            self.elon_gs_inputs[3].insert(tk.END,'200')

        if self.elon_gs_inputs[4].get() == '':
            self.elon_gs_inputs[4].insert(tk.END,'50')

        #self.ss_rates = [gl_e,eer_e,eir_e]

        self.insert_entry(self.ss_rates_elon[0],str(self.sms.POI.total_length))

        if self.ss_rates_elon[1].get() == '':
            self.ss_rates_elon[1].insert(tk.END,'10')

        if self.ss_rates_elon[2].get() == '':
            self.ss_rates_elon[2].insert(tk.END,'0.03')
        # Keep this handle alive, or else figure will disappear

    def __update_codon_frame(self,event):
        '''
        If there is a protien of intrest update the dialog boxes on the the stochastic sim frame
        '''
        try:
            self.sms.POI.name
        except:
            return

        try:
            self.insert_entry(self.co_gs_inputs[0],self.gb_rec.id)
        except:
            pass


        #self.ss_gs_inputs = [gbid_e,sfe,sste,repe]

        self.insert_entry(self.co_gs_inputs[1],self.seq_name)
        if self.co_gs_inputs[2].get() == '':
            self.co_gs_inputs[2].insert(tk.END,'1000')

        if self.co_gs_inputs[3].get() == '':
            self.co_gs_inputs[3].insert(tk.END,'50')

        #self.ss_rates = [gl_e,eer_e,eir_e]

        self.insert_entry(self.ss_rates_codon[0],str(self.sms.POI.total_length))

        if self.ss_rates_codon[1].get() == '':
            self.ss_rates_codon[1].insert(tk.END,'10')

        if self.ss_rates_codon[2].get() == '':
            self.ss_rates_codon[2].insert(tk.END,'0.03')
        # Keep this handle alive, or else figure will disappear


        cmnorm = self.sms.POI.codon_sensitivity/np.max(self.sms.POI.codon_sensitivity).flatten().tolist()
        #colors = self.__cmap_map(lambda x: .9*x +.2,cm.Wistia)(cmnorm[0])
        colors = cm.coolwarm(cmnorm[0])


        for i in range(len(self.aa_entries)):

            self.aa_entries[i].config(state='normal')
            self.aa_entries[i].delete(0,tk.END)
            self.aa_entries[i].insert(tk.END,str(self.sms.POI.aa_seq.count(self.sms.aa_keys[i])))
            self.aa_entries[i].config(state='readonly')
            self.aa_entries[i].config(readonlybackground='#FFFFFF')


            self.sensitivity_entries[i].config(state='normal')
            self.sensitivity_entries[i].delete(0,tk.END)
            self.sensitivity_entries[i].insert(tk.END,str(self.sms.POI.codon_sensitivity[0,i]))
            self.sensitivity_entries[i].config(state='readonly')
            color = (colors[i][0:3])*255


            dc=[0,0,0]

            dc[0] = hex(int(color[0]))[2:]
            dc[1] = hex(int(color[1]))[2:]
            dc[2] = hex(int(color[2]))[2:]
            if len(dc[0]) < 2:
                dc[0] = '0'+dc[0]
            if len(dc[1]) < 2:
                dc[1] = '0'+dc[1]
            if len(dc[2]) < 2:
                dc[2] = '0'+dc[2]
            color = '#' + dc[0] + dc[1] + dc[2]

            self.sensitivity_entries[i].config(readonlybackground=color)



class ssa():
    '''
    SSA container class

    holds intensity / ribosome data as well as the propensities used
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




    def save_txt(self,filename):

        if '.txt' in filename:
            f = open(filename,'a')
            for key in self.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    f.write((key + '\r\n'))
                    np.savetxt(f,np.atleast_2d(self.__dict__[key]),delimiter=',',fmt='%s')
                    f.write(('\r\n'))

        else:
            filename = filename + '.txt'
            f = open(filename,'a')
            for key in self.__dict__.keys():

                if key != 'rib_vec' and key != 'ribosome_means':
                    f.write((key + '\r\n'))
                    np.savetxt(f,np.atleast_2d(self.__dict__[key]),delimiter=',',fmt='%s')
                    f.write(('\r\n'))
        f.close()

    def load_from_txt(self,filename):
        if '.txt' in filename:
            ssa_obj = np.loadtxt(filename,dtype=str)


            for i in range(   0,len(ssa_obj)-1  ):
                label = ssa_obj[i]
                if label in ['rib_means','rib_vec','n_traj','start_time','k','time_vec_fixed','dwelltime','intensity_vec','mean_autocorr','no_rib_per_mrna','ke_sim','autocorr_vec','ribosome_means','error_autocorr','rib_density']:
                    if label == 'intensity_vec':

                        tvec = self.time_vec_fixed[np.where(self.time_vec_fixed >= self.start_time)]
                        i_vec = np.zeros((self.n_traj,len(tvec)))

                        for j in range(self.n_traj):
                            array = np.fromstring(ssa_obj[i+j+1],dtype=float,sep=',')
                            i_vec[j] = array


                        exec(('self.'+label+ '=i_vec'))

                    elif label in ['start_time','no_rib_per_mrna','ke_sim','dwelltime']:
                        array = np.fromstring(ssa_obj[i+1],dtype=float,sep=',')[0]
                        exec(('self.'+label+ '=array'))
                    elif label in ['n_traj']:
                        array = int(np.fromstring(ssa_obj[i+1],dtype=float,sep=',')[0])
                        exec(('self.'+label+ '=array'))
                    else:
                        array = np.fromstring(ssa_obj[i+1],dtype=float,sep=',')
                        exec(('self.'+label+ '=array'))



    def save_from_json(self,filename):

        if '.json' in filename:

            ssadict = {}
            for key in self.ssa_harr.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa_harr.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa_harr.__dict__[key]

            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

        else:
            filename =  filename + '.json'

            ssadict = {}
            for key in self.ssa_harr.__dict__.keys():
                if key != 'rib_vec' and key != 'ribosome_means':
                    try:
                        ssadict[key] = self.ssa_harr.__dict__[key].tolist()
                    except:
                        ssadict[key] = self.ssa_harr.__dict__[key]

            json.dump(ssadict, codecs.open(filename, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)


    def load_json(self,filename):
        if '.json' in filename:

            obj_text = codecs.open(filename, 'r', encoding='utf-8').read()
            ssadict = json.loads(obj_text)

            for key in ssadict.keys():
                if key in ['rib_means','time_vec_fixed','intensity_vec','mean_autocorr','autocorr_vec','error_autocorr','rib_density']:

                    self.__dict__[key] = np.array(ssadict[key])
                else:
                    self.__dict__[key] = ssadict[key]


class GenericSSA():
    def __init__(self,type='linear'):

        self.time_variant = False
        self.xi=np.array([])
        self.ti= None
        self.tf=None
        self.S=np.array([])
        self.type=type
        self.ptimes=100
        self.params={}
        if type=='linear':
            #self.fast_rxn = 0.5
            self.W0=np.array([])
            self.W1=np.array([])
        if type == 'nonlinear':
            #self.fast_rxn = 0.5
            self.P=lambda x,t:None



    def gettvec(self):
        return np.linspace(self.ti,self.tf,self.ptimes)

    def updateparams(self):
        for k, v in self.params.iteritems():
            setattr(self, k, v)

    def _run_trajectory(self):#!!!!!!!!!!!!!!!!!renamed run to solve(big deal)
        x=self.xi
        t=self.ti
        __n=len(x)
        self.time=self.gettvec()
        data=np.zeros((len(self.xi),self.ptimes))
        ip=0

        if self.type=='linear':
            if self.time_variant == False:
                while t<self.tf:
                    rate=np.atleast_2d(np.dot(self.W1,x))+self.W0
                    rate=np.cumsum(rate)
                    with np.errstate(divide='ignore', invalid='ignore'):
                        t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()
                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            data[:,ip]=x.reshape(__n)
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i
                            break
                    x=x+np.atleast_2d(self.S[:,event]).T

            else:


                x = np.concatenate((np.array([0]),self.xi.flatten()))

                #x = self.xi
                __n=len(x)
                self.time=self.gettvec()
                data=np.zeros((len(self.xi),self.ptimes))
                a,b = self.S.shape
                S = np.vstack((np.zeros(b),self.S))
                S = np.hstack((np.zeros((a+1,1)),S))
                while t<self.tf:
                    __n=len(x)
                    self.time=self.gettvec()
                    data=np.zeros((len(self.xi),self.ptimes))
                    a,b = self.S.shape
                    S = np.vstack((np.zeros(b),self.S))
                    S = np.hstack((np.zeros((a+1,1)),S))
                    while t<self.tf:
                        trate=self.get_P(x[1:],t)

                        rate = np.concatenate((np.array([self.fast_rxn]),trate))
                        rate=np.cumsum(rate)


                        t=(t-np.log(np.random.rand(1))/rate[-1])
                        ro=rate[-1]*np.random.rand()

                        while t>self.time[ip]:
                            if t>self.tf:
                                b = len(self.time[ip:])
                                fill = np.repeat(x[1:],b)
                                data[:,ip:]=fill.reshape(__n-1,b)
                                return data
                            else:
                                #data[:,ip]=x.reshape(__n)
                                data[:,ip]=x[1:]
                                ip=ip+1
                        for i in range(len(rate)):
                            if rate[i]>=ro:
                                event=i

                                break

                        x=x+S[:,event].ravel()
                    '''
                    rate=np.atleast_2d(np.dot(self.W1(t),x))+self.W0(t)

                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    print(t)
                    ro=rate[-1]*np.random.rand()
                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            data[:,ip]=x.reshape(__n)
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i
                            break
                    x=x+np.atleast_2d(self.S[:,event]).T

                    '''



        elif self.type=='nonlinear':
            if self.time_variant == True:  #if time variant use fast reaction
                x = np.concatenate((np.array([0]),self.xi.flatten()))

                #x = self.xi
                __n=len(x)
                self.time=self.gettvec()
                data=np.zeros((len(self.xi),self.ptimes))
                a,b = self.S.shape
                S = np.vstack((np.zeros(b),self.S))
                S = np.hstack((np.zeros((a+1,1)),S))
                while t<self.tf:
                    trate=self.get_P(x[1:],t)
                    rate = np.concatenate((np.array([self.fast_rxn]),trate))
                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()

                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x[1:],b)
                            data[:,ip:]=fill.reshape(__n-1,b)
                            return data
                        else:
                            #data[:,ip]=x.reshape(__n)
                            data[:,ip]=x[1:]
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i

                            break

                    x=x+S[:,event].ravel()

            else:   #if not time variant ignore fast reaction

                x = self.xi.flatten()

                #x = self.xi
                __n=len(x)
                self.time=self.gettvec()

                while t<self.tf:
                    rate=self.get_P(x,t)

                    rate=np.cumsum(rate)

                    t=(t-np.log(np.random.rand(1))/rate[-1])
                    ro=rate[-1]*np.random.rand()

                    while t>self.time[ip]:
                        if t>self.tf:
                            b = len(self.time[ip:])
                            fill = np.repeat(x,b)
                            data[:,ip:]=fill.reshape(__n,b)
                            return data
                        else:
                            #data[:,ip]=x.reshape(__n)
                            data[:,ip]=x
                            ip=ip+1
                    for i in range(len(rate)):
                        if rate[i]>=ro:
                            event=i

                            break

                    x=x+self.S[:,event].ravel()



        else:
            'Error'
        self.data=data
        return data

    def _solve(self,n):
        __data=np.zeros((len(self.xi),self.ptimes,n))

        for i in range(n):
            __d=self._run_trajectory()
            __data[:,:,i]=__d
        self.data = __data

        return __data

    def setpar(self,key,val):
        self.params[key]=val

    def get_dist(self,specID=0):
        '''
        build distribution (non-normalized and pdf)
        of rna for the model)
        '''
        n_specs, n_times, n_traj = self.data.shape
        n_specs = int(n_specs)
        n_times = int(n_times)
        n_traj = int(n_traj)
        specID = int(specID)
        max_rna = int(np.max(self.data[specID,:,:]))
        self.pdf = np.zeros((n_times,max_rna+1))
        self.fdist = np.zeros((n_times,max_rna+1))
        for i in range(n_times):
            ind = int(i)
            for j in range(n_traj):
                jnd = int(j)
                self.fdist[ind,int(self.data[int(specID),ind,jnd])] +=1
            self.pdf[ind,:] = self.fdist[ind,:] / np.sum(self.fdist[ind,:])

    def get_traj(self,specID=0,ntraj='all'):
        n_specs, n_times, n_traj = self.data.shape
        n_specs = int(n_specs)
        n_times = int(n_times)




        if isinstance(specID,int):
            if ntraj == 'all':
                return self.data[specID],ntraj
            else:
                try:
                    ntraj = ntraj.flatten().astype(int).tolist()
                except:
                    ntraj = int(ntraj)
                    pass

                return self.data[specID][:,ntraj],ntraj
        else:
            if specID == 'all':
                if ntraj == 'all':
                    return self.data,ntraj
                else:

                    try:
                        ntraj = ntraj.flatten().astype(int).tolist()
                    except:
                        pass

                    return self.data,ntraj

    def get_means(self,specID=0):
        '''
        get the first moment.
        '''
        n_specs, n_times, n_traj = self.data.shape
        max_rna = int(np.max(self.data[int(specID),:,:])+1)
        self.means = np.zeros(n_times)
        for i in range(n_times):
            self.means[i] = np.sum(np.arange(max_rna)*self.pdf[i,:])

    def get_variances(self,specID=0):
        '''
        get the second moment.
        '''
        self.get_dist()
        n_specs, n_times, n_traj = self.data.shape
        max_rna = int(np.max(self.data[specID,:,:])+1)
        self.variances = np.zeros(n_times)
        self.covariances = np.zeros((n_specs,n_specs,n_times))
        for i in range(n_times):
            self.variances[i] = np.sum((np.arange(max_rna)**2)*self.pdf[i,:])-(np.sum(np.arange(max_rna)*self.pdf[i,:])**2)
            self.covariances[:,:,i] = np.cov(self.data[:,i,:])

    def return_all_var(self):
        all_members = self.__dict__.keys()

        return [ (item, self.__dict__[item]) for item in all_members if not item.startswith("_")]

    def return_names(self):
        all_members = self.__dict__.keys()
        return [ item for item in all_members if not item.startswith("_")]



def main():

    root = tk.Tk() #create GUI app

    root.geometry("+10+10") #define GUI shape
    icondata = '''
                R0lGODlhIwAhAHAAACH5BAEAAPwALAAAAAAjACEAhwAAAAAAMwAAZgAAmQAAzAAA/
                wArAAArMwArZgArmQArzAAr/wBVAABVMwBVZgBVmQBVzABV/wCAAACAMwCAZgCAmQ
                CAzACA/wCqAACqMwCqZgCqmQCqzACq/wDVAADVMwDVZgDVmQDVzADV/wD/AAD/MwD
                /ZgD/mQD/zAD//zMAADMAMzMAZjMAmTMAzDMA/zMrADMrMzMrZjMrmTMrzDMr/zNV
                ADNVMzNVZjNVmTNVzDNV/zOAADOAMzOAZjOAmTOAzDOA/zOqADOqMzOqZjOqmTOqz
                DOq/zPVADPVMzPVZjPVmTPVzDPV/zP/ADP/MzP/ZjP/mTP/zDP//2YAAGYAM2YAZm
                YAmWYAzGYA/2YrAGYrM2YrZmYrmWYrzGYr/2ZVAGZVM2ZVZmZVmWZVzGZV/2aAAGa
                AM2aAZmaAmWaAzGaA/2aqAGaqM2aqZmaqmWaqzGaq/2bVAGbVM2bVZmbVmWbVzGbV
                /2b/AGb/M2b/Zmb/mWb/zGb//5kAAJkAM5kAZpkAmZkAzJkA/5krAJkrM5krZpkrm
                ZkrzJkr/5lVAJlVM5lVZplVmZlVzJlV/5mAAJmAM5mAZpmAmZmAzJmA/5mqAJmqM5
                mqZpmqmZmqzJmq/5nVAJnVM5nVZpnVmZnVzJnV/5n/AJn/M5n/Zpn/mZn/zJn//8wA
                AMwAM8wAZswAmcwAzMwA/8wrAMwrM8wrZswrmcwrzMwr/8xVAMxVM8xVZsxVmcxVzMx
                V/8yAAMyAM8yAZsyAmcyAzMyA/8yqAMyqM8yqZsyqmcyqzMyq/8zVAMzVM8zVZszVmcz
                VzMzV/8z/AMz/M8z/Zsz/mcz/zMz///8AAP8AM/8AZv8Amf8AzP8A//8rAP8rM/8rZv8
                rmf8rzP8r//9VAP9VM/9VZv9Vmf9VzP9V//+AAP+AM/+AZv+Amf+AzP+A//+qAP+qM/+
                qZv+qmf+qzP+q///VAP/VM//VZv/Vmf/VzP/V////AP//M///Zv//mf//zP///wAAAAA
                AAAAAAAAAAAj/ANUIHEiwoMGDCBMqXMiwocOHEAXu2xex4cSJAComvLgPgMeMGg1i/Ph
                xI0eKJkl6DBVqI8mJCEeuDAWK5UGZH1EWlMmyJs2WIlV6vEhQZs2joUyF0imx40qVJ53
                OPIUUFMymJGsKJQkPHgCfYFkCmJgGZ89QW8vBUxevK7yfoKiiHYsTwE+xJMmVU3tR3du
                wacuRewsXr9py6jh2PQtKqF548fz+vbu3nL6TfgvP9bh3rWKWcvWSo3cSsrpQcu+GQtz
                1M9LK+fIlu7y2K1LGaiNfnni6sF6/ktctPvtTaee18/T5vb26MjzhXU/fDvu7reSkd0E
                JPi6cJfaaSn02ble7GKxczs4Ze6+Kdvu39YxnFg57V+lH+OxXTj9FM3XVmeD1B5Z+crE
                XoE/3sVfcXEj5R9xd8imlmVZ2zQcXVbfpZyFSHjnYn2pxUUhcTQVuaKJ4KJqI2m0Shjd
                idqm5aOGMIcJooI31ZadjKAEBADs='''

    iconimage = tk.PhotoImage(data=icondata)

    #img = tk.Image("photo", file=os.path.join("icons","icon.gif"))



    root.tk.call('wm','iconphoto',root._w,iconimage)




    app = GUI(root)#def app


    root.mainloop()  #do the GUI stuff


    try:
        root.destroy() #destroy the GUI if the user hits file exit
    except:
        pass # dont if they hit the x button because its already gone




if __name__ == '__main__': #do the main
    main()


