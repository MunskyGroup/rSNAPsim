# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 19:42:52 2018

@author: willi
"""

from . import CodonDictionaries
cdict = CodonDictionaries.CodonDictionaries()

from . import AuxCodonDicts
aux_cdict = AuxCodonDicts.AuxCodonDicts()

from . import DiffusionRateCalc
diffcalc = DiffusionRateCalc.DiffusionRateCalc()


from . import SequenceManipMethods 
seqmanip = SequenceManipMethods.SequenceManipMethods()

from . import IntensityAnalyses
inta = IntensityAnalyses.IntensityAnalyses()

from . import IntensityAnalysesRagged
inta_r = IntensityAnalysesRagged.IntensityAnalysesRagged()

from . import expv
expv = expv.expv

from . import poi
poi = poi.poi

from . import FileParser
fp = FileParser.FileParser()

from . import ProbeVectorFactory
probef = ProbeVectorFactory.ProbeVectorFactory()

from . import PropensityFactory
propf = PropensityFactory.PropensityFactory()

from . import ModelBuilder
model_builder = ModelBuilder.ModelBuilder()

from . import TranslationSolvers 
solver = TranslationSolvers.TranslationSolvers()

from . import SSA_Soln
solution_obj = SSA_Soln

from . import ODE_Soln
ode_solution_obj = ODE_Soln

from . import TranslationOptimization
optimizer = TranslationOptimization
from ._version import __version__
from . import GenericMetaData



# from . import Core
# from . import Seq
# from .Core import *
# from .IO import *
# from .Seq import *
# from .Opt import *

# from .IntA import *
# from .Pfact import *
# from .Solvers import *
