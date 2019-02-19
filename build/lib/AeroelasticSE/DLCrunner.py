"""
Demonstration of setting up an OpenMDAO 1.x problem using the FST8Workflow component
(in FST8_aeroelasticsolver), which executes the FST8 reader, writer, and wrapper and assigns
all variables in the FAST outlist to OpenMDAO 'Unknowns'. It also
implements an "input config" function which allows the user to put all variables that they
wish to explicitly define into a dictionary. The input config function assigns these
variables to the correct locations in the variable tree.
"""
# Hacky way of doing relative imports
import numpy as np
import os, sys
sys.path.insert(0, os.path.abspath(".."))
import matplotlib.pyplot as plt
from openmdao.api import Group, Problem, Component, IndepVarComp, ParallelGroup
from openmdao.api import SqliteRecorder
from AeroelasticSE.FAST_mdao.FST8_aeroelasticsolver import FST8Workflow
from AeroelasticSE.Turbsim_mdao.turbsim_openmdao import turbsimGroup