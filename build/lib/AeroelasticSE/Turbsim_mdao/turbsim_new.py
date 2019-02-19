from openmdao.api import Group, Problem, Component, IndepVarComp, ParallelGroup
from AeroelasticSE.Turbsim_mdao.turbsim_writer import TurbsimBuilder
from AeroelasticSE.Turbsim_mdao.turbsim_wrapper import Turbsim_wrapper
from AeroelasticSE.Turbsim_mdao.turbsim_reader import turbsimReader
import sys
import os
import numpy as np


# ===================== OpenMDAO Components and Groups =====================

class TurbsimWorkflow (Component):
    """ An OpenMDAO Component for running the FST (FAST) workflow

    This class is a basic workflow utlizing the FST reader, writer, and wrapper
    in OpenMDAO. User is expected to subclass or modify
    this class to enable OpenMDAO-style connectivity to FAST variables of interest."""

    def __init__(self, config, case):
        """ A FAST component that executes the workflow. Takes a single dictionary, config, as well as a case
        name as the input. Executes FAST workflow in the associated case running directory.
        """
        super(TurbsimWorkflow, self).__init__()

        # Initialize objects
        self.reader = turbsimReader()
        self.writer = TurbsimBuilder()
        self.wrapper = Turbsim_wrapper()
        
        for key in config:
            if key == 'Turbsim_masterfile':
                self.reader.Turbsim_masterfile = config[key]
                #print self.reader.Turbsim_masterfile
            elif key == 'Turbsim_masterdir':
                self.reader.Turbsim_masterdir = config[key]
                #print self.reader.Turbsim_masterdir
            elif key == 'Turbsim_runfile':
                self.writer.Turbsim_runfile = config[key]
                #print self.writer.Turbsim_runfile 
            elif key == 'Turbsim_rundir':
                self.writer.Turbsim_rundir = config[key]
                #print self.writer.Turbsim_rundir
            elif key == 'Turbsim_exe':
                self.wrapper.Turbsim_exe = config[key]
                #print self.wrapper.Turbsim_exe
                
        self.reader.execute()   #Read/populate vartrees
        self.writer.turbsim_vt = self.reader.turbsim_vt   #Pass vartrees from reader to write
        self.writer.InputConfig(**config)   #Edit vartrees according to keys in config dictionary
        
        # Pass file and directory names from writer to wrapper
        self.wrapper.Turbsim_runfile = self.writer.Turbsim_runfile
        self.wrapper.Turbsim_rundir = self.writer.Turbsim_rundir
        
        #print 'check wrapper rundir'
        #print self.wrapper.Turbsim_rundir
        
        #print 'writer rundir'
        #print self.writer.Turbsim_rundir
        
        #print 'check wrapper runfile'
        #print self.wrapper.Turbsim_runfile
        
        #print 'writer runfile'
        #print self.writer.Turbsim_runfile
        
    def solve_nonlinear(self, params, unknowns, resids):

        # Create running directory if it doesn't exist
        if not os.path.isdir(self.writer.Turbsim_rundir):
            os.makedirs(self.writer.Turbsim_rundir)

        # Write new analysis files
        self.writer.execute()

        # Execute analysis
        self.wrapper.execute()



class TurbSim_Solver(Group):
    """
    OpenMDAO group to execute Turbsim components in parallel.

    Here, 'configs' is a dictionary of dictionaries, unlike in FSTWorkflow
    where 'config' is merely a dictionary. 'caseids' is similarly an array of
    caseid strings.
    """
    def __init__(self, configs, caseids):
        super(TurbSim_Solver, self).__init__()

        #self._check_config(configs, caseids) #could write function to check setup

        pgT = self.add('pgT', ParallelGroup()) # add parallel group
        
        #Loop over cases, add them to the parallel group
        case_num = len(caseids)
        for i in range(case_num):
            pgT.add(caseids[i], TurbsimWorkflow(configs[caseids[i]], caseids[i]))


if __name__=="__main__":
    pass
