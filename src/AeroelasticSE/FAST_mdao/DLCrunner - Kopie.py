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
import string
import os, sys
sys.path.insert(0, os.path.abspath(".."))
import matplotlib.pyplot as plt

from openmdao.api import Group, Problem, Component, IndepVarComp
from openmdao.api import ParallelGroup, ParallelFDGroup
from openmdao.api import SqliteRecorder
from AeroelasticSE.FAST_mdao.FST8_aeroelasticsolver import FST8Workflow
from AeroelasticSE.Turbsim_mdao.turbsim_openmdao import turbsimGroup
from openmdao.core.mpi_wrap import MPI
if MPI:
    from openmdao.core.petsc_impl import PetscImpl as impl
else:
    from openmdao.core.basic_impl import BasicImpl as impl
from FST8_aeroelasticsolver import FST8Workflow, FST8AeroElasticSolver

# Initial OpenMDAO problem setup for parallel group
top = Problem(impl=impl, root=ParallelFDGroup(1))
root = top.root

# Setup genral dlc config dictionary.
dlc_cfg = {} #dictionary describing what we do
dlc_cfg['DLC'] = 1.2
dlc_cfg['RunTurbsim'] = False #shall we run turbsim first to generate turbulent wind files?
dlc_cfg['RunFAST'] = True #shall we run FAST?
dlc_cfg['WSBin1'] = 3.0 #First wind speed to consider (generally equal to cut-in)
dlc_cfg['WSBin2'] = 24.0 #Last wind speed to consider (generally equal to cut-out)
dlc_cfg['WSBinWidth'] = 2 #Bin width
dlc_cfg['NoSeeds'] = 6 #Number of seeds per wind speed bin
dlc_cfg['NoTIBins'] = 6 #Is the turbulence intensity weibull distributed and should be binned?
dlc_cfg['DiscreteYaw'] = True #Is the turbulence intensity weibull distributed and should be binned?
dlc_cfg['YawAngles'] = [-8.0, 0.0, 8.0] #Yaw angles to be used




dlc_cfg['WSBins'] = np.arange(dlc_cfg['WSBin1'], dlc_cfg['WSBin2'], dlc_cfg['WSBinWidth'])
if dlc_cfg['WSBins'][-1] != dlc_cfg['WSBin2']:
    dlc_cfg['WSBins'] = np.append(dlc_cfg['WSBins'], dlc_cfg['WSBin2'])#
dlc_cfg['SeedList'] = list(string.ascii_lowercase)[0:dlc_cfg['NoSeeds']]
dlc_cfg['TIBins'] = np.arange(dlc_cfg['NoTIBins']) + 1

print dlc_cfg['WSBins']
print dlc_cfg['SeedList']
print dlc_cfg['TIBins']


# Setup input config dictionary of dictionaries.
cfg_master = {} #master config dictionary (dictionary of dictionaries)

#generate case ids
runvar = 0
caseids = [];
short_caseids = []
for i in range(len(dlc_cfg['WSBins'])):
    for j in range(len(dlc_cfg['SeedList'])):   
            for k in range(len(dlc_cfg['TIBins'])):
                for l in range(len(dlc_cfg['YawAngles'])):
                    caseid = 'case%02.0f_%s_%i' % (dlc_cfg['WSBins'][i], dlc_cfg['SeedList'][j], dlc_cfg['TIBins'][k])
                    short_caseid = '%02.0f_%s_%i' % (dlc_cfg['WSBins'][i], dlc_cfg['SeedList'][j], dlc_cfg['TIBins'][k]) 
    
                    caseids.append(caseid)
                    short_caseids.append(short_caseid)
                    cfg = {}
                    cfg['fst_runfile'] = '{0}.fst'.format(short_caseids[runvar])
                    cfg['fst_rundir'] = os.path.join('./ExampleCase/02_Run/dlc12/', short_caseids[runvar])
    
                    # These parameters the same for all cases
                    cfg['fst_masterfile'] = 'Test01.fst' 
                    cfg['fst_masterdir']= './ExampleCase/02_Run/model/'
                    cfg['fst_exe'] = 'FAST_gwin64'
                    cfg['ad_file_type'] = 1 
    
                    # Put dictionary into master dictionary, keyed by caseid
                    cfg_master[caseids[runvar]] = cfg
                
                    runvar = runvar + 1
            
                       
print cfg_master[caseids[i]]
print cfg['fst_runfile']
            
# Add parallel group to omdao problem, pass in master config file
root.add('ParallelFASTCases', FST8AeroElasticSolver(cfg_master, caseids))

top.setup()
top.run()

top.cleanup()   #Good practice, especially when using recorder



