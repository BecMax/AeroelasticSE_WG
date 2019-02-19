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
import random
import ast
sys.path.insert(0, os.path.abspath(".."))
import matplotlib.pyplot as plt
import pandas as pd

from openmdao.api import Group, Problem, Component, IndepVarComp
from openmdao.api import ParallelGroup, ParallelFDGroup
from openmdao.api import SqliteRecorder
#from AeroelasticSE.FAST_mdao.FST8_aeroelasticsolver import FST8Workflow
#from AeroelasticSE.Turbsim_mdao.turbsim_new import TurbsimWorkflow
from openmdao.core.mpi_wrap import MPI
if MPI:
    from openmdao.core.petsc_impl import PetscImpl as impl
else:
    from openmdao.core.basic_impl import BasicImpl as impl
from FST8_aeroelasticsolver import FST8Workflow, FST8AeroElasticSolver
from AeroelasticSE.Turbsim_mdao.turbsim_new import TurbsimWorkflow, TurbSim_Solver

from openmdao.api import view_tree

# Initial OpenMDAO problem setup for parallel group
top = Problem(impl=impl, root=ParallelFDGroup(1))
root = top.root

# ===================== Input Section =====================
# Setup genral dlc config dictionary
dlc_cfg = {} #dictionary describing what we do
dlc_cfg['SettingsFromExcel'] = True
dlc_cfg['ExcelFile'] = 'dlc12_settings.xlsx'
dlc_cfg['DLC'] = [1.2]
dlc_cfg['RunTurbsim'] = False #shall we run turbsim first to generate turbulent wind files?
dlc_cfg['RunFAST'] = False #shall we run FAST?

# Wind speed bins
dlc_cfg['CalcWSBins'] = False #Shall the Wind Speed Bins be calculated with a fixed bin width; alternatively they can be set
dlc_cfg['WSBinLow'] = [3.0] #First Wind Speed to consider (generally equal to cut-in); only used if CalcWSBins is True
dlc_cfg['WSBinHigh'] = [24.0] #Last Wind Speed to consider (generally equal to cut-out); only used if CalcWSBins is True
dlc_cfg['WSBinWidth'] = [2] #Bin Width; only used if CalcWSBins is True
dlc_cfg['WSBins'] = [12.0, 15.0]  #WSBins; only used if CalcWSBins is False

# Seeds, TI Bins, Shear
dlc_cfg['NoSeeds'] = [2] #Number of Seeds per wind speed bin
dlc_cfg['NoTIBins'] = [6] #Is the Turbulence Intensity weibull distributed and should be binned?
dlc_cfg['ShearExp'] = [0.2, 0.16] # Shear Exponents

# Turbine Angles
dlc_cfg['WindDir'] = [-8.0, 0.0, 8.0] #Yaw Angles to be used; sets automatically to [0.0] if DiscreteYaw is False
dlc_cfg['RotorAngles'] = [0.0, 90.0] #Rotor Position Angles to be used; sets automatically to [0.0] if DiscreteRotorPos' is False
dlc_cfg['PitchAngles'] = [0.0, 90.0] #Pitch Angles to be used; sets automatically to [0.0] if DiscreteRotorPos' is False


# ===================== Some Parameters from Inputs =====================

if dlc_cfg['SettingsFromExcel'] == True:
    df = pd.read_excel(os.path.join('./ExampleCaseOF/02_Run/dlc12', dlc_cfg['ExcelFile']), comment='#' ) 

    df = df.set_index('Parameters').T
    for key in dlc_cfg:
        for col in df.columns:
            if key == col:
                if col == 'RunTurbsim' or col == 'RunFAST' or col == 'CalcWSBins':
                    dlc_cfg[key] = ast.literal_eval(df[col].dropna().to_numpy()[0])
                    print 'The setting for %s is read from Excel file it is %r' % (col, dlc_cfg[key])
                else:
                    dlc_cfg[key] = df[col].dropna().to_numpy()
                    print 'The setting for %s is read from Excel file it first value equals %10.4f' % (col, dlc_cfg[key][0])

print dlc_cfg['WSBinLow']
print dlc_cfg['WSBinHigh']
         
if dlc_cfg['CalcWSBins'] == True:
    dlc_cfg['WSBins'] = np.arange(dlc_cfg['WSBinLow'][0], dlc_cfg['WSBinHigh'][0], dlc_cfg['WSBinWidth'][0])
    if dlc_cfg['WSBins'][-1] != dlc_cfg['WSBinHigh']:
        dlc_cfg['WSBins'] = np.append(dlc_cfg['WSBins'], dlc_cfg['WSBinHigh'])
dlc_cfg['NoWSBins'] = len(dlc_cfg['WSBins'])

dlc_cfg['NoSeeds'] = dlc_cfg['NoSeeds'][0]
dlc_cfg['NoTIBins'] = dlc_cfg['NoTIBins'][0]

dlc_cfg['SeedList'] = list(string.ascii_lowercase)[0:dlc_cfg['NoSeeds']]
dlc_cfg['TIBins'] = np.arange(dlc_cfg['NoTIBins']) + 1

dlc_cfg['NoShearExp'] = len(dlc_cfg['ShearExp'])
dlc_cfg['NoWindDir'] = len(dlc_cfg['WindDir'])
dlc_cfg['NoRotorAngles'] = len(dlc_cfg['RotorAngles'])
dlc_cfg['NoPitchAngles'] = len(dlc_cfg['PitchAngles'])


if dlc_cfg['SettingsFromExcel'] == True:
    if dlc_cfg['NoTIBins'] > 1:
        dlc_cfg['TIBins'] = np.zeros(shape=(dlc_cfg['NoWSBins'], dlc_cfg['NoTIBins']))
        for col in df.columns:
            if col[0:2] == 'TI':
                if col[2:6] == 'Bins':	
                    data =  df[col].dropna().to_numpy()
                    idxspeed = data[0]
                    TIdata = data[1:]
                    idx = np.where(dlc_cfg['WSBins']==idxspeed)
                    dlc_cfg['TIBins'][idx[0][0],:] = TIdata
                
print dlc_cfg['TIBins']

# ===================== Coding =====================

# Setup input config dictionary of dictionaries.
cfg_master_wind = {} #master config dictionary (dictionary of dictionaries)
cfg_master = {} #master config dictionary (dictionary of dictionaries)

#allocate case ids
casestr = 'case_'
caseids_wind = []
caseids = []
short_caseids = []

for i in range(dlc_cfg['NoWSBins']):
    for j in range(dlc_cfg['NoSeeds']):
        for k in range(dlc_cfg['NoTIBins']):
            for l in range(dlc_cfg['NoShearExp']):

                windid = '%02.0f'  % dlc_cfg['WSBins'][i]
                if dlc_cfg['NoSeeds'] > 1:
                    windid = windid + '_%s'  % dlc_cfg['SeedList'][j]
                if dlc_cfg['NoTIBins'] > 1:
                    windid = windid + '_I%01.0f' % dlc_cfg['TIBins'][i,k]                   
                if dlc_cfg['NoShearExp'] > 1:
                    shearstr = 100.0*dlc_cfg['ShearExp'][l]
                    windid = windid + '_x%01.0f' % shearstr
                    
                caseid_wind = casestr + windid
                caseids_wind.append(caseid_wind)
      
                cfg_wind = {}
                cfg_wind['Turbsim_runfile'] ='{0}.inp'.format(windid)
                cfg_wind['Turbsim_rundir'] = './ExampleCaseOF/04_Wind/'
      
                # These parameters the same for all turbsim cases
                cfg_wind['Turbsim_masterfile'] = 'Test18.inp' 
                cfg_wind['Turbsim_masterdir']= './ExampleCaseOF/04_Wind/model/'
                cfg_wind['Turbsim_exe']= 'TurbSim_x64'
       
                # change wind speed
                cfg_wind['URef'] = dlc_cfg['WSBins'][i]
                # random seed in range according to turbsim
                cfg_wind['RandSeed1'] = random.randrange(-2147483648,2147483647,1)
  
                # Put dictionary into master dictionary, keyed by caseid
                cfg_master_wind[caseid_wind] = cfg_wind
        
                for m in range(dlc_cfg['NoWindDir']):
                    for n in range(dlc_cfg['NoRotorAngles']):
                        for o in range(dlc_cfg['NoPitchAngles']):

                            short_caseid = windid
                            if dlc_cfg['NoWindDir'] > 1:
                                short_caseid = short_caseid + '_%01.0f' % dlc_cfg['WindDir'][m]
                            if dlc_cfg['NoRotorAngles'] > 1:
                                short_caseid = short_caseid + '_A%01.0f' % dlc_cfg['RotorAngles'][n]
                            if dlc_cfg['NoPitchAngles'] > 1:
                                short_caseid = short_caseid + '_P%01.0f' % dlc_cfg['PitchAngles'][o]
                        
                        caseid = casestr + short_caseid
                        caseid = caseid.replace("-", "m") # replace minus with m as minus string causes problems in OpenMDAO group definition
                        print caseid
                        
                    
                        caseids.append(caseid)
                        short_caseids.append(short_caseid)
                        
                        
                        cfg = {}
                        cfg['fst_runfile'] = '{0}.fst'.format(short_caseid)
                        cfg['fst_rundir'] = os.path.join('./ExampleCaseOF/02_Run/dlc12/', short_caseid)
                        
                        cfg['Filename'] = os.path.join('../../../04_Wind/', '{0}.bts'.format(windid))
        
                        # These parameters the same for all cases
                        cfg['fst_masterfile'] = 'dlc12.fst' 
                        cfg['fst_masterdir']= './ExampleCaseOF/02_Run/dlc12/case/'
                        cfg['fst_exe'] = 'openfast' #'FAST_gwin64' #'openfast'
                        cfg['pass2Numpy'] = True

                        
                       
                        cfg['writeElasto'] = True
                        cfg['writeBladeStruc'] = False 
                        cfg['writeTower'] = False 
                        cfg['writeInflow'] = True
                        cfg['writeAero'] = False
                        cfg['writeServo'] = True
                        
                        cfg['copyBeamDyn'] = False
                        cfg['copyAirfoils'] = False
                        cfg['copyAeroBlade'] = False
                        cfg['copyDLL'] = False
                        cfg['copyDLLinfile'] = False
                        cfg['copyTMDamp'] = False
        
                        cfg['PropagationDir'] = dlc_cfg['WindDir'][m]
                        cfg['Azimuth'] = dlc_cfg['RotorAngles'][n]
        
                        # Put dictionary into master dictionary, keyed by caseid
                        cfg_master[caseid] = cfg
            
                       
print cfg_master[caseids[i]]
            
# Add parallel group to omdao problem, pass in master config file
# strangely Turbsim needs to come second to be run first
if dlc_cfg['RunFAST'] == True:
    root.add('ParallelFASTCases', FST8AeroElasticSolver(cfg_master, caseids))
if dlc_cfg['RunTurbsim'] == True:   
    root.add('ParallelTurbsimCases', TurbSim_Solver(cfg_master_wind, caseids_wind))


top.setup()




top.check_setup()
#view_model(top)


top.run()

top.cleanup()   #Good practice, especially when using recorder



