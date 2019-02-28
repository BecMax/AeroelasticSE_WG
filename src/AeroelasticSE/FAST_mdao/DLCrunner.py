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
from os.path import relpath

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

import argparse

# https://stackoverflow.com/questions/40807862/python-argparse-with-mandatory-input-file-argument
#Create the parser
parser = argparse.ArgumentParser(description='Does some stuff with an input file.')

#add the argument
parser.add_argument('-i', '--input', dest='infile', required=True, metavar='EXCEL_FILE', help='The excel input file.')
parser.add_argument('-s', '--string', dest='insheet', type=str, required=True, metavar='EXCEL_SHEET', help='The excel sheet as a string.')

#parse and assign to the variable
args = parser.parse_args()
infile=args.infile
insheet=args.insheet

# Initial OpenMDAO problem setup for parallel group
top = Problem(impl=impl, root=ParallelFDGroup(1))
root = top.root

# ===================== Input Section =====================

rundir = '02_Run'
winddir = '04_Wind'
postdir = '03_Postpro'

# Initialize general dlc config dictionary
dlc_cfg = {} #dictionary describing what we do
dlc_cfg['SettingsFromExcel'] = True
dlc_cfg['Standard'] = '' #IEC, DIBT, DNVGL
dlc_cfg['DLC'] = ''
dlc_cfg['WindModel'] = ''
dlc_cfg['WindFileFormat'] = '' # Either Bladed (.wnd) or Turbsim format (.bts)
dlc_cfg['Generate2ndSeed'] = ''
dlc_cfg['RunTurbsim'] = False #shall we run turbsim first to generate turbulent wind files?
dlc_cfg['OutFileFmt'] = [0]
dlc_cfg['RunFAST'] = False #shall we run FAST?

# Wind speed bins
dlc_cfg['CalcWSBins'] = False #Shall the Wind Speed Bins be calculated with a fixed bin width; alternatively they can be set
dlc_cfg['WSBinLow'] = [0.0] #First Wind Speed to consider (generally equal to cut-in); only used if CalcWSBins is True
dlc_cfg['WSBinHigh'] = [0.0] #Last Wind Speed to consider (generally equal to cut-out); only used if CalcWSBins is True
dlc_cfg['WSBinWidth'] = [0] #Bin Width; only used if CalcWSBins is True
dlc_cfg['WSBins'] = [0.0, 0.0]  #WSBins; only used if CalcWSBins is False
dlc_cfg['WSRated'] = [0.0]  #Rated Windspeed required to calculate difference to Bins from Rated for IECWind generated files

# Seeds, TI Bins, Shear
dlc_cfg['NoSeeds'] = [0] #Number of Seeds per wind speed bin
dlc_cfg['NoTIBins'] = [0] #Is the Turbulence Intensity weibull distributed and should be binned?
dlc_cfg['ShearExp'] = [0.0, 0.0] # Shear Exponents
dlc_cfg['SignChanges'] = ['',''] # SignChanges negative and positive

# Turbine Angles
dlc_cfg['WindAsYawMisalignment'] = ''
dlc_cfg['WindDir'] = [0.0, 0.0] #Yaw Angles to be used; sets automatically to [0.0] if DiscreteYaw is False
dlc_cfg['RotorAngles'] = [0.0, 0.0] #Rotor Position Angles to be used; sets automatically to [0.0] if DiscreteRotorPos' is False
dlc_cfg['BladesToPitch'] = [0]
dlc_cfg['PitchAngles'] = [0.0, 0.0] #Pitch Angles to be used; sets automatically to [0.0] if DiscreteRotorPos' is False


# ===================== Some Parameters from Inputs =====================

if dlc_cfg['SettingsFromExcel'] == True:
    df = pd.read_excel(os.path.join(os.getcwd(), infile), comment='#', sheet_name= insheet) 

    df = df.set_index('Parameters').T
    for key in dlc_cfg:
        for col in df.columns:
            if key == col:
                # read as string
                if col == 'DLC':
                    dlc_cfg[key] = str(df[col].dropna().to_numpy()[0])
                # read strings
                elif col == 'Standard' or col == 'WindModel' or col == 'WindFileFormat' or col == 'Generate2ndSeed':
                    dlc_cfg[key] = df[col].dropna().to_numpy()[0]
                    print 'The setting for %s is read from Excel file it is %s' % (col, dlc_cfg[key])
                # read list of strings
                elif col == 'SignChanges':
                    dlc_cfg[key] = df[col].dropna().to_numpy()
                    print 'The setting for %s is read from Excel file its first value equals %s' % (col, dlc_cfg[key][0])
                # read bools
                elif col == 'RunTurbsim' or col == 'RunFAST' or col == 'CalcWSBins' or col == 'WindAsYawMisalignment':
                    dlc_cfg[key] = ast.literal_eval(df[col].dropna().to_numpy()[0])
                    print 'The setting for %s is read from Excel file it is %r' % (col, dlc_cfg[key])
                # read floats into array
                else:
                    dlc_cfg[key] = df[col].dropna().to_numpy()
                    print 'The setting for %s is read from Excel file its first value equals %10.4f' % (col, dlc_cfg[key][0])

if dlc_cfg['Standard'] == 'IEC':
    if dlc_cfg['DLC'][:3] == '1.2' or dlc_cfg['DLC'][:3] == '2.4'  or dlc_cfg['DLC'][:3] == '3.1' or dlc_cfg['DLC'][:3] == '4.1' or dlc_cfg['DLC'][:3] == '6.4':
        rundir = os.path.join(rundir, 'FLS')
        postdir = os.path.join(postdir, 'FLS')
    else:
        rundir = os.path.join(rundir, 'ULS')
        postdir = os.path.join(postdir, 'ULS')
elif dlc_cfg['Standard'] == 'DIBt':
    if dlc_cfg['DLC'][:3] == 'D.4' or dlc_cfg['DLC'][:3] == '1.2'  or dlc_cfg['DLC'][:3] == '2.4' or dlc_cfg['DLC'][:3] == '3.1' or dlc_cfg['DLC'][:3] == '4.1' or dlc_cfg['DLC'][:3] == '6.4':
        rundir = os.path.join(rundir, 'FLS')
        postdir = os.path.join(postdir, 'FLS')
    else:
        rundir = os.path.join(rundir, 'ULS')
        postdir = os.path.join(postdir, 'ULS')
#elif dlc_cfg['Standard'] == 'DNVGL':        
         
if dlc_cfg['CalcWSBins'] == True:
    dlc_cfg['WSBins'] = np.arange(dlc_cfg['WSBinLow'][0], dlc_cfg['WSBinHigh'][0], dlc_cfg['WSBinWidth'][0])
    if dlc_cfg['WSBins'][-1] != dlc_cfg['WSBinHigh']:
        dlc_cfg['WSBins'] = np.append(dlc_cfg['WSBins'], dlc_cfg['WSBinHigh'])
dlc_cfg['NoWSBins'] = len(dlc_cfg['WSBins'])

dlc_cfg['NoSeeds'] = dlc_cfg['NoSeeds'][0]
dlc_cfg['NoTIBins'] = dlc_cfg['NoTIBins'][0]

dlc_cfg['SeedList'] = list(string.ascii_lowercase)[0:dlc_cfg['NoSeeds']]
#dlc_cfg['TIBins'] = np.arange(dlc_cfg['NoTIBins']) + 1

dlc_cfg['NoShearExp'] = len(dlc_cfg['ShearExp'])
dlc_cfg['NoSignChanges'] = len(dlc_cfg['SignChanges'])

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
                


# ===================== Coding =====================

# Setup input config dictionary of dictionaries.
cfg_master_wind = {} #master config dictionary (dictionary of dictionaries)
cfg_master = {} #master config dictionary (dictionary of dictionaries)

#allocate case ids
casestr = 'case_'
caseids_wind = []
caseids = []
short_caseids = []
counter = 0

for i in range(dlc_cfg['NoWSBins']):
    for j in range(dlc_cfg['NoSeeds']):
        for k in range(dlc_cfg['NoTIBins']):
            for l in range(dlc_cfg['NoShearExp']):
                for m in range(dlc_cfg['NoSignChanges']):

                    windid = '%02.0f'  % dlc_cfg['WSBins'][i]
                    if dlc_cfg['NoSeeds'] > 1:
                        windid = windid + '_%s'  % dlc_cfg['SeedList'][j]
                    if dlc_cfg['NoTIBins'] > 1:
                        windid = windid + '_I%01.0f' % dlc_cfg['TIBins'][i,k]                   
                    if dlc_cfg['NoShearExp'] > 1:
                        shearstr = 100.0*dlc_cfg['ShearExp'][l]
                        windid = windid + '_x%01.0f' % shearstr
                    if dlc_cfg['NoSignChanges'] > 1:
                        print str(dlc_cfg['SignChanges'][m])[4:7]
                        if str(dlc_cfg['SignChanges'][m])[0:3] == 'Pos':
                            print 'in the if'
                            #positive change:
                            dirchangestr1 = '+'
                            dirchangestr11 = 'Pos' #OpenMDAO cant deal with +/- signs in case ids
                        elif str(dlc_cfg['SignChanges'][m])[0:3] == 'Neg':
                            #negative change:
                            dirchangestr1 = '-'
                            dirchangestr11 = 'Neg' #OpenMDAO cant deal with +/- signs in case ids
                        else:
                            print 'WindDirChange needs to be either positive or negative sweep'
                        windid = windid + '_%s' % dirchangestr11
      
                    cfg_wind = {}
              
                    # Turbsim case windfiles
                    if dlc_cfg['WindModel'] == 'NTM':
                        cfg_wind['Turbsim_masterfile'] = 'NTM.inp' 
                        cfg_wind['Turbsim_masterdir']= winddir
                        windinput = windid
                    elif dlc_cfg['WindModel'] == 'NWP':
                        # Just NTM without any turbulence, thus TI = 0
                        cfg_wind['Turbsim_masterfile'] = 'NWP.inp' 
                        cfg_wind['Turbsim_masterdir']= winddir
                        windinput = windid + '_NWP'
                        pass
                    elif dlc_cfg['WindModel'] == 'ETM':
                        cfg_wind['Turbsim_masterfile'] = 'ETM.inp' 
                        cfg_wind['Turbsim_masterdir']= winddir
                        windinput = windid + '_ETM'
                    
                    # These parameters are the same for all turbsim cases    
                    if dlc_cfg['WindModel'] == 'NTM' or dlc_cfg['WindModel'] == 'NWP' or dlc_cfg['WindModel'] == 'ETM':
                        cfg_wind['Turbsim_runfile'] ='{0}.inp'.format(windinput)
                        cfg_wind['Turbsim_rundir'] = winddir
                        cfg_wind['Turbsim_exe']= 'TurbSim_gwin64'
                    else: 
                        difftorated = dlc_cfg['WSBins'][i] - dlc_cfg['WSRated'][0]
                        difftocutin = dlc_cfg['WSBins'][i] - dlc_cfg['WSBinLow'][0]
                        difftocutout = dlc_cfg['WSBins'][i] - dlc_cfg['WSBinHigh'][0]


                    if  dlc_cfg['WindModel'] == 'ECD':
                        if difftorated == 0.0:
                            windinput = 'ECD%sR' % dirchangestr1
                        else:
                            windinput = 'ECD%sR%+01.1f' % (dirchangestr1, difftorated)
                    if  dlc_cfg['WindModel'] == 'EWS':
                        if str(dlc_cfg['SignChanges'][m])[4:7] == 'Ver':
                            windinput = 'EWSV%+02.1f' % dlc_cfg['WSBins'][i]
                            windid = windid + '_%s' % 'Ver'
                        elif str(dlc_cfg['SignChanges'][m])[4:7] == 'Hor':
                            windinput = 'EWSH%+02.1f' % dlc_cfg['WSBins'][i]
                            windid = windid + '_%s' % 'Hor'    
                    if  dlc_cfg['WindModel'] == 'EOG':
                        if difftorated == 0.0:
                            windinput = 'EOGR' 
                        elif difftocutin == 0.0:
                            windinput = 'EOGI' 
                        elif difftocutout == 0.0:
                            windinput = 'EOGO' 
                        else:
                            windinput = 'EOGR%+01.1f' % difftorated
                    if  dlc_cfg['WindModel'] == 'EDC':
                        if difftorated == 0.0:
                            windinput = 'EDC%sR' % dirchangestr1
                        elif difftocutin == 0.0:
                            windinput = 'EDC%sI' % dirchangestr1
                        elif difftocutout == 0.0:
                            windinput = 'EDC%sO' % dirchangestr1
                        else:
                            windinput = 'EDC%sR%+01.1f' % (dirchangestr1, difftorated)

                
                    if dlc_cfg['WindFileFormat'] == '.wnd':   
                        windinputfile = '{0}.wnd'.format(windinput)
                    elif dlc_cfg['WindFileFormat'] == '.bts':
                        windinputfile = '{0}.bts'.format(windinput)

                    # Only run Turbsim if input wind files are not there
                    if os.path.isfile(os.path.join(winddir, windinputfile)) == False:
                        if dlc_cfg['WindModel'] == 'NTM' or dlc_cfg['WindModel'] == 'NWP' or dlc_cfg['WindModel'] == 'ETM':
                            print 'Wind input file %s does not exist running Turbsim to generate it' % windinputfile
                   
                            caseid_wind = casestr + windid
                            caseids_wind.append(caseid_wind)
                            
                            # change wind speed
                            cfg_wind['URef'] = dlc_cfg['WSBins'][i]
                            # random seed in range according to turbsim
                            cfg_wind['RandSeed1'] = random.randrange(-2147483648,2147483647,1)
                            if dlc_cfg['Generate2ndSeed'] == 'RandInt':
                                cfg_wind['RandSeed2'] = random.randrange(-2147483648,2147483647,1)
                            # TI
                            if dlc_cfg['NoTIBins'] > 1:
                                cfg_wind['IECturbc'] = dlc_cfg['TIBins'][i][k] #check this
                            # Shear exp
                            if dlc_cfg['NoShearExp'] > 1:
                                cfg_wind['PLExp'] = dlc_cfg['ShearExp'][l]

     
                            dlc_cfg['RunTurbsim'] = True
                   
                            # Put dictionary into master dictionary, keyed by caseid
                            cfg_master_wind[caseid_wind] = cfg_wind
                        else:
                            print 'Wind input file %s does not exist need to run IECWind to generate it' % windinputfile
                  
                    else:
                        print 'Wind input file %s does exist do not need to run Turbsim or IECWind' % windinputfile
                        pass
        
                    for n in range(dlc_cfg['NoWindDir']):
                        for o in range(dlc_cfg['NoRotorAngles']):
                            for p in range(dlc_cfg['NoPitchAngles']):

                                short_caseid = windid
                                if dlc_cfg['NoWindDir'] > 1:
                                    short_caseid = short_caseid + '_%01.0f' % dlc_cfg['WindDir'][n]
                                if dlc_cfg['NoRotorAngles'] > 1:
                                    short_caseid = short_caseid + '_A%01.0f' % dlc_cfg['RotorAngles'][o]
                                if dlc_cfg['NoPitchAngles'] > 1:
                                    short_caseid = short_caseid + '_P%01.0f' % dlc_cfg['PitchAngles'][p]
                        
                            caseid = casestr + short_caseid
                            caseid = caseid.replace("-", "m") # replace minus with m as minus string causes problems in OpenMDAO group definition
                            print caseid
                        
                    
                            caseids.append(caseid)
                            short_caseids.append(short_caseid)
                        
                        
                            cfg = {}
                            
                            dlcstring = dlc_cfg['DLC']
                            
                            cfg['fst_runfile'] = '{0}.fst'.format(short_caseid)
                            cfg['fst_rundir'] = os.path.join(rundir, 'DLC%s', short_caseid) % dlcstring
                        
                            # Filename of inputfile for Wind in InflowWind
                            reldir = relpath(winddir, cfg['fst_rundir'])
                            cfg['Filename'] = os.path.join(reldir, windinputfile)
                            cfg['FilenameRoot'] = cfg['Filename'][:-4]
                            
                            # Output file format
                            cfg['OutFileFmt'] = dlc_cfg['OutFileFmt'][0]
                            if dlc_cfg['OutFileFmt'][0] == 1:
                                cfg['fst_outfile'] = '{0}.out'.format(short_caseid)
                            else:
                                cfg['fst_outfile'] = '{0}.outb'.format(short_caseid)
        
                            # These parameters the same for all cases
                            cfg['fst_masterfile'] = 'DLC%s.fst' % dlcstring
                            cfg['fst_masterdir']= os.path.join(rundir, 'DLC%s', 'case') % dlcstring
                            cfg['fst_exe'] = 'openfast' #'FAST_gwin64'
                            cfg['passToNumpy'] = True

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

                            # Set ElastoDyn inputs (if 'WindAsYawMisalignment' false set in InflowWind)
                            if dlc_cfg['WindAsYawMisalignment'] == True:
                                cfg['NacYaw'] = dlc_cfg['WindDir'][n]
                            else:
                                cfg['PropagationDir'] = dlc_cfg['WindDir'][n]
                                
                            cfg['Azimuth'] = dlc_cfg['RotorAngles'][o]
                            
                            if dlc_cfg['BladesToPitch'][0] == 1.:
                                cfg['BlPitch1'] = dlc_cfg['PitchAngles'][p]
                            elif dlc_cfg['BladesToPitch'][0] == 2.:
                                cfg['BlPitch1'] = dlc_cfg['PitchAngles'][p]
                                cfg['BlPitch2'] = dlc_cfg['PitchAngles'][p]
                            elif dlc_cfg['BladesToPitch'][0] == 3.:
                                cfg['BlPitch1'] = dlc_cfg['PitchAngles'][p]
                                cfg['BlPitch2'] = dlc_cfg['PitchAngles'][p]
                                cfg['BlPitch3'] = dlc_cfg['PitchAngles'][p]
        
                            # Put dictionary into master dictionary, keyed by caseid
                            cfg_master[caseid] = cfg
                            
                            counter = counter + 1
                        
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

# ===================== Write Output files into txt =====================

postdir = os.path.join(postdir, 'DLC%s') % dlcstring
if not os.path.exists(postdir):
    os.makedirs(postdir)
print postdir
outfiles_file = os.path.join(postdir, 'Outputfiles.txt')
f = open(outfiles_file, 'w')
f.write('-----  Input Files  ------------------------------------------------------------\n')
f.write('{:3}                 NumFiles          The number of input files to read.\n'.format(counter))
#for id in short_caseids:
for key, value in sorted(cfg_master.iteritems()):
    outpath = relpath(cfg_master[key]['fst_rundir'], postdir)
    print outpath
    outfile = os.path.join(outpath, cfg_master[key]['fst_outfile'])
    print outfile
    f.write('"{:}"\n'.format(outfile))
f.write('==EOF==                             DO NOT REMOVE OR CHANGE.  MUST COME JUST AFTER LAST LINE OF VALID INPUT.')	


