from numpy import zeros, array
from FST_vartrees_out import FstOutput

# This variable tree contains all parameters required to create a FAST model
# for FAST versions 7 and 8.

# .fst Simulation Control
class FstSimCtrl(object):
	def __init__(self):
		self.Echo = False
		self.AbortLevel = ''
		self.TMax = 0.0
		self.DT = 0.0
		self.InterpOrder = 0
		self.NumCrctn = 0
		self.DT_UJac = 0.0
		self.UJacSclFact = 0.0
		self.ADAMSPrep = 0   #FAST7 only
		self.AnalMode = 0   #FAST7 only

# Feature Switches and Flags
class FtrSwtchsFlgs(object):
	def __init__(self):
		self.CompElast = 0
		self.CompInflow = 0
		self.CompAero = 0
		self.CompServo = 0
		self.CompHydro = 0
		self.CompSub = 0
		self.CompMooring = 0
		self.CompIce = 0
		self.CompNoise = 0   #FAST7 only

# Input Files
class InputFiles(object):
	def __init__(self):
		self.EDFile = ''
		self.BDBldFile1 = ''
		self.BDBldFile2 = ''
		self.BDBldFile3 = ''
		self.InflowFile = ''
		self.AeroFile = ''
		self.ServoFile = ''
		self.HydroFile = ''
		self.SubFile = ''
		self.MooringFile = ''
		self.IceFile = ''
		self.ADFile = ''   #FAST7 only
		self.NoiseFile = ''   #FAST7 only
		self.ADAMSFile = ''   #FAST7 only
		self.LinFile = ''   #FAST7 only

# FAST Output Parameters
class FstOutputParams(object):
	def __init__(self):
		self.SumPrint   = False
		self.SttsTime   = 0.0
		self.ChkptTime  = 0.0
		self.DT_Out     = ''
		self.OutFileFmt = 0
		self.TabDelim   = False
		self.OutFmt     = ''

# Linearization
class Linearization(object):
	def __init__(self):
		self.Linearize = False
		self.NLinTimes = 2
		self.LinTimes = (30, 60)
		self.LinInputs = 1
		self.LinOutputs = 1
		self.LinOutJac = False
		self.LinOutMod = False

# Visualization
class Visualization(object):
	def __init__(self):
		self.WrVTK = 0
		self.VTK_type = 0
		self.VTK_fields = False
		self.VTK_fps = 0

# ElastoDyn Simulation Control
class EdSimCtrl(object):
	def __init__(self):
		self.Echo = False
		self.Method = 0
		self.DT = 0.0

# Environmental Condition
class EnvirCond(object):
	def __init__(self):
		self.Gravity = 0.0

# Degrees of Freedom
class DOF(object):
	def __init__(self):
		self.FlapDOF1 = False
		self.FlapDOF2 = False
		self.EdgeDOF = False
		self.TeetDOF = False
		self.DrTrDOF = False
		self.GenDOF = False
		self.YawDOF = False
		self.TwFADOF1 = False
		self.TwFADOF2 = False
		self.TwSSDOF1 = False
		self.TwSSDOF2 = False
		self.PtfmSgDOF = False
		self.PtfmSwDOF = False
		self.PtfmHvDOF = False
		self.PtfmRDOF = False
		self.PtfmPDOF = False
		self.PtfmYDOF = False

# Initial Conditions
class InitConds(object):
	def __init__(self):
		self.OoPDefl    = 0.0
		self.IPDefl     = 0.0
		self.BlPitch1   = 0.0
		self.BlPitch2   = 0.0
		self.BlPitch3   = 0.0
		self.TeetDefl   = 0.0
		self.Azimuth    = 0.0
		self.RotSpeed   = 0.0
		self.NacYaw     = 0.0
		self.TTDspFA    = 0.0
		self.TTDspSS    = 0.0
		self.PtfmSurge  = 0.0
		self.PtfmSway   = 0.0
		self.PtfmHeave  = 0.0
		self.PtfmRoll   = 0.0
		self.PtfmPitch  = 0.0
		self.PtfmYaw    = 0.0

# Turbine Configuration
class TurbConfig(object):
	def __init__(self):
		self.NumBl      = 0
		self.TipRad     = 0.0
		self.HubRad     = 0.0
		self.PreCone1   = 0.0
		self.PreCone2   = 0.0
		self.PreCone3   = 0.0
		self.HubCM      = 0.0
		self.UndSling   = 0.0
		self.Delta3     = 0.0
		self.AzimB1Up   = 0.0
		self.OverHang   = 0.0
		self.ShftGagL   = 0.0
		self.ShftTilt   = 0.0
		self.NacCMxn    = 0.0
		self.NacCMyn    = 0.0
		self.NacCMzn    = 0.0
		self.NcIMUxn    = 0.0
		self.NcIMUyn    = 0.0
		self.NcIMUzn    = 0.0
		self.Twr2Shft   = 0.0
		self.TowerHt    = 0.0
		self.TowerBsHt  = 0.0
		self.PtfmCMxt   = 0.0
		self.PtfmCMyt   = 0.0
		self.PtfmCMzt   = 0.0
		self.PtfmRefzt  = 0.0

# Mass and Inertia
class MassInertia(object):
	def __init__(self):
		self.TipMass1   = 0.0
		self.TipMass2   = 0.0
		self.TipMass3   = 0.0
		self.HubMass    = 0.0
		self.HubIner    = 0.0
		self.GenIner    = 0.0
		self.NacMass    = 0.0
		self.NacYIner   = 0.0
		self.YawBrMass  = 0.0
		self.PtfmMass   = 0.0
		self.PtfmRIner  = 0.0
		self.PtfmPIner  = 0.0
		self.PtfmYIner  = 0.0

# ED Blade (Structure)
class BladeStruc(object):
	def __init__(self):
		self.BldNodes = 0
		self.BldFile1 = ''
		self.BldFile2 = ''
		self.BldFile3 = ''

		# Including the blade files and properties in the same object,
		# as is done here, implies that the properties are done for all
		# blades (assumed for now)

		# General Model Inputs
		self.NBlInpSt = 0   #Number of blade input stations (-)
		self.CalcBMode = False   #Calculate blade mode shapes internally {T: ignore mode shapes from below, F: use mode shapes from below} [CURRENTLY IGNORED] (flag)
		self.BldFlDmp1 = 0.0   #Blade flap mode #1 structural damping in percent of critical (%)
		self.BldFlDmp2 = 0.0   #Blade flap mode #2 structural damping in percent of critical (%)
		self.BldEdDmp1 = 0.0   #Blade edge mode #1 structural damping in percent of critical (%)
		self.FlStTunr1 = 0.0   #Blade flapwise modal stiffness tuner, 1st mode (-)
		self.FlStTunr2 = 0.0   #Blade flapwise modal stiffness tuner, 2nd mode (-)
		self.AdjBlMs = 0.0   #Factor to adjust blade mass density (-)
		self.AdjFlSt = 0.0   #Factor to adjust blade flap stiffness (-)
		self.AdjEdSt = 0.0   #Factor to adjust blade edge stiffness (-)
		
		# Distributed Blade Properties
		self.BlFract = zeros([1])
		self.AeroCent = zeros([1])
		self.PitchAxis = zeros([1])
		self.StrcTwst = zeros([1])
		self.BMassDen = zeros([1])
		self.FlpStff = zeros([1])
		self.EdgStff = zeros([1])
		self.GJStff = zeros([1])
		self.EAStff = zeros([1])
		self.Alpha = zeros([1])
		self.FlpIner = zeros([1])
		self.EdgIner = zeros([1])
		self.PrecrvRef = zeros([1])
		self.PreswpRef = zeros([1]) #[AH] Added during openmdao1 update
		self.FlpcgOf = zeros([1])
		self.Edgcgof = zeros([1])
		self.FlpEAOf = zeros([1])
		self.EdgEAOf = zeros([1])
		
		# Blade Mode Shapes
		self.BldFl1Sh = zeros([1])
		self.BldFl2Sh = zeros([1])
		self.BldEdgSh = zeros([1])

# Rotor-Teeter
class RotorTeeter(object):
	def __init__(self):
		self.TeetMod  = 0
		self.TeetDmpP = 0.0
		self.TeetDmp  = 0.0
		self.TeetCDmp = 0.0
		self.TeetSStP = 0.0
		self.TeetHStP = 0.0
		self.TeetSSSp = 0.0
		self.TeetHSSp = 0.0

class DriveTrain(object):
	def __init__(self):
		self.GBoxEff  = 0.0
		self.GBRatio  = 0.0
		self.DTTorSpr = 0.0
		self.DTTorDmp = 0.0
		self.GBRevers = 0.0   #FAST7 only
		self.DynBrkFi = 0.0   #FAST7 only

class Furling(object):
	def __init__(self):
		self.Furling = False
		self.FurlFile = ''

class Tower(object):
	def __init__(self):
		self.TwrNodes = 0
		self.TwrFile = ''

		# General Tower Parameters
		self.NTwInptSt = 0   #Number of input stations to specify tower geometry
		self.CalcTMode = False   #calculate tower mode shapes internally {T: ignore mode shapes from below, F: use mode shapes from below} [CURRENTLY IGNORED] (flag)
		self.TwrFADmp1 = 0.0   #Tower 1st fore-aft mode structural damping ratio (%)
		self.TwrFADmp2 = 0.0   #Tower 2nd fore-aft mode structural damping ratio (%)
		self.TwrSSDmp1 = 0.0   #Tower 1st side-to-side mode structural damping ratio (%)
		self.TwrSSDmp2 = 0.0   #Tower 2nd side-to-side mode structural damping ratio (%)

		# Tower Adjustment Factors
		self.FAStTunr1 = 0.0   #Tower fore-aft modal stiffness tuner, 1st mode (-)
		self.FAStTunr2 = 0.0   #Tower fore-aft modal stiffness tuner, 2nd mode (-)
		self.SSStTunr1 = 0.0   #Tower side-to-side stiffness tuner, 1st mode (-)
		self.SSStTunr2 = 0.0   #Tower side-to-side stiffness tuner, 2nd mode (-)
		self.AdjTwMa = 0.0   #Factor to adjust tower mass density (-)
		self.AdjFASt = 0.0   #Factor to adjust tower fore-aft stiffness (-)
		self.AdjSSSt = 0.0   #Factor to adjust tower side-to-side stiffness (-)
	 
		# Distributed Tower Properties   
		self.HtFract = zeros([1])
		self.TMassDen = zeros([1])
		self.TwFAStif = zeros([1])
		self.TwSSStif = zeros([1])
		self.TwGJStif = zeros([1])
		self.TwEAStif = zeros([1])
		self.TwFAIner = zeros([1])
		self.TwSSIner = zeros([1])
		self.TwFAcgOf = zeros([1])
		self.TwSScgOf = zeros([1])
		
		# Tower Mode Shapes
		self.TwFAM1Sh = zeros([1])   #Tower Fore-Aft Mode 1 Shape Coefficients x^2, x^3, x^4, x^5, x^6
		self.TwFAM2Sh = zeros([1])   #Tower Fore-Aft Mode 2 Shape Coefficients x^2, x^3, x^4, x^5, x^6
		self.TwSSM1Sh = zeros([1])   #Tower Side-to-Side Mode 1 Shape Coefficients x^2, x^3, x^4, x^5, x^6
		self.TwSSM2Sh = zeros([1])   #Tower Side-to-Side Mode 2 Shape Coefficients x^2, x^3, x^4, x^5, x^6

# Currently FAST7 only
class Platform(object):
	def __init__(self):
		self.PtfmModel = 0
		self.PtfmFile = ''

class EdOutParams(object):
	def __init__(self):
		self.SumPrint = False
		self.OutFile  = 0
		self.TabDelim = False
		self.OutFmt   = ''
		self.TStart   = 0.0
		self.DecFact  = 0.0
		self.NTwGages = 0
		self.TwrGagNd = []
		self.NBlGages = 0
		self.BldGagNd = []

# Inflow Wind General Parameters
class InflowWind(object):
	def __init__(self):
		self.Echo = False
		self.WindType = 0
		self.PropagationDir = 0.0
		self.NWindVel = 0
		self.WindVxiList = 0.0
		self.WindVyiList = 0.0
		self.WindVziList = 0.0

# Parameters for Steady Wind Conditions [used only for WindType = 1]
class SteadyWindParams(object):
	def __init__(self):
		self.HWindSpeed = 0.0
		self.RefHt = 0.0
		self.PLexp = 0.0

# Parameters for Uniform wind file   [used only for WindType = 2]
class UniformWindParams(object):
	def __init__(self):
		self.Filename = ''
		self.RefHt = 0.0
		self.RefLength = 0.0

# Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3]
class TurbSimWindParams(object):
	def __init__(self):
		self.Filename = ''

# Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4]
class BladedWindParams(object):
	def __init__(self):
		self.FilenameRoot = ''
		self.TowerFile = False

# Parameters for HAWC-format binary files  [Only used with WindType = 5]
class HAWCWindParams(object):
	def __init__(self):
		self.FileName_u  = ''
		self.FileName_v  = ''
		self.FileName_w  = ''
		self.nx          = 0
		self.ny          = 0
		self.nz          = 0
		self.dx          = 0.0
		self.dy          = 0.0
		self.dz          = 0.0
		self.RefHt       = 0.0
		self.ScaleMethod = 0
		self.SFx         = 0.0
		self.SFy         = 0.0
		self.SFz         = 0.0
		self.SigmaFx     = 0.0
		self.SigmaFy     = 0.0
		self.SigmaFz     = 0.0
		self.URef        = 0.0
		self.WindProfile = 0
		self.PLExp       = 0.0
		self.Z0          = 0.0

# Inflow Wind Output Parameters (actual OutList included in master OutList)
class InflowOutParams(object):
	def __init__(self):
		self.SumPrint = False

# Wnd Wind File Parameters
class WndWind(object):
	def __init__(self):
		self.TimeSteps = 0   #number of time steps
		self.Time = zeros([1])   #time steps
		self.HorSpd = zeros([1])   #horizontal wind speed
		self.WindDir = zeros([1])   #wind direction
		self.VerSpd = zeros([1])   #vertical wind speed
		self.HorShr = zeros([1])   #horizontal shear
		self.VerShr = zeros([1])   #vertical power-law shear
		self.LnVShr = zeros([1])   #vertical linear shear
		self.GstSpd = zeros([1])   #gust speed not sheared by Aerodyn

# AeroDyn Parameters
class AeroDyn(object):
	def __init__(self):
		# General Model Inputs
		self.SysUnits = 'SI'   #Enum('SI', ('SI','ENG, desc='System of units for used for input and output [must be SI for FAST] (unquoted string)
		self.StallMod = 'BEDOES'   #Enum('BEDDOES', ('BEDDOES', 'STEADY, desc = 'Dynamic stall included [BEDDOES or STEADY] (unquoted string)
		self.UseCM = 'NO_CM'   #Enum('NO_CM', ('NO_CM', 'USE_CM, desc = 'Use aerodynamic pitching moment model? [USE_CM or NO_CM] (unquoted string)
		self.InfModel = 'EQUIL'   #Enum('EQUIL', ('EQUIL', 'DYNIN, desc = 'Inflow model [DYNIN or EQUIL] (unquoted string)
		self.IndModel = 'SWIRL'   #Enum('SWIRL', ('NONE', 'WAKE', 'SWIRL, desc = 'Induction-factor model [NONE or WAKE or SWIRL] (unquoted string)
		self.AToler = 0.0   #Induction-factor tolerance (convergence criteria) (-)
		self.TLModel = 'PRANDtl'   #Enum('PRANDtl', ('PRANDtl', 'GTECH', 'NONE, desc = 'Tip-loss model (EQUIL only) [PRANDtl, GTECH, or NONE] (unquoted string)
		self.HLModel = 'PRANDTl'   #Enum('PRANDtl', ('PRANDtl', 'GTECH', 'NONE, desc = 'Hub-loss model (EQUIL only) [PRANdtl or NONE] (unquoted string)

		# Turbine Inputs
		self.WindFile = ''   #Initial wind file from template import (quoted string)
		self.wind_file_type = 'hh'   #Enum('hh', ('hh', 'bts', 'wnd, desc='type of wind file
		self.HH = 0.0   #units='m', desc= 'Wind reference (hub) height [TowerHt+Twr2Shft+OverHang*SIN(ShftTilt)] (m)
		self.TwrShad = 0.0   #Tower-shadow velocity deficit (-)
		self.ShadHWid = 9999.9   #units='m', desc='Tower-shadow half width (m)'
		self.T_Shad_Refpt = 9999.9   #units='m', desc='Tower-shadow reference point (m)

		# Wind Aero Inputs
		self.AirDens = 1.225   #units='kg / (m**3)', desc='Air density (kg/m^3)
		self.KinVisc = 1.464e-5   #units='m**2 / s', desc='Kinematic air viscosity [CURRENTLY IGNORED] (m^2/sec)
		self.DTAero = 0.02479   #units='s', desc = 'Time interval for aerodynamic calculations (sec)
		

# AeroDyn Blade
class AeroDynBlade(object):
	def __init__(self):
		self.NumFoil = 0   #Number of airfoil files (-)
		self.FoilNm = zeros([1])   #Names of the airfoil files [NumFoil lines] (quoted strings)

		self.af_data = []   #list of airfoild data sets

		self.BldNodes = 0   #Number of blade nodes used for analysis
		self.RNodes = zeros([1])   #Distance from blade root to center of element
		self.AeroTwst = zeros([1])   #units='deg',desc='Twist at RNodes locations
		self.DRNodes = zeros([1])   #Span-wise width of the blade element, measured along the span of the blade
		self.Chord = zeros([1])   #units='m',desc='Chord length at node locations; planform area = chord*DRNodes
		self.NFoil = zeros([1])   #Airfoil ID Number
		self.PrnElm = zeros([1])   #Flag for printing element ouput for blade section

# AeroDyn Airfoil Polar
class ADAirfoilPolar(object):
	def __init__(self):
	    self.IDParam = 0.0   #Table ID Parameter (Typically Reynolds number)
	    self.StallAngle = 0.0   #Stall angle (deg)
	    self.ZeroCn = 0.0   #Zero lift angle of attack (deg)
	    self.CnSlope = 0.0   #Cn slope for zero lift (dimensionless)
	    self.CnPosStall = 0.0   #Cn at stall value for positive angle of attack
	    self.CnNegStall = 0.0   #Cn at stall value for negative angle of attack
	    self.alphaCdMin = 0.0   #Angle of attack for minimum CD (deg)
	    self.CdMin = 0.0   #Minimum Cd Value

	    self.alpha = zeros([1])   #angle of attack
	    self.cl = zeros([1])   #coefficient of lift
	    self.cd = zeros([1])   #coefficient of drag
	    self.cm = zeros([1])   #coefficient of the pitching moment

# AeroDyn airfoil
class ADAirfoil(object):
	def __init__(self):
	    self.description = ''   #description of airfoil
	    self.number_tables = 0   #number of airfoil polars
	    self.af_tables = []   #list of airfoil polars
	    
# AeroDyn15 Parameters
class AeroDyn15(object):
	def __init__(self):
		# General Options
		self.Echo = False #Echo the input to "<rootname>.AD.ech"?  (flag)
		self.DTAero = ''    #Time interval for aerodynamic calculations {or "default"} (s)
		self.WakeMod = 0   #Type of wake/induction model (switch) {0=none, 1=BEMT}
		self.AFAeroMod = 0   #Type of blade airfoil aerodynamics model (switch) {1=steady model, 2=Beddoes-Leishman unsteady model}
		self.TwrPotent = 0   #Type tower influence on wind based on potential flow around the tower (switch) {0=none, 1=baseline potential flow, 2=potential flow with Bak correction}
		self.TwrShadow = False  #Calculate tower influence on wind based on downstream tower shadow? (flag)
		self.TwrAero = False   #Calculate tower aerodynamic loads? (flag)
		self.FrozenWake = False   #Assume frozen wake during linearization? (flag) [used only when WakeMod=1 and when linearizing]
		# openfast
		self.CavitCheck = False   #Perform cavitation check? (flag) TRUE will turn off unsteady aerodynamics
       
		# Environmental Conditions 
		self.AirDens = 0.0   #Air density (kg/m^3)
		self.KinVisc = 0.0   #Kinematic air viscosity (m^2/s)
		self.SpdSound = 0.0   #Speed of sound (m/s)
		# openfast
		self.Patm = 0.0   # Atmospheric pressure (Pa) [used only when CavitCheck=True]	
		self.Pvap = 0.0   # Vapour pressure of fluid (Pa) [used only when CavitCheck=True]
		self.FluidDepth = 0.0   # Water depth above mid-hub height (m) [used only when CavitCheck=True]		
		
		# Blade-Element/Momentum Theory Options
		self.SkewMod = 0			# Type of skewed-wake correction model (switch) {1=uncoupled, 2=Pitt/Peters, 3=coupled} [used only when WakeMod=1]
		
		#openfast
		self.SkewModFactor = ''       # Constant used in Pitt/Peters skewed wake model {or "default" is 15/32*pi} (-) [used only when SkewMod=2; unused when WakeMod=0]
		
		self.TipLoss = False			# Use the Prandtl tip-loss model? (flag) [used only when WakeMod=1]
		self.HubLoss = False			# Use the Prandtl hub-loss model? (flag) [used only when WakeMod=1]
		self.TanInd = False			 # Include tangential induction in BEMT calculations? (flag) [used only when WakeMod=1]
		self.AIDrag = False			 # Include the drag term in the axial-induction calculation? (flag) [used only when WakeMod=1]
		self.TIDrag = False			 # Include the drag term in the tangential-induction calculation? (flag) [used only when WakeMod=1 and TanInd=TRUE]
		self.IndToler = ''		   # Convergence tolerance for BEMT nonlinear solve residual equation {or "default"} (-) [used only when WakeMod=1]
		self.MaxIter = 0			# Maximum number of iteration steps (-) [used only when WakeMod=1]
		# openfast - Dynamic Blade-Element/Momentum Theory Options
		self.DBEMT_Mod = 0         # Type of dynamic BEMT (DBEMT) model {1=constant tau1, 2=time-dependent tau1} (-) [used only when WakeMod=2]
		self.tau1_const = 0.0      # Time constant for DBEMT (s) [used only when WakeMod=2 and DBEMT_Mod=1]
		# Beddoes-Leishman Unsteady Airfoil Aerodynamics Options
		self.UAmod = 0             # Unsteady Aero Model Switch (switch) {1=Baseline model (Original), 2=Gonzalez variant (changes in Cn,Cc,Cm), 3=Minemma/Pierce variant (changes in Cc and Cm)} [used only when AFAeroMod=2]
		self.FLookup = False    # Flag to indicate whether a lookup for f' will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files (flag) [used only when AFAeroMod=2]
		
		# Airfoil Information
		self.InCol_Alfa	= 0	 # The column in the airfoil tables that contains the angle of attack (-)
		self.InCol_Cl = 0 		   # The column in the airfoil tables that contains the lift coefficient (-)
		self.InCol_Cd = 0		   # The column in the airfoil tables that contains the drag coefficient (-)
		self.InCol_Cm = 0		   # The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)
		self.InCol_Cpmin = 0		# The column in the airfoil tables that contains the Cpmin coefficient; use zero if there is no Cpmin column (-)

		# Rotor/Blade Properties
		self.UseBlCm = False   # Include aerodynamic pitching moment in calculations?  (flag)
		self.ADBlFile1 = ''   # Name of file containing distributed aerodynamic properties for Blade #1 (-)
		self.ADBlFile2 = ''   # Name of file containing distributed aerodynamic properties for Blade #2 (-) [unused if NumBl < 2]
		self.ADBlFile3 = ''   # Name of file containing distributed aerodynamic properties for Blade #3 (-) [unused if NumBl < 3]
		
# AeroDyn15 Blade
class AeroDyn15Blade(object):
	def __init__(self):
		self.NumAFfiles	= 0		 # Number of airfoil files used (-)
		self.AFNames = zeros([1])   # Airfoil file names (NumAFfiles lines) (quoted strings)

		self.af_data = []   #list of airfoild data sets
		
		
class AeroDyn15Tower(object):
	def __init__(self):
		self.NumTwrNds  = 0   # Number of tower nodes used in the analysis  (-) [used only when TwrPotent/=0, TwrShadow=True, or TwrAero=True]
		self.TwrElev = zeros([1])   #
		self.TwrDiam = zeros([1])   #
		self.TwrCd = zeros([1])   #
		

# AeroDyn airfoil
class AD15Airfoil(object):
	def __init__(self):
	    self.description = ''   #description of airfoil
	    self.InterpOrd = ''	# Interpolation order to use for quasi-steady table lookup {1=linear; 3=cubic spline; "default"} [default=3]						
	    self.NonDimArea = 0.0	# The non-dimensional area of the airfoil (area/chord^2) (set to 1.0 if unsure or unneeded)						
	    self.NumCoords = 0         # The number of coordinates in the airfoil shape file.  Set to zero if coordinates not included.						
	    self.NumTabs = 0	# Number of airfoil tables in this file.  Each table must have lines for Re and Ctrl.
	    self.af_tables = []   #list of airfoil polars

	    
# AeroDyn Airfoil Polar
class AD15AirfoilPolar(object):
	def __init__(self):			
		self.Re = 0.0	# Reynolds number in millions						
		self.Ctrl = 0	# Control setting (must be 0 for current AirfoilInfo)						
		self.InclUAdata = False	# Is unsteady aerodynamics data included in this table? If TRUE, then include 30 UA coefficients below this line											
		self.alpha0 = 0.0 # 0-lift angle of attack, depends on airfoil.						
		self.alpha1 = 0.0 	# Angle of attack at f=0.7, (approximately the stall angle) for AOA>alpha0. (deg)						
		self.alpha2 = 0.0 	# Angle of attack at f=0.7, (approximately the stall angle) for AOA<alpha0. (deg)						
		self.eta_e = 0.0 	# Recovery factor in the range [0.85 - 0.95] used only for UAMOD=1, it is set to 1 in the code when flookup=True.	(-)					
		self.C_nalpha = 0.0 	# Slope of the 2D normal force coefficient curve. (1/rad)						
		self.T_f0 = 0.0 	# Initial value of the time constant associated with Df in the expression of Df and f''. [default = 3]						
		self.T_V0 = 0.0 	# Initial value of the time constant associated with the vortex lift decay process; it is used in the expression	of C	vn.	It	depends on	Re,	M, and a
		self.T_p = 0.0 	# Boundary-layer,leading edge pressure gradient time constant in the expression of Dp. It should be tuned based o	n ai	rfoi	l e	xperimenta	l da	ta. [def
		self.T_VL = 0.0 	# Initial value of the time constant associated with the vortex advection process; it represents the non-dimensio	nal	time	in	semi-chor	ds,	needed f
		self.b1 = ''	# Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils,	but	may	be	different	for	turbine
		self.b2 = ''	# Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils,	but	may	be	different	for	turbine
		self.b5 = 0.0	# Constant in the expression of K'''_q,Cm_q^nc, and k_m,q.  [from  experimental results, defaults to 5]						
		self.A1 = ''	# Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils,	but	may	be	different	for	turbine
		self.A2 = ''	# Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils,	but	may	be	different	for	turbine
		self.A5 = ''	# Constant in the expression of K'''_q,Cm_q^nc, and k_m,q. [from experimental results, defaults to 1]						
		self.S1 = 0.0	# Constant in the f curve best-fit for alpha0<=AOA<=alpha1; by definition it depends on the airfoil. [ignored if	UAMo	d<>1	]			
		self.S2 = 0.0	# Constant in the f curve best-fit for         AOA> alpha1; by definition it depends on the airfoil. [ignored if	UAMo	d<>1	]			
		self.S3 = 0.0	# Constant in the f curve best-fit for alpha2<=AOA< alpha0; by definition it depends on the airfoil. [ignored if	UAMo	d<>1	]			
		self.S4 = 0.0	# Constant in the f curve best-fit for         AOA< alpha2; by definition it depends on the airfoil. [ignored if	UAMo	d<>1	]			
		self.Cn1 = 0.0	# Critical value of C0n at leading edge separation. It should be extracted from airfoil data at a given Mach and	Reyn	olds	nu	mber. It c	an b	e calcul
		self.Cn2 = 0.0	# As Cn1 for negative AOAs.						
		self.St_sh = 0.0 	# Strouhal's shedding frequency constant.  [default = 0.19]						
		self.Cd0 = 0.0	# 2D drag coefficient value at 0-lift.						
		self.Cm0 = 0.0	# 2D pitching moment coefficient about 1/4-chord location, at 0-lift, positive if nose up. [If the aerodynamics c	oeff	icie	nts	table doe	s no	t includ
		self.k0 = 0.0	# Constant in the \hat(x)_cp curve best-fit; = (\hat(x)_AC-0.25).  [ignored if UAMod<>1]						
		self.k1 = 0.0	# Constant in the \hat(x)_cp curve best-fit.  [ignored if UAMod<>1]						
		self.k2 = 0.0	# Constant in the \hat(x)_cp curve best-fit.  [ignored if UAMod<>1]						
		self.k3 = 0.0	# Constant in the \hat(x)_cp curve best-fit.  [ignored if UAMod<>1]						
		self.k1_hat = 0.0	# Constant in the expression of Cc due to leading edge vortex effects.  [ignored if UAMod<>1]						
		self.x_cp_bar = ''	# Constant in the expression of \hat(x)_cp^v. [ignored if UAMod<>1, default = 0.2]						
		self.UACutout = ''	# Angle of attack above which unsteady aerodynamics are disabled (deg). [Specifying the string "Default" sets UAC	utou	t to	45	degrees]		
		self.filtCutOff = ''	# Cut-off frequency (-3 dB corner frequency) for low-pass filtering the AoA input to UA, as well as the 1st and 2	nd d	eriv	ati	ves (Hz) [	defa	ult = 20
		self.Alpha = zeros([1])   #angle of attack
		self.Cl = zeros([1])   #coefficient of lift
		self.Cd = zeros([1])   #coefficient of dragself.Cm = zeros([1])   #coefficient of the pitching moment	 
	    
class AD15OutParams(object):
	def __init__(self):
		self.SumPrint = False
		self.NBlOuts  = 0
		self.BlOutNd  = []
		self.NTwOuts  = 0
		self.TwOutNd  = []


# ServoDyn Simulation Control
class SdSimCtrl(object):
	def __init__(self):
		self.Echo = False
		self.DT = 0.0

# Pitch Control
class PitchCtrl(object):
	def __init__(self):
		self.PCMode       = 0
		self.TPCOn        = 0.0
		self.TPitManS1    = 0.0
		self.TPitManS2    = 0.0
		self.TPitManS3    = 0.0
		self.TPitManE1    = 0.0   #FAST7 only
		self.TPitManE2    = 0.0   #FAST7 only
		self.TPitManE3    = 0.0   #FAST7 only
		self.PitManRat1   = 0.0
		self.PitManRat2   = 0.0
		self.PitManRat3   = 0.0
		self.BlPitchF1    = 0.0
		self.BlPitchF2    = 0.0
		self.BlPitchF3    = 0.0
		self.BlPitch1     = 0.0   #FAST7 only
		self.BlPitch2     = 0.0   #FAST7 only
		self.BlPitch3     = 0.0   #FAST7 only

# Generator and Torque Control
class GenTorqCtrl(object):
	def __init__(self):
		self.VSContrl = 0
		self.GenModel = 0
		self.GenEff   = 0.0
		self.GenTiStr = False
		self.GenTiStp = False
		self.SpdGenOn = 0.0
		self.TimGenOn = 0.0
		self.TimGenOf = 0.0

# Simple Variable-Speed Torque Control
class VarSpeedTorqCtrl(object):
	def __init__(self):
		self.VS_RtGnSp = 0.0
		self.VS_RtTq   = 0.0
		self.VS_Rgn2K  = 0.0
		self.VS_SlPc   = 0.0

# Simple Induction Generator
class InductGen(object):
	def __init__(self):
		self.SIG_SlPc = 0.0
		self.SIG_SySp = 0.0
		self.SIG_RtTq = 0.0
		self.SIG_PORt = 0.0

# Thevenin-Equivalent Induction Generator
class ThevEqInductGen(object):
	def __init__(self):
		self.TEC_Freq = 0.0
		self.TEC_NPol = 0
		self.TEC_SRes = 0.0
		self.TEC_RRes = 0.0
		self.TEC_VLL  = 0.0
		self.TEC_SLR  = 0.0
		self.TEC_RLR  = 0.0
		self.TEC_MR   = 0.0

# High-Speed Shaft Brake
class ShaftBrake(object):
	def __init__(self):
		self.HSSBrMode = 0
		self.THSSBrDp  = 0.0
		self.HSSBrDT   = 0.0
		self.HSSBrTqF  = 0.0

# Nacelle-Yaw Control
class NacYawCtrl(object):
	def __init__(self):
		self.YCMode    = 0
		self.TYCOn     = 0.0
		self.YawNeut   = 0.0
		self.YawSpr    = 0.0
		self.YawDamp   = 0.0
		self.TYawManS  = 0.0
		self.YawManRat = 0.0
		self.NacYawF   = 0.0

# Tip Brake (used in FAST7 only)
class TipBrake(object):
	def __init__(self):
		self.TiDynBrk  = 0.0
		self.TTpBrDp1  = 0.0
		self.TTpBrDp2  = 0.0
		self.TTpBrDp3  = 0.0
		self.TBDepISp1 = 0.0
		self.TBDepISp2 = 0.0
		self.TBDepISp3 = 0.0
		self.TBDrConN = 0.0
		self.TBDrConD = 0.0
		self.TpBrDT   = 0.0

# Tuned Mass Damper
class TunedMassDamp(object):
	def __init__(self):
		self.CompNTMD = False
		self.NTMDfile = ''
		self.CompTTMD = False
		self.TTMDfile = ''

# Bladed Interface
class BladedInterface(object):
	def __init__(self):
		self.DLL_FileName = ''
		self.DLL_InFile   = ''
		self.DLL_ProcName = ''
		self.DLL_DT       = ''
		self.DLL_Ramp     = False
		self.BPCutoff     = 0.0
		self.NacYaw_North = 0.0
		self.Ptch_Cntrl   = 0.0
		self.Ptch_SetPnt  = 0.0
		self.Ptch_Min     = 0.0
		self.Ptch_Max     = 0.0
		self.PtchRate_Min = 0.0
		self.PtchRate_Max = 0.0
		self.Gain_OM      = 0.0
		self.GenSpd_MinOM = 0.0
		self.GenSpd_MaxOM = 0.0
		self.GenSpd_Dem   = 0.0
		self.GenTrq_Dem   = 0.0
		self.GenPwr_Dem   = 0.0
		self.DLL_NumTrq = 0.0
		self.GenSpd_TLU = zeros([0])
		self.GenTrq_TLU = zeros([0])

# ServoDyn Output Params
class SdOutParams(object):
	def __init__(self):
		self.SumPrint = False
		self.OutFile  = 0
		self.TabDelim = False
		self.OutFmt   = ''
		self.TStart   = 0.0


# ====== INITIALIZE FAST MODEL BY INITIALIZING ALL VARIABLE TREES ======

class FstModel(object):
	def __init__(self, ):

		# Description
		self.description = ''

		# Fst file vartrees
		self.fst_sim_ctrl = FstSimCtrl()
		self.ftr_swtchs_flgs = FtrSwtchsFlgs()
		self.input_files = InputFiles()
		self.fst_output_params = FstOutputParams()
		self.visualization = Visualization()
		self.linearization = Linearization()
		self.fst_out_params = FstOutputParams()
		
		# Elastodyn vartrees
		self.ed_sim_ctrl = EdSimCtrl()
		self.envir_cond = EnvirCond()
		self.dof = DOF()
		self.init_conds = InitConds()
		self.turb_config = TurbConfig()
		self.mass_inertia = MassInertia()
		self.blade_struc = BladeStruc()
		self.rotor_teeter = RotorTeeter()
		self.drivetrain = DriveTrain()
		self.furling = Furling()
		self.platform = Platform()
		self.tower = Tower()
		self.ed_out_params = EdOutParams()

		# Wind vartrees
		self.inflow_wind = InflowWind()
		self.steady_wind_params = SteadyWindParams()
		self.uniform_wind_params = UniformWindParams()
		self.turbsim_wind_params = TurbSimWindParams()
		self.bladed_wind_params = BladedWindParams()
		self.hawc_wind_params = HAWCWindParams()
		self.inflow_out_params = InflowOutParams()
		self.wnd_wind = WndWind()

		# AeroDyn vartrees
		self.aerodyn = AeroDyn()
		self.blade_aero = AeroDynBlade()
		self.aerodyn15 = AeroDyn15()
		self.blade_aero15 = AeroDyn15Blade()
		self.tower_aero15 = AeroDyn15Tower()
		self.ad15_out_params = AD15OutParams()
			
		# ServoDyn vartrees
		self.sd_sim_ctrl = SdSimCtrl()
		self.pitch_ctrl = PitchCtrl()
		self.gen_torq_ctrl = GenTorqCtrl()
		self.var_speed_torq_ctrl = VarSpeedTorqCtrl()
		self.induct_gen = InductGen()
		self.theveq_induct_gen = ThevEqInductGen()
		self.shaft_brake = ShaftBrake()
		self.nac_yaw_ctrl = NacYawCtrl()
		self.tip_brake = TipBrake()
		self.tuned_mass_damper = TunedMassDamp()
		self.bladed_interface = BladedInterface()
		self.sd_out_params = SdOutParams()
		
		# List of Outputs (all input files -- FST, ED, SD)
		# TODO: Update FstOutput for a few new outputs in FAST8
		self.outlist = FstOutput()   # 











