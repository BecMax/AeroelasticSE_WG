from numpy import zeros, array
import numpy as np

# This variable tree contains all parameters required to create a FAST model
# for FAST versions 7 and 8.


##########################
### Main input file ######
##########################
# Runtime Options
class runtime_options(object):
    def __init__(self):
         self.echo = False
         self.RandSeed1 = 0
         self.RandSeed2 = ''
         self.WrBHHTP = False
         self.WrFHHTP = False
         self.WrADHH = False
         self.WrADFF = False
         self.WrBLFF = False
         self.WrADTWR = False
         self.WrFMTFF = False
         self.WrACT = False
         self.Clockwise = False
         self.ScaleIEC = 0

# Turbine/Model Specifications
class tmspecs(object):
    def __init__(self):
         self.NumGrid_Z = 0
         self.NumGrid_Y = 0
         self.TimeStep = 0.0
         self.AnalysisTime = 0.0
         self.UsableTime = ''
         self.HubHt = 0.0
         self.GridHeight = 0.0
         self.GridWidth = 0.0
         self.VFlowAng = 0.0
         self.HFlowAng = 0.0

# Meteorological Boundary Conditions
class metboundconds(object):
    def __init__(self):
         self.TurbModel = ''
         self.UserFile = ''
         self.IECstandard = ''
         self.IECturbc = ''
         self.IEC_WindType = ''
         self.ETMc = ''
         self.WindProfileType = ''
         self.ProfileFile = ''
         self.RefHt = 0.0
         self.URef = 0.0
         self.ZJetMax = ''
         self.PLExp = ''
         self.Z0 = ''

# Non-IEC Meteorological Boundary Conditions
class noniecboundconds(object):
    def __init__(self):
         self.Latitude = ''
         self.RICH_NO = 0.0
         self.UStar = ''
         self.ZI = ''
         self.PC_UW = 0.0
         self.PC_UV = 0.0
         self.PC_VW = 0.0

# Spatial Coherence Parameters
class spatialcoherance(object):
    def __init__(self):
         self.SCMod1 = ''
         self.SCMod2 = ''
         self.SCMod3 = ''
         self.InCDec1 = [0.0, 0.0]
         self.InCDec2 = [0.0, 0.0]
         self.InCDec3 = [0.0, 0.0]
         self.CohExp = 0.0

# Coherent Turbulence Scaling Parameters
class coherentTurbulence(object):
    def __init__(self):
         self.CTEventPath = ''
         self.CTEventFile = ''
         self.Randomize = False
         self.DistScl = 0.0
         self.CTLy = 0.0
         self.CTLz = 0.0
         self.CTStartTime = 0.0

class turbsiminputs(object):
    def __init__(self):
         self.runtime_options = runtime_options()
         self.tmspecs = tmspecs()
         self.metboundconds = metboundconds()
         self.noniecboundconds = noniecboundconds()
         self.spatialcoherance = spatialcoherance()
         self.coherentTurbulence = coherentTurbulence()

##########################
### Turbulance input file ######
##########################

