from turbsim_vartrees import turbsiminputs
from turbulence_spectrum import turb_specs
from wind_profile_writer import write_wind
import os
import numpy as np
class TurbsimBuilder(turbsiminputs):
    def __init__(self):
        self.turbsim_vt = turbsiminputs()
        self.Turbsim_runfile = ''
        self.Turbsim_rundir = ''
         #self.tsim_turbulence_file = 'turbulence_default.in'
         #self.tsim_profile_file = 'default.shear'
 
         # Turbulence file parameters
         #self.wind_speed = 8.
         #self.L_u = 2.54e+02
         #self.L_v=1.635e+02
         #self.L_w=4.7e+01
         ##self.sigma_u=1.325
         #self.sigma_v=0.9
         #self.sigma_w=0.7625
         #self.turbulence_file_name = 'tsim_user_turbulence_default.inp'
         #self.turbulence_template_file = 'TurbsimInputFiles/turbulence_user.inp'
         
         # profile file parameters
         #self.profile_template = 'TurbsimInputFiles/shear.profile'
         #self.shear_exponent = 0.7
         #self.veer = 12.5
         #self.turbsim_vt.metboundconds.ProfileFile = 'default.profile'

         #self.run_dir = 'run%d'%np.random.uniform(0,1e10)
         
    def InputConfig(self, **kwargs):
        # for k, w in kwargs.iteritems():
            # try:
            #     success = False
            #     if hasattr(self, k):
            #         setattr(self,k,w)
            #         success = True
            #     else:
            #         # Currently only can assign variables at the first level (fst_vt.subtree.variable)
            #         # This could be re-written to do this recursively
            #         data = k.split('.')   #split input into keys
            #         var_tree = data[-2]
            #         variable = data[-1]
            #         for key, val in self.fst_vt.__dict__.iteritems():
            #             if key == var_tree:
            #                 setattr(self.fst_vt.__dict__[key],variable,w)
            #                 success = True
            #     if not success:
            #         print "Unable to assign attribute '{0}'.".format(k)
            # except:
            #     pass
            #     # print "Error: Could not assign attribute '{0}'".format(k)


        """
        The approach below assigns values from input config dictionary to variables
        in sub-variable trees with a name that matches the value's key. Aside from checking
        the first-level members (i.e. fst_infile), it only navigates
        "down" a prescribed number of object members (2), so things like af_data (which are a 
        part of the sub-variable tree blade_aero) would have to be changed directly (or
        another approach would have to be used). Another drawback is that if any to channels
        have the same name (i.e. "Echo") they will both be assigned the value in the config
        dictionary.
        """
        for k, w in kwargs.iteritems():                
            try:
                success = False
                if hasattr(self, k):
                    setattr(self,k,w)
                    success = True
                else:
                    for key in self.turbsim_vt.__dict__:
                        subvartree = self.turbsim_vt.__dict__[key]
                        if hasattr(subvartree, k):
                            setattr(subvartree,k,w)
                            success = True
                if not success:
                    # These items are specific to FSTWorkflow and are assigned elsewhere
                    if k not in ['Turbsim_masterfile','Turbsim_masterdir','Turbsim_runfile',\
                        'Turbsim_rundir','Turbsim_exe']:
                        print "Could not find attribute '{0}'.".format(k)
            except:
                print "Something went wrong with assignment of attribute '{0}'.".format(k)
                pass

    def execute(self):
         if not os.path.exists(self.Turbsim_rundir): os.makedirs(self.Turbsim_rundir)
         #self.turbsim_vt.metboundconds.UserFile = os.sep.join([self.run_dir, self.turbulence_file_name])
         #self.turbsim_vt.metboundconds.ProfileFile = os.sep.join([self.run_dir, self.turbsim_vt.metboundconds.ProfileFile])

         # Write turbulence file
         #turb_specs(V_ref=float(self.wind_speed), L_u=float(self.L_u), L_v=float(self.L_v), L_w=float(self.L_w), sigma_u=float(self.sigma_u), sigma_v=float(self.sigma_v), sigma_w=float(self.sigma_w), filename=self.turbsim_vt.metboundconds.UserFile, template_file=self.turbulence_template_file)

         # Write profile file
         #write_wind(V_ref=float(self.wind_speed), alpha=float(self.shear_exponent), Beta=float(self.veer), Z_hub=float(self.turbsim_vt.tmspecs.HubHt), filename=self.turbsim_vt.metboundconds.ProfileFile, template_file=self.profile_template)

         #self.turbsim_vt.metboundconds.UserFile = os.sep.join(['..', self.run_dir, self.turbulence_file_name])
         #self.turbsim_vt.metboundconds.ProfileFile = os.sep.join(['..', self.turbsim_vt.metboundconds.ProfileFile])
         #self.turbsim_vt.metboundconds.ProfileFile = os.sep.join(['..', self.run_dir, self.turbsim_vt.metboundconds.ProfileFile])

         tsinp = open(os.sep.join([self.Turbsim_rundir, self.Turbsim_runfile]), 'w')
         tsinp.write("-----\n")
         tsinp.write("-----\n")
         tsinp.write("-----\n")

         # runtime options
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.echo))
         tsinp.write("{}\n".format(int(self.turbsim_vt.runtime_options.RandSeed1)))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.RandSeed2))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.WrBHHTP))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.WrFHHTP))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.WrADHH))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.WrADFF))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.WrBLFF))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.WrADTWR))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.WrFMTFF))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.WrACT))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.Clockwise))
         tsinp.write("{}\n".format(self.turbsim_vt.runtime_options.ScaleIEC))

         # Turbine/Model Specifications
         tsinp.write("\n")
         tsinp.write("----\n")
         tsinp.write("{}\n".format(self.turbsim_vt.tmspecs.NumGrid_Z))
         tsinp.write("{}\n".format(self.turbsim_vt.tmspecs.NumGrid_Y))
         tsinp.write("{}\n".format(self.turbsim_vt.tmspecs.TimeStep))
         tsinp.write("{}\n".format(self.turbsim_vt.tmspecs.AnalysisTime))
         tsinp.write("{}\n".format(self.turbsim_vt.tmspecs.UsableTime))
         tsinp.write("{}\n".format(self.turbsim_vt.tmspecs.HubHt))
         tsinp.write("{}\n".format(self.turbsim_vt.tmspecs.GridHeight))
         tsinp.write("{}\n".format(self.turbsim_vt.tmspecs.GridWidth))
         tsinp.write("{}\n".format(self.turbsim_vt.tmspecs.VFlowAng))
         tsinp.write("{}\n".format(self.turbsim_vt.tmspecs.HFlowAng))

         # Meteorological Boundary Conditions
         tsinp.write("\n")
         tsinp.write("----\n")
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.TurbModel))
         tsinp.write('"{}"\n'.format(self.turbsim_vt.metboundconds.UserFile))
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.IECstandard))
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.IECturbc))
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.IEC_WindType))
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.ETMc))
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.WindProfileType))
         tsinp.write('"{}"\n'.format(self.turbsim_vt.metboundconds.ProfileFile))
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.RefHt))
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.URef))
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.ZJetMax))
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.PLExp))
         tsinp.write("{}\n".format(self.turbsim_vt.metboundconds.Z0))

         # Non-IEC Meteorological Boundary Conditions
         tsinp.write("\n")
         tsinp.write("----\n")
         tsinp.write("{}\n".format(self.turbsim_vt.noniecboundconds.Latitude))
         tsinp.write("{}\n".format(self.turbsim_vt.noniecboundconds.RICH_NO))
         tsinp.write("{}\n".format(self.turbsim_vt.noniecboundconds.UStar))
         tsinp.write("{}\n".format(self.turbsim_vt.noniecboundconds.ZI))
         tsinp.write("{}\n".format(self.turbsim_vt.noniecboundconds.PC_UW))
         tsinp.write("{}\n".format(self.turbsim_vt.noniecboundconds.PC_UV))
         tsinp.write("{}\n".format(self.turbsim_vt.noniecboundconds.PC_VW))

         # Spatial Coherence Parameters
         tsinp.write("\n")
         tsinp.write("----\n")
         tsinp.write("{}\n".format(self.turbsim_vt.spatialcoherance.SCMod1))
         tsinp.write("{}\n".format(self.turbsim_vt.spatialcoherance.SCMod2))
         tsinp.write("{}\n".format(self.turbsim_vt.spatialcoherance.SCMod3))
         try:
             tsinp.write('"%f %f"\n'%(float(self.turbsim_vt.spatialcoherance.InCDec1[0]), float(self.turbsim_vt.spatialcoherance.InCDec1[1])))
         except:
             tsinp.write('"{}\n'.format(self.turbsim_vt.spatialcoherance.InCDec1[0]))
         try:    
             tsinp.write('"%f %f"\n'%(float(self.turbsim_vt.spatialcoherance.InCDec2[0]), float(self.turbsim_vt.spatialcoherance.InCDec2[1])))
         except:
         	 tsinp.write('"{}\n'.format(self.turbsim_vt.spatialcoherance.InCDec2[0]))
         try:
             tsinp.write('"%f %f"\n'%(float(self.turbsim_vt.spatialcoherance.InCDec3[0]), float(self.turbsim_vt.spatialcoherance.InCDec3[1])))
         except:
		 	 tsinp.write('"{}\n'.format(self.turbsim_vt.spatialcoherance.InCDec3[0]))
         tsinp.write("{}\n".format(self.turbsim_vt.spatialcoherance.CohExp))

         # Coherent Turbulence Scaling Parameters
         tsinp.write("\n")
         tsinp.write("----\n")
         tsinp.write("{}\n".format(self.turbsim_vt.coherentTurbulence.CTEventPath))
         tsinp.write("{}\n".format(self.turbsim_vt.coherentTurbulence.CTEventFile))
         tsinp.write("{}\n".format(self.turbsim_vt.coherentTurbulence.Randomize))
         tsinp.write("{}\n".format(self.turbsim_vt.coherentTurbulence.DistScl))
         tsinp.write("{}\n".format(self.turbsim_vt.coherentTurbulence.CTLy))
         tsinp.write("{}\n".format(self.turbsim_vt.coherentTurbulence.CTLz))
         tsinp.write("{}\n".format(self.turbsim_vt.coherentTurbulence.CTStartTime))



if __name__=='__main__':
    
    pass

