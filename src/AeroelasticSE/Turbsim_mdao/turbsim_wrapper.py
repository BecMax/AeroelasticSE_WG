import os
import subprocess

from AeroelasticSE.Turbsim_mdao.turbsim_writer import TurbsimBuilder
from AeroelasticSE.Turbsim_mdao.turbsim_reader import turbsimReader
from turbsim_vartrees import turbsiminputs

class TurbsimExternalCode(object):

    pass

class Turbsim_wrapper(TurbsimExternalCode):
    def __init__(self):
        super(Turbsim_wrapper, self).__init__()
        
        self.Turbsim_exe = None
        self.Turbsim_runfile = None
        self.Turbsim_rundir = None
         
    def execute(self):
        
        print "Executing Turbsim v2"
        self.input_file = os.path.join(self.Turbsim_rundir, self.Turbsim_runfile)
        
        print "Calling ", self.Turbsim_exe
        print "Input file = ", self.input_file

        # Get only tail of input_file (we are changing running directory)
        (head,tail) = os.path.split(self.input_file)

        # Construct absolute path of executable
        Turbsimexe_abs = self.Turbsim_exe
        print "Turbsimexe_new: ", Turbsimexe_abs

        exec_str = []
        exec_str.append(Turbsimexe_abs)
        exec_str.append(tail)

        olddir = os.getcwd()
        os.chdir(self.Turbsim_rundir)
        print "current directory: ", os.getcwd()
        print "exec_str: ", exec_str
        subprocess.call(exec_str)#, stdin=None, stdout=None, stderr=None, shell=False)
        os.chdir(olddir)

         
if __name__=='__main__':
    
    pass

