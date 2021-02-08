## @ingroup Methods-Aerodynamics
# SU2_inviscid.py
#
# Created:  Sep 2016, E. Botero
# Modified: Jan 2017, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# SUAVE imports
import SUAVE
from SUAVE.Core import Data, Units
from subprocess import call

# Local imports
from SUAVE.Analyses.Aerodynamics import Aerodynamics
from SUAVE.Input_Output.SU2.call_SU2_CFD import call_SU2_CFD
from SUAVE.Input_Output.SU2.write_SU2_cfg import write_SU2_cfg
from SUAVE.Input_Output.OpenVSP import vsp_write_Tecnam
from SUAVE.Input_Output.OpenVSP import vspaero_Tecnam

# Package imports
import sys
import numpy as np
import time
import pylab as plt
import sklearn
from sklearn import gaussian_process
from sklearn import neighbors
from sklearn import svm


# ----------------------------------------------------------------------
#  Class
# ----------------------------------------------------------------------
## @ingroup Analyses-Aerodynamics
class VSP_inviscid_No_Surrogates_Tecnam(Aerodynamics):
    """This builds a surrogate and computes lift and drag using SU2

    Assumptions:
    Inviscid, subsonic

    Source:
    None
    """   
    def __defaults__(self):
        """This sets the default values and methods for the analysis.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """ 
        self.tag = 'SU2_inviscid'

        self.geometry = Data()
        self.settings = Data()
        
        self.iters = 5 #OpenVSP number of iterations


    def evaluate(self,state,settings,geometry):
        """Evaluates lift and drag using available surrogates.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        state.conditions.
          mach_number      [-]
          angle_of_attack  [radians]

        Outputs:
        inviscid_lift      [-] CL
        inviscid_drag      [-] CD

        Properties Used:
        self.surrogates.
          lift_coefficient [-] CL
          drag_coefficient [-] CD
        """  
        # Unpack      
        conditions = state.conditions
        
        mach = conditions.freestream.mach_number
        AoA  = conditions.aerodynamics.angle_of_attack
        engines_number_tot = geometry.propulsors.propulsor.number_of_engines
        rpm = conditions.propulsion.rpm
        
        #print conditions.freestream
        
        rho=conditions.freestream.density
        vel_sound=conditions.freestream.speed_of_sound
    
        
        Cp=conditions.propulsion.propeller_power_coefficient
        
        Ct=conditions.propulsion.propeller_thrust_coefficient
        
        # Build OpenVSP geometry file
        tag  = geometry.tag
        vsp_write_Tecnam(geometry, tag)
        
        # Inviscid lift
        data_len = len(AoA)
        inviscid_lift = np.zeros([data_len,1])
        CD = np.zeros([data_len,1])
        time0 = time.time()
        
        initial_time = geometry.fuselages.fuselage.time
        
        flaps_angle = geometry.wings['main_wing'].control_surfaces.flap.deflection / Units.deg
        
        f = open("flapsangle.txt", "a")
        f.write(str(flaps_angle)+" ")
        f.close()
        
        for ii,_ in enumerate(AoA):
            inviscid_lift[ii], CD[ii] = vspaero_Tecnam(vel_sound[ii][0],tag + ".vsp3", rho[ii][0], AoA[ii][0], mach[ii][0], self.iters, rpm[ii][0], engines_number_tot, Cp[ii][0], Ct[ii][0], flaps_angle)
            print ('CL='+str(inviscid_lift[ii]))
            print ('CD='+str(CD[ii]))
            
        time1 = time.time()
            
        print ('The total elapsed time to run VSPAERO: '+ str(time1-time0) + '  seconds')
        print ('Current Simulation time is :' + str(time1-initial_time) + '  seconds')
            
        conditions.aerodynamics.lift_breakdown.inviscid_wings_lift       = Data()
        conditions.aerodynamics.lift_breakdown.inviscid_wings_lift.total = inviscid_lift
        conditions.aerodynamics.lift_coefficient                        = inviscid_lift
        state.conditions.aerodynamics.lift_coefficient                   = inviscid_lift
        #state.conditions.aerodynamics.lift_breakdown.compressible_wings  = inviscid_lift
        
        # Inviscid drag, zeros are a placeholder for possible future implementation
        #inviscid_drag                                              = np.zeros([data_len,1])
        conditions.aerodynamics.drag_breakdown.induced.inviscid   = CD
        state.conditions.aerodynamics.inviscid_drag_coefficient    = CD
        conditions.aerodynamics.drag_breakdown.induced.total       = CD
        
        # Now calculate the vehicle level oswald efficiency
        
        CL      = conditions.aerodynamics.lift_coefficient
        
        AR = geometry.wings.main_wing.aspect_ratio
        
        e_osw = (CL**2)/(np.pi*AR*CD)
        
        
        conditions.aerodynamics.drag_breakdown.induced.oswald_efficiency_factor = e_osw
        
        return inviscid_lift