## @ingroup Analyses-Aerodynamics
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

# Local imports
from SUAVE.Analyses.Aerodynamics import Aerodynamics
from SUAVE.Input_Output.SU2.call_SU2_CFD import call_SU2_CFD
from SUAVE.Input_Output.SU2.write_SU2_cfg import write_SU2_cfg
from SUAVE.Input_Output.OpenVSP import vsp_write_x57
from SUAVE.Input_Output.OpenVSP import vsp_write
from SUAVE.Input_Output.OpenVSP import vspaero

# Package imports
import numpy as np
import time
import pylab as plt
import sklearn
import pickle
import pandas as pd
from sklearn import gaussian_process
from sklearn import neighbors
from sklearn import svm
import statsmodels.formula.api as smf
import pandas
import math
import os

# ----------------------------------------------------------------------
#  Class
# ----------------------------------------------------------------------
## @ingroup Analyses-Aerodynamics
class VSP_Analysis_MLP_Tecnam(Aerodynamics):
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
        #self.settings.half_mesh_flag     = True
        #self.settings.parallel           = False
        #self.settings.processors         = 1
        #self.settings.maximum_iterations = 1500

        self.training_file             = None
        
        # Surrogate model
        self.vel_sound_range = [325.0, 345.0] #(m/s)
        self.rho_range       = [0.85, 1.3]
        self.alpha_range     = [0.0, 20.0] #º
        self.mach_range      = [0.0, 0.3] 
        self.rpm_range       = [0.0, 6000.0]
        self.cp_range        = [0.0, 0.8]
        self.ct_range        = [0.0, 0.5]
        self.flap_range      = [0.0, 30.0] #º
        self.filenameCL = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CL.sav'
        self.filenameCD = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CD.sav'
        self.surrogates = Data()


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
        
        flaps_angle = geometry.wings['main_wing'].control_surfaces.flap.deflection / Units.deg
        
        # Build x_data vector
        
        data_len = len(AoA)
        
        df = pd.DataFrame(columns=['vel_sound','rho', 'alpha', 'mach', 'rpm', 'cp', 'ct', 'flap'])
        
        for ii,_ in enumerate(AoA):
            if math.isnan(AoA[ii])==True:
                AoA[ii]=0.0
            if math.isnan(mach[ii])==True:
                mach[ii]=0.0
            if math.isnan(rpm[ii])==True:
                rpm[ii]=0.001
            if math.isnan(Cp[ii])==True:
                Cp[ii]=0.0
            if math.isnan(Ct[ii])==True:
                Ct[ii]=0.001
                
            df = df.append({'vel_sound': vel_sound[ii][0],'rho': rho[ii][0], 'alpha': AoA[ii][0], 'mach': mach[ii][0], 'rpm': rpm[ii][0], 'cp': Cp[ii][0], 'ct': Ct[ii][0], 'flap': flaps_angle}, ignore_index=True)
            
        x_data=df
        
        #Normalize data between 0 and 1

        x_data.values[:,0]=x_data.values[:,0]/self.vel_sound_range[1]
        x_data.values[:,1]=x_data.values[:,1]/self.rho_range[1]
        x_data.values[:,2]=x_data.values[:,2]/self.alpha_range[1]
        #x_data.values[:,3]=x_data.values[:,3]/mach_range[1]
        x_data.values[:,4]=x_data.values[:,4]/self.rpm_range[1]
        #x_data.values[:,5]=x_data.values[:,5]/cp_range[1]
        #x_data.values[:,6]=x_data.values[:,6]/ct_range[1]
        x_data.values[:,7]=x_data.values[:,7]/self.flap_range[1]
                
        
        lift_model = pickle.load(open(self.filenameCL, 'rb'))
        drag_model = pickle.load(open(self.filenameCD, 'rb'))
        
        # Inviscid lift and Drag
        inviscid_lift = lift_model.predict(x_data.values)
        CD            = drag_model.predict(x_data.values)
        
        #for ii,_ in enumerate(inviscid_lift):
        #    print (inviscid_lift[ii])
        #    inviscid_lift[ii]=[inviscid_lift[ii]]
        #    CD[ii]=[CD[ii]]
        
        inviscid_lift = inviscid_lift.reshape((len(inviscid_lift), 1))
        CD = CD.reshape((len(CD), 1))
        
        #print ('CL='+str(inviscid_lift))
        #print ('CD='+str(CD))
        
            
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
