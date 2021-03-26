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
class VSP_Analysis_MLP_Tecnam_Updated(Aerodynamics):
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
        self.filenameCDi_HS = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CDi_HS.sav'
        self.filenameCDi_VS = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CDi_VS.sav'
        self.filenameCDi_Wing = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CDi_Wing.sav'
        self.filenameCL_HS = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CL_HS.sav'
        self.filenameCL_VS = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CL_VS.sav'
        self.filenameCL_Wing = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CL_Wing.sav'
        self.filenameCDi = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CDi.sav'
        self.filenameCD = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CD.sav'
        self.filenameCL = '/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/2ndMethology/Surrogate/NN_Surrogate_Model_MLPRegressor_CL.sav'
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
        AoA  = np.degrees(conditions.aerodynamics.angle_of_attack)
        #print(AoA)
        engines_number_tot = geometry.propulsors.propulsor.number_of_engines
        rpm = conditions.propulsion.rpm
        
        #print conditions.freestream
        
        
        # Create Result Data Structures 
        
        conditions.aerodynamics.drag_breakdown.induced                     = Data()
        conditions.aerodynamics.drag_breakdown.induced.inviscid_wings      = Data()
        conditions.aerodynamics.lift_breakdown                             = Data()
        conditions.aerodynamics.lift_breakdown.inviscid_wings              = Data()
        conditions.aerodynamics.lift_breakdown.compressible_wings          = Data()
        conditions.aerodynamics.drag_breakdown.compressible                = Data()
        
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
        #print(x_data.values[:,2])
        #x_data.values[:,3]=x_data.values[:,3]/mach_range[1]
        x_data.values[:,4]=x_data.values[:,4]/self.rpm_range[1]
        #x_data.values[:,5]=x_data.values[:,5]/cp_range[1]
        #x_data.values[:,6]=x_data.values[:,6]/ct_range[1]
        x_data.values[:,7]=x_data.values[:,7]/self.flap_range[1]
         
        #Load files
        
        CDi_HS_model = pickle.load(open(self.filenameCDi_HS, 'rb'))
        CDi_VS_model = pickle.load(open(self.filenameCDi_VS, 'rb'))
        CDi_Wing_model = pickle.load(open(self.filenameCDi_Wing, 'rb'))
        CL_HS_model = pickle.load(open(self.filenameCL_HS, 'rb'))
        CL_VS_model = pickle.load(open(self.filenameCL_VS, 'rb'))
        CL_Wing_model = pickle.load(open(self.filenameCL_Wing, 'rb'))
        CDi_model = pickle.load(open(self.filenameCDi, 'rb'))
        CD_model = pickle.load(open(self.filenameCD, 'rb'))
        CL_model = pickle.load(open(self.filenameCL, 'rb'))
        
        # Induced Drag, Drag and Inviscid lift
        CDi     = CDi_model.predict(x_data.values)
        CD      = CD_model.predict(x_data.values)
        CL      = CL_model.predict(x_data.values)
        
        CDi = CDi.reshape((len(CDi), 1))
        CD = CD.reshape((len(CD), 1))
        CL = CL.reshape((len(CL), 1))
        
        #print(CL)
        
        #print ('CL='+str(inviscid_lift))
        #print ('CD='+str(CD))
        
        # Pack
        conditions.aerodynamics.lift_coefficient                = CL
        conditions.aerodynamics.lift_breakdown.total            = CL
        conditions.aerodynamics.drag_breakdown.induced.inviscid = CDi
        
        
        for wing in geometry.wings.keys():
            #print(wing)
         
            if wing == 'horizontal_stabilizer':
                #print('Entering HS')
                CDi_Wing     = CDi_HS_model.predict(x_data.values)
                CL_Wing      = CL_HS_model.predict(x_data.values)
                
            elif wing == 'vertical_stabilizer':
                #print('Entering VS')
                CDi_Wing     = CDi_VS_model.predict(x_data.values)
                CL_Wing      = CL_VS_model.predict(x_data.values)
                
            elif wing == 'main_wing':
                #print('Entering Wing')
                CDi_Wing     = CDi_Wing_model.predict(x_data.values)
                CL_Wing      = CL_Wing_model.predict(x_data.values)
                
            CDi_Wing = CDi_Wing.reshape((len(CDi_Wing), 1))
            CL_Wing = CL_Wing.reshape((len(CL_Wing), 1))
                
            # Pack 
            conditions.aerodynamics.lift_breakdown.inviscid_wings[wing]         = CL_Wing
            conditions.aerodynamics.lift_breakdown.compressible_wings[wing]     = CL_Wing
            conditions.aerodynamics.drag_breakdown.induced.inviscid_wings[wing] = CDi_Wing
        
        
        return 
