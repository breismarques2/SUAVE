#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 22:59:23 2017
@author: root
"""

import SUAVE
from SUAVE.Core import Units, Data
from subprocess import call

import os
import math

try:
    import sys
    sys.path.insert(0, '/Users/Bruno/OpenVSP-python3/OpenVSP/build/python_api')
    import vsp_g as vsp
    import vsp as vsp1

except ImportError:
    # This allows SUAVE to build without OpenVSP
    pass
import numpy as np


## @ingroup Input_Output-OpenVSP
def vspaero_Tecnam(vel_sound,tag,rho,AoA,MachNumber,NumberIterations, rpm, engines_number_tot, Cp, Ct, flap_deflection):
    
    if 1==1:

        try:
            
            vsp.ClearVSPModel()
            
        except NameError:
            print ('VSP import failed - vspaero')
            return -1
        
        print("Testing Setup - Start\n")
            
        vsp.VSPCheckSetup()
        vsp.VSPRenew()
        
        print("Testing Setup - End")


        #open the file created in vsp_write
    
        vsp.ReadVSPFile(tag)
        
        #//==== Execute Mass Properties Analysis ====//
        #// Set defaults
        vsp.SetAnalysisInputDefaults( "MassProp" );
        #PrintAnalysisInputs( "MassProp" );
        ridmp = vsp.ExecAnalysis( "MassProp" );
        #PrintResults( ridmp );
        #csvname='Mass_Prop_'+tag[:-5]
        #vsp.WriteResultsCSVFile( ridmp, csvname+'.csv');
           
        data = []
        with open(tag[:-5]+"_MassProps.txt") as f:
            for line in f:
                data.append([word for word in line.split(" ") if word])
        f.close()
        
        line_CenterGravity = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'Center']
        
        Xcg=float(data[line_CenterGravity[0]][0])
        Ycg=float(data[line_CenterGravity[0]][1])
        Zcg=float(data[line_CenterGravity[0]][2])

        vsp.DeleteGeomVec( vsp.GetStringResults( ridmp, "Mesh_GeomID" ) );
        
        #==== Analysis: VSPAero Compute Geometry to Create Vortex Lattice DegenGeom File ====//
        
        compgeom_name = "VSPAEROComputeGeometry";

        # Set defaults
        vsp.SetAnalysisInputDefaults(compgeom_name);

        # list inputs, type, and current values
        #vsp.PrintAnalysisInputs(compgeom_name);

        # Execute
        compgeom_resid=vsp.ExecAnalysis(compgeom_name);

        # Get & Display Results
        #vsp.PrintResults( compgeom_resid );
        
    
        #==== Analysis: VSPAero Compute Geometry ====//
    
        analysis_name="VSPAEROSweep"
    
        #Set defaults
    
        vsp.SetAnalysisInputDefaults(analysis_name)
    
        #Change some input values
        #    Analysis method
    
        analysis_method = [vsp.VORTEX_LATTICE]
    
    
        vsp.SetIntAnalysisInput( analysis_name, "AnalysisMethod", analysis_method, 0 )
        
        #Reference geometry set
        
        geom_set=[0]
        vsp.SetIntAnalysisInput( analysis_name, "GeomSet", geom_set );
                               
        #Reference areas, lengths
        
        wing_id=vsp.FindGeomsWithName("main_wing")
        
        #print wing_id[0]
        
        #sref=[float(0)]*1
        #bref=[float(0)]*1
        #cref=[float(0)]*1
        
        
        #sref[0]=(vsp.GetParmVal(wing_id[0],"TotalArea", "WingGeom"))
        #bref[0]=(vsp.GetParmVal(wing_id[0],"TotalSpan", "WingGeom"))
        #cref[0]=(vsp.GetParmVal(wing_id[0],"TotalChord", "WingGeom"))
        
        #print sref
        #print bref
        #print cref
        
        ref_flag=[1]
        #ref_wing=[1]
        
        
        #vsp.SetDoubleAnalysisInput( analysis_name, 'Sref', sref )
        #vsp.SetDoubleAnalysisInput( analysis_name, 'bref', bref )
        #vsp.SetDoubleAnalysisInput( analysis_name, 'cref', cref )
        vsp.SetIntAnalysisInput( analysis_name, "RefFlag", ref_flag )
        vsp.SetStringAnalysisInput( analysis_name, 'WingID', wing_id )
        
        
        #Center of Gravity
        
        vsp.SetDoubleAnalysisInput( analysis_name, "Xcg", [Xcg] )
        vsp.SetDoubleAnalysisInput( analysis_name, "Ycg", [Ycg] )
        vsp.SetDoubleAnalysisInput( analysis_name, "Zcg", [Zcg] )

        
                                 
        #Freestream parameters
        #Alpha
        
        alpha_start=[float(math.degrees(AoA))]
        alpha_end=[float(math.degrees(AoA))]
        alpha_npts=[1]
        
        print ("Alpha")
        print (alpha_start)
        
        vsp.SetDoubleAnalysisInput( analysis_name, "AlphaStart", alpha_start )
        vsp.SetDoubleAnalysisInput( analysis_name, "AlphaEnd", alpha_end )
        vsp.SetIntAnalysisInput( analysis_name, "AlphaNpts", alpha_npts )
        
        #Beta
        beta_start=[float(0)]
        beta_end=[float(0)]
        beta_npts=[1]
        
        
        vsp.SetDoubleAnalysisInput( analysis_name, "BetaStart", beta_start )
        vsp.SetDoubleAnalysisInput( analysis_name, "BetaEnd", beta_end )
        vsp.SetIntAnalysisInput( analysis_name, "BetaNpts", beta_npts );
                               
        #Mach
        
        mach_start=[float(MachNumber)]
        mach_end=[float(MachNumber)]
        mach_npts=[1]
        
        print ("Mach")
        print (mach_start)
        
        vsp.SetDoubleAnalysisInput( analysis_name, "MachStart", mach_start )
        vsp.SetDoubleAnalysisInput( analysis_name, "MachEnd", mach_end )
        vsp.SetIntAnalysisInput( analysis_name, "MachNpts", mach_npts )
                                
        vsp.Update()
        
        #Actuator Disk Flag
        
        vsp.SetIntAnalysisInput( analysis_name, 'ActuatorDiskFlag',[1])
        
        #vsp.Update()
        
        
        #Vinf
        
        vel=MachNumber * vel_sound
        
        vel_rotor=float(vel)
        print ("Vel Rotor")
        print (vel_rotor)
        
        vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
        vel_id = vsp.FindParm( vspaero_settings_container_id, 'Vinf', 'VSPAERO')
        vsp.SetParmVal( vel_id, vel_rotor)
        
        #rho
        
        rho_rotor=float(rho)
        print ("Rho")
        print (rho_rotor)
        
        vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
        rho_id = vsp.FindParm( vspaero_settings_container_id, 'Rho', 'VSPAERO')
        vsp.SetParmVal( rho_id, rho_rotor)
        
        
        #rpm

        number_eng = int(engines_number_tot)
        print ("Number of Engines")
        print (number_eng)
        rpm_float=float(rpm)
        print ("RPM")
        print (rpm_float)
        
        g=int(0)
        
        
        while g<number_eng:
            
            if rpm_float==0.0:
                if g == 0:
                    aux_rpm_1=0.001
                else:
                    aux_rpm_1=-0.001
                #vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                #rpm_id = vsp.FindParm( vspaero_settings_container_id, 'RotorRPM', 'Rotor_'+str(g))
                #vsp.SetParmVal( rpm_id, aux_rpm_1)  
                disk_id = vsp.FindActuatorDisk(g)
                vsp.SetParmVal( vsp.FindParm( disk_id, 'RotorRPM', 'Rotor' ), aux_rpm_1)
            elif math.isnan(rpm_float)==True:
                if g == 0:
                    aux_rpm_1=0.001
                else:
                    aux_rpm_1=-0.001
                #vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                #rpm_id = vsp.FindParm( vspaero_settings_container_id, 'RotorRPM', 'Rotor_'+str(g))
                #vsp.SetParmVal( rpm_id, aux_rpm_1)
                disk_id = vsp.FindActuatorDisk(g)
                vsp.SetParmVal( vsp.FindParm( disk_id, 'RotorRPM', 'Rotor' ), aux_rpm_1)

            else:
                if g == 0:
                    aux_rpm_1=rpm_float
                else:
                    aux_rpm_1=-rpm_float
                #vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                #rpm_id = vsp.FindParm( vspaero_settings_container_id, 'RotorRPM', 'Rotor_'+str(g))
                #vsp.SetParmVal( rpm_id, aux_rpm_1)
                disk_id = vsp.FindActuatorDisk(g)
                vsp.SetParmVal( vsp.FindParm( disk_id, 'RotorRPM', 'Rotor' ), aux_rpm_1 );

            vsp.Update()
            g=g+1
            
            

        #Power Coefficient
        
        Cp_float=float(Cp)
        print ("Cp")
        print (Cp_float)
    
        
        
        g=int(0)
        
        while g<number_eng:

                
            if Cp_float<-1000.00000:
                aux_cp_1=-1000.00000
                vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                cp_id = vsp.FindParm( vspaero_settings_container_id, 'RotorCP', 'Rotor_'+str(g))
                vsp.SetParmVal( cp_id, aux_cp_1)
                    
            elif Cp_float>1000.00000:
                aux_cp_1=1000.00000
                vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                cp_id = vsp.FindParm( vspaero_settings_container_id, 'RotorCP', 'Rotor_'+str(g))
                vsp.SetParmVal( cp_id, aux_cp_1)
            elif math.isnan(Cp_float)==True:
                aux_cp_1=0.0
                vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                cp_id = vsp.FindParm( vspaero_settings_container_id, 'RotorCP', 'Rotor_'+str(g))
                vsp.SetParmVal( cp_id, aux_cp_1)

                    
            else:
                aux_cp_1=Cp_float
                vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                cp_id = vsp.FindParm( vspaero_settings_container_id, 'RotorCP', 'Rotor_'+str(g))
                vsp.SetParmVal( cp_id, aux_cp_1)
                    

            g=g+1
            
        
            
        
        #Thrust Coefficient
        
        Ct_float=float(Ct)
        print ("Ct")
        print (Ct_float)
        
        
        g=int(0)
        
        while g<number_eng:
                
            if 0.001 < Ct_float< 1000.00000:
                aux_ct_1=Ct_float
                vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                ct_id = vsp.FindParm( vspaero_settings_container_id, 'RotorCT', 'Rotor_'+str(g))
                vsp.SetParmVal( ct_id, aux_ct_1)
                    
            elif Ct_float > 1000.00000:
                aux_ct_1=1000.00000
                vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                ct_id = vsp.FindParm( vspaero_settings_container_id, 'RotorCT', 'Rotor_'+str(g))
                vsp.SetParmVal( ct_id, aux_ct_1)
                    
            elif math.isnan(Ct_float)==True:
                aux_ct_1=0.001
                vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                ct_id = vsp.FindParm( vspaero_settings_container_id, 'RotorCT', 'Rotor_'+str(g))
                vsp.SetParmVal( ct_id, aux_ct_1)

                    
            else:
                aux_ct_1=0.001
                vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
                ct_id = vsp.FindParm( vspaero_settings_container_id, 'RotorCT', 'Rotor_'+str(g))
                vsp.SetParmVal( ct_id, aux_ct_1)
                              

            g=g+1
            
            
        
        #vsp.SetDoubleAnalysisInput( analysis_name, "MachABFDS", mach_start )
        
        #vspaero_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
        #rpm_cruise_id = vsp.FindParm( vspaero_settings_container_id, 'RotorRPM', 'Rotor_0')
        #vsp.SetParmValUpdate( rpm_cruise_id, 1000.0)
                                
        #vsp.Update()
                  
        ##Case Setup
        
        wakeNumIter=[int(NumberIterations)]
        wakeSkipUntilIter=[int(NumberIterations+1)]
        batch_mode_flag=[1]
        
        vsp.SetIntAnalysisInput( analysis_name, "WakeNumIter", wakeNumIter )
        #vsp.SetIntAnalysisInput( analysis_name, "WakeSkipUntilIter", wakeSkipUntilIter )
        vsp.SetIntAnalysisInput( analysis_name, "BatchModeFlag", batch_mode_flag )
                                
        vsp.Update()
        
        # Control Surfaces
        
        group_index = vsp.CreateVSPAEROControlSurfaceGroup()
        
        vsp.SetVSPAEROControlGroupName("main_wing_flaps_group", group_index)
        
        vsp.AddAllToVSPAEROControlSurfaceGroup( group_index );
        
        vsp.Update()
        
        control_group_settings_container_id = vsp.FindContainer( "VSPAEROSettings", 0 )
        
        subsur_ids = vsp.GetAllSubSurfIDs()	
        
        deflection_gain_id = vsp.FindParm(control_group_settings_container_id,'Surf_' + subsur_ids[0] + '_1_Gain', "ControlSurfaceGroup" )
        
        vsp.SetParmVal(deflection_gain_id, -1.000)
        
        deflection_angle_id = vsp.FindParm( control_group_settings_container_id, "DeflectionAngle", "ControlSurfaceGroup_0" )
        
        vsp.SetParmVal(deflection_angle_id, flap_deflection)
        
    
        #list inputs, type, and current values
    
        #vsp.PrintAnalysisInputs(analysis_name)
        
    
        #Execute
    
        rid = vsp.ExecAnalysis(analysis_name)
            
    
        #Get & Display Results
    
        #vsp.PrintResults(rid)
                        
        #Write in CSV
        
        csvname='Simulation_'+tag[:-5]
        
        vsp.WriteResultsCSVFile(rid,csvname+'.csv')
        
        #Close File
        
        vsp.ClearVSPModel();
    
        # Check for errors

        #errorMgr = vsp.ErrorMgrSingleton_getInstance()
        #num_err = errorMgr.GetNumTotalErrors()
        #for i in range(0, num_err):
        #    err = errorMgr.PopLastError()
        #    print("error = ", err.m_ErrorString)
            
        print ('FINISHED VSPAERO SIMULATION')
        
        #data = []
        #with open(csvname+'.csv') as f:
        #    for line in f:
        #        data.append([word for word in line.split(",") if word])
        #f.close()
        
        #try:
        #    lift=data[18][NumberIterations]
        #    drag=data[14][NumberIterations]
            
        #    CL = float(lift[:-1])
        #    CD = float(drag[:-1])
        
        #except IndexError:
            
        #os.system('/anaconda/bin/vspaero -fs '+str(mach_start[0])+'  END '+str(alpha_start[0])+'  END '+str(beta_start[0])+'  END -omp 4 -nowake 4 '+tag[:-5]+'_DegenGeom')
        
        data = []
        try:
            
            # CL ,CDi and CD0 Total
            
            with open(tag[:-5]+'_DegenGeom.fem') as f:
                for line in f:
                    data.append([word for word in line.split(" ") if word])
            f.close()
            
            line_TotalForces = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'Total']
            CL=float(data[line_TotalForces[1]][2])
            CDi=float(data[line_TotalForces[2]][2])
            
            data = []
            
            with open(tag[:-5]+'_DegenGeom.polar') as f:
                for line in f:
                    data.append([word for word in line.split(" ") if word])
            f.close()
            
            CD0 = float(data[1][5])
            
            # Skin Friction Drag Break Out
            
            data = []
            
            
            with open(tag[:-5]+'_DegenGeom.history') as f:
                for line in f:
                    data.append([word for word in line.split(" ") if word])
            f.close()
            
            Line_main_wing_cdo = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'main_wing']
            Line_horizontal_stabilizer_cdo = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'horizontal_stabilizer']
            Line_vertical_stabilizer_cdo = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'vertical_stabilizer']
            Line_fuselage_cdo = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'fuselage']
            
            main_wing_cdo = float(data[Line_main_wing_cdo[0]][1])+float(data[Line_main_wing_cdo[1]][1])
            horizontal_stabilizer_cdo = float(data[Line_horizontal_stabilizer_cdo[0]][1])+float(data[Line_horizontal_stabilizer_cdo[1]][1])
            vertical_stabilizer_cdo = float(data[Line_vertical_stabilizer_cdo[0]][1])
            fuselage_cdo = float(data[Line_fuselage_cdo[0]][1])+float(data[Line_fuselage_cdo[1]][1])+float(data[Line_fuselage_cdo[2]][1])+float(data[Line_fuselage_cdo[3]][1])
            
            # Lift and Induced Drag Break Out
            
            data = []
            
            
            with open(tag[:-5]+'_DegenGeom.lod') as f:
                for line in f:
                    data.append([word for word in line.split(" ") if word])
            f.close()
            
            Line_main_wing = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'main_wing']
            Line_horizontal_stabilizer = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'horizontal_stabilizer']
            Line_vertical_stabilizer = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'vertical_stabilizer']
            Line_fuselage = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'fuselage']
            
            main_wing_cl = float(data[Line_main_wing[0]][5])+float(data[Line_main_wing[1]][5])
            horizontal_stabilizer_cl = float(data[Line_horizontal_stabilizer[0]][5])+float(data[Line_horizontal_stabilizer[1]][5])
            vertical_stabilizer_cl = float(data[Line_vertical_stabilizer[0]][5])
            fuselage_cl = float(data[Line_fuselage[0]][5])+float(data[Line_fuselage[1]][5])+float(data[Line_fuselage[2]][5])+float(data[Line_fuselage[3]][5])
            
            main_wing_cdi = float(data[Line_main_wing[0]][6])+float(data[Line_main_wing[1]][6])
            horizontal_stabilizer_cdi = float(data[Line_horizontal_stabilizer[0]][6])+float(data[Line_horizontal_stabilizer[1]][6])
            vertical_stabilizer_cdi = float(data[Line_vertical_stabilizer[0]][6])
            fuselage_cdi = float(data[Line_fuselage[0]][6])+float(data[Line_fuselage[1]][6])+float(data[Line_fuselage[2]][6])+float(data[Line_fuselage[3]][6])
            
            # Lift and Drag Main Wing Distribution
            
            data = []
                       
            with open(tag[:-5]+'_DegenGeom.lod') as f:
                for line in f:
                    data.append([word for word in line.split(" ") if word])
            f.close()
            
            Line_Wing = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'Wing']
            
            wing_lift_distribution = []
            wing_drag_distribution = []
            
            current_line = Line_Wing[0]+1
            
            while data[current_line][0]=='1':
                
                wing_lift_distribution.append(data[current_line][5])
                wing_drag_distribution.append(data[current_line][6])
                current_line = current_line+1 
            
            
        except IOError:
            
            print ("Other part of code")
            os.system('/anaconda/envs/py3/bin/vspaero -fs '+str(mach_start[0])+'  END '+str(alpha_start[0])+'  END '+str(beta_start[0])+'  END -omp 4 -nowake 4 '+tag[:-5]+'_DegenGeom')
            with open(tag[:-5]+'_DegenGeom.fem') as f:
                for line in f:
                    data.append([word for word in line.split(" ") if word])
            f.close()
            
            line_TotalForces = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'Total']
        
            CL=float(data[line_TotalForces[1]][2])
            CDi=float(data[line_TotalForces[2]][2])
            
            data = []
            
            with open(tag[:-5]+'_DegenGeom.polar') as f:
                for line in f:
                    data.append([word for word in line.split(" ") if word])
            f.close()
            
            CD0 = float(data[1][5])
            
            # Skin Friction Drag Break Out
            
            data = []
            
            
            with open(tag[:-5]+'_DegenGeom.history') as f:
                for line in f:
                    data.append([word for word in line.split(" ") if word])
            f.close()
            
            Line_main_wing_cdo = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'main_wing']
            Line_horizontal_stabilizer_cdo = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'horizontal_stabilizer']
            Line_vertical_stabilizer_cdo = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'vertical_stabilizer']
            Line_fuselage_cdo = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'fuselage']
            
            main_wing_cdo = float(data[Line_main_wing_cdo[0]][1])+float(data[Line_main_wing_cdo[1]][1])
            horizontal_stabilizer_cdo = float(data[Line_horizontal_stabilizer_cdo[0]][1])+float(data[Line_horizontal_stabilizer_cdo[1]][1])
            vertical_stabilizer_cdo = float(data[Line_vertical_stabilizer_cdo[0]][1])
            fuselage_cdo = float(data[Line_fuselage_cdo[0]][1])+float(data[Line_fuselage_cdo[1]][1])+float(data[Line_fuselage_cdo[2]][1])+float(data[Line_fuselage_cdo[3]][1])
            
            # Lift and Induced Drag Break Out
            
            data = []
            
            
            with open(tag[:-5]+'_DegenGeom.lod') as f:
                for line in f:
                    data.append([word for word in line.split(" ") if word])
            f.close()
            
            Line_main_wing = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'main_wing']
            Line_horizontal_stabilizer = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'horizontal_stabilizer']
            Line_vertical_stabilizer = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'vertical_stabilizer']
            Line_fuselage = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'fuselage']
            
            main_wing_cl = float(data[Line_main_wing[0]][5])+float(data[Line_main_wing[1]][5])
            horizontal_stabilizer_cl = float(data[Line_horizontal_stabilizer[0]][5])+float(data[Line_horizontal_stabilizer[1]][5])
            vertical_stabilizer_cl = float(data[Line_vertical_stabilizer[0]][5])
            fuselage_cl = float(data[Line_fuselage[0]][5])+float(data[Line_fuselage[1]][5])+float(data[Line_fuselage[2]][5])+float(data[Line_fuselage[3]][5])
            
            main_wing_cdi = float(data[Line_main_wing[0]][6])+float(data[Line_main_wing[1]][6])
            horizontal_stabilizer_cdi = float(data[Line_horizontal_stabilizer[0]][6])+float(data[Line_horizontal_stabilizer[1]][6])
            vertical_stabilizer_cdi = float(data[Line_vertical_stabilizer[0]][6])
            fuselage_cdi = float(data[Line_fuselage[0]][6])+float(data[Line_fuselage[1]][6])+float(data[Line_fuselage[2]][6])+float(data[Line_fuselage[3]][6])
            
            # Lift and Drag Main Wing Distribution
            
            data = []
                       
            with open(tag[:-5]+'_DegenGeom.lod') as f:
                for line in f:
                    data.append([word for word in line.split(" ") if word])
            f.close()
            
            Line_Wing = [ix for ix, row in enumerate(data) for iy, i in enumerate(row) if i == 'Wing']
            
            wing_lift_distribution = []
            wing_drag_distribution = []
            
            current_line = Line_Wing[0]+1
            
            while data[current_line][0]=='1':
                
                wing_lift_distribution.append(data[current_line][5])
                wing_drag_distribution.append(data[current_line][6])
                current_line = current_line+1 
            
        
        #try:
            
        #drag=data[14][NumberIterations]
            
        #except IndexError:
            
        #    print 'INDEX ERROR'
        #    drag='-10.0000'
        
        #myfile="/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/Code/Bruno_Aircraft/Optimization_Lo_Fid/"+csvname+'.csv'

        ## If file exists, delete it ##
        #if os.path.isfile(myfile):
        #    os.remove(myfile)
        #else:    ## Show an error ##
        #    print("Error: %s file not found" % myfile)
        
        cwd = os.getcwd()
            
        myfile=cwd+'/'+tag[:-5]+'_DegenGeom.lod'

        ## If file exists, delete it ##
        if os.path.isfile(myfile):
            os.remove(myfile)
        else:    ## Show an error ##
            print("Error: %s file not found" % myfile)
            
        myfile=cwd+'/'+tag[:-5]+'_DegenGeom.adb.cases'

        ## If file exists, delete it ##
        if os.path.isfile(myfile):
            os.remove(myfile)
        else:    ## Show an error ##
            print("Error: %s file not found" % myfile)
            
        myfile=cwd+'/'+tag[:-5]+'_DegenGeom.adb'

        ## If file exists, delete it ##
        if os.path.isfile(myfile):
            os.remove(myfile)
        else:    ## Show an error ##
            print("Error: %s file not found" % myfile)
            
        myfile=cwd+'/'+tag[:-5]+'_DegenGeom.history'

        ## If file exists, delete it ##
        if os.path.isfile(myfile):
            os.remove(myfile)
        else:    ## Show an error ##
            print("Error: %s file not found" % myfile)
            
        myfile=cwd+'/'+tag[:-5]+'_DegenGeom.polar'

        ## If file exists, delete it ##
        if os.path.isfile(myfile):
            os.remove(myfile)
        else:    ## Show an error ##
            print("Error: %s file not found" % myfile)
            
        myfile=cwd+'/'+tag[:-5]+'_DegenGeom.fem'

        ## If file exists, delete it ##
        if os.path.isfile(myfile):
            os.remove(myfile)
        else:    ## Show an error ##
            print("Error: %s file not found" % myfile)
            
        myfile=cwd+'/'+tag[:-5]+'_DegenGeom.vspaero'

        ## If file exists, delete it ##
        if os.path.isfile(myfile):
            os.remove(myfile)
        else:    ## Show an error ##
            print("Error: %s file not found" % myfile)
            
        myfile=cwd+'/'+tag[:-5]+'_DegenGeom.csv'

        ## If file exists, delete it ##
        if os.path.isfile(myfile):
            os.remove(myfile)
        else:    ## Show an error ##
            print("Error: %s file not found" % myfile)
            
            
        myfile=cwd+'/'+tag[:-5]+"_MassProps.txt"

        # If file exists, delete it ##
        if os.path.isfile(myfile):
            os.remove(myfile)
        else:    ## Show an error ##
            print("Error: %s file not found" % myfile)
        
        
        
        #print float(lift[:-1])
        #print float(drag[:-1])
        
            
                
        #wb = openpyxl.Workbook()
        #sheet = wb.active
        #for row_index in range(len(data)):
        #    for col_index, letter in zip(range(len(data[row_index])), string.ascii_uppercase):
        #        sheet[letter+str(row_index+1)]= data[row_index][col_index]

        #wb.save(csvname+'.xlsx')
        
        #xlsx_filename=csvname+'.xlsx'
        
        #with open(xlsx_filename, "rb") as f:
        #    in_mem_file = io.BytesIO(f.read())

        #book = openpyxl.load_workbook(in_mem_file, read_only=True)
        
        #book = openpyxl.load_workbook(csvname+'.xlsx')
        
        #sheet = book.active
        
        #CL = sheet[chr(NumberIterations+65)+'19'].value
        #CD = sheet[chr(NumberIterations+65)+'15'].value
                   
        #CL = float(lift[:-1])
        #CD = float(drag[:-1])
                   
        #book.save(csvname+'.xlsx')
                   
        #call(["cd","/Users/Bruno/Documents/Delft/Courses/2016-2017/Thesis/Code/Bruno_Aircraft","&&","echo","iforgot","|","sudo","-S","rm","Simulation_takeoff.xlsx"])

        #call(["echo","iforgot","|","sudo","-S","rm","Simulatio.xlsx"])    
        
        
    
    return CL, CDi, CD0, main_wing_cl, horizontal_stabilizer_cl, vertical_stabilizer_cl, fuselage_cl, main_wing_cdi, horizontal_stabilizer_cdi, vertical_stabilizer_cdi, fuselage_cdi, main_wing_cdo, horizontal_stabilizer_cdo, vertical_stabilizer_cdo, fuselage_cdo, wing_lift_distribution, wing_drag_distribution

