## @ingroup Analyses-Aerodynamics
# SU2_Euler.py
#
# Created:  Sep 2016, E. Botero
# Modified: Jan 2017, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Data
from .Markup import Markup
from SUAVE.Analyses import Process
import numpy as np

from SUAVE.Input_Output.OpenVSP.write_vsp_mesh import write_vsp_mesh
from SUAVE.Input_Output.GMSH.write_geo_file import write_geo_file
from SUAVE.Input_Output.GMSH.mesh_geo_file import mesh_geo_file

# Default aero Results
#from .Results import Results

# The aero methods
from SUAVE.Methods.Aerodynamics import Fidelity_Zero as Methods
from SUAVE.Methods.Aerodynamics.Common import Fidelity_Zero as Common
from .Process_Geometry import Process_Geometry
from SUAVE.Methods.Aerodynamics import VSP_Analysis_MLP_Tecnam_Updated

# ----------------------------------------------------------------------
#  Analysis
# ----------------------------------------------------------------------
## @ingroup Analyses-Aerodynamics
class Open_VSP_Analysis_Tecnam_Updated(Markup):
    """This uses SU2 to compute lift.

    Assumptions:
    Subsonic

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
        self.tag    = 'fidelity_zero_markup'
    
        # correction factors
        settings = self.settings
        settings.fuselage_lift_correction           = 1.14
        settings.trim_drag_correction_factor        = 1.02
        settings.wing_parasite_drag_form_factor     = 1.1
        settings.fuselage_parasite_drag_form_factor = 2.3
        settings.oswald_efficiency_factor           = None
        settings.span_efficiency                    = None
        settings.viscous_lift_dependent_drag_factor = 0.38
        settings.drag_coefficient_increment         = 0.0000
        settings.spoiler_drag_increment             = 0.00 
        settings.maximum_lift_coefficient           = np.inf
        settings.number_panels_spanwise             = None 
        settings.number_panels_chordwise            = None 
        settings.use_surrogate                      = False 
        settings.include_slipstream_effect          = False 
        settings.plot_vortex_distribution           = False
        
        # build the evaluation process
        compute = self.process.compute
        
        compute.lift = Process()

        # Run OpenVSP to determine lift
        compute.lift.inviscid_wings                = VSP_Analysis_MLP_Tecnam_Updated.VSP_Analysis_MLP_Tecnam_Updated()
        compute.lift.vortex                        = SUAVE.Methods.skip
        compute.lift.fuselage                      = Common.Lift.fuselage_correction
        compute.lift.total                         = Common.Lift.aircraft_total
        
        compute.drag = Process()
        compute.drag.parasite                      = Process()
        compute.drag.parasite.wings                = Process_Geometry('wings')
        compute.drag.parasite.wings.wing           = Common.Drag.parasite_drag_wing 
        compute.drag.parasite.fuselages            = Process_Geometry('fuselages')
        compute.drag.parasite.fuselages.fuselage   = Common.Drag.parasite_drag_fuselage
        compute.drag.parasite.propulsors           = Process_Geometry('propulsors')
        compute.drag.parasite.propulsors.propulsor = Common.Drag.parasite_drag_propulsor
        compute.drag.parasite.pylons               = Common.Drag.parasite_drag_pylon
        compute.drag.parasite.total                = Common.Drag.parasite_total
        compute.drag.induced                       = Common.Drag.induced_drag_aircraft
        compute.drag.compressibility               = Process()
        compute.drag.compressibility.wings         = Process_Geometry('wings')
        compute.drag.compressibility.wings.wing    = Common.Drag.compressibility_drag_wing
        compute.drag.compressibility.total         = Common.Drag.compressibility_drag_wing_total
        compute.drag.miscellaneous                 = Common.Drag.miscellaneous_drag_aircraft_ESDU
        compute.drag.untrimmed                     = Common.Drag.untrimmed
        compute.drag.trim                          = Common.Drag.trim
        compute.drag.spoiler                       = Common.Drag.spoiler_drag
        compute.drag.total                         = Common.Drag.total_aircraft
        
        
    def initialize(self):
        pass
        """Initializes the surrogate needed for lift calculation.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        self.geometry
        """                  
        super(Open_VSP_Analysis_Tecnam_Updated, self).initialize()
        
        use_surrogate             = self.settings.use_surrogate
        include_slipstream_effect = self.settings.include_slipstream_effect 
        vortex_distribution_flag  = self.settings.plot_vortex_distribution 
        n_sw                      = self.settings.number_panels_spanwise    
        n_cw                      = self.settings.number_panels_chordwise  
                                  
        self.process.compute.lift.inviscid_wings.geometry = self.geometry 
        self.process.compute.lift.inviscid_wings.initialize(use_surrogate , vortex_distribution_flag , n_sw ,  n_cw ,include_slipstream_effect )          
                                                            
    finalize = initialize                                          