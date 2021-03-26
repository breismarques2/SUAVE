## @defgroup Methods-Aerodynamics Aerodynamics
# Aerodynamic methods contain the functions for the aerodynamic analyses.
# @ingroup Methods

from . import AVL
from . import AERODAS
from . import Fidelity_Zero
from . import Supersonic_Zero
from . import Common
from . import Lifting_Line
from .VSP_inviscid_No_Surrogates import VSP_inviscid_No_Surrogates
from .VSP_inviscid_No_Surrogates_Tecnam import VSP_inviscid_No_Surrogates_Tecnam
from .VSP_Analysis_OLS import VSP_Analysis_OLS
from .VSP_Analysis_MLP_Tecnam import VSP_Analysis_MLP_Tecnam
from .VSP_Analysis_MLP_Tecnam_Updated import VSP_Analysis_MLP_Tecnam_Updated