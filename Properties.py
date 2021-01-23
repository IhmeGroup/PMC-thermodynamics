# -*- coding: utf-8 -*-
"""
PMC properties classes

Copyright 2020, Danyal Mohaddes, All rights reserved.

Refer to D. Mohaddes, C.T. Chang, M. Ihme, "Thermodynamic cycle analysis of superadiabatic matrix-stabilized combustion for gas turbine engines," Energy (207) 2020.
"""

# Interface-stabilized PMC with dimensional parameters based on Ref. 12
class InterfaceStabilized:
    def __init__(self):
        self.lambda_eff = 8.7     # Matrix thermal conductivity          [W/m.K]
        self.L = 3.0*0.0254       # Total burner length                  [m]
        self.massFlux = 2.0       # Mass flux through burner             [kg/m**2.s]
        self.voidFrac = 0.805     # Average void fraction                [-]
        self.dp_upstr = 0.492E-3  # Pore diameter (65PPI)                [m]
        self.dp_dwnstr = 3.192E-3 # Pore diameter (10PPI)                [m]
        self.emissivity = 0.9     # Thermal emissivity                   [-]
        self.DeltaT_PH = 20.      # Temperature difference (Eq. 24d)     [K]
        self.DeltaT_RC = 20.      # Temperature difference (Eq. 24e)     [K]

# Interface-stabilized PMC with non-dimensional parameters
class InterfaceStabilizedNondim:
    def __init__(self):
        self.iota = 0.05      # See Eq. 26,  [-]
        self.delta = 1.38     # See Eq. 26,  [-]
        self.mu = 0.0         # See Eq. 26,  [-]
        self.DeltaT_PH = 20.0 # See Eq. 24d, [K]
        self.DeltaT_RC = 20.0 # See Eq. 24e, [K]

