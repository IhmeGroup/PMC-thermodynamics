"""
Driver code for dimensional version of porous media combustor (PMC)

Copyright 2020, Danyal Mohaddes, All rights reserved.

Refer to D. Mohaddes, C.T. Chang, M. Ihme, "Thermodynamic cycle analysis of superadiabatic matrix-stabilized combustion for gas turbine engines," Energy (207) 2020.
"""

from PMC import PMC
from Properties import InterfaceStabilized

from matplotlib import pyplot as plt
import numpy as np

# For pretty plots
from matplotlib import rc
rc('text', usetex=True)
plt.rcParams.update({'font.size': 14})

# Gas turbine engine properties
f_st = 1./17.16 # Stoichiometric fuel/air ratio   [-]
LHV = 50.06e6   # Fuel lower heating value        [J/kg]
T2 = 300.       # Compressor inlet temperature    [K]
T4 = 1400.      # Turbine inlet temperature       [K]
gamma = 1.4     # Specific heat ratio             [-]
cp = 1200.      # Constant gaseous heat capacity  [J/kg.K]
eta_c = 0.85    # Compressor adiabatic efficiency [-]
eta_t = 0.9     # Turbine adiabatic efficiency    [-]
phi = 0.3       # Global equivalence ratio        [-]
beta = 0.102    # Combustor dilution ratio        [-]

# Function to get pressure ratio (Eq. 6b)
def piFunc(phi,eta_c):
    f = phi*f_st
    return (1. + eta_c*(1./T2*(T4-f/(1.+f)*LHV/cp) - 1.))**(gamma/(gamma-1.))

# Create PMC properties object
myProps = InterfaceStabilized()
# Can change properties here
myProps.emissivity = 0
# Create PMC object
myPMC = PMC(myProps)

# Get global fuel/air ratio
f = phi*f_st
# Get combustor inlet temperature (Eq. 6a), assuming T3a = T3
T3a = T4 - (f/(f+1.))*(LHV/cp)
# Get combustor local equivalence ratio (Eq. 3)
f_loc = f/(1.0 - beta)
# Run PMC model
myPMC.run(cp,T2,T3a,T4,f_loc,f_st,LHV)
# Get PMC temperatures
T_vec = [T3a,myPMC.T3x,myPMC.T3y,myPMC.T3b]
    
# Plotting
# Tick locations as per Fig. 8
PMC_xticks = [r'$\Phi_{3a}$',r'$\Phi_{3x}$',r'$\Phi_{3y}$',r'$\Phi_{3b}$']
# Plotting locations are non-physical, only for visualization purposes
xlocations = [0,1.4,1.65,3]

plt.plot(xlocations,T_vec, marker='o', linestyle='--', linewidth=1)

plt.xticks(xlocations,PMC_xticks)
plt.ylabel(r'$T_{ti} \ $[K]')

# Saving
plt.savefig('./PMC_dim_temperature_profile.eps',bbox_inches='tight')
