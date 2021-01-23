# -*- coding: utf-8 -*-
"""
Porous Media Combustor (PMC) Class

Copyright 2020, Danyal Mohaddes, All rights reserved.

Refer to D. Mohaddes, C.T. Chang, M. Ihme, "Thermodynamic cycle analysis of superadiabatic matrix-stabilized combustion for gas turbine engines," Energy (207) 2020.
"""

import numpy as np
from sympy import Symbol, nsolve

# Class for solving steady-state PMC problem with dimensional parameters
class PMC:
    def __init__(self,PMCData):
        self.lambda_eff = PMCData.lambda_eff
        self.voidFrac = PMCData.voidFrac
        self.massFlux = PMCData.massFlux
        self.L = PMCData.L
        self.dp_upstr = PMCData.dp_upstr
        self.dp_dwnstr = PMCData.dp_dwnstr
        self.emissivity = PMCData.emissivity
        self.DeltaT_PH = PMCData.DeltaT_PH
        self.DeltaT_RC = PMCData.DeltaT_RC
       
    # Solve the PMC problem 
    def run(self,cp,T0,T3a,T4,f,f_st,LHV):
        
        th_3a = T3a/T0
        th_4 = T4/T0
        
        K_R_1 = 3*(1.-self.voidFrac)/(self.dp_upstr*1000.) #from Hsu & Powell,
        K_R_2 = 3*(1.-self.voidFrac)/(self.dp_dwnstr*1000.) # pore diameter in mm and K_R
        K_R = np.average([K_R_1,K_R_2])                     # in m^-1
        
        iota = self.lambda_eff/(self.massFlux*cp*self.L)
        delta = (16.*ct.stefan_boltzmann*T0**3)/(3.*K_R*self.lambda_eff)
        mu = (3.*self.emissivity*K_R*(1.-self.voidFrac)*self.L)/16.
        
        print('iota =',iota,'delta =',delta,'mu =',mu)
        
        Delta1 = self.DeltaT_PH/T0
        Delta2 = self.DeltaT_RC/T0
            
        if f < f_st:
            LHV = LHV #lean
        else:
            LHV = LHV/(f/f_st) #rich
        
        th_3x = Symbol('THETA_3x')
        th_3y = Symbol('THETA_3y')
        th_3b = Symbol('THETA_3b')
        th_S1 = Symbol('THETA_S1')
        th_S2 = Symbol('THETA_S2')
        
        f1 = iota*((1.+delta*((th_S1+th_S2)/2.)**3)*(th_S2-th_S1) - 
                   delta*mu*(th_S1**4-th_3a**4)) - (th_3x-th_3a)
        f2 = f/(f+1.)*(LHV/(cp*T0)) - (th_3y - th_3x) #local eq ratio (DON'T use beta here!!)
        f3 = iota*((1.+delta*((th_S1+th_S2)/2.)**3)*(th_S2-th_S1) + 
                   delta*mu*(th_S2**4-th_4**4)) - (th_3y-th_3b)
        f4 = Delta1 - (th_S1 - th_3x)
        f5 = Delta2 - (th_3b - th_S2)
        
        SOL = nsolve((f1,f2,f3,f4,f5),(th_3x,th_3y,th_3b,th_S1,th_S2),\
                     ((T3a+1000)/T3a,(T3a+1500)/T3a,\
                      (T3a+1000)/T3a,(T3a+1000)/T3a,\
                      (T3a+1000)/T3a),tol=1e-6)
        SOL *= T0 # return to dimensional units
        self.T3x,self.T3y,self.T3b = SOL[0],SOL[1],SOL[2] #gas temperatures
        self.TS1,self.TS2 = SOL[3],SOL[4] #solid temperatures
        return

# Class for solving steady-state PMC problem with non-dimensional parameters
class PMCNondim(PMC):

    def __init__(self,PMCData):
        self.iota = PMCData.iota
        self.delta = PMCData.delta
        self.mu = PMCData.mu
        self.DeltaT_PH = PMCData.DeltaT_PH
        self.DeltaT_RC = PMCData.DeltaT_RC

    # Solve the PMC problem 
    def run(self,cp,T0,T3a,T4,f,f_st,LHV):
        th_3a = T3a/T0
        th_4 = T4/T0
        
        Delta1 = self.DeltaT_PH/T0
        Delta2 = self.DeltaT_RC/T0
        
        if f < f_st:
            LHV = LHV #lean
        else:
            LHV = LHV/(f/f_st) #rich
        
        th_3x = Symbol('THETA_3x')
        th_3y = Symbol('THETA_3y')
        th_3b = Symbol('THETA_3b')
        th_S1 = Symbol('THETA_S1')
        th_S2 = Symbol('THETA_S2')
        
        iota = self.iota
        delta = self.delta
        mu = self.mu

        f1 = iota*((1.+delta*((th_S1+th_S2)/2.)**3)*(th_S2-th_S1) - 
                   delta*mu*(th_S1**4-th_3a**4)) - (th_3x-th_3a)
        f2 = f/(f+1.)*(LHV/(cp*T0)) - (th_3y - th_3x) #local eq ratio (DON'T use beta here!!)
        f3 = iota*((1.+delta*((th_S1+th_S2)/2.)**3)*(th_S2-th_S1) + 
                   delta*mu*(th_S2**4-th_4**4)) - (th_3y-th_3b)
        f4 = Delta1 - (th_S1 - th_3x)
        f5 = Delta2 - (th_3b - th_S2)
        
        SOL = nsolve((f1,f2,f3,f4,f5),(th_3x,th_3y,th_3b,th_S1,th_S2),\
                     ((T3a+1000)/T3a,(T3a+1500)/T3a,\
                      (T3a+1000)/T3a,(T3a+1000)/T3a,\
                      (T3a+1000)/T3a),tol=1e-6)

        SOL *= T0 # return to dimensional units
        self.T3x,self.T3y,self.T3b = SOL[0],SOL[1],SOL[2] #gas temperatures
        self.TS1,self.TS2 = SOL[3],SOL[4] #solid temperatures
        return
