# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from math import *
import cmath as cmath
#from Cit_par import *

hp0    =  2     	      # pressure altitude in the stationary flight condition [m]
V0     =   100          # true airspeed in the stationary flight condition [m/sec]
alpha0 =    2         # angle of attack in the stationary flight condition [rad]
th0    =     2        # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      =   2000          # mass [kg]

# aerodynamic properties
e      =   0.75          # Oswald factor [ ]
CD0    =   2          # Zero lift drag coefficient [ ]
CLa    =   2          # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma    =   -0.56          # longitudinal stabilty [ ]
Cmde   =  2         # elevator effectiveness [ ]

# Aircraft geometry

S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	          # mean aerodynamic cord [m]
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	          # [ ]
ih     = -2 * np.pi / 180   # stabiliser angle of incidence [rad]

# Constant values concerning atmosphere and gravity

rho0   = 1.2250          # air density at sea level [kg/m^3]
Lambda = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)

# air density [kg/m^3]
#rho    = rho0 ** ( ((1+(Lambda * hp0 / Temp0))), (-((g / (Lambda*R)) + 1)))
rho =1.225
W      = m * g            # [N]       (aircraft weight)

# Constant values concerning aircraft inertia

muc    = m / (rho * S * c)
mub    = m / (rho * S * b)
KX2    = 0.019**2
KZ2    = 0.042**2
KXZ    = 0.002
KY2    = 1.25 * 1.114

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * np.pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Lift and drag coefficient

CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e) # Drag coefficient [ ]

# Stabiblity derivatives

CX0    = W * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu    = -0.02792
CXa    = +0.47966		# Positive! (has been erroneously negative since 1993)
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZ0    = -W * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu    = -0.37616
CZa    = -5.74340
CZadot = -0.00350
CZq    = -5.66290
CZde   = -0.69612

Cmu    = +0.06990
Cmadot = +0.17800
Cmq    = -8.79415

CYb    = -0.7500
CYbdot =  0
CYp    = -0.0304
CYr    = +0.8495
CYda   = -0.0400
CYdr   = +0.2300

Clb    = -0.10260
Clp    = -0.71085
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.1348
Cnbdot =   0
Cnp    =  -0.0602
Cnr    =  -0.2061
Cnda   =  -0.0120
Cndr   =  -0.0939



#Transformations moeten nog!!!

######### Short period motion #############
def short_period():
    # A_spm = muc ** 2 * KY2      #KY2 is already the squared one
    # B_spm = -2 * muc * (KY2 * CZa + Cmadot + Cmq)
    # C_spm = CZa * Cmq - 2 * muc * Cma

    #Refined approximation
    A_spm = 4 * muc ** 2 * KY2      #KY2 is already the squared one
    B_spm = -2 * muc * (KY2 * CZa + Cmadot + Cmq)
    C_spm = CZa * Cmq - 2 * muc * Cma

    lambda_c_spm1 = (-B_spm + 1j*cmath.sqrt(4 * A_spm * C_spm - B_spm **2)) / (2 * A_spm) * V0/c
    lambda_c_spm2 = (-B_spm - 1j*cmath.sqrt(4 * A_spm * C_spm - B_spm **2)) / (2 * A_spm) * V0/c
    period_spm = 2*np.pi/abs(lambda_c_spm1.imag) #* c/V0
    t_half_spm = np.log(0.5)/lambda_c_spm1.imag #* c/V0
    zeta_spm  = -lambda_c_spm1.real / np.sqrt((lambda_c_spm1.imag)**2+(lambda_c_spm1.real)**2)
    omega_n_spm = np.sqrt((lambda_c_spm1.imag)**2+(lambda_c_spm1.real)**2) * V0/c * np.sqrt(1-zeta_spm**2)

    # #Coarse approximation
    # A_spm_2 = -2 * muc * KY2_squared
    # B_spm_2 = Cmadot + Cmq
    # C_spm_2 = Cma
    # spm_char1 = np.poly1d([A_spm_1, B_spm_1, C_spm_1])
    # lambda_c_spm1 = np.roots(spm_char1)
    # spm_char2 = np.poly1d([A_spm_2, B_spm_2, C_spm_2])
    # lambda_c_spm2 = np.roots(spm_char2)

    return lambda_c_spm1, lambda_c_spm2, period_spm, t_half_spm, zeta_spm, omega_n_spm

######### Phugoid motion #############
def phugoid():
    # A_phug = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
    # B_phug = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq*(CZu * CXa - CXu * CZa)
    # C_phug = CZ0 * (Cmu * CZa - CZu * Cma)

    #Refined approximation
    A_phug_1 = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
    B_phug_1 = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq*(CZu * CXa - CXu * CZa)
    C_phug_1 = CZ0 * (Cmu * CZa - CZu * Cma)

    lambda_c_phug1 = (-B_phug_1 + cmath.sqrt(4 * A_phug_1 * C_phug_1 - B_phug_1 **2)) / (2 * A_phug_1) * V0/c
    lambda_c_phug2 = (-B_phug_1 - cmath.sqrt(4 * A_phug_1 * C_phug_1 - B_phug_1 **2)) / (2 * A_phug_1) * V0/c

    period_phug = 2*np.pi/abs(lambda_c_phug1.imag) * c/V0
    t_half_phug = np.log(0.5)/lambda_c_phug1.imag * c/V0
    zeta_phug  = -lambda_c_phug1.real / np.sqrt((lambda_c_phug1.imag)**2+(lambda_c_phug1.real)**2)
    omega_n_phug = np.sqrt((lambda_c_phug1.imag)**2+(lambda_c_phug1.real)**2) * V0/c * np.sqrt(1-zeta_phug**2)

    # #Coarse approximation
    # A_phug_2 = -4 * muc ** 2
    # B_phug_2 = 2 * muc * CXu
    # C_phug_2 = -CZu * CZ0
    # phug_char2 = np.poly1d([A_phug_2, B_phug_2, C_phug_2])
    # lambda_c_phug2 = np.roots(phug_char2)
    #
    # lambda_c_phug1 = (-B_phug + 1j*np.sqrt(4 * A_phug_1 * C_phug_1 - B_phug_1 **2)) / (2 * A_phug_1)
    # lambda_c_phug2 = (-B_phug - 1j*np.sqrt(4 * A_phug_1 * C_phug_1 - B_phug_1 **2)) / (2 * A_phug_1)
    # t_half = np.log(0.5) * c / (lambda_c_phug1 * V0)
    # time_cst = -c / (lambda_c_phug1 * V0)

    return lambda_c_phug1, lambda_c_phug2, period_phug, t_half_phug, zeta_phug, omega_n_phug

######### A-periodic roll motion #############
def aperiodic_roll():
    lambda_b_apr = Clp / (4 * mub * KX2)
    t_half = np.log(0.5) * b / (lambda_b_apr * V0)
    time_cst = -c / (lambda_b_apr * V0)
    return lambda_b_apr, t_half, time_cst


######### Spiral motion #############
def spiral():
    lambda_b_spir = (2 * CL * (Clb * Cnr - Cnb * Clr)) / ((Clp * (CYb * Cnr + 4 * mub * Cnb )) - Cnp * (CYb * Clr + 4 * mub * Clb))
    t_half = np.log(0.5) * c / (lambda_b_spir * V0)
    time_cst = -c / (lambda_b_spir * V0)
    return lambda_b_spir, t_half, time_cst


######### Dutch roll motion  case 1 #############
def dutchroll():
    A_dutch = 8 * mub ** 2 * KZ2
    B_dutch = -2 * mub * (Cnr + 2 * KZ2 * CYb)
    C_dutch = 4 * mub * Cnb + CYb * Cnr

    dutch_char = np.poly1d([A_dutch, B_dutch, C_dutch])
    lambda_b_dutch = np.roots(dutch_char) #misschien V0/b toevoegen eraan
    period_dutch = 2 * np.pi * b / (lambda_b_dutch[0].imag * V0)
    t_half_dutch = np.log(0.5) * b / (lambda_b_dutch[0].real * V0)
    zeta_dutch = -lambda_b_dutch[0].real / np.sqrt((lambda_b_dutch[0].imag)**2+(lambda_b_dutch[0].real)**2)
    omega_n_dutch = np.sqrt((lambda_b_dutch[0].imag)**2+(lambda_b_dutch[0].real)**2) * V0/b * np.sqrt(1-zeta_dutch**2)

    return lambda_b_dutch[0], lambda_b_dutch[1], period_dutch, t_half_dutch, zeta_dutch, omega_n_dutch


######### Dutch roll motion  case 2 #############
def dutchroll_2():
    A_dutch2 = -2 * mub * KZ2
    B_dutch2 = 1/2 * Cnr
    C_dutch2 = -Cnb

    lambda_b_dutch2_1 = (-B_dutch2 + np.sqrt(4 * A_dutch2 * C_dutch2 - B_dutch2 **2)) / (2 * A_dutch2)
    lambda_b_dutch2_2 = (-B_dutch2 - np.sqrt(4 * A_dutch2 * C_dutch2 - B_dutch2 **2)) / (2 * A_dutch2)
    return lambda_b_dutch2_1, lambda_b_dutch2_2

lambda_c_spm1, lambda_c_spm2, period_spm, t_half_spm, zeta_spm, omega_n_spm = short_period()
lambda_c_phug1, lambda_c_phug2, period_phug, t_half_phug, zeta_phug, omega_n_phug = phugoid()

lambda_b_apr, t_half_apr, time_cst_apr = aperiodic_roll()
lambda_b_spir, t_half_spir, time_cst_spir = spiral()
lambda_b_dutch1, lambda_b_dutch2, period_dutch, t_half_dutch, zeta_dutch, omega_n_dutch = dutchroll()

print("short_period")
print("-----------------------")
print("eigenvalues:", lambda_c_spm1, lambda_c_spm2)
print("t_half:", t_half_spm)
print("zeta:", zeta_spm)
print("omega_n:", omega_n_spm)
print()
print("phugoid")
print("-----------------------")
print("eigenvalues:", lambda_c_phug1, lambda_c_phug2)
print("t_half:", t_half_phug)
print("zeta:", zeta_phug)
print("omega_n:", omega_n_phug)
print()
print("aperiodic roll")
print("-----------------------")
print("eigenvalue:", lambda_b_apr)
print("t_half:", t_half_apr)
print("time constant", time_cst_apr)
print()
print("spiral")
print("-----------------------")
print("eigenvalue:", lambda_b_spir)
print("t_half:", t_half_spir)
print("time constant", time_cst_spir)
print()
print("dutchroll")
print("-----------------------")
print("eigenvalues:", lambda_c_phug1, lambda_c_phug2)
print("t_half:", t_half_phug)
print("zeta:", zeta_phug)
print("omega_n:", omega_n_phug)
