# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from math import *
import cmath as cmath
#from Cit_par import *

# A i r c r a f t & Atmospheric p r op e r t i e s
g = 9.81 #m/ s^2
m = 6100 #kg
rho = 0.89 #kg /m^3
S = 30 #m^2
c = 2.0569 #m
b = 15.911 #m
K_yy = 1.3925**( 0.5 )
KY2 = 1.3925
K_xx = 0.019**( 1/ 2 )
KX2 = 0.019
K_zz = 0.042**( 1/ 2 )
KZ2 = 0.042
K_xz = 0.002
muc = m / ( rho * S * c)
mub = m / ( rho * S * b )
V0 = 88 #m/ s
theta = 4.5 #deg
theta = theta *np.pi /180 #rad
CL = 2*m*g / ( rho * V0 ** 2 * S )
#Symmetric f l i g h t p r op e r t i e s
CXu = -0.0279
CXa = -0.4797
CXadot = 0.0833
CXq = -0.2817
CX0 = m * g * np.sin(theta) / ( 0.5 *rho * V0 ** 2. * S )
CXdelta = -0.0373
CZ0 = - m * g * np.cos(theta) / ( 0.5 *rho * V0 ** 2. * S )
CZa = -5.7434*1.15 #was −5.7437
CZadot = -0.0035
CZu = -0.37616*1.5 #was −0.37616
CZq = -5.6629
CZdelta = -0.6961
Cmu = 0.0699
Cma = -0.7317 #was −0.5626 −−> rece i ved value from s t a t . measurements
Cmadot = 0.1780
Cmq = -8.7941
Cmdelta = -1.75631 #was −1.1642 −−> rece i ved value from s t a t . measurements
#Asymmetric p r op e r t i e s
CYb = -0.75
CYbeta_dot = 0
CYp = -0.0304
CYr = 0.8495
C_Y_delta_a = 0.04 #−0.04
C_Y_delta_r = -0.23 # 0.23
Clb = -0.1026
Clp = -0.7108
Clr = 0.2376
C_l_delta_a = 0.2309 #−0.2309
C_l_delta_r = -0.0344 #0.0344
Cnb = 0.1348
C_n_beta_dot = 0
Cnp = 0.1502 #was −0.0602
Cnr = -0.2061
C_n_delta_a = 0.012 #−0.012
C_n_delta_r = 0.0939 #−0.0939



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

    lambda_c_phug1 = (-B_phug_1 + 1j*cmath.sqrt(4 * A_phug_1 * C_phug_1 - B_phug_1 **2)) / (2 * A_phug_1) * V0/c
    lambda_c_phug2 = (-B_phug_1 - 1j*cmath.sqrt(4 * A_phug_1 * C_phug_1 - B_phug_1 **2)) / (2 * A_phug_1) * V0/c

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
