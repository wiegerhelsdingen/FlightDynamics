# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from math import *
from Cit_par_test import *



######### Short period motion #############
def short_period():
    A_spm = muc ** 2 * KY2      #KY2 is already the squared one
    B_spm = -2 * muc * (KY2 * CZa + Cmadot + Cmq)
    C_spm = CZa * Cmq - 2 * muc * Cma

    lambda_c_spm1 = (-B_spm + np.sqrt(4 * A_spm * C_spm - B_spm **2)) / (2 * A_spm)
    lambda_c_spm2 = (-B_spm - np.sqrt(4 * A_spm * C_spm - B_spm **2)) / (2 * A_spm)
    t_half = np.log(0.5) * c / (lambda_c_spm1 * V0)
    time_cst = -c / (lambda_c_spm1 * V0)
    return lambda_c_spm1, lambda_c_spm2, t_half, time_cst

######### Phugoid motion #############
def phugoid():
    A_phug = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
    B_phug = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq*(CZu * CXa - CXu * CZa)
    C_phug = CZ0 * (Cmu * CZa - CZu * Cma)

    lambda_c_phug1 = (-B_phug + np.sqrt(4 * A_phug * C_phug - B_phug **2)) / (2 * A_phug)
    lambda_c_phug2 = (-B_phug - np.sqrt(4 * A_phug * C_phug - B_phug **2)) / (2 * A_phug)
    t_half = np.log(0.5) * c / (lambda_c_phug1 * V0)
    time_cst = -c / (lambda_c_phug1 * V0)
    return lambda_c_phug1, lambda_c_phug2, t_half, time_cst

######### A-periodic roll motion #############
def aperiodic_roll():
    lambda_b_apr = Clp / (4 * mub * KX2)
    t_half = np.log(0.5) * b / (lambda_b_apr * V0)
    time_cst = -c / (lambda_B_apr * V0)
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
    t_half = np.log(0.5) * b / (lambda_b_dutch[0].real * V0)
    time_cst = -c / (lambda_b_dutch[0].real * V0)
    p = 2 * np.pi * b / (lambda_b_dutch[0].imag * V0)
    return lambda_b_dutch, t_half, time_cst, p


######### Dutch roll motion  case 2 #############
def dutchroll_2():
    A_dutch2 = -2 * mub * KZ2
    B_dutch2 = 1/2 * Cnr
    C_dutch2 = -Cnb
       
    lambda_b_dutch2_1 = (-B_dutch2 + np.sqrt(4 * A_dutch2 * C_dutch2 - B_dutch2 **2)) / (2 * A_dutch2)
    lambda_b_dutch2_2 = (-B_dutch2 - np.sqrt(4 * A_dutch2 * C_dutch2 - B_dutch2 **2)) / (2 * A_dutch2)
    return lambda_b_dutch2_1, lambda_b_dutch2_2