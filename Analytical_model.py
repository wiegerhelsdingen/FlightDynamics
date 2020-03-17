# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from math import *
from Cit_par import *



######### Short period motion #############
def short_period():

    #Refined approximation
    A_spm_1 = 4 * muc ** 2 * KY2_squared      #KY2 is already the squared one
    B_spm_1 = -2 * muc * (KY2_squared * CZa + Cmadot + Cmq)
    C_spm_1 = CZa * Cmq - 2 * muc * Cma

    #Coarse approximation
    A_spm_2 = -2 * muc * KY2_squared
    B_spm_2 = Cmadot + Cmq
    C_spm_2 = Cma

    spm_char1 = np.poly1d([A_spm_1, B_spm_1, C_spm_1])
    lambda_c_spm1 = np.roots(spm_char1)

    spm_char2 = np.poly1d([A_spm_2, B_spm_2, C_spm_2])
    lambda_c_spm2 = np.roots(spm_char2)

    return lambda_c_spm1, lambda_c_spm2
print(short_period())

######### Phugoid motion #############
def phugoid():

    #Refined approximation
    A_phug_1 = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
    B_phug_1 = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq*(CZu * CXa - CXu * CZa)
    C_phug_1 = CZ0 * (Cmu * CZa - CZu * Cma)

    A_phug_2 = -4 * muc ** 2
    B_phug_2 = 2 * muc * CXu
    C_phug_2 = -CZu * CZ0

    phug_char1 = np.poly1d([A_phug_1, B_phug_1, C_phug_1])
    lambda_c_phug1 = np.roots(phug_char1)

    phug_char2 = np.poly1d([A_phug_2, B_phug_2, C_phug_2])
    lambda_c_phug2 = np.roots(phug_char2)

    return lambda_c_phug1, lambda_c_phug2

######### A-periodic roll motion #############
def aperiodic_roll():
    lambda_b_apr = Clp / (4 * mub * KX2)
    return lambda_b_apr


######### Spiral motion #############
def spiral():
    lambda_b_spir = (2 * CL * (Clb * Cnr - Cnb * Clr)) / ((Clp * (CYb * Cnr + 4 * mub * Cnb )) - Cnp * (CYb * Clr + 4 * mub * Clb))
    return lambda_b_spir


######### Dutch roll motion  case 1 #############
def dutchroll():
    A_dutch = 8 * mub ** 2 * KZ2
    B_dutch = -2 * mub * (Cnr + 2 * KZ2 * CYb)
    C_dutch = 4 * mub * Cnb + CYb * Cnr


    lambda_b_dutch1 = (-B_dutch + np.sqrt(4 * A_dutch * C_dutch - B_dutch **2)) / (2 * A_dutch)
    lambda_b_dutch2 = (-B_dutch - np.sqrt(4 * A_dutch * C_dutch - B_dutch **2)) / (2 * A_dutch)

    dutch_char = np.poly1d([A_dutch, B_dutch, C_dutch])
    lambda_b_dutch11 = np.roots(dutch_char)
    return lambda_b_dutch1, lambda_b_dutch2
######### Dutch roll motion  case 2 #############
def dutchroll_2():
    A_dutch2 = -2 * mub * KZ2
    C_dutch2 = -Cnb
    #B_ dutch = 1/2 * Cnr


    lambda_b_dutch2_1 = (-B_dutch2 + np.sqrt(4 * A_dutch2 * C_dutch2 - B_dutch2 **2)) / (2 * A_dutch2)
    lambda_b_dutch2_2 = (-B_dutch2 - np.sqrt(4 * A_dutch2 * C_dutch2 - B_dutch2 **2)) / (2 * A_dutch2)
    return lambda_b_dutch2_1, lambda_b_dutch2_2
