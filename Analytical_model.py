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
    A_spm = muc ** 2 * KY2      #KY2 is already the squared one
    B_spm = -2 * muc * (KY2 * CZa + Cmadot + Cmq)
    C_spm = CZa * Cmq - 2 * muc * Cma

    lambda_c_spm1 = (-B_spm + np.sqrt(4 * A_spm * C_spm - B_spm **2)) / (2 * A_spm)
    lambda_c_spm2 = (-B_spm - np.sqrt(4 * A_spm * C_spm - B_spm **2)) / (2 * A_spm)
    return lambda_c_spm1, lambda_c_spm2

######### Phugoid motion #############
def phugoid():
    A_phug = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
    B_phug = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq*(CZu * CXa - CXu * CZa)
    C_phug = CZ0 * (Cmu * CZa - CZu * Cma)

    lambda_c_phug1 = (-B_phug + np.sqrt(4 * A_phug * C_phug - B_phug **2)) / (2 * A_phug)
    lambda_c_phug2 = (-B_phug - np.sqrt(4 * A_phug * C_phug - B_phug **2)) / (2 * A_phug)
    return lambda_c_phug1, lambda_c_phug2

######### A-periodic roll motion #############
def aperiodic_roll():
    lambda_b_apr = Clp / (4 * mub * KX2)
    return lambda_b_apr


######### Spiral motion #############
def lambda_b_spir():
    lambda_b_spir = (2 * CL * (Clb * Cnr - Cnb * Clr)) / ((Clp * (CYb * Cnr + 4 * mub * Cnb )) - Cnp * (CYb * Clr + 4 * mub * Clb))
    return lambda_b_spir


######### Dutch roll motion  case 1 #############
def dutchroll():
    A_dutch = 8 * mub ** 2 * KZ2
    B_dutch = -2 * mub * (Cnr + 2 * KZ2 * CYb)
    C_dutch = 4 * mub * Cnb + CYb * Cnr


    lambda_b_dutch1 = (-B_dutch + np.sqrt(4 * A_dutch * C_dutch - B_dutch **2)) / (2 * A_dutch)
    lambda_b_dutch2 = (-B_dutch - np.sqrt(4 * A_dutch * C_dutch - B_dutch **2)) / (2 * A_dutch)
    return lambda_b_dutch1, lambda_b_dutch2, lambda_d
print(dutchroll())

######### Dutch roll motion  case 2 #############
# A_dutch2 =
# B_ dutch2 =
# C_dutch2 =
