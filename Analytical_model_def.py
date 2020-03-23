# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from math import *
import cmath as cmath
from parameters import *
# from Cit_par import *

######### Short period motion #############
def short_period(V0, m, rho, muc, mub, CL, CD, CX0, CZ0):
    # A_spm = muc ** 2 * KY2      #KY2 is already the squared one
    # B_spm = -2 * muc * (KY2 * CZa + Cmadot + Cmq)
    # C_spm = CZa * Cmq - 2 * muc * Cma

    #Refined approximation
    A_spm = 4 * muc ** 2 * KY2    #KY2 is already the squared one
    B_spm = -2 * muc * (KY2 * CZa + Cmadot + Cmq)
    C_spm = CZa * Cmq - 2 * muc * Cma

    lambda_c_spm1 = (-B_spm + 1j*cmath.sqrt(4 * A_spm * C_spm - B_spm **2)) / (2 * A_spm)
    lambda_c_spm2 = (-B_spm - 1j*cmath.sqrt(4 * A_spm * C_spm - B_spm **2)) / (2 * A_spm)
    period_spm = 2*np.pi/abs(lambda_c_spm1.imag) * c/V0
    t_half_spm = np.log(0.5)/lambda_c_spm1.real * c/V0
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
def phugoid(V0, m, rho, muc, mub, CL, CD, CX0, CZ0):
    print(m)
    # A_phug = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
    # B_phug = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq*(CZu * CXa - CXu * CZa)
    # C_phug = CZ0 * (Cmu * CZa - CZu * Cma)

    #Refined approximation
    A_phug_1 = 2 * muc * (CZa * Cmq - 2 * muc * Cma)
    B_phug_1 = 2 * muc * (CXu * Cma - Cmu * CXa) + Cmq*(CZu * CXa - CXu * CZa)
    C_phug_1 = CZ0 * (Cmu * CZa - CZu * Cma)

    lambda_c_phug1 = (-B_phug_1 + 1j*cmath.sqrt(4 * A_phug_1 * C_phug_1 - B_phug_1 **2)) / (2 * A_phug_1)
    lambda_c_phug2 = (-B_phug_1 - 1j*cmath.sqrt(4 * A_phug_1 * C_phug_1 - B_phug_1 **2)) / (2 * A_phug_1)

    period_phug = 2*np.pi/abs(lambda_c_phug1.imag) * c/V0
    t_half_phug = np.log(0.5)/lambda_c_phug1.real * c/V0
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
def aperiodic_roll(V0, m, rho, muc, mub, CL, CD, CX0, CZ0):
    lambda_b_apr = Clp / (4 * mub * KX2)
    t_half = np.log(0.5) * b / (lambda_b_apr * V0)
    time_cst = -b / (lambda_b_apr * V0)
    return lambda_b_apr, t_half, time_cst


######### Spiral motion #############
def spiral(V0, m, rho, muc, mub, CL, CD, CX0, CZ0):
    lambda_b_spir = (2 * CL * (Clb * Cnr - Cnb * Clr)) / ((Clp * (CYb * Cnr + 4 * mub * Cnb )) - Cnp * (CYb * Clr + 4 * mub * Clb))
    t_half = -np.log(0.5) * b / (lambda_b_spir * V0)
    time_cst = -b / (lambda_b_spir * V0)
    return lambda_b_spir, t_half, time_cst


######### Dutch roll motion  case 1 #############
def dutchroll(V0, m, rho, muc, mub, CL, CD, CX0, CZ0):
    A_dutch = 8 * mub ** 2 * KZ2
    B_dutch = -2 * mub * (Cnr + 2 * KZ2 * CYb)
    C_dutch = 4 * mub * Cnb + CYb * Cnr

    dutch_char = np.poly1d([A_dutch, B_dutch, C_dutch])
    lambda_b_dutch = np.roots(dutch_char) #misschien V0/b toevoegen eraan
    period_dutch = 2 * np.pi / abs(lambda_b_dutch[0].imag) * b/V0
    t_half_dutch = np.log(0.5) / lambda_b_dutch[0].real * b/V0
    zeta_dutch = -lambda_b_dutch[0].real / np.sqrt((lambda_b_dutch[0].imag)**2+(lambda_b_dutch[0].real)**2)
    omega_n_dutch = np.sqrt((lambda_b_dutch[0].imag)**2+(lambda_b_dutch[0].real)**2) * V0/b * np.sqrt(1-zeta_dutch**2)

    return lambda_b_dutch[0], lambda_b_dutch[1], period_dutch, t_half_dutch, zeta_dutch, omega_n_dutch


######### Dutch roll motion  case 2 #############
def dutchroll_2(V0, m, rho, muc, mub, CL, CD, CX0, CZ0):
    A_dutch2 = -2 * mub * KZ2
    B_dutch2 = 1/2 * Cnr
    C_dutch2 = -Cnb

    lambda_b_dutch2_1 = (-B_dutch2 + np.sqrt(4 * A_dutch2 * C_dutch2 - B_dutch2 **2)) / (2 * A_dutch2)
    lambda_b_dutch2_2 = (-B_dutch2 - np.sqrt(4 * A_dutch2 * C_dutch2 - B_dutch2 **2)) / (2 * A_dutch2)
    return lambda_b_dutch2_1, lambda_b_dutch2_2

lambda_c_spm1, lambda_c_spm2, period_spm, t_half_spm, zeta_spm, omega_n_spm = short_period(V0_spm,m_spm,rho_spm,muc_spm,mub_spm,CL_spm,CD_spm,CX0_spm,CZ0_spm)
lambda_c_phug1, lambda_c_phug2, period_phug, t_half_phug, zeta_phug, omega_n_phug = phugoid(V0_spm,m_spm,rho_spm,muc_spm,mub_spm,CL_spm,CD_spm,CX0_spm,CZ0_spm)
lambda_b_apr, t_half_apr, time_cst_apr = aperiodic_roll(V0_spm,m_spm,rho_spm,muc_spm,mub_spm,CL_spm,CD_spm,CX0_spm,CZ0_spm)
lambda_b_spir, t_half_spir, time_cst_spir = spiral(V0_spm,m_spm,rho_spm,muc_spm,mub_spm,CL_spm,CD_spm,CX0_spm,CZ0_spm)
lambda_b_dutch1, lambda_b_dutch2, period_dutch, t_half_dutch, zeta_dutch, omega_n_dutch = dutchroll(V0_spm,m_spm,rho_spm,muc_spm,mub_spm,CL_spm,CD_spm,CX0_spm,CZ0_spm)

# lambda_c_spm1, lambda_c_spm2, period_spm, t_half_spm, zeta_spm, omega_n_spm = short_period(V0_spm,m_spm,rho_spm,muc_spm,mub_spm,CL_spm,CD_spm,CX0_spm,CZ0_spm)
# lambda_c_phug1, lambda_c_phug2, period_phug, t_half_phug, zeta_phug, omega_n_phug = phugoid(V0_phug,m_phug,rho_phug,muc_phug,mub_phug,CL_phug,CD_phug,CX0_phug,CZ0_phug)
# lambda_b_apr, t_half_apr, time_cst_apr = aperiodic_roll(V0_apr,m_apr,rho_apr,muc_apr,mub_apr,CL_apr,CD_apr,CX0_apr,CZ0_apr)
# lambda_b_spir, t_half_spir, time_cst_spir = spiral(V0_spir,m_spir,rho_spir,muc_spir,mub_spir,CL_spir,CD_spir,CX0_spir,CZ0_spir)
# lambda_b_dutch1, lambda_b_dutch2, period_dutch, t_half_dutch, zeta_dutch, omega_n_dutch = dutchroll(V0_dutch1,m_dutch1,rho_dutch1,muc_dutch1,mub_dutch1,CL_dutch1,CD_dutch1,CX0_dutch1,CZ0_dutch1)

print("short_period")
print("-----------------------")
print("eigenvalues:", lambda_c_spm1, lambda_c_spm2)
print("period:",period_spm)
print("t_half:", t_half_spm)
print("zeta:", zeta_spm)
print("omega_n:", omega_n_spm)
print()
print("phugoid")
print("-----------------------")
print("eigenvalues:", lambda_c_phug1, lambda_c_phug2)
print("period",period_phug)
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
print("eigenvalues:", lambda_b_dutch1, lambda_b_dutch2)
print("period",period_dutch)
print("t_half:", t_half_dutch)
print("zeta:", zeta_dutch)
print("omega_n:", omega_n_dutch)
