# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:08:35 2020

@author: thijs
"""
import numpy as np
import math
import control as ctrl
import matplotlib.pyplot as plt
import control as cntrl
from parameters import *


def numres(V0,m,rho,muc,mub,CL,CD,CX0,CZ0):

    #C1,2,3-matrices for symmetric case

    C_1s = np.matrix([[-2.*muc*c/V0, 0, 0, 0],
            [0, (CZadot-2.*muc)*c/V0, 0, 0],
            [0, 0, -c/V0, 0],
            [0, Cmadot*c/V0, 0, -2.*muc*KY2*(c/V0)]])

    C_2s = np.matrix([[CXu, CXa, CZ0, CXq],
            [CZu, CZa, -CX0, (CZq + 2.*muc)],
            [0, 0, 0, 1.],
            [Cmu, Cma, 0, Cmq]])

    C_3s = np.matrix([[CXde],
            [CZde],
            [0],
            [Cmde]])

    #A and B-matrices
    C_1s_inv = np.linalg.inv(C_1s)

    A_s = np.matmul(-C_1s_inv, C_2s)
    B_s = np.matmul(-C_1s_inv, C_3s)
    C_s =           np.matrix([[1, 0, 0, 0],
                    [0, 1, 0, 1],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]])
    D_s = np.matrix([[0] ,
            [0] ,
            [0] ,
            [0]])


    sys_s = ctrl.ss(A_s, B_s, C_s, D_s)

    eigenvals_s = np.linalg.eigvals(A_s)
    period_s1 = 2*np.pi/eigenvals_s[0].imag
    t_half_s1 = np.log(0.5)/eigenvals_s[0].real
    zeta_s1 = -eigenvals_s[0].real/np.sqrt((eigenvals_s[0].imag)**2+(eigenvals_s[0].real)**2)
    omega_s1= np.sqrt((eigenvals_s[0].imag)**2+(eigenvals_s[0].real)**2)*V0/c*np.sqrt(1-zeta_s1**2)

    period_s2 = 2*np.pi/eigenvals_s[2].imag
    t_half_s2 = np.log(0.5)/eigenvals_s[2].real
    zeta_s2 = -eigenvals_s[2].real/np.sqrt((eigenvals_s[2].imag)**2+(eigenvals_s[2].real)**2)
    omega_s2= np.sqrt((eigenvals_s[2].imag)**2+(eigenvals_s[2].real)**2)*V0/c*np.sqrt(1-zeta_s2**2)

    num_sym_par = []
    num_sym_par.append([eigenvals_s[:2],period_s1,t_half_s1,zeta_s1,omega_s1])
    num_sym_par.append([eigenvals_s[2:],period_s2,t_half_s2,zeta_s2,omega_s2])
    #C1,2,3-matrices for Asymmetric case

    C_1a = np.matrix([[(CYbdot-2*mub)*b/V0, 0, 0, 0],
             [0, -0.5*(b/V0), 0, 0 ],
             [0, 0, -4*mub*KX2*(b/V0), 4*mub*KXZ*(b/V0)],
             [Cnbdot*b/V0, 0, 4*mub*KXZ*(b/V0), -4*mub*KZ2*(b/V0) ]])

    C_2a = np.matrix([[CYb, CL, CYp , (CYr - 4*mub)],
            [0        , 0   , 1      , 0],
            [Clb, 0   , Clp , Clr],
            [Cnb, 0   , Cnp , Cnr]])

    C_3a = -np.matrix([[CYda, CYdr],
            [0, 0 ],
            [Clda, Cldr],
            [Cnda, Cndr]])

    #A and B-matrices
    C_1a_inv = np.linalg.inv(C_1a)

    A_a = np.matmul(-C_1a_inv, C_2a)
    B_a = np.matmul(-C_1a_inv,C_3a)
    C_a = np.matrix([[1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]])
    D_a = np.matrix([ [ 0 , 0 ] ,
            [0 , 0] ,
            [0 , 0] ,
            [0 , 0]])

    sys_a = ctrl.ss(A_a, B_a, C_a, D_a)
    eigenvals_a = np.linalg.eigvals(A_a)
    period_a1 = 2*np.pi/eigenvals_a[1].imag
    t_half_a1 = np.log(0.5)/eigenvals_a[1].real
    zeta_a1 = -eigenvals_s[1].real/np.sqrt((eigenvals_a[1].imag)**2+(eigenvals_a[1].real)**2)
    omega_a1= np.sqrt((eigenvals_a[1].imag)**2+(eigenvals_a[1].real)**2)*V0/b*np.sqrt(1-zeta_a1**2)
    t_half_a2 = np.log(0.5)/(eigenvals_a[0].real)
    t_cst_a2 = -1/(eigenvals_a[0].real)
    t_half_a3 = -np.log(0.5)/(eigenvals_a[3].real)
    t_cst_a3 = 1/(eigenvals_a[3].real)

    num_asym_par = []
    num_asym_par.append([eigenvals_a[1:3],period_a1,t_half_a1,zeta_a1,omega_a1])
    num_asym_par.append([eigenvals_a[0],t_half_a2,t_cst_a2])
    num_asym_par.append([eigenvals_a[3],t_half_a3,t_cst_a3])
    return sys_s, num_sym_par, sys_a, num_asym_par
