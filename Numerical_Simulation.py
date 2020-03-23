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


def numres(eigenmotion,V0,m,rho,muc,mub,CL,CD,CX0,CZ0):

    #C1,2,3-matrices for symmetric case
    u_transform = 1/V0
    q_transform = c/V0

    C_1s = np.array([[-2*muc*c/V0*u_transform, 0, 0, 0],
            [0, (CZadot-2*muc)*c/V0, 0, 0],
            [0, 0, -c/V0, 0],
            [0, Cmadot*c/V0, 0, -2*muc*KY2**2*(c/V0)*q_transform]], dtype='float')

    C_2s = - np.array([[CXu*u_transform, CXa, CZ0, CXq*q_transform],
            [CZu*u_transform, CZa, -CX0, (CZq + 2*muc)*q_transform],
            [0, 0, 0, 1*q_transform],
            [Cmu*u_transform, Cma, 0, Cmq*q_transform]], dtype='float')

    C_3s = np.array([[CXde],
            [CZde],
            [0],
            [Cmde]], dtype='float')

    #A and B-matrices
    C_1s_inv = np.linalg.inv(C_1s)

    A_s = np.matmul(-C_1s_inv, C_2s)
    B_s = np.matmul(-C_1s_inv, C_3s)
    C_s = np.array([[1, 0, 0, 0],\
                    [0, 1, 0, 0],\
                    [0, 0, 1, 0],\
                    [0, 0, 0, 1]])
    D_s = np.zeros((4,1))

    sys_s = ctrl.ss(A_s, B_s, C_s, D_s)
    eigenvals_s = np.linalg.eigvals(A_s)
    period_s1 = 2*np.pi/abs(eigenvals_s[1])*c/V0
    t_half_s1 = np.log(0.5)/abs(eigenvals_s[1])*c/V0
    zeta_s1 = -eigenvals_s[0]/np.sqrt((eigenvals_s[1])**2+(eigenvals_s[0])**2)
    omega_s1= np.sqrt((eigenvals_s[1])**2+(eigenvals_s[0])**2)*V0/c*np.sqrt(1-zeta_s1**2)
    t_half_s2 = np.log(0.5)/abs(eigenvals_s[3])*c/V0
    period_s2 = 2*np.pi/abs(eigenvals_s[3])*c/V0
    zeta_s2 = -eigenvals_s[2]/np.sqrt((eigenvals_s[3])**2+(eigenvals_s[2])**2)
    omega_s2= np.sqrt((eigenvals_s[3])**2+(eigenvals_s[2])**2)*V0/c*np.sqrt(1-zeta_s2**2)
    print('Symmetric flight:')
    print('Eigenvalue 1: ', eigenvals_s[:2])
    print('Period, halftime, zeta, omega: ', period_s1,t_half_s1,zeta_s1,omega_s1)
    print('Eigenvalue 2: ', eigenvals_s[2:])
    print('Period, halftime, zeta, omega: ', period_s2,t_half_s2,zeta_s2,omega_s2)


    #C1,2,3-matrices for Asymmetric case
    p_transform = b/(2*V0)
    r_transform = b/(2*V0)

    C_1a = np.array([[(CYbdot-2*mub)*b/V0, 0, 0, 0],
             [0, -0.5*(b/V0), 0, 0 ],
             [0, 0, -4*mub*KX2**2*(b/V0)*p_transform, 4*mub*KXZ*(b/V0)*r_transform],
             [Cnb, 0, 4*mub*KXZ*(b/V0)*p_transform, -4*mub*KZ2**2*(b/V0)*r_transform ]], dtype='float')

    C_2a = np.array([[CYb, CL, CYp*p_transform , (CYr - 4*mub)*r_transform],
            [0        , 0   , 1*p_transform      , 0],
            [Clb, 0   , Clp*p_transform , Clr*r_transform],
            [Cnb, 0   , Cnp*p_transform , Cnr*r_transform]], dtype='float')

    C_3a = np.array([[CYda, CYdr],
            [0, 0 ],
            [Clda, Cldr],
            [Cnda, Cndr]], dtype='float')

    #A and B-matrices
    C_1a_inv = np.linalg.inv(C_1a)

    A_a = np.matmul(-C_1a_inv, C_2a)
    B_a = np.matmul(-C_1a_inv,C_3a)
    C_a = np.array([[1, 0, 0, 0],\
                    [0, 1, 0, 1],\
                    [0, 0, 1, 0],\
                    [0, 0, 0, 1]])
    D_a = np.zeros((4,2))

    sys_a = ctrl.ss(A_a, B_a, C_a, D_a)
    eigenvals_a = np.linalg.eigvals(A_a)
    period_a1 = 2*np.pi/abs(eigenvals_a[1])*b/V0
    t_half_a1 = np.log(0.5)/abs(eigenvals_a[1])*b/V0
    zeta_a1 = -eigenvals_s[0]/np.sqrt((eigenvals_a[1])**2+(eigenvals_a[0])**2)
    omega_a1= np.sqrt((eigenvals_a[1])**2+(eigenvals_a[0])**2)*V0/b*np.sqrt(1-zeta_a1**2)
    t_half_a2 = np.log(0.5)*b/(eigenvals_a[2]*V0)
    t_cst_a2 = -b/(eigenvals_a[2]*V0)
    t_half_a3 = np.log(0.5)*b/(eigenvals_a[3]*V0)
    t_cst_a3 = -b/(eigenvals_a[3]*V0)
    print('Asymmetric flight:')
    print('Eigenvalue 1: ' , eigenvals_a[:2])
    print('Period, halftime, zeta, omega: ', period_a1,t_half_a1,zeta_a1,omega_a1)
    print('Eigenvalue 2: ', eigenvals_a[2])
    print('Period, halftime, zeta, omega: ', t_half_a2,t_cst_a2)
    print('Eigenvalue 3: ', eigenvals_a[3])
    print('Period, halftime, zeta, omega: ', t_half_a3,t_cst_a3)

    return sys_s, eigenvals_s, sys_a, eigenvals_a

num_spm   = numres('spm',V0_spm,m_spm,rho_spm,muc_spm,mub_spm,CL_spm,CD_spm,CX0_spm,CZ0_spm)
# num_phug  = numres('phug',V0_phug,m_phug,rho_phug,muc_phug,mub_phug,CL_phug,CD_phug,CX0_phug,CZ0_phug)
# num_apr   = numres('apr',V0_apr,m_apr,rho_apr,muc_apr,mub_apr,CL_apr,CD_apr,CX0_apr,CZ0_apr)
# num_spir  = numres('spir',V0_spir,m_spir,rho_spir,muc_spir,mub_spir,CL_spir,CD_spir,CX0_spir,CZ0_spir)
# num_dutch1= numres('dutch1',V0_dutch1,m_dutch1,rho_dutch1,muc_dutch1,mub_dutch1,CL_dutch1,CD_dutch1,CX0_dutch1,CZ0_dutch1)
