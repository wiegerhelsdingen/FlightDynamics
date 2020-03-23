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
    # u_transform = 1
    # q_transform = 1

    u_transform = 1/V0
    q_transform = c/V0

    C_1s = [[-2.*muc*c/V0*u_transform, 0, 0, 0],
            [0, (CZadot-2.*muc)*c/V0, 0, 0],
            [0, 0, -c/V0, 0],
            [0, Cmadot*c/V0, 0, -2.*muc*KY2*(c/V0)*q_transform]]

    C_2s = [[CXu*u_transform, CXa, CZ0, CXq*q_transform],
            [CZu*u_transform, CZa, -CX0, (CZq + 2.*muc)*q_transform],
            [0, 0, 0, 1.*q_transform],
            [Cmu*u_transform, Cma, 0, Cmq*q_transform]]

    C_3s = [[CXde],
            [CZde],
            [0],
            [Cmde]]

    #A and B-matrices
    C_1s_inv = np.linalg.inv(C_1s)

    A_s = np.matmul(-C_1s_inv, C_2s)
    B_s = np.matmul(-C_1s_inv, C_3s)
    C_s =           [[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]]
    D_s = [[0] ,
            [0] ,
            [0] ,
            [0]]

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


    print(eigenvals_s)
    print("-----------------------")
    print('Symmetric flight:')
    print("-----------------------")
    print('Eigenvalue 1: ', eigenvals_s[:2])
    print('Period: ', period_s1)
    print('halftime: ', t_half_s1 )
    print('zeta: ', zeta_s1)
    print('omega: ', omega_s1)
    print('')

    print('Eigenvalue 2: ', eigenvals_s[2:])
    print('Period: ', period_s2)
    print('halftime: ', t_half_s2 )
    print('zeta: ', zeta_s2)
    print('omega: ', omega_s2)
    #C1,2,3-matrices for Asymmetric case
    # p_transform = 1
    # r_transform = 1
    p_transform = b/(2*V0)
    r_transform = b/(2*V0)

    C_1a = [[(CYbdot-2*mub)*b/V0, 0, 0, 0],
             [0, -0.5*(b/V0), 0, 0 ],
             [0, 0, -4*mub*KX2*(b/V0)*p_transform, 4*mub*KXZ*(b/V0)*r_transform],
             [Cnbdot*b/V0, 0, 4*mub*KXZ*(b/V0)*p_transform, -4*mub*KZ2*(b/V0)*r_transform ]]

    C_2a = [[CYb, CL, CYp*p_transform , (CYr - 4*mub)*r_transform],
            [0        , 0   , 1*p_transform      , 0],
            [Clb, 0   , Clp*p_transform , Clr*r_transform],
            [Cnb, 0   , Cnp*p_transform , Cnr*r_transform]]

    C_3a = [[CYda, CYdr],
            [0, 0 ],
            [Clda, Cldr],
            [Cnda, Cndr]]

    #A and B-matrices
    C_1a_inv = np.linalg.inv(C_1a)

    A_a = np.matmul(-C_1a_inv, C_2a)
    B_a = np.matmul(-C_1a_inv,C_3a)
    C_a = [[1, 0, 0, 0],\
            [0, 1, 0, 1],\
            [0, 0, 1, 0],\
            [0, 0, 0, 1]]
    D_a = [ [ 0 , 0 ] ,
            [0 , 0] ,
            [0 , 0] ,
            [0 , 0]]

    sys_a = ctrl.ss(A_a, B_a, C_a, D_a)
    eigenvals_a = np.linalg.eigvals(A_a)
    period_a1 = 2*np.pi/abs(eigenvals_a[1])*b/V0
    t_half_a1 = np.log(0.5)/abs(eigenvals_a[1])*b/V0
    zeta_a1 = -eigenvals_s[0]/np.sqrt((eigenvals_a[1])**2+(eigenvals_a[0])**2)
    omega_a1= np.sqrt((eigenvals_a[1])**2+(eigenvals_a[0])**2)*V0/b*np.sqrt(1-zeta_a1**2)
    t_half_a2 = np.log(0.5)*b/(abs(eigenvals_a[2])*V0)
    t_cst_a2 = -b/(abs(eigenvals_a[2])*V0)
    t_half_a3 = np.log(0.5)*b/(eigenvals_a[3]*V0)
    t_cst_a3 = -b/(eigenvals_a[3]*V0)
    print(eigenvals_a)
    print("-----------------------")
    print('Asymmetric flight:')
    print("-----------------------")
    print('Eigenvalue 1: ', eigenvals_a[1:3])
    print('Period: ', period_a1)
    print('halftime: ', t_half_a1 )
    print('zeta: ', zeta_a1)
    print('omega: ', omega_a1)
    print('')

    print('Eigenvalue 2: ', eigenvals_a[0])
    print('half time :', t_half_a2)
    print('time constant: ', t_cst_a2)
    print('')

    print('Eigenvalue 3: ', eigenvals_a[3])
    print('half time :', t_half_a3)
    print('time constant: ', t_cst_a3)

    return sys_s, eigenvals_s, sys_a, eigenvals_a
