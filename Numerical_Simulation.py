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


def numres(hp0,V0,alpha0,th0):
    # Aircraft mass
    m      =      1000       # mass [kg]

    # aerodynamic properties
    e      = 0.8         # Oswald factor [ ]
    CD0    = 0.04        # Zero lift drag coefficient [ ]
    CLa    = 5.084       # Slope of CL-alpha curve [ ]

    # Longitudinal stability
    Cma    = -0.5626     # longitudinal stabilty [ ]
    Cmde   = -1.1642     # elevator effectiveness [ ]

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
    a = -0.0065         # temperature gradient in ISA [K/m]
    Temp0  = 288.15          # temperature at sea level in ISA [K]
    R      = 287.05          # specific gas constant [m^2/sec^2K]
    g      = 9.81            # [m/sec^2] (gravity constant)

    # air density [kg/m^3]
    rho    = rho0 * ((1+(a * hp0 / Temp0)))**(-((g / (a*R)) + 1))
    W      = m * g            # [N]       (aircraft weight)

    # Constant values concerning aircraft inertia
    muc    = m / (rho * S * c)
    mub    = m / (rho * S * b)
    KX2    = 0.019
    KZ2    = 0.042
    KXZ    = 0.002
    KY2    = 1.25 * 1.114
    KX2_squared = KX2**2
    KZ2_squared = KZ2**2
    KXZ_squared = KXZ**2
    KY2_squared = KY2**2

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

    Cm0    =  0.02970
    Cmu    = +0.06990
    Cmadot = +0.17800
    Cmq    = -8.79415
    CmTc   = -0.00064

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
    # print('symmetric eigenvalues' , np.linalg.eigvals(A_s))


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

    return sys_s, eigenvals_s, sys_a, eigenvals_a

print('Asymmetric eigenvalues' , np.linalg.eig(A_a)[0])
