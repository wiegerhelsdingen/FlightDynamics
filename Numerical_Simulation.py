# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:08:35 2020

@author: thijs
"""
import numpy as np
import math
import control as ctrl
import matplotlib.pyplot as plt

# Citation 550 - Linear simulation

# xcg = 0.25 * c

# Stationary flight condition

hp0    =  2     	      # pressure altitude in the stationary flight condition [m]
V0     =   2          # true airspeed in the stationary flight condition [m/sec]
alpha0 =    2         # angle of attack in the stationary flight condition [rad]
th0    =     2        # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      =   2          # mass [kg]

# aerodynamic properties
e      =   2          # Oswald factor [ ]
CD0    =   2          # Zero lift drag coefficient [ ]
CLa    =   2          # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma    =   2          # longitudinal stabilty [ ]
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
ih     = -2 * math.pi / 180   # stabiliser angle of incidence [rad]

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
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.25 * 1.114

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * math.pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Lift and drag coefficient

CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (math.pi * A * e) # Drag coefficient [ ]

# Stabiblity derivatives

CX0    = W * math.sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu    = -0.02792
CXa    = +0.47966		# Positive! (has been erroneously negative since 1993)
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZ0    = -W * math.cos(th0) / (0.5 * rho * V0 ** 2 * S)
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


#C1,2,3-matrices for symmetric case

C_1s = np.array([[-2*muc*c/V0**2, 0, 0, 0],
        [0, (CZadot-2*muc)*c/V0, 0, 0],
        [0, 0, -c/V0, 0],
        [0, Cmadot*c/V0, 0, -2*muc*KY2**2*(c**2/V0**2)]])

C_2s = - np.array([[CXu, CXa, CZ0, CXq],
        [CZu/V0, CZa, -CX0, c/V0*(CZq + 2*muc)],
        [0, 0, 0, c/V0],
        [Cmu/V0, Cma, 0, c/V0*Cmq]])

C_3s = np.array([[CXde],
        [CZde],
        [0],
        [Cmde]])

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

# print('symmetric eigenvalues' , np.linalg.eigvals(A_s))


#C1,2,3-matrices for Asymmetric case

C_1a = [[(CYbdot-2*mub)*b/V0, 0, 0, 0],
         [0, -0.5*(b/V0), 0, 0 ],
         [0, 0, -4*mub*KX2**2*(b/V0)*(b/(2*V0)), 4*mub*KXZ*(b/V0)*(b/(2*V0))],
         [Cnb, 0, 4*mub*KXZ*(b/V0)*(b/(2*V0)), -4*mub*KZ2**2*(b/V0)*(b/(2*V0)) ]]

C_2a = [[CYb, CL, CYp*(b/(2*V0)) , (CYr - 4*mub)*(b/(2*V0))],
        [0        , 0   , (b/(2*V0))       , 0],
        [Clb, 0   , Clp*(b/(2*V0)) , Clr*(b/(2*V0))],
        [Cnb, 0   , Cnp*(b/(2*V0)) , Cnr*(b/(2*V0))]]


#Geen tweede kolom??
C_3a = [[CYda, CYdr],
        [0, 0 ],
        [Clda, Cldr],
        [Cnda, Cndr]]

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

# print('Asymmetric eigenvalues' , np.linalg.eig(A_a)[0])
