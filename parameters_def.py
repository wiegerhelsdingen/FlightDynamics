import numpy as np
from math import *
from Data_extraction import *
import sys


# ----------------------------------------------------------------------------
#Flight independent parameters
# ----------------------------------------------------------------------------
# mass
W_empty = 9165*0.453592 #kg
bf_kg = 4100*0.453592  #kg
masspas = np.array([90,102,80,83,94,84,74,79,103]) #kg
# aerodynamic properties
e      = 0.8         # Oswald factor [ ]
CD0    = 0.04        # Zero lift drag coefficient [ ]
CLa    = 4.68 #5.084       # Slope of CL-alpha curve [ ]
# Longitudinal stability
Cma    = -0.58256477 #-0.5626     # longitudinal stabilty [ ]
Cmde   = -1.21166653 #-1.1642     # elevator effectiveness [ ]
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

# Constant values concerning aircraft inertia
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.3925
# Aerodynamic constants
Cmac   = 0                      # Moment coefficient about the aerodynamicf centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * np.pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

CXu    = -0.085
CXa    = +0.47966		# Positive! (has been erroneously negative since 1993)
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZu    = -0.61616
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

Clb    = -0.2260
Clp    = -0.81085
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.0968
Cnbdot =   0
Cnp    =  -0.0602
Cnr    =  -0.3061
Cnda   =  -0.0120
Cndr   =  -0.0939

# ----------------------------------------------------------------------------
#Flight dependent parameters
# ----------------------------------------------------------------------------

time = parameters[0][2]

# Eigenmotions times
# spm = 43:46
# phug = 44:37
# apr  = 48:21
# dutch1 = 50:23
# dutch2 = 51:30
# spir = 55:01

#Time interval of eigenmotions
spmt0    = 43*60+47.5
spm_t    = 10
# time =
phugt0   = 44*60+37
phug_t   = 180
aprt0    = 48*60+21
apr_t    = 100
dutch1t0 = 50*60+23
dutch1_t  = 40
dutch2t0 = 51*60+30
dutch2_t = 40
spirt0   = 55*60+1
spir_t   = 100

# Stationary parameters
hp = parameters[38][2] * 0.3048         #m
FF1 = parameters[5][2] * 0.000125997881 #kg/s
FF2 = parameters[5][2] * 0.000125997881 #kg/s
AoA = parameters[1][2]  * np.pi/180     #rad
Vt  = parameters[43][2] * 0.514444444   #m/s
th  = parameters[23][2] * np.pi/180     #rad

def dependent_parameters(hp0,V0,alpha0,th0,fuel):
    m      =      W_empty + bf_kg + np.sum(masspas) - fuel    # mass [kg]
    rho0   = 1.2250          # air density at sea level [kg/m^3]
    a = -0.0065         # temperature gradient in ISA [K/m]
    Temp0  = 288.15          # temperature at sea level in ISA [K]
    R      = 287.05          # specific gas constant [m^2/sec^2K]
    g      = 9.81            # [m/sec^2] (gravity constant)
    rho    = rho0 * ((1+(a * hp0 / Temp0)))**(-((g / (a*R)) + 1))
    W      = m * g            # [N]       (aircraft weight)
    muc    = m / (rho * S * c)
    mub    = m / (rho * S * b)
    CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
    CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e) # Drag coefficient [ ]
    CX0    = W * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
    CZ0    = -W * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
    return m,rho,muc,mub,CL,CD,CX0,CZ0

def eigenmotion_parameters(time,t0):
    dt = 0.1
    i_t0 = int(np.where(time==t0)[0])
    #Stationary flight Parameters
    hp0 = hp[i_t0][0]
    alpha0 = AoA[i_t0][0]
    th0 = th[i_t0][0]
    V0 = Vt[i_t0][0]
    fuel = np.sum(FF1[(time<t0)]*dt) + np.sum(FF2[(time<t0)]*dt)
    m,rho,muc,mub,CL,CD,CX0,CZ0 = dependent_parameters(hp0,V0,alpha0,th0,fuel)
    return V0,m,rho,muc,mub,CL,CD,CX0,CZ0
