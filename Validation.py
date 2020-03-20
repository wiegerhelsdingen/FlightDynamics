import numpy as np
from math import *
import matplotlib.pyplot as plt
from Data_extraction import *


time = parameters[0][2]
time_hrs = time/3600.
dt = 0.1

# Eigenmotions times
# spm = 43:46
# phug = 44:37
# apr  = 48:21
# dutch1 = 50:23
# dutch2 = 51:30
# spir = 55:01

#Symmetric parameters
AoA = parameters[1][2]
Vt  = parameters[43][2]
Vc  = parameters[42][2]
th  = parameters[23][2]
q   = parameters[28][2]


#Asymmetric Parameters
yaw     = parameters[24][2]
roll    = parameters[22][2]
rolldot = parameters[27][2]
yawdot  = parameters[29][2]

#Time interval of eigenmotions
spmt0    = 43*60+46
spm_t    = 3
phugt0   = 44*60+37
phug_t   = 100
aprt0    = 48*60+21
apr_t    = 100
dutch1t0 = 50*60+23
dutch1_t  = 40
dutch2t0 = 51*60+30
dutch2_t = 40
spirt0   = 55*60+1
spir_t   = 100


def symplot(time,t0,t1,AoA,Vt,Vc,th,q):
    interval0 = (time > t0)
    interval1 = (time < (t0+t1))
    interval = interval0
    for i in range(len(interval)):
        if interval0[i] == False or interval1[i] == False:
            interval[i] = False
        else:
            interval[i] = True

    time = np.arange(0,t1-dt,dt)
    AoA = AoA[interval]
    Vt = Vt[interval]
    Vc = Vc[interval]
    th = th[interval]
    q = q[interval]

    plt.subplot(411)
    plt.plot(time, th, lw = 3, c = 'red', label='$\theta$ recorded')
    plt.ylabel('deg [-]')
    plt.grid()

    plt.subplot(412)
    plt.plot(time, AoA, lw=3, c='red', label='$\alpha$ recorded')
    plt.ylabel('deg [-]')
    plt.grid()

    plt.subplot(413)
    plt.plot(time, Vt, lw=3, c='red', label='True airspeed recorded')
    plt.ylabel('m/s')
    plt.grid()

    plt.subplot(414)
    plt.plot(time, q, lw=3, c='red', label='pitch rate recorded')
    plt.ylabel('deg/sec [1/s]')
    plt.grid()

    plt.show()
    return time, AoA, Vt, Vc, th, q

def asymplot(time,t0,t1,yaw,roll,rolldot,yawdot):
    interval0 = (time > t0)
    interval1 = (time < (t0+t1))
    interval = interval0
    for i in range(len(interval)):
        if interval0[i] == False or interval1[i] == False:
            interval[i] = False
        else:
            interval[i] = True

    time = np.arange(0,t1-dt,dt)
    yaw = yaw[interval]
    roll = roll[interval]
    rolldot = rolldot[interval]
    yawdot = yawdot[interval]

    plt.subplot(411)
    plt.plot(time, yaw, lw = 3, c = 'red', label='$\beta$ recorded')
    plt.ylabel('deg [-]')
    plt.grid()

    plt.subplot(412)
    plt.plot(time, roll, lw=3, c='red', label='$\phi$ recorded')
    plt.ylabel('deg [-]')
    plt.grid()

    plt.subplot(413)
    plt.plot(time, rolldot, lw=3, c='red', label='roll rate recorded')
    plt.ylabel('m/s')
    plt.grid()

    plt.subplot(414)
    plt.plot(time, yawdot, lw=3, c='red', label='yaw rate recorded')
    plt.ylabel('deg/sec [1/s]')
    plt.grid()

    plt.show()

    return time, yaw, roll, rolldot, yawdot


#Short period motion
time_spm, AoA_spm, Vt_spm, Vc_spm, th_spm, q_spm = symplot(time,spmt0,spm_t,AoA,Vt,Vc,th,q)
#Phugoid motion
time_phug, AoA_phug, Vt_phug, Vc_phug, th_phug, q_phug = symplot(time,phugt0,phug_t,AoA,Vt,Vc,th,q)
#Aperiodic roll motion
time_apr, yaw_apr, roll_apr, rolldot_apr, yawdot_apr = asymplot(time,aprt0,apr_t,yaw,roll,rolldot,yawdot)
#Spiral motion
time_spir, yaw_spir, roll_spir, rolldot_spir, yawdot_spir = asymplot(time,spirt0,spir_t,yaw,roll,rolldot,yawdot)
#Dutch roll motion
time_dutch1, yaw_dutch1, roll_dutch1, rolldot_dutch1, yawdot_dutch1 = asymplot(time,dutch1t0,dutch1_t,yaw,roll,rolldot,yawdot)
time_dutch2, yaw_dutch2, roll_dutch2, rolldot_dutch2, yawdot_dutch2 = asymplot(time,dutch2t0,dutch2_t,yaw,roll,rolldot,yawdot)
