import numpy as np
from math import *
import matplotlib.pyplot as plt
from Data_extraction import *
import control as cntrl
from Numerical_Simulation import *
import control as cntrl
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

# Stationary parameters
hp = parameters[38][2] * 0.3048
#Symmetric parameters
AoA = parameters[1][2]  * np.pi/180
Vt  = parameters[43][2] * 0.514444444
Vc  = parameters[42][2] * 0.514444444
th  = parameters[23][2] * np.pi/180
q   = parameters[28][2] * np.pi/180
#Asymmetric parameters
yaw     = parameters[24][2] * np.pi/180
roll    = parameters[22][2] * np.pi/180
rolldot = parameters[27][2] * np.pi/180
yawdot  = parameters[29][2] * np.pi/180
#Force parameters
de = parameters[18][2] * np.pi/180
dr = parameters[19][2] * np.pi/180
da = parameters[17][2] * np.pi/180

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


def symplot(time,t0,t1,AoA,Vt,Vc,th,q,force):
    #Specific time interval for eigenmotion
    interval0 = (time > t0)
    interval1 = (time < (t0+t1))
    interval = interval0
    for i in range(len(interval)):
        if interval0[i] == False or interval1[i] == False:
            interval[i] = False
        else:
            interval[i] = True
    #Time
    i_t0 = int(np.where(time==t0)[0])
    time = np.arange(0,t1-dt,dt)
    #Stationary Parameters
    hp0 = hp[i_t0]
    alpha0 = AoA[i_t0]
    th0 = th[i_t0]
    V0 = Vt[i_t0]
    q0 = q[i_t0]
    print(V0)
    #Experimental data
    AoA = AoA[interval]
    Vt = Vt[interval]
    Vc = Vc[interval]
    th = th[interval]
    q = q[interval]
    #force
    force = de[interval]
    #Numerical Model
    numres = ctrl.forced_response(sys_s, time, force)
    numres_V = numres[1][0] + V0
    numres_alpha = numres[1][1] + alpha0
    numres_th = numres[1][2] + th0
    numres_q = numres[1][3] + q0
    #Plots
    # plt.subplot(411)
    # plt.plot(time, th, lw = 2, c = 'red', label='$\theta$ recorded')
    # plt.plot(time, numres_th, lw = 2, c = 'blue', label='$\theta$ numerical')
    # plt.ylabel('deg [-]')
    # plt.grid()
    # plt.subplot(412)
    # plt.plot(time, AoA, lw=2, c='red', label='$\alpha$ recorded')
    # plt.plot(time, numres_alpha, lw = 2, c = 'blue', label='$\alpha$ numerical')
    # plt.ylabel('deg [-]')
    # plt.grid()
    # plt.subplot(413)
    # plt.plot(time, Vt, lw=3, c='red', label='True airspeed recorded')
    # plt.plot(time, numres_V, lw = 2, c = 'blue', label='True airspeed numerical')
    # plt.ylabel('m/s')
    # plt.grid()
    # plt.subplot(414)
    # plt.plot(time, q, lw=3, c='red', label='pitch rate recorded')
    # plt.plot(time, numres_q, lw = 2, c = 'blue', label='pitch rate numerical')
    # plt.ylabel('deg/sec [1/s]')
    # plt.grid()
    # plt.show()
    return time, AoA, Vt, Vc, th, q, force

def asymplot(time,t0,t1,yaw,roll,rolldot,yawdot,force1,force2):

    #Specific time interval for eigenmotion
    interval0 = (time > t0)
    interval1 = (time < (t0+t1))
    interval = interval0
    for i in range(len(interval)):
        if interval0[i] == False or interval1[i] == False:
            interval[i] = False
        else:
            interval[i] = True
    #Time
    i_t0 = int(np.where(time==t0)[0])
    time = np.arange(0,t1-dt,dt)
    #Stationary parameters
    hp0 = hp[i_t0]
    alpha0 = AoA[i_t0]
    th0 = th[i_t0]
    V0 = Vt[i_t0]
    yaw0 = yaw[i_t0]
    roll0 = roll[i_t0]
    rolldot0 = rolldot[i_t0]
    yawdot0 = yawdot[i_t0]
    print(V0)
    #Experimental Data
    yaw = yaw[interval]
    roll = roll[interval]
    rolldot = rolldot[interval]
    yawdot = yawdot[interval]
    #Forces
    force1 = force1[interval]
    force2 = force2[interval]
    force = np.array([force1, force2])
    #Numerical Model
    numres = ctrl.forced_response(sys_a, time, force)
    numres_yaw = numres[1][0] + yaw0
    numres_roll = numres[1][1] + roll0
    numres_rolldot = numres[1][2] + rolldot0
    numres_yawdot = numres[1][3] + yawdot0
    # Plots
    # plt.subplot(411)
    # plt.plot(time, yaw, lw = 2, c = 'red', label='$\beta$ recorded')
    # plt.plot(time, numres_yaw, lw = 2, c = 'blue', label='$\beta$ numerical')
    # plt.ylabel('deg [-]')
    # plt.grid()
    # plt.subplot(412)
    # plt.plot(time, roll, lw=3, c='red', label='$\phi$ recorded')
    # plt.plot(time, numres_roll, lw = 2, c = 'blue', label='$\phi$ numerical')
    # plt.ylabel('deg [-]')
    # plt.grid()
    # plt.subplot(413)
    # plt.plot(time, rolldot, lw=3, c='red', label='roll rate recorded')
    # plt.plot(time, numres_rolldot, lw = 2, c = 'blue', label='roll rate numerical')
    # plt.ylabel('m/s')
    # plt.grid()
    # plt.subplot(414)
    # plt.plot(time, yawdot, lw=3, c='red', label='yaw rate recorded')
    # plt.plot(time, numres_yawdot, lw = 2, c = 'blue', label='yaw rate numerical')
    # plt.ylabel('deg/sec [1/s]')
    # plt.grid()
    # plt.show()
    return time, yaw, roll, rolldot, yawdot, force1, force2

#Short period motion
time_spm, AoA_spm, Vt_spm, Vc_spm, th_spm, q_spm, de_spm = symplot(time,spmt0,spm_t,AoA,Vt,Vc,th,q,de)
#Phugoid motion
time_phug, AoA_phug, Vt_phug, Vc_phug, th_phug, q_phug, de_phug = symplot(time,phugt0,phug_t,AoA,Vt,Vc,th,q,de)
#Aperiodic roll motion
time_apr, yaw_apr, roll_apr, rolldot_apr, yawdot_apr, da_apr, dr_apr = asymplot(time,aprt0,apr_t,yaw,roll,rolldot,yawdot,da,dr)
#Spiral motion
time_spir, yaw_spir, roll_spir, rolldot_spir, yawdot_spir, da_spir, dr_spir = asymplot(time,spirt0,spir_t,yaw,roll,rolldot,yawdot,da,dr)
#Dutch roll motion
time_dutch1, yaw_dutch1, roll_dutch1, rolldot_dutch1, yawdot_dutch1, da_dutch1, dr_dutch1 = asymplot(time,dutch1t0,dutch1_t,yaw,roll,rolldot,yawdot,da,dr)
time_dutch2, yaw_dutch2, roll_dutch2, rolldot_dutch2, yawdot_dutch2, da_dutch2, dr_dutch2 = asymplot(time,dutch2t0,dutch2_t,yaw,roll,rolldot,yawdot,da,dr)
