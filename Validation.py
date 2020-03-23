import numpy as np
from math import *
import matplotlib.pyplot as plt
from Data_extraction import *
import control as cntrl
from Numerical_Simulation import *
import control as cntrl

#Recorded time
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
hp = parameters[38][2] * 0.3048         #m
FF1 = parameters[5][2] * 0.000125997881 #kg/s
FF2 = parameters[5][2] * 0.000125997881 #kg/s
#Symmetric parameters
AoA = parameters[1][2]  * np.pi/180     #rad
Vt  = parameters[43][2] * 0.514444444   #m/s
Vc  = parameters[42][2] * 0.514444444   #m/s
th  = parameters[23][2] * np.pi/180     #rad
q   = parameters[28][2] * np.pi/180     #rad/s
#Asymmetric parameters
yaw     = parameters[24][2] * np.pi/180 #rad
roll    = parameters[22][2] * np.pi/180 #rad
rolldot = parameters[27][2] * np.pi/180 #rad/s
yawdot  = parameters[29][2] * np.pi/180 #rad/s
#Force parameters
de = parameters[18][2] * np.pi/180      #rad
dr = parameters[19][2] * np.pi/180      #rad
da = parameters[17][2] * np.pi/180      #rad
#Time interval of eigenmotions
spmt0    = 43*60+47.5
spm_t    = 10
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


def symplot(eigenmotion, time,t0,t1,AoA,Vt,Vc,th,q,force):
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
    fuel = np.sum(FF1[(time<t0)]*dt) + np.sum(FF2[(time<t0)]*dt)
    # print(hp0,V0,alpha0,th0)
    #Experimental data
    #Experimental data input
    AoA = AoA[interval]
    Vt = Vt[interval]
    Vc = Vc[interval]
    th = th[interval]
    q = q[interval]
    #Experimental Data Period, half time and eigenvalues
    filter = (abs(th-th0) < 10E-4)
    P = time[filter][1]
    # realeigen_exp =
    #force
    force = de[interval]
    #Numerical Model
    num_solution = numres(hp0,V0,alpha0,th0,fuel)
    sys = num_solution[0]
    num_response = ctrl.forced_response(sys, time, force)
    num_eigenval = num_solution[1]
    print(num_eigenval)
    num_response_V = num_response[1][0] +V0 #transformation
    num_response_alpha = num_response[1][1] +alpha0 #transformation
    num_response_th = num_response[1][2] +th0 #transformation
    num_response_q = num_response[1][3] +q0 #transformation
    # Plots
    fig, axs = plt.subplots(2, 2, constrained_layout=True)
    fig.suptitle(eigenmotion,fontsize=16)
    axs[0][0].plot(time, th, lw = 2, c = 'red', label='theta recorded')
    # axs[0][0].plot(time, num_response_th, lw = 2, c = 'blue', label='theta numerical')
    axs[0][0].set_ylabel('pitch angle recorded [rad]')
    axs[0][0].grid()
    axs[0][1].plot(time, AoA, lw=2, c='red', label='Ao recorded')
    # axs[0][1].plot(time, num_response_alpha, lw = 2, c = 'blue', label='$\alpha$ numerical')
    axs[0][1].set_ylabel('AoA recorded [rad]')
    axs[0][1].grid()
    axs[1][0].plot(time, Vt, lw=3, c='red', label='True airspeed recorded')
    # axs[1][0].plot(time, num_response_V, lw = 2, c = 'blue', label='True airspeed numerical')
    axs[1][0].set_ylabel('True airspeed recorded [m/s]')
    axs[1][0].grid()
    axs[1][1].plot(time, q, lw=3, c='red', label='pitch rate recorded')
    # axs[1][1].plot(time, num_response_q, lw = 2, c = 'blue', label='pitch rate numerical')
    axs[1][1].set_ylabel('pitch rate recorded [rad/s]')
    axs[1][1].grid()
    plt.show()
    return num_response, num_eigenval

def asymplot(eigenmotion,time,t0,t1,yaw,roll,rolldot,yawdot,force1,force2):
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

    #Experimental Data
    #Experimental Data input
    yaw = yaw[interval]
    roll = roll[interval]
    rolldot = rolldot[interval]
    yawdot = yawdot[interval]
    #Experimental Data Period, half time and eigenvalues
    if eigenmotion == 'aperiodic roll' or eigenmotion == 'spiral':
        filter = (yaw > (yaw0/np.exp(1)))
        t_const = time[filter][-1]
    # if eigenmotion == 'dutch roll 1' or eigenmotion == 'dutch roll 2':

    #Forces
    force1 = force1[interval]
    force2 = force2[interval]
    force = np.array([force1, force2])

    #Numerical rodel
    #Numerical response and eigenvalues
    num_solution = numres(hp0,V0,alpha0,th0,m)
    sys = num_solution[2]
    num_response = ctrl.forced_response(sys, time, force)
    num_eigenval = num_solution[3]
    print(num_eigenval)
    num_response_yaw = num_response[1][0] + yaw0
    num_response_roll = num_response[1][1] + roll0
    num_response_rolldot = num_response[1][2] + rolldot0
    numr_response_yawdot = num_response[1][3] + yawdot0
    # Plots
    # plt.subplot(411)
    # plt.plot(time, yaw, lw = 2, c = 'red', label='$\beta$ recorded')
    # plt.plot(time, num_response_yaw, lw = 2, c = 'blue', label='$beta$ numerical')
    # plt.ylabel('$\beta$ recorded [rad]')
    # plt.grid()
    # plt.subplot(412)
    # plt.plot(time, roll, lw=3, c='red', label='$\phi$ recorded')
    # plt.plot(time, num_response_roll, lw = 2, c = 'blue', label='$phi$ numerical')
    # plt.ylabel('$\phi$ recorded [rad]')
    # plt.grid()
    # plt.subplot(413)
    # plt.plot(time, rolldot, lw=3, c='red', label='roll rate recorded')
    # plt.plot(time, num_response_rolldot, lw = 2, c = 'blue', label='roll rate numerical')
    # plt.ylabel('roll rate recorded [m/s]')
    # plt.grid()
    # plt.subplot(414)
    # plt.plot(time, yawdot, lw=3, c='red', label='yaw rate recorded')
    # plt.plot(time, numr_response_yawdot, lw = 2, c = 'blue', label='yaw rate numerical')
    # plt.ylabel('yaw rate recorded [1/s]')
    # plt.grid()
    # plt.show()
    return num_response, num_eigenval

#Short period motion
numres_spm, numeigen_spm = symplot('short period',time,spmt0,spm_t,AoA,Vt,Vc,th,q,de)
#Phugoid motion
numres_phug, numeigen_spm = symplot('phugoid',time,phugt0,phug_t,AoA,Vt,Vc,th,q,de)
#Aperiodic roll motion
numres_apr, numeigen_apr = asymplot('aperiodic roll',time,aprt0,apr_t,yaw,roll,rolldot,yawdot,da,dr)
#Spiral motion
numres_spir, numeigen_spir = asymplot('spiral',time,spirt0,spir_t,yaw,roll,rolldot,yawdot,da,dr)
#Dutch roll motion
numres_dutch1, numeigen_dutch1 = asymplot('dutch roll 1',time,dutch1t0,dutch1_t,yaw,roll,rolldot,yawdot,da,dr)
numres_dutch2, numeigen_dutch2 = asymplot('dutch roll 2',time,dutch2t0,dutch2_t,yaw,roll,rolldot,yawdot,da,dr)
