import numpy as np
from math import *
import matplotlib.pyplot as plt
from Data_extraction import *
import control as cntrl
from Numerical_Simulation import *
import control as cntrl
from parameters_def import *

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
spm_time = np.arange(0,spm_t-dt,dt)
phugt0   = 44*60+37
phug_t   = 150
phug_time= np.arange(0,phug_t-dt,dt)
aprt0    = 48*60+21
apr_t    = 20
apr_time = np.arange(0,apr_t-dt,dt)
dutch1t0 = 50*60+23
dutch1_t  = 20
dutch1_time = np.arange(0,dutch1_t-dt,dt)
dutch2t0 = 51*60+30
dutch2_t = 20
dutch2_time = np.arange(0,dutch2_t-dt,dt)
spirt0   = 55*60+1
spir_t   = 100
spir_time = np.arange(0,spir_t-dt,dt)


def symplot(eigenmotion, time,t0,t1,AoA,Vt,Vc,th,q,force):
    #------------------------------------
    #Specific time interval for eigenmotion
    #------------------------------------
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
    time_eigen = np.arange(0,t1-dt,dt)
    #------------------------------------
    #Experimental Data
    #------------------------------------
    #Stationary flight Parameters
    alpha0 = AoA[i_t0][0]
    th0 = th[i_t0][0]
    V0 = Vt[i_t0][0]
    q0 = q[i_t0][0]
    #force
    force = de[interval]

    #------------------------------------
    #Experimental response
    #------------------------------------
    AoA = AoA[interval]
    Vt = Vt[interval]
    Vc = Vc[interval]
    th = th[interval]
    q = q[interval]

    #------------------------------------
    #Numerical response
    #------------------------------------
    #Flight Dependent parameters
    V0,m,rho,muc,mub,CL,CD,CX0,CZ0 = eigenmotion_parameters(time,t0)
    #Numerical solution
    num_solution = numres(V0,m,rho,muc,mub,CL,CD,CX0,CZ0)
    sys = num_solution[0]
    #Response
    num_response = ctrl.forced_response(sys, time_eigen, force)
    #Eigenvalues
    num_eigenval = num_solution[1]
    num_response_V = num_response[1][0]*V0 + V0 #transformation
    num_response_alpha = num_response[1][1] + alpha0 #transformation
    num_response_th = num_response[1][2] + th0 #transformation
    num_response_q = num_response[1][3] *V0/c #transformation
    # Plots
    fig, axs = plt.subplots(2, 2, constrained_layout=True)
    fig.suptitle(eigenmotion,fontsize=16)
    axs[0][0].plot(time_eigen, th, lw = 1, c = 'orange', label='ecorded')
    axs[0][0].plot(time_eigen, num_response_th, lw = 1, c = 'deepskyblue', label='numerical')
    axs[0][0].set_ylabel('pitch angle [rad]')
    axs[0][0].grid()
    axs[0][1].plot(time_eigen, AoA, lw=1, c='orange', label='recorded')
    axs[0][1].plot(time_eigen, num_response_alpha, lw = 1, c = 'deepskyblue', label='numerical')
    axs[0][1].set_ylabel('AoA [rad]')
    axs[0][1].grid()
    axs[1][0].plot(time_eigen, Vt, lw=1, c='orange', label='recorded')
    axs[1][0].plot(time_eigen, num_response_V, lw = 1, c = 'deepskyblue', label='numerical')
    axs[1][0].set_ylabel('True airspeed [m/s]')
    axs[1][0].grid()
    axs[1][1].plot(time_eigen, q, lw=1, c='orange', label='recorded')
    axs[1][1].plot(time_eigen, num_response_q, lw = 1, c = 'deepskyblue', label='numerical')
    axs[1][1].set_ylabel('pitch rate [rad/s]')
    axs[1][1].grid()
    axs[0][1].legend(bbox_to_anchor=( -0.016 , 1 ),loc=3,fontsize='medium' )
    plt.show()
    return num_response, num_eigenval

def asymplot(eigenmotion,time,t0,t1,roll,rolldot,yawdot,force1,force2):
    #------------------------------------
    #Specific time interval for eigenmotion
    #------------------------------------
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
    time_eigen = np.arange(0,t1-dt,dt)

    #------------------------------------
    #Experimental Data
    #------------------------------------
    #Stationary parameters
    hp0 = hp[i_t0]
    alpha0 = AoA[i_t0]
    th0 = th[i_t0]
    V0 = Vt[i_t0]
    roll0 = roll[i_t0]
    rolldot0 = rolldot[i_t0]
    yawdot0 = yawdot[i_t0]
    #Forces
    force1 = force1[interval]
    force2 = force2[interval]
    force = np.array([force1, force2])

    #------------------------------------
    #Experimental response
    #------------------------------------
    roll = roll[interval]
    rolldot = rolldot[interval]
    yawdot = yawdot[interval]

    #------------------------------------
    #Numerical response and eigenvalues
    #------------------------------------
    #Flight dependent variables
    V0,m,rho,muc,mub,CL,CD,CX0,CZ0 = eigenmotion_parameters(time,t0)
    #Numerical solution
    num_solution = numres(V0,m,rho,muc,mub,CL,CD,CX0,CZ0)
    sys = num_solution[2]
    #Response
    num_response = ctrl.forced_response(sys, time_eigen, force)
    #Eigenvalues
    num_eigenval = num_solution[3]
    num_response_yaw = num_response[1][0]
    num_response_roll = num_response[1][1] + roll0
    num_response_rolldot = num_response[1][2] * (2*V0)/b + rolldot0
    numr_response_yawdot = num_response[1][3] * (2*V0)/b + yawdot0

    #------------------------------------
    # Plots
    #------------------------------------
    fig, axs = plt.subplots(2, 2, constrained_layout=True)
    fig.suptitle(eigenmotion,fontsize=16)
    # axs[0][0].plot(time_eigen, yaw, lw = 2, c = 'red', label='$\beta$ recorded')
    axs[0][0].plot(time_eigen, num_response_yaw, lw = 1, c = 'deepskyblue', label=' numerical')
    axs[0][0].set_ylabel('yaw angle [rad]')
    axs[0][0].grid()
    axs[0][1].plot(time_eigen, roll, lw=1, c='orange', label='recorded')
    axs[0][1].plot(time_eigen, num_response_roll, lw = 1, c = 'deepskyblue', label='numerical')
    axs[0][1].set_ylabel('$roll angle [rad]')
    axs[0][1].grid()
    axs[1][0].plot(time_eigen, rolldot, lw=1, c='orange', label=' recorded')
    axs[1][0].plot(time_eigen, num_response_rolldot, lw = 1, c = 'deepskyblue', label=' numerical')
    axs[1][0].set_ylabel('roll rate [m/s]')
    axs[1][0].grid()
    axs[1][1].plot(time_eigen, yawdot, lw=1, c='orange', label='recorded')
    axs[1][1].plot(time_eigen, numr_response_yawdot, lw = 1, c = 'deepskyblue', label='numerical')
    axs[1][1].set_ylabel('yaw rate [1/s]')
    axs[1][1].grid()
    axs[0][1].legend(bbox_to_anchor=( -0.016 , 1 ),loc=3,fontsize='medium' )
    plt.show()
    return num_response, num_eigenval

# #Short period motion
numres_spm, numeigen_spm = symplot('short period',time,spmt0,spm_t,AoA,Vt,Vc,th,q,de)
# #Phugoid motion
numres_phug, numeigen_spm = symplot('phugoid',time,phugt0,phug_t,AoA,Vt,Vc,th,q,de)
# #Aperiodic roll motion
numres_apr, numeigen_apr = asymplot('aperiodic roll',time,aprt0,apr_t,roll,rolldot,yawdot,da,dr)
# #Spiral motion
numres_spir, numeigen_spir = asymplot('spiral',time,spirt0,spir_t,roll,rolldot,yawdot,da,dr)
# #Dutch roll motion
# numres_dutch1, numeigen_dutch1 = asymplot('dutch roll 1',time,dutch1t0,dutch1_t,roll,rolldot,yawdot,da,dr)
numres_dutch2, numeigen_dutch2 = asymplot('dutch roll 2',time,dutch2t0,dutch2_t,roll,rolldot,yawdot,da,dr)
