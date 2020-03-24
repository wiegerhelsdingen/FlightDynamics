import numpy as np
import math
import control as ctrl
import matplotlib.pyplot as plt
import control as cntrl
from parameters import *
from Numerical_Simulation import *

V0,m,rho,muc,mub,CL,CD,CX0,CZ0 = eigenmotion_parameters(time,spmt0)
num_solution = numres(V0,m,rho,muc,mub,CL,CD,CX0,CZ0)
sys_s = num_solution[0]
sys_a = num_solution[2]

def initial_value_s(title,sys,time,X0):
    num_response = ctrl.initial_response(sys,time,X0)
    num_response = num_response[1]

    fig, axs = plt.subplots(4, 1, constrained_layout=True)
    fig.suptitle(title,fontsize=16)
    axs[0].plot(time, num_response[0], lw = 1, c = 'blue')
    axs[0].set_ylabel('True airspeed [m/s]')
    axs[0].grid()
    axs[1].plot(time, num_response[1], lw = 1, c = 'orange')
    axs[1].set_ylabel('pitch angle [rad]')
    axs[1].grid()
    axs[2].plot(time, num_response[2], lw = 1, c = 'red')
    axs[2].set_ylabel('AoA [rad]')
    axs[2].grid()
    axs[3].plot(time[:200], num_response[3][:200], lw = 1, c = 'green')
    axs[3].set_ylabel('pitch rate [rad/s]')
    axs[3].grid()
    plt.show()
    return

def initial_value_a(title,sys,time,X0):
    sys_a = num_solution[2]
    num_response = ctrl.initial_response(sys,time,X0)
    num_response = num_response[1]

    fig, axs = plt.subplots(4, 1, constrained_layout=True)
    fig.suptitle(title,fontsize=16)
    axs[0].plot(time, num_response[0], lw = 1, c = 'blue')
    axs[0].set_ylabel('yaw angle [rad]')
    axs[0].grid()
    axs[1].plot(time, num_response[1], lw = 1, c = 'orange')
    axs[1].set_ylabel('roll angle  [rad]')
    axs[1].grid()
    axs[2].plot(time, num_response[2], lw = 1, c = 'red')
    axs[2].set_ylabel('yaw rate [rad/s]')
    axs[2].grid()
    axs[3].plot(time, num_response[3], lw = 1, c = 'green')
    axs[3].set_ylabel('roll rate [rad/s]')
    axs[3].grid()
    plt.show()

time_s = np.arange(0,80,0.1)

# num_response_V = initial_value_s('Response to initial V', sys_s,time_s,[1,0,0,0])
# num_response_th = initial_value_s('Response to initial pitch angle',sys_s,time_s,[0,1,0,0])
# num_response_alpha = initial_value_s('Response to initial AoA',sys_s,time_s,[0,0,1,0])
# num_response_q = initial_value_s('Response to initial pitch rate',sys_s,time_s,[0,0,0,1])

time_a = np.arange(0,15,0.1)

num_response_b = initial_value_a('Response to initial yaw angle',sys_a,time_a,[1,0,0,0])
num_response_phi = initial_value_a('Response to initial roll angle',sys_a,time_a,[0,1,0,0])
# num_response_p = initial_value_a('Response to initial yaw rate',sys_a,time_a,[0,0,1,0])
# num_response_r = initial_value_a('Response to initial roll rate',sys_a,time_a,[0,0,0,1])
