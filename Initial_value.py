import numpy as np
import math
import control as ctrl
import matplotlib.pyplot as plt
import control as cntrl
from parameters import *
from Numerical_Simulation import *

V0,m,rho,muc,mub,CL,CD,CX0,CZ0 = eigenmotion_parameters(time,spmt0)
time = np.arange(0,100,0.1)
num_solution = numres(V0,m,rho,muc,mub,CL,CD,CX0,CZ0)
sys_s = num_solution[0]

num_response = ctrl.impulse_response(sys_s,time)
num_response = num_response[1]

plt.plot(time, num_response[0], lw = 1, c = 'blue')
plt.ylabel('True airspeed [m/s]')
plt.grid()
plt.show()
plt.plot(time, num_response[1], lw = 1, c = 'blue')
plt.ylabel('pitch angle [rad]')
plt.grid()
plt.show()
plt.plot(time, num_response[2], lw = 1, c = 'blue')
plt.ylabel('AoA [rad]')
plt.grid()
plt.show()
plt.plot(time[:200], num_response[3][:200], lw = 1, c = 'blue')
plt.ylabel('pitch rate [rad/s]')
plt.grid()
plt.show()

time = np.arange(0,60,0.1)
sys_a = num_solution[2]
num_response = ctrl.impulse_response(sys_a,time)
num_response = num_response[1]

plt.plot(time, num_response[0], lw = 1, c = 'blue')
plt.ylabel('yaw angle [rad]')
plt.grid()
plt.show()
plt.plot(time, num_response[1], lw = 1, c = 'blue')
plt.ylabel('roll angle  [rad]')
plt.grid()
plt.show()
plt.plot(time, num_response[2], lw = 1, c = 'blue')
plt.ylabel('yaw rate [rad/s]')
plt.grid()
plt.show()
plt.plot(time, num_response[3], lw = 1, c = 'blue')
plt.ylabel('roll rate [rad/s]')
plt.grid()
plt.show()
