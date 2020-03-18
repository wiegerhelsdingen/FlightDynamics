# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 17:47:47 2020

@author: tomvancranenburgh, reinier vos
"""

import numpy as np
import scipy.io as sc
import math
import matplotlib.pyplot as plt

#cl-cd series 1
mat1 = np.matrix([[7000,248,13.8,745,803,367,1.7],
                     [7000,221,11.8,641,687,400,2.5],
                     [6980,188,9.2,548,593,430,3.7],
                     [7000,162,7.8,456,502,470,5.5],
                     [7000,140,6.8,438,472,497,7.7],
                     [6980,120,5.8,456,507,515,10.6]])
#elevator trim curve
mat2 = np.matrix([[7120,162,5.5,-0.2,2.6,0,472,513,580,8.5],
                [7420,152,6.3,-0.7,2.6,-11,465,506,602,7.8],
                [7730,142,7.4,-1.2,	2.6,-29,462,500,626,6.5],
                [8020,133,8.6,-1.7,2.6,-37,456,494,644,5.8],
                [7390,172,4.6,0.1,2.6,17,471,512,671,8.2],
                [7130,182,4.1,0.5,2.6,36,476,516,681,9.2],
                [6780,192,3.5,0.7,2.6,58,483,523,700,10.5]])

#cg shift 
mat3 = np.matrix([[7240,165,5.2,0,2.8,0,471,511,735,8.0],
                    [7290,167,5.2,-0.7,2.8,-17,469,	511,773,8.2]])

# paramters
W_empty = 9165 #kg
blockfuel = 4100 #lbs
masspas = np.array([90,102,80,83,94,84,74,79,103]) #kg

#%% data unit conversion
bf_kg = blockfuel*0.453592  #kg
#mat1 
h_mat1 = mat1[:,0]*0.3048           # m
IAS_mat1 = mat1[:,1]*0.514444       # m/s
TAT_mat1 = mat1[:,2]+273.15         #temperature
FFL_mat1 = mat1[:,3]* (1/7936.64)   # kg/s
FFR_mat1 = mat1[:,4]* (1/7936.64)   #kg/s
WF_mat1 = mat1[:,5]* 0.453592       #kg
AOA_mat1 = mat1[:,6]                #degree  

#mat2
h_mat2 = mat2[:,0]*0.3048           # m
IAS_mat2 = mat2[:,1]*0.514444       # m/s
TAT_mat2 = mat2[:,2]+273.15         #KELVIN
DE_mat2 = mat2[:,3]                #degrees
DETR_mat2 = mat2[:,4]                #degrees
Fe_mat2 = mat2[:,5]                #N
FFL_mat2 = mat2[:,6]* (1/7936.64)   # kg/s
FFR_mat2 = mat2[:,7]* (1/7936.64)   #kg/s
WF_mat2 = mat2[:,8]* 0.453592       #kg
AOA_mat2 = mat2[:,9]                #degrees 

#mat3
h_mat3 = mat3[:,0]*0.3048           # m
IAS_mat3 = mat3[:,1]*0.514444       # m/s
TAT_mat3 = mat3[:,2] +273.15        #KELVIN
DE_mat3 = mat3[:,3]                #degrees
DETR_mat3 = mat3[:,4]                #degrees
Fe_mat3 = mat3[:,5]                #N
FFL_mat3 = mat3[:,6]* (1/7936.64)   # kg/s
FFR_mat3 = mat3[:,7]* (1/7936.64)   #kg/s
WF_mat3 = mat3[:,8]* 0.453592       #kg
AOA_mat3 = mat3[:,9]                #DEGREE  

#constants 
g0=9.81 #not needed for x cg calculation
S=30.00  #m^2
mf_s=0.048 #kg/sec, standard engine fuel flow per engine
c= 2.0569 #m, average chord
b= 15.911 #m, wing span

rho0=1.225   #kg/m^3 
p0=101325    #N/m^2 = Pa
T0=288.15    #K
R=287.05     #gas constant, [m^2 / K*sec^2]
lamb=-0.0065 #lambda for ISA pressure calculations
gamma=1.4    #ratio specific heats

'''
Ws=60500
Cm_0=
Cn_alpha
Cm_deltaf
deltaf
Cm_Tc=

'''

#%% Calibration

#convert indicated air speed to calibrated air speed, appendix A
def ias_cas(ias):
    cas = ias - 2
    return cas

#convert inidcated mach number to calibrated mach number WHERE DO WE USE THIS??
def im_cm(im): 
    if im>0.4 and im<0.705:
        cm=im-0.007
    else: 
        im=0
    return cm

#get pressure at altitude 
def pressure(hp):
    p=p0 * (1 + (lamb*hp)/T0) ** (- g0 / (lamb * hp))
    return p 

#calculated air density from perfect gas law
def density(p, T):
    rho = p / (R * T)
    return rho 

def mach(Vc, p):
    M=math.sqrt( (2/(gamma-1)) * ((1+ (p0/p)* ((1 + (gamma-1) / (2*gamma) * (rho0/p0 ) * Vc**2 )**(gamma/(gamma-1)) -1 ))**((gamma-1)/gamma)-1))
    return M

#corrected static air temperature for ram rise
def temperature(TAT, M):
    T = TAT / (1 + (gamma-1)/2 * M**2)    
    return T

#calculate true airspeed
def Vtrue(M, T):
    Vt=M*math.sqrt(gamma * R * T)
    return Vt
    
#equivalent airspeed
def Vequivalent(rho):
    Ve=Vt * math.sqrt(rho/rho0)
    return Ve

#drag curve
#Tp =   #DEFINE Tp HERE
def drag(Tp, AOA):
    D =  Tp * math.cos(AOA)   #CHECK IF ALPHA IS AOA
    return D



#%% Center of gravity 

#unit conversions to SI
inc_m=0.0254
m_inc=1/0.0254
pound_kg=0.453592 
ft_m=0.3048

#change everything so that unit change xcg only happens at the end

masspas = np.array([90,102,80,83,94,84,74,79,103]) #kg

#passenger weights [kg]
WP1=masspas[0]*m_inc
WP2=masspas[1]*m_inc
WL1=masspas[3]*m_inc
WR1=masspas[4]*m_inc
WL2=masspas[5]*m_inc
WR2=masspas[6]*m_inc
WL3=masspas[7]*m_inc
WR3=masspas[8]*m_inc
WCO=masspas[2]*m_inc

x0=131
x1=214
x2=251
x3=288
xC=170

#empty weight
ew_arm=291.65
ew_moment=W_empty*ew_arm
    
#fuel
fuel=blockfuel #-(ENTER FUNCTION FOR FUEL FLOW HERE)
fuel_moment=fuel*11705.5/100
    
#payload
payload_moment=(WP1+WP2)*x0+(WL1+WR1)*x1+(WL2+WR2)*x2+(WL3+WR3)*x3+WCO*xC
payload_weight=WP1+WP2+WL1+WR1+WL2+WR2+WL3+WR3+WCO
'''
def centergravity():
  
    #cg calculation
    xcg_nonsi=(ew_moment+fuel_moment+payload_moment)/(W_empty+fuel+payload_weight)
    xcg=xcg_nonsi*inc_m
    
    return xcg_nonsi, xcg
'''

#%%Measurement set 1 

def weight(WF):
    W=(W_empty+bf_kg+payload_weight- WF)*g0
    return W

# Lift and drag coefficient
Cl_mat1_list=[]
Cd_mat1_list=[]
alpha=[]

for i in range(0,len(IAS_mat1)):
    hp=float(h_mat1[i])
    ias=float(IAS_mat1[i])
    W=float(weight(WF_mat1[i]))
    cas=ias_cas(ias)
    Vc=cas
    p=pressure(hp)
    M=mach(Vc,p)
    TAT=float(TAT_mat1[i])
    T=temperature(TAT, M)
    rho=density(p, T)
    Vt=Vtrue(M, T)
    a1 = float(AOA_mat1[i])
    
    #D=drag(a1[i])                                #FIX THIS, ADD DEFINITION FOR D, Maybe take Tp out of brackets here
    
    alpha.append(a1)
    Cl = (2 * W) / (rho * Vt ** 2 * S)              # Lift coefficient [-]
    #Cd = (2 * D) / (rho * Vt **2 * S)
    Cl_mat1_list.append(Cl)
    #Cd_mat1_list.append(Cd)
   
#%% CL-alpha curve

plt.scatter(alpha,Cl_mat1_list)
plt.xlabel('angle of attack [degrees]')
plt.ylabel('lift coefficient [-]')
plt.grid()

z=np.polyfit(alpha,Cl_mat1_list,1)
t=np.poly1d(z)
plt.plot(alpha,t(alpha),"r-")
plt.title('Lift coefficient versus angle of attack')
print("y=%.6fx+%.6f"%(z[0],z[1])) 

plt.show()

'''
#%% Cd-alpha curve
plt.plot(alpha,Cd_mat1_list)
plt.xlabel('angle of attack [degrees]')
plt.ylabel('drag coefficient [-]')
plt.grid()

zz=np.polyfit(alpha,Cd_mat1_list,1)
tt=np.poly1d(zz)
plt.plot(alpha,tt(alpha),"r-")
plt.title('Drag coefficient versus angle of attack')

plt.show()



#%% Cl-Cd curve, Cl^2-Cd plot

plt.plot(Cl_mat1_list,Cd_mat1_list)
plt.xlabel('lift coefficient [-]')
plt.ylabel('drag coefficient [-]')
plt.grid()

zzz=np.polyfit(alpha,Cd_mat1_list,1)
ttt=np.poly1d(zzz)
plt.plot(alpha,ttt(alpha),"r-")
plt.title('Drag coefficient versus lift coefficient')

plt.show()


plt.plot(Cl_mat1_list**2,Cd_mat1_list)
plt.xlabel('lift coefficient [-]')
plt.ylabel('drag coefficient [-]')
plt.grid()

zzz=np.polyfit(alpha,Cd_mat1_list,1)
ttt=np.poly1d(zzz)
plt.plot(alpha,ttt(alpha),"r-")
plt.title('Drag coefficient versus lift coefficient')

plt.show()



#%% Oswald efficiency factor

#CD = CD0 + (CLa * alpha0) ** 2 / (math.pi * A * e) # Drag coefficient [-]
#dCldCd=
#e=    


'''











'''


#%% Elevator trim curve 

#reduction airspeed
W= 
Ve= 
def reductionairspeed(Tm, W, hp, ias):
    Ve=ias_cas(ias)
    Vr=Ve*sqrt(Ws/W) * #HERE ALSO CONVERT KTS TO CORRECT UNITS

#reduction non-standard engine thrust



#calculate everything for measurements 1
for i in range(len(IAS_mat1)):
    


#thrust coefficients
Tps =    #standard engine thrust, output of thrust.exe 
D_e =     #engine inlet diameter

Tcs=Tps / ( 0.5 * rho * V_e_reduced**2 * 2 * D_e**2 )


#elevator trim
# taken from model parameters
Cl_alpha=5.084    #will be compared
Cm0=0.0297
Cm_alpha=-0.5626
Cm_Tc=-0.0064

CN_alpha=Cl_alpha #take this from the graph



delta_e=- 1/(C_m_delta) * (CM0 + (Cm_alpha/CN_alpha) * (W/(0.5*rho*V_e_reduced**2*S) ) + Cm_deltaf*deltaf + Cm_Tc* Tc_s )



#%% Elevator effectiveness
xcg2= xc( )   #make cg a definition with inputs

delta_xcg=xcg2-xcg1 
delta_delta_e= delta_e2-delta_e1
Cm_delta= - (1/ delta_delta_e) * (W / (0.5*rho*Vt**2*S)) *(delta_xcg / c)


'''





