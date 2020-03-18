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
W_empty = 9165*0.453592 #kg
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
WF_mat1_lbs = mat1[:,5]              #lbs needed for cg 
AOA_mat1 = mat1[:,6]                #degree  

#mat2
h_mat2 = mat2[:,0]*0.3048           # m
IAS_mat2 = mat2[:,1]*0.514444       # m/s
TAT_mat2 = mat2[:,2]+273.15         #KELVIN
DE_mat2 = mat2[:,3]                 #degrees
DETR_mat2 = mat2[:,4]               #degrees
Fe_mat2 = mat2[:,5]                 #N
FFL_mat2 = mat2[:,6]* (1/7936.64)   # kg/s
FFR_mat2 = mat2[:,7]* (1/7936.64)   #kg/s
WF_mat2 = mat2[:,8]* 0.453592       #kg
WF_mat2_lbs = mat2[:,8]              #lbs needed for cg 
AOA_mat2 = mat2[:,9]                #degrees 

#mat3
h_mat3 = mat3[:,0]*0.3048           # m
IAS_mat3 = mat3[:,1]*0.514444       # m/s
TAT_mat3 = mat3[:,2] +273.15        #KELVIN
DE_mat3 = mat3[:,3]                 #degrees
DETR_mat3 = mat3[:,4]               #degrees
Fe_mat3 = mat3[:,5]                 #N
FFL_mat3 = mat3[:,6]* (1/7936.64)   # kg/s
FFR_mat3 = mat3[:,7]* (1/7936.64)   #kg/s
WF_mat3 = mat3[:,8]* 0.453592       #kg
WF_mat3_lbs = mat3[:,8]             #lbs needed for cg 
AOA_mat3 = mat3[:,9]                #DEGREE  

#constants 
g0=9.81 #not needed for x cg calculation
S=30.00  #m^2
mf_s=0.048 #kg/sec, standard engine fuel flow per engine
c= 2.0569 #m, average chord
b= 15.911 #m, wing span
A= b**2 / S

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

#%% OBTAIN FUEL MOMENT (FM) POLYNOMIAL FOR CG CALCULATIONS (NOTE: entire section is based on table E2 and is in lbs and inches)
FM_MOMENTS = [298.16, 591.18,879.08,1165.42,1448.40,1732.53,2014.80,2298.84,2581.92,2866.30,3150.18,3434.52,3718.52,4003.23,4287.76,4572.24,
               4856.56,5141.16,5425.64,5709.90,5994.04,6278.47,6562.82,6846.96,7131.00,7415.33,7699.60,7984.34,8269.06,8554.05,8839.04,9124.80,
               9410.62,9696.97,9983.40,10270.08,10556.84,10843.87,11131.00,11418.20,11705.50,11993.31,12281.18,12569.04,12856.86,13144.73,
               13432.48,13720.56,14008.46,14320.34] # ALL FUEL MOMENTS GIVEN
FM_MOMENTS = [i * 100 for i in FM_MOMENTS] # ALL FUEL MOMENTS * 100
FM_WEIGHTS = [100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,
              2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5008,] #ALL FUEL MASSES
### OBTAINING POLYNOMIAL -> COPY THIS TO CG FUEL CALC
"""
plt.scatter(FM_WEIGHTS,FM_MOMENTS)
plt.xlabel('FM_WEIGHTS')
plt.ylabel('FM_MOMENTS')
plt.grid()
FMPOLY =np.polyfit(FM_WEIGHTS,FM_MOMENTS,1)
FMPOLY_t=np.poly1d(FMPOLY)
plt.plot(FM_WEIGHTS,FMPOLY_t(FM_WEIGHTS),"r-")
plt.title('FM_MOMENTS VS FM_WEIGHTS')
print("FM_y=%.6fx+%.6f"%(FMPOLY[0],FMPOLY[1])) 
##POLYNOMIAL SEEMS LINEAR ! CHECK THIS, OTHERWISE CHANGE ORDER
"""
#%% Center of gravity 
#unit conversions to SI
inc_m=0.0254
m_inc=1/0.0254
pound_kg=0.453592 
ft_m=0.3048
lbsin_kgm = 1/86.796166214519 # lbs inches to kg m

#change everything so that unit change xcg only happens at the end
masspas = np.array([90,102,80,83,94,84,74,79,103]) #kg
#passenger weights [kg]
WP1=masspas[0]
WP2=masspas[1]
WL1=masspas[3]
WR1=masspas[4]
WL2=masspas[5]
WR2=masspas[6]
WL3=masspas[7]
WR3=masspas[8]
WCO=masspas[2]
# passenger distances from tip [m]
x0=131*inc_m
x1=214*inc_m
x2=251*inc_m
x3=288*inc_m
xC=170*inc_m
# baggage distances (only use if necessary) [m] from tip
x_nose_b = 74*inc_m
x_aft1_b = 321*inc_m
x_aft1_b = 338*inc_m

W_payload =WP1+WP2+WL1+WR1+WL2+WR2+WL3+WR3+WCO  #[kg]

def centerofgravity(x3R, WF):
    #IMPORTANT READ BELOW:
    #ensure that all calculations are done for kg and meters not inches and lbs!!!
    #WF has to be an entry of the WF_MAT_LBS matrices,  not the [kg] ones!
    #empty weight
    ew_arm=291.65*inc_m # [m] mass balance report
    ew_moment=W_empty*ew_arm  # [kg m]
    xcg_BEM = 291.65*inc_m
    #payload 
    payload_moment=(WP1+WP2)*x0+(WL1+WR1)*x1+(WL2+WR2)*x2+(WL3)*x3+(WR3)*x3R+WCO*xC # [kg m]
    #ZFM CG
    xcg_ZFM = (payload_moment + ew_moment)/(W_payload + W_empty) # [m ]
    #fuel
    W_fuel_lbs= float(blockfuel - WF)            # [lbs]!!!
    W_fuel= float(blockfuel - WF)*pound_kg      #[kg]
    fuel_moment_lbs = 285.2562*W_fuel_lbs +989.5738 # polynomial obtained earlier (has to be changed manually)
    fuel_moment= fuel_moment_lbs*lbsin_kgm   #from lbs in to kg m , current unit [kg  m]
    #RAMP MASS(RM) CG
    xcg_RM = ((payload_moment + ew_moment+ fuel_moment)/(W_payload + W_empty+ W_fuel))
    return xcg_BEM, xcg_ZFM, xcg_RM
# plot of xcg rm change due to wf
xcgRM = []
WF_RM = []
for i in range(0,len(WF_mat1_lbs)):
    
    WF = float(WF_mat1_lbs[i])
    xcgs = centerofgravity(x3, WF)
    xcgRM1 = xcgs[2] 
    xcgRM.append(xcgRM1)
    WF_RM.append(WF)

plt.figure
plt.scatter(xcgRM, WF_RM)
plt.ylabel('WF [kg]')
plt.xlabel('xcg_RM from tip [m]')
plt.title('XCG_RM vs WF')
plt.grid()
#%% Measurement set 1 

def weight(WF):
    W=(W_empty+bf_kg+W_payload- WF )* g0
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
    Cl = (2 * W) / (rho * (Vt ** 2) * S)              # Lift coefficient [-]
    #Cd = (2 * D) / (rho * (Vt **2) * S)
    Cl_mat1_list.append(Cl)
    #Cd_mat1_list.append(Cd)
   
#%% CL-alpha curve
#change alpha0 to root location when you know trendline, following includes alpha_CL=0 point
alphacl0 = -0.87444612
rootcl0 = 0
alpha.insert(0,alpha0)
Cl_mat1_list.insert(0,rootcl)
plt.scatter(alpha,Cl_mat1_list)
plt.xlabel('angle of attack [radians]')
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
plt.scatter(alpha,Cd_mat1_list)
plt.xlabel('angle of attack [degrees]')
plt.ylabel('drag coefficient [-]')
plt.grid()

zz=np.polyfit(alpha,Cd_mat1_list,1)
tt=np.poly1d(zz)
plt.plot(alpha,tt(alpha),"r-")

plt.show()



#%% Cl-Cd curve, Cl^2-Cd plot

plt.scatter(Cl_mat1_list,Cd_mat1_list)
plt.xlabel('lift coefficient [-]')
plt.ylabel('drag coefficient [-]')
plt.grid()

zzz=np.polyfit(alpha,Cd_mat1_list,1)
ttt=np.poly1d(zzz)
plt.plot(alpha,ttt(alpha),"r-")

plt.show()

Cl2_mat1_list=[]
for i in range(len(Cl_mat1_list)):
    Cl2[i]=(Cl_mat1_list[i])**2
    Cl2_mat1_list.append(Cl2)
    i=i+1

plt.scatter(Cl2_mat1_list,Cd_mat1_list)
plt.xlabel('lift coefficient [-]')
plt.ylabel('drag coefficient [-]')
plt.grid()

zzz=np.polyfit(alpha,Cd_mat1_list,1)
ttt=np.poly1d(zzz)
plt.plot(alpha,ttt(alpha),"r-")

plt.show()



#%% Oswald efficiency factor

#CD = CD0 + (CLa * alpha0) ** 2 / (math.pi * A * e) # Drag coefficient [-]
dCddCl2 = 
e = 1 / (math.pi * A dCLdCd2 )
print('oswald efficiency factor e =', e)

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





