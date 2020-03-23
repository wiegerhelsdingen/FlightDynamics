
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

#thrust computations from thrust.exe for first measurement series, in N
thrust = np.matrix([[3349.31, 3756.94],
                    [2821.05, 3133.59],
                    [2384.87, 2699.04],
                    [1876.38, 2208],
                    [1858.28, 2111.42],
                    [2104.39, 2495.74]])

#thrust computations from thrust.exe for second measurement series, in N
# non standardized
thrust2 = np.matrix([[1990.94, 2286.60],
                     [2005.42, 2306.69],
                     [2053.98, 2337.53],
                     [2073.09, 2358.38],
                     [1945.13, 2237.69],
                     [1917.62, 2198.03],
                     [1898.25, 2174.97]])

#thrust computations from thrust.exe for second measurement series, in N
# standardized -> left and right thrust must match!
thrust2_s =np.matrix([[1346.94, 1346.94],
                      [1401.15, 1401.15],
                      [1460.50, 1460.50],
                      [1514.92, 1514.92],
                      [1318.04, 1318.04],
                      [1266.70, 1266.70],
                      [1211.88, 1211.88]])

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
TAT_mat2 = mat2[:,9]+273.15         #KELVIN
DE_mat2 = mat2[:,3]                 #degrees
DETR_mat2 = mat2[:,4]               #degrees
Fe_mat2 = mat2[:,5]                 #N
FFL_mat2 = mat2[:,6]* (1/7936.64)   # kg/s
FFR_mat2 = mat2[:,7]* (1/7936.64)   #kg/s
WF_mat2 = mat2[:,8]* 0.453592       #kg
WF_mat2_lbs = mat2[:,8]              #lbs needed for cg 
AOA_mat2 = mat2[:,2]                #degrees 

#mat3
h_mat3 = mat3[:,0]*0.3048           # m
IAS_mat3 = mat3[:,1]*0.514444       # m/s
TAT_mat3 = mat3[:,9] +273.15        #KELVIN
DE_mat3 = mat3[:,3]                 #degrees
DETR_mat3 = mat3[:,4]               #degrees
Fe_mat3 = mat3[:,5]                 #N
FFL_mat3 = mat3[:,6]* (1/7936.64)   # kg/s
FFR_mat3 = mat3[:,7]* (1/7936.64)   #kg/s
WF_mat3 = mat3[:,8]* 0.453592       #kg
WF_mat3_lbs = mat3[:,8]             #lbs needed for cg 
AOA_mat3 = mat3[:,2]                #DEGREE  

#thrust
Tleft = thrust[:,0]                  # N
Tright = thrust[:,1]                 # N
Tleft2 = thrust2[:,0]                # N
Tright2 = thrust2[:,1]               # N
Tleft2_s = thrust2_s[:,0]            # N
Tright2_s = thrust2_s[:,1]           # N
#constants 
g0=9.81 #not needed for x cg calculation
S=30.00  #m^2
mf_s=0.048 #kg/sec, standard engine fuel flow per engine
c= 2.0569 #m, average chord
b= 15.911 #m, wing span
A= b**2 / S
u0 = 1.802*(10**-5)   # dynamic viscosity [kg/m-s]
rho0=1.225   #kg/m^3 
p0=101325    #N/m^2 = Pa
T0=288.15    #K
R=287.05     #gas constant, [m^2 / K*sec^2]
lamb=-0.0065 #lambda for ISA pressure calculations
gamma=1.4    #ratio specific heats
R_engine = 0.69/2    #radius engine [m] average of range given in report
A_engine = math.pi*(R_engine**2)     #area engine [m^2]
Ws = 60500      #N, needed for  reduced airspeed & standardization
#%% Calibration

#convert indicated air speed to calibrated air speed, appendix A
def ias_cas(ias):
    cas = ias - 2*0.514444
    return cas

#convert inidcated mach number to calibrated mach number WHERE DO WE USE THIS??
def im_cm(im): 
    if im>0.4 and im<0.705:
        cm=im-0.007
    else: 
        cm=0
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
def Vequivalent(rho, Vt):
    Ve=Vt * math.sqrt(rho/rho0)
    return Ve

#drag curve
def drag(Tp, AOA):
    D =  Tp * math.cos(math.radians(AOA))   #CHECK IF ALPHA IS AOA -> must be radians!
    return D
# reduced velocity
def Vreduction(W, Ve):
    Vred = Ve * math.sqrt(Ws / W)  #Reduced airspeed
    return (Vred)
def Dviscosity(T):
    u = u0*((T/T0)**(3/2))*((T0+S)/(T+S))
    return (u)
def reynolds(rho,V,u):
    Re = rho*V*c/u
    return Re
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
### plot of xcg rm change due to wf

xcgRM = []
WF_RM = []
for i in range(0,len(WF_mat1_lbs)):
    
    WF = float(WF_mat1_lbs[i])
    xcgs = centerofgravity(x3, WF)
    xcgRM1 = xcgs[2] 
    xcgRM.append(xcgRM1)
    WF_RM.append(WF)
"""
plt.figure(9)
plt.scatter(xcgRM, WF_RM)
plt.ylabel('WF [kg]')
plt.xlabel('xcg_RM from tip [m]')
plt.title('XCG_RM vs WF')
plt.grid()

### check whether ZFM xcg matches curve for fuel in wing
xcg1 = []
FUELinwing1 = []
for i in range(0,len(FM_WEIGHTS)):
    WF = float(FM_WEIGHTS[i])
    FUELinwing = float(blockfuel - WF)*pound_kg
    xcgs = centerofgravity(x3, WF)
    xcgRM1 = xcgs[2] 
    xcg1.append(xcgRM1)
    FUELinwing1.append(FUELinwing)
    
plt.figure(10)
plt.scatter(xcg1, FUELinwing1)
plt.scatter(xcgs[1], 0)
plt.ylabel('FUELinwing [kg]')
plt.xlabel('xcg_RM from tip [m]')
plt.title('XCG_RM vs FUELinwing')
plt.grid()
"""
# conclusion: slight discrepancy due to table E2: wing can never be totally empty!

#%% Measurement set 1 

def weight(WF):
    W=(W_empty+bf_kg+W_payload- WF )* g0
    return W

# Lift and drag coefficient
Cl_mat1_list=[]
Cd_mat1_list=[]
alpha=[]
Tp=[]

for i in range(0,len(Tleft)):
    Tprop=float(Tleft[i])+float(Tright[i])
    Tp.append(Tprop)
    
machlist1 = []
reynoldslist1 = []
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
    u1 = Dviscosity(T)
    Re1 = reynolds(rho,Vt,u1)
    machlist1.append(M)
    reynoldslist1.append(Re1)
    a1 = float(AOA_mat1[i])
    Tprop = float(Tp[i])
    D=drag(Tprop, a1)             
    
    alpha.append(a1)
    Cl = (2 * W) / (rho * (Vt ** 2) * S)              # Lift coefficient [-]
    Cd = (2 * D) / (rho * (Vt **2) * S)
    Cl_mat1_list.append(Cl)
    Cd_mat1_list.append(Cd)

#Mach & Reynolds number range series 1
print("Mach range series 1 =", min(machlist1),"to", max(machlist1))
print("Reynolds range series 1= ", min(reynoldslist1),"to",max(reynoldslist1))

#%% CL-alpha curve
#FIRST ORDER
#ROOT INSERTED
z=np.polyfit(alpha,Cl_mat1_list,1)
t=np.poly1d(z)
alphacl0 = -z[1]/z[0]
rootcl0 = 0
alphalist = np.linspace(alphacl0,float(max(AOA_mat1)),100)
CLA_CL = Cl_mat1_list.copy()
CLA_ALPHA = alpha.copy()
CLA_ALPHA.insert(0,alphacl0)
CLA_CL.insert(0,rootcl0)

plt.figure(1)
plt.scatter(CLA_ALPHA,CLA_CL, label='Measured data')
plt.xlabel('angle of attack [radians]')
plt.ylabel('lift coefficient [-]')
plt.grid()
plt.plot(CLA_ALPHA,t(CLA_ALPHA),"r-", label='Regression line')
print("y=%.6fx+%.6f"%(z[0],z[1])) 
plt.legend(loc ="upper left")
plt.show()
CLA_GRAD = z[0]

#%% Cd-alpha curve
#SECOND ORDER
plt.figure(2)
plt.scatter(alpha,Cd_mat1_list, label='Measured data')
plt.xlabel('angle of attack [degrees]')
plt.ylabel('drag coefficient [-]')
plt.title('DRAG - ALPHA')
plt.grid()
zz=np.polyfit(alpha,Cd_mat1_list,2)
tt=np.poly1d(zz)    
plt.plot(alphalist,tt(alphalist),"r-", label='Regression line')
plt.legend(loc ="upper left")
plt.show()

#%% Cl-Cd curve
#(HIGHER ORDER) 
#THE ROOT (CD0) IS FOUND BY FITTING A 4TH ORDER POLYNOMIAL 
#THEN THE ROOT IS ADDED AND A 3RD ORDER POLY IS FITTED ON THE RESULT
ROOTPOLY=np.polyfit(Cd_mat1_list[0:5],Cl_mat1_list[0:5],4)
ROOTt2=np.poly1d(ROOTPOLY)
t2root = np.real((np.roots(ROOTt2))[3])
CLCD_CL = Cl_mat1_list.copy()
CLCD_CD = Cd_mat1_list.copy()
CLCD_CL.insert(0,0)
CLCD_CD.insert(0,t2root)

plt.figure(3)
plt.scatter(CLCD_CD, CLCD_CL, label='Measured data')
plt.xlabel('drag coefficient [-]')
plt.ylabel('lift coefficient [-]')
plt.title('LIFT - DRAG')
plt.grid()
CDCL=np.polyfit(CLCD_CD[0:3],CLCD_CL[0:3],3)
t2=np.poly1d(CDCL)
cdlist = np.linspace(t2root,CLCD_CD[2],100)
CDCL_rest = np.polyfit(CLCD_CD[2:7],CLCD_CL[2:7],3)
t2_rest =np.poly1d(CDCL_rest)
cdlist_rest = np.linspace(CLCD_CD[2],CLCD_CD[6],100)
plt.plot(cdlist_rest,t2_rest(cdlist_rest),"r-", label='Regression line')
plt.plot(cdlist,t2(cdlist),"r-")
plt.legend(loc ="upper left")
plt.show()

print("CD0=", t2root)

#%% Cl^2-Cd plot
#FIRST ORDER
Cl2_mat1_list=[]

for i in range(0, len(Cl_mat1_list)):
    Cl2=(float(Cl_mat1_list[i]))**2
    Cl2_mat1_list.append(Cl2)

#ROOT INSERTED
CL2CD=np.polyfit(Cl2_mat1_list, Cd_mat1_list,1)
t3=np.poly1d(CL2CD)
CL2CD_CL0 = -t3[0]/t3[1]
rootCL2CD0 = 0
CL2CD_CL2 = Cl2_mat1_list.copy()
CL2CD_CD = Cd_mat1_list.copy()
CL2CD_CL2.insert(0,CL2CD_CL0)
CL2CD_CD.insert(0,rootCL2CD0)

plt.figure(4)
plt.scatter(CL2CD_CL2, CL2CD_CD, label='Measured data')
plt.xlabel('lift coefficient^2 [-]')
plt.ylabel('drag coefficient [-]')
plt.title('LIFT^2 - DRAG')
plt.grid()
plt.plot(CL2CD_CL2,t3(CL2CD_CL2),"r-", label='Regression line')
plt.legend(loc ="upper left")
plt.show()

print('CL^2/CD line gradient =',t3[1])
CL2CDGRAD = t3[1]

#%% Oswald efficiency factor
#CD = CD0 + (CLa * alpha0) ** 2 / (math.pi * A * e) # Drag coefficient [-]

e = 1 / (math.pi * A * CL2CDGRAD)
print('oswald efficiency factor e =', e)

CD0 = t2root 
CLVAL = np.linspace(0,1,100)
CDVAL = CD0 + (CLVAL**2)/(math.pi*e*A)
plt.figure(3)
plt.plot(CDVAL,CLVAL,'g', label='CD standard formula')
plt.legend(loc ="upper left")
plt.show()
#%% Code for center gravity shit

x3R1=x3          #m
x3R2=x0
WF1=WF_mat3[0]   #kg
WF2=WF_mat3[1]   #kg

sta_ref=261.45*inc_m     #to define center of gravity with respect to imaginary foward end MAC
xcg1=centerofgravity(x3R1, WF1)[2]-sta_ref
xcg2=centerofgravity(x3R2, WF2)[2]-sta_ref

delta_xcg=xcg2-xcg1   #this turns out negative, define xcg from MAC for it to become positive??  
print('center gravity 1', xcg1, 'm', xcg1*m_inc, 'inch')
print('center gravity 2', xcg2, 'm', xcg2*m_inc, 'inch')
print('change center gravity', delta_xcg, 'm')   


#%% Elevator effectiveness

def elevatoreffectiveness():
    delta_e1 = DE_mat3[0]
    delta_e2 = DE_mat3[1]
    delta_e = delta_e2-delta_e1
    delta_e_rad = np.radians(delta_e)
    W3=float(weight(WF_mat3[0]))    # WHICH MEASUREMENT OF THE TWO FOR USED FUEL TO USE HERE?
    
    hp3=float(h_mat3[0])            # TAKE AVERAGE OF BOTH VALUES OR ASSUME ONE?
    ias3=float(IAS_mat3[0])
    cas3=ias_cas(ias3)
    Vc3=cas3
    p3=pressure(hp3)
    M3=mach(Vc3,p3)
    TAT3=float(TAT_mat3[0])
    T3=temperature(TAT3, M3)
    rho3=density(p3, T3)
    Vt3=Vtrue(M3, T3)
    
    Cm_delta= - (1/ delta_e_rad) * (W3 / (0.5*rho3*(Vt3**2)*S)) *(delta_xcg / c)
    
    return Cm_delta, delta_e_rad

elevator_effectiveness=elevatoreffectiveness()[0]
delta=elevatoreffectiveness()[1]
print(delta)
print('elevator effectivenenss [-] = ', elevator_effectiveness)

#%% cm alpha

Cmdelta = elevator_effectiveness
alpha2 = []
De2 = []
for i in range(0,len(AOA_mat2)):
    alpha2.append(float(AOA_mat2[i]))
    De2.append(float(DE_mat2[i]))

DE_A = np.polyfit(alpha2,De2, 1)
u = np.poly1d(DE_A)
de_da = u[1]
"""
plt.figure(6)
plt.scatter(alpha2,De2)
plt.show()
"""
Cmalpha = - Cmdelta * de_da
print('Cm alpha = ', Cmalpha)
Cm0 = 0.0297 # obtained from appendix B

#%% Cm elevator equal
# C_N = C_L (Assummed to be equal)
# 
C_N = CLA_GRAD
CM_TC = -0.0064  #obtained from appendix B
Tp2=[]
for i in range(0,len(Tleft2)):
    Tprop2=float(Tleft2[i])+float(Tright2[i])
    Tp2.append(Tprop2)
    
Tp2_s=[]
for i in range(0,len(Tleft2_s)):
    Tprop2_s=float(Tleft2_s[i])+float(Tright2_s[i])
    Tp2_s.append(Tprop2_s)
machlist2 = []
reynoldslist2 = []


delta_redlist = []
Vred2list = []
Fstlist = []
for i in range(0,len(IAS_mat2)):
    hp2=float(h_mat2[i])
    ias2=float(IAS_mat2[i])
    W2=float(weight(WF_mat2[i]))
    cas2=ias_cas(ias2)
    Vc2=cas2
    p2=pressure(hp2)
    M2=mach(Vc2,p2)
    machlist2.append(M2)
    TAT2=float(TAT_mat2[i])
    T2=temperature(TAT2, M2)
    rho2=density(p2, T2)
    Vt2=Vtrue(M2, T2)
    Veq2 = Vequivalent(rho2, Vt2)
    Vred2 = Vreduction(W2, Veq2)
    Tp22 = float(Tp2[i])
    Tp22_s = float(Tp2_s[i])
    T_c = (2*Tp22)/(rho0*(R_engine**2)*(Veq2**2))
    T_cs = (2*Tp22_s)/(rho2*(R_engine**2)*(Vred2**2))
    delta_eq_meas = float(DE_mat2[i])
    delta_red = delta_eq_meas - (1/Cmdelta)*CM_TC*(T_cs-T_c)
    delta_red = math.radians(delta_red)
    delta_redlist.append(float(delta_red))
    Vred2list.append(float(Vred2))
    Fe = float(Fe_mat2[i])
    Fst = Fe * (Ws/W2)
    Fstlist.append(float(Fst))
    u2 = Dviscosity(T2)
    Re2 = reynolds(rho2,Vt2,u2)
    reynoldslist2.append(Re2)
"""
#DONT USE THIS,  MUST ASK TAS
delta_redlist = []
Vred2list = []
Fstlist = []
for i in range(0,len(IAS_mat2)):
    hp2=float(h_mat2[i])
    ias2=float(IAS_mat2[i])
    W2=float(weight(WF_mat2[i]))
    cas2=ias_cas(ias2)
    Vc2=cas2
    p2=pressure(hp2)
    M2=mach(Vc2,p2)
    TAT2=float(TAT_mat2[i])
    T2=temperature(TAT2, M2)
    rho2=density(p2, T2)
    Vt2=Vtrue(M2, T2)
    Veq2 = Vequivalent(rho2, Vt2)
    Vred2 = Vreduction(W2, Veq2)
    Tp22 = float(Tp2[i])
    Tp22_s = float(Tp2_s[i])
    T_c = (2*Tp22)/(rho0*A_engine*(Veq2**2))
    T_cs = (2*Tp22_s)/(rho2*A_engine*(Vred2**2))
    delta_eq_meas = -(1/Cmdelta)*(Cm0 + (Cmalpha/C_N)*(2*W2/(rho2*(Vred2**2)*S))+CM_TC*T_c)
    delta_red = delta_eq_meas - (1/Cmdelta)*CM_TC*(T_cs-T_c)
    delta_red = math.radians(delta_red)
    delta_redlist.append(float(delta_red))
    Vred2list.append(float(Vred2))
    Fe = float(Fe_mat2[i])
    Fst = Fe * (Ws/W2)
    Fstlist.append(float(Fst))
#DONT USE THIS,  MUST ASK TAS
"""
deltaveq = np.polyfit(Vred2list,delta_redlist,2)
dveq = np.poly1d(deltaveq)
FeVeq =np.polyfit(Vred2list,Fstlist,2)
Fveq = np.poly1d(FeVeq)
vlist = np.linspace(min(Vred2list),max(Vred2list),100)
plt.figure(7)
plt.scatter(Vred2list, delta_redlist, label = 'Measured data')
plt.grid()
plt.xlabel('Reduced velocity [m/s]')
plt.ylabel('Reduced elevator deflection [m]')
plt.plot(vlist,dveq(vlist),"r-",label = 'Regression line')
plt.gca().invert_yaxis()
plt.title('Reduced velocity [m/s] vs. Reduced elevator deflection [m]')
plt.legend(loc = 'upper right')
plt.show()

plt.figure(8)
plt.scatter(Vred2list, Fstlist, label = 'Measured data')
plt.xlabel('Reduced velocity [m/s]')
plt.ylabel('Elevator force [N]')
plt.gca().invert_yaxis()
plt.plot(vlist,Fveq(vlist),"r-",label = 'Regression line')
plt.title('Reduced velocity [m/s] vs. Elevator force [N]')
plt.grid()
plt.legend(loc = 'upper right')
plt.show()

print("Mach range series 2 =", min(machlist2),"to", max(machlist2))
print("Reynolds range series 2 (Vt used)= ", min(reynoldslist2),"to",max(reynoldslist2))
