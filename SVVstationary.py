import numpy as np
import scipy.io as sc

# paramters
W_empty = 9165 #kg
blockfuel = 4100 #lbs
masspas = np.array([90,102,80,83,94,84,74,79,103]) #kg
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

#%% data unit conversion
bf_kg = blockfuel*0.453592  #kg
#mat1 
h_mat1 = mat1[:,0]*0.3048           # m
IAS_mat1 = mat1[:,1]*0.514444       # m/s
AOA_mat1 = mat1[:,2]                #degrees
FFL_mat1 = mat1[:,3]* (1/7936.64)   # kg/s
FFR_mat1 = mat1[:,4]* (1/7936.64)   #kg/s
WF_mat1 = mat1[:,5]* 0.453592       #kg
TAT_mat1 = mat1[:,6]+273.15         #kelvin    

#mat2
h_mat2 = mat2[:,0]*0.3048           # m
IAS_mat2 = mat2[:,1]*0.514444       # m/s
AOA_mat2 = mat2[:,2]                #degrees
DE_mat2 = mat2[:,3]                #degrees
DETR_mat2 = mat2[:,4]                #degrees
Fe_mat2 = mat2[:,5]                #N
FFL_mat2 = mat2[:,6]* (1/7936.64)   # kg/s
FFR_mat2 = mat2[:,7]* (1/7936.64)   #kg/s
WF_mat2 = mat2[:,8]* 0.453592       #kg
TAT_mat2 = mat2[:,9]+273.15         #kelvin  

#mat3
h_mat3 = mat3[:,0]*0.3048           # m
IAS_mat3 = mat3[:,1]*0.514444       # m/s
AOA_mat3 = mat3[:,2]                #degrees
DE_mat3 = mat3[:,3]                #degrees
DETR_mat3 = mat3[:,4]                #degrees
Fe_mat3 = mat3[:,5]                #N
FFL_mat3 = mat3[:,6]* (1/7936.64)   # kg/s
FFR_mat3 = mat3[:,7]* (1/7936.64)   #kg/s
WF_mat3 = mat3[:,8]* 0.453592       #kg
TAT_mat3 = mat3[:,9]+273.15         #kelvin  
#%% CL-alpha curve
CL = np.zeros(6)
p = np.zeros(6)
rho = np.zeros(6)
L = np.zeros(6)
V = np.zeros(6)
M = np.zeros(6)
T = np.zeros(6)
for i in range (0,6):
    p[i] = p0*(1+(lda*(float(h_mat1[i])))/(T0))**(-g0/(R*lda))
    M[i] = m.sqrt((2/(y-1))*(((1+p0/float(p[i]))*(1+((y-1)/(2*y))*(rho0/p0)*float(IAS_mat1[i])**2)**(y/(y-1))-1)**((y-1)/y)-1))
    rho[i] = float(p[i])/(R*float(TAT_mat1[i]))
    L[i] = (M_total - float(WF_mat1[i]))*g0
    T[i] = float(TAT_mat1[i])/(1+((y-1)/2*float(M[i])**2))
    V[i] = float(M[i])*m.sqrt(y*R*float(T[i]))
    CL[i] = 2*float(L[i])/(float(rho[i])*(float(V[i])**2*S))
#%%
plt.figure
plt.plot(np.radians(AOA_mat1),np.flip(CL))
plt.show()











