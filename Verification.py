import numpy as np
from math import *
from parameters_def import *
from Analytical_model_def import *
from Numerical_Simulation import *

#Flight dependent paramters
V0,m,rho,muc,mub,CL,CD,CX0,CZ0 = eigenmotion_parameters(time,spmt0)



#-------------------------------------------------------------------------------------------------------------
#Analytical model
#-------------------------------------------------------------------------------------------------------------

lambda_c_spm1, lambda_c_spm2, period_spm, t_half_spm, zeta_spm, omega_n_spm = short_period(V0,m,rho,muc,mub,CL,CD,CX0,CZ0)
lambda_c_phug1, lambda_c_phug2, period_phug, t_half_phug, zeta_phug, omega_n_phug = phugoid(V0,m,rho,muc,mub,CL,CD,CX0,CZ0)
lambda_b_apr, t_half_apr, time_cst_apr = aperiodic_roll(V0,m,rho,muc,mub,CL,CD,CX0,CZ0)
lambda_b_spir, t_half_spir, time_cst_spir = spiral(V0,m,rho,muc,mub,CL,CD,CX0,CZ0)
lambda_b_dutch1, lambda_b_dutch2, period_dutch, t_half_dutch, zeta_dutch, omega_n_dutch = dutchroll(V0,m,rho,muc,mub,CL,CD,CX0,CZ0)


print("short_period")
print("-----------------------")
print("eigenvalues:", lambda_c_spm1, lambda_c_spm2)
print("period:",period_spm)
print("t_half:", t_half_spm)
print("zeta:", zeta_spm)
print("omega_n:", omega_n_spm)
print()
print("phugoid")
print("-----------------------")
print("eigenvalues:", lambda_c_phug1, lambda_c_phug2)
print("period",period_phug)
print("t_half:", t_half_phug)
print("zeta:", zeta_phug)
print("omega_n:", omega_n_phug)
print()
print("aperiodic roll")
print("-----------------------")
print("eigenvalue:", lambda_b_apr)
print("t_half:", t_half_apr)
print("time constant", time_cst_apr)
print()
print("spiral")
print("-----------------------")
print("eigenvalue:", lambda_b_spir)
print("t_half:", t_half_spir)
print("time constant", time_cst_spir)
print()
print("dutchroll")
print("-----------------------")
print("eigenvalues:", lambda_b_dutch1, lambda_b_dutch2)
print("period",period_dutch)
print("t_half:", t_half_dutch)
print("zeta:", zeta_dutch)
print("omega_n:", omega_n_dutch)

#-------------------------------------------------------------------------------------------------------------
#Numerical model
#-------------------------------------------------------------------------------------------------------------
sys_s, num_sym_par, sys_a, num_asym_par = numres(V0,m,rho,muc,mub,CL,CD,CX0,CZ0)
# print(num_sym_par)
# print(num_asym_par)

print("-----------------------")
print('Symmetric flight:')
print("-----------------------")
print('Eigenvalue 1: ', num_sym_par[0][0])
print('Period: ', num_sym_par[0][1])
print('halftime: ', num_sym_par[0][2] )
print('zeta: ', num_sym_par[0][3])
print('omega: ', num_sym_par[0][4])
print('')

print('Eigenvalue 2: ', num_sym_par[1][0])
print('Period: ', num_sym_par[1][1])
print('halftime: ', num_sym_par[1][2] )
print('zeta: ', num_sym_par[1][3])
print('omega: ', num_sym_par[1][4])

print("-----------------------")
print('Asymmetric flight:')
print("-----------------------")
print('Eigenvalue 1: ', num_asym_par[0][0])
print('Period: ', num_asym_par[0][1])
print('halftime: ', num_asym_par[0][2] )
print('zeta: ', num_asym_par[0][3])
print('omega: ', num_asym_par[0][4])
print('')

print('Eigenvalue 2: ', num_asym_par[1][0])
print('half time :', num_asym_par[1][1])
print('time constant: ', num_asym_par[1][2])
print('')

print('Eigenvalue 3: ', num_asym_par[2][0])
print('half time :', num_asym_par[2][1])
print('time constant: ', num_asym_par[2][2])
