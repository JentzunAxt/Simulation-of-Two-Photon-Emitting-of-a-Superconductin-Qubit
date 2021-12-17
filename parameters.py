
'''
Call these functions to load parameters
'''


import numpy as np


'''
Objects have a dimention of ([N,N],[2,2])
'''
'''
scaling parameters (natural unit system)
Here we use GHz
'''

hbar=1/(2*np.pi)
scaling =1e9 #GHz

N =10 # The number of Fock states of the system
#omega_g=3.3
omega_g=3.3# The coupling energy between the quibit and the resonator
omega_r = 2.48 # Fundamental frequency of the resonator
Delta =6.391 # level separation

'''
Dissipation factors
'''
Gamma_0=0.125e-4# The spontaneous emission rate
Gamma_phi=0.013e-4 # Dephasing rate
T=1.3 #temperature of the resonator
k=1.7e-6 #resonator loss rate
'''
planck constants and scale factor (Ghz recommended)
'''
def scaling_param():
    scaling_args={'hbar': hbar, 'scaling':scaling}
    return scaling_args
  

'''
physical parameters of the system,no omega_d and epsilon
'''
def get_param(Omega_R0,omega_d,epsilon) :
    args = {'N': N, 'Omega_R0': Omega_R0,'omega_g': omega_g,
            'omega_r': omega_r,'omega_d': omega_d,
            'epsilon': epsilon,'Delta': Delta}

    return args

'''
parameters that dictates the dissipation terms
'''
def get_dissipation_param() :
    args = {'Gamma_0': Gamma_0, 'Gamma_phi': Gamma_phi,
            'T': T,'k':k}

    return args


'''
parameters in the Hamiltonian
'''
def get_Hamiltonian_param(Omega_R0,omega_d,epsilon) :
    omega_q = np.sqrt(epsilon ** 2 + Delta ** 2)  # eigen frequency of qubit system
    delta_omega = omega_d - omega_q
    eta = np.arctan(epsilon / Delta)
    Omega_R = np.sqrt(Omega_R0 ** 2 * np.cos(eta) ** 2 + delta_omega ** 2)
    beta = np.arcsin(delta_omega / Omega_R )

    Nthermal = 1 / (np.exp(omega_r / T) - 1)  # thermal distribution for the resonator

    args = {'omega_q': omega_q,
            'delta_omega':delta_omega ,
            'eta': eta,
            'beta': beta,
            'Omega_R': Omega_R,
            'Nthermal': Nthermal,
            }
    return args
