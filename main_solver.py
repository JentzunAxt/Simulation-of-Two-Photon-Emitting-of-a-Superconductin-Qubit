from __future__ import division, print_function
from re import I
import matplotlib.pyplot as plt
import numpy as np
from qutip import *



def intialization(N,omega_r,T):
    '''
    the resonator at thermal equilibrium state
    '''
    Nthermal=1 / (np.exp(omega_r / T) - 1)
    ar = destroy(N)
    adr = ar.dag()
    H_r = -omega_r * (adr * ar) / T
    rho_thermal = H_r.expm()
    Z = rho_thermal.tr()

    rho0 = tensor(rho_thermal/Z, fock_dm(2, 0))  # density matrix
    #print(fock_dm(2, 0))
    return rho0



def evo_solver(N,timelist,args,Hamiltonian_args,dissipation_args,opts):
    '''
    load parameters
    '''

    omega_r=args['omega_r']
    omega_g = args['omega_g']
    omega_d= args['omega_d']
    epsilon = args['epsilon']

    eta = Hamiltonian_args['eta']
    beta = Hamiltonian_args['beta']
    omega_q = Hamiltonian_args['omega_q']
    Omega_R = Hamiltonian_args['Omega_R']
    Nthermal = Hamiltonian_args['Nthermal']

    Gamma_0=dissipation_args['Gamma_0']
    Gamma_phi=dissipation_args['Gamma_phi']
    T=dissipation_args['T']
    k=dissipation_args['k']



    '''
    excitation,relaxation,dephasing rates 
    in rotating representation
    '''
    Gamma_u = (Gamma_0 / 4.0) * np.cos(eta) ** 2 * (1 +
              np.sin(beta)) ** 2 + (Gamma_phi / 2) * np.sin(eta) ** 2 * np.cos(
              beta) ** 2
    Gamma_d = (Gamma_0 / 4.0) * np.cos(eta) ** 2 * (1 -
              np.sin(beta)) ** 2 + (Gamma_phi / 2) * np.sin(eta) ** 2 * np.cos(
              beta) ** 2

    Gamma_phibar = (Gamma_0 / 2) * np.cos(eta) ** 2 * np.cos(beta) ** 2 \
                   + Gamma_phi * np.sin(eta) ** 2 * np.sin(beta) ** 2

    '''
    Define the Creation and Annihilation Operators:
    '''
    a = tensor(destroy(N), qeye(2))
    ad = a.dag()

    '''
    Define the Creation and Annihilation Operators:
    '''
    sx = tensor(qeye(N), sigmax())
    sz = tensor(qeye(N), sigmaz())
    sm = tensor(destroy(N), qeye(2))
    sp = sm.dag()


    '''
    solve master equation
    '''

    H = omega_r * ad * a + 1 / 2 * Omega_R * sz + omega_g * np.sin(eta) * (
            np.sin(beta) * sz - np.cos(beta) * sx) * (a + ad) - (omega_g ** 2 / omega_q) * \
        np.cos(eta) ** 2 * (np.sin(beta) * sz - np.cos(beta) * sx) * (a + ad) ** 2

    c_ops = [np.sqrt(Gamma_d) * sm, np.sqrt(
        Gamma_u) * sp, np.sqrt(Gamma_phibar) * sz, np.sqrt(k * (Nthermal + 1)) * a,
             np.sqrt(k * Nthermal) * ad]

    rho0=intialization(N, omega_r, T) # initial states
    Number = a * ad  # photon number
    e_ops = [Number, a]  # The observables whose expectations to be calculated

     # mesolver option, store all the density matrises.
    result = mesolve(H, rho0, timelist, c_ops, e_ops, options=opts)  # solve me.

    return result


def final_solver(N,args,Hamiltonian_args,dissipation_args,method='direct'):
    '''
    load parameters
    '''

    omega_r=args['omega_r']
    omega_g = args['omega_g']
    omega_d= args['omega_d']
    epsilon = args['epsilon']

    eta = Hamiltonian_args['eta']
    beta = Hamiltonian_args['beta']
    omega_q = Hamiltonian_args['omega_q']
    Omega_R = Hamiltonian_args['Omega_R']
    Nthermal = Hamiltonian_args['Nthermal']

    Gamma_0=dissipation_args['Gamma_0']
    Gamma_phi=dissipation_args['Gamma_phi']
    T=dissipation_args['T']
    k=dissipation_args['k']



    '''
    excitation,relaxation,dephasing rates 
    in rotating representation
    '''
    Gamma_u = (Gamma_0 / 4.0) * np.cos(eta) ** 2 * (1 +
              np.sin(beta)) ** 2 + (Gamma_phi / 2) * np.sin(eta) ** 2 * np.cos(
              beta) ** 2
    Gamma_d = (Gamma_0 / 4.0) * np.cos(eta) ** 2 * (1 -
              np.sin(beta)) ** 2 + (Gamma_phi / 2) * np.sin(eta) ** 2 * np.cos(
              beta) ** 2

    Gamma_phibar = (Gamma_0 / 2) * np.cos(eta) ** 2 * np.cos(beta) ** 2 \
                   + Gamma_phi * np.sin(eta) ** 2 * np.sin(beta) ** 2

    '''
    Define the Creation and Annihilation Operators:
    '''
    a = tensor(destroy(N), qeye(2))
    ad = a.dag()

    '''
    Define the Creation and Annihilation Operators:
    '''
    sx = tensor(qeye(N), sigmax())
    sz = tensor(qeye(N), sigmaz())
    sm = tensor(destroy(N), qeye(2))
    sp = sm.dag()
    '''
    solve master equation
    '''

    H = omega_r * ad * a + 1 / 2 * Omega_R * sz + omega_g * np.sin(eta) * (
            np.sin(beta) * sz - np.cos(beta) * sx) * (a + ad) - (omega_g ** 2 / omega_q) * \
        np.cos(eta) ** 2 * (np.sin(beta) * sz - np.cos(beta) * sx) * (a + ad) ** 2

    c_ops = [np.sqrt(Gamma_d) * sm, np.sqrt(
        Gamma_u) * sp, np.sqrt(Gamma_phibar) * sz, np.sqrt(k * (Nthermal + 1)) * a,
             np.sqrt(k * Nthermal) * ad]

    Number = ad * a  # photon number
     # solve steady state
    rhot = steadystate(H, c_ops, method=method, use_rcm=True)
    N=qutip.expect(Number, rhot)
    T=np.real(qutip.expect(a, rhot)*qutip.expect(ad, rhot))

    return N,T










