from __future__ import division, print_function
from time import time
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from parameters import get_param, get_dissipation_param, get_Hamiltonian_param
from  main_solver import evo_solver,final_solver


def epsilon_Omega_R0_sweep(N,omega_d, epsilon_list, Omega_R0_list, dissipation_args):
    '''
    one order ,set N1,T1 matrix
    '''
    Number = np.zeros((len(epsilon_list), len(Omega_R0_list)))
    Trans = np.zeros((len(epsilon_list), len(Omega_R0_list)))
    i = 0
    j = 0

    for epsilon in epsilon_list:
        for Omega_R0 in Omega_R0_list:
            args = get_param(Omega_R0, omega_d, epsilon)
            '''
            calculate evolution of red and blue detune
            '''
            Hamiltonian_args = get_Hamiltonian_param(Omega_R0, omega_d, epsilon)  # (Omega_R0,omega_d,epsilon)
            method = 'direct'
            Number_expectation, t = final_solver(N, args, Hamiltonian_args, dissipation_args, method)

            '''
            extract expectations 
            '''

             # photon number expectation and transmittance of the resonator
            t_conj = np.real(t) - np.imag(t) * 1j  # t conjugate
            T = np.real(t_conj * t)  # Transmittivity T=|t|**2, take real part only

            '''
            save to value matrix
            '''
            Number[i, j] = Number_expectation
            Trans[i, j] = T

            print('%d,%d:\t \t number:\t\t %.8f \t\t T:\t\t%.8f' % (i, j, Number[i, j], Trans[i, j]))
            j = j + 1
            if (j == len(Omega_R0_list)): j = 0  # reset

        i = i + 1
        if ( i == len(epsilon_list) ): break

    return Number, Trans

def find12summit(N,epsilon,omega_d_list,Omega_R0_list,dissipation_args):
    '''
    one order ,set N1,T1 matrix
    '''
    Number = np.zeros((len(Omega_R0_list), len(omega_d_list)))
    Trans = np.zeros((len(Omega_R0_list), len(omega_d_list)))
    i = 0
    j = 0

    for Omega_R0 in Omega_R0_list:
        for omega_d in omega_d_list:
            args = get_param(Omega_R0, omega_d, epsilon)
            '''
            calculate evolution of red and blue detune
            '''
            Hamiltonian_args = get_Hamiltonian_param(Omega_R0, omega_d, epsilon)  # (Omega_R0,omega_d,epsilon)
            method = 'direct'
            Number_expectation, t = final_solver(N,args,Hamiltonian_args,dissipation_args,method)

            '''
            extract expectations 
            '''

            t_conj = np.real(t) - np.imag(t) * 1j  # t conjugate
            T = np.real(t_conj * t)  # Transmittivity T=|t|**2, take real part only

            '''
            save to value matrix
            '''
            Number[i, j] = Number_expectation
            Trans[i, j] = T

            print('%d,%d:\t \t number:\t\t %.8f \t\t T:\t\t%.8f' % (i, j, Number[i, j], Trans[i, j]))
            j = j + 1
            if (j == len(omega_d_list)):
                j = 0  # reset
                print('%d,%d\n' % (i, j))

        i = i + 1
        if (i == len(Omega_R0_list)): break

    return Number,Trans

def omega_d_sweep(N,epsilon,Omega_R0,omega_d_list,dissipation_args):
    '''
    one order ,set N1,T1 matrix
    '''
    Number = np.zeros(len(omega_d_list))
    Trans = np.zeros(len(omega_d_list))
    i = 0

    for i in np.arange(len(omega_d_list)):
            args = get_param(Omega_R0, omega_d_list[i], epsilon)
            '''
            calculate evolution of red and blue detune
            '''
            Hamiltonian_args = get_Hamiltonian_param(Omega_R0, omega_d_list[i], epsilon)  # (Omega_R0,omega_d,epsilon)
            method = 'svd'
            Number_expectation, t = final_solver(N,args,Hamiltonian_args,dissipation_args,method)

            '''
            extract expectations 
            '''

            #t_conj = np.real(t) - np.imag(t) * 1j  # t conjugate
           # T = np.real(t_conj * t)  # Transmittivity T=|t|**2, take real part only

            '''
            save to value matrix
            '''
            print(Number_expectation)
            Number[i] = Number_expectation
            Trans[i] = t

            print('%d:\t \t number:\t\t %.8f \t\t T:\t\t%.8f' % (i, Number[i], Trans[i]))

    return Number,Trans

def Omega_R0_sweep(N,epsilon,omega_d,Omega_R0_list,dissipation_args):
    '''
    one order ,set N1,T1 matrix
    '''
    Number = np.zeros(len(Omega_R0_list))
    Trans = np.zeros(len(Omega_R0_list))
    i = 0

    for i in np.arange(len(Omega_R0_list)):
            args = get_param(Omega_R0_list[i], omega_d, epsilon)
            '''
            calculate evolution of red and blue detune
            '''
            Hamiltonian_args = get_Hamiltonian_param(Omega_R0_list[i], omega_d, epsilon)  # (Omega_R0,omega_d,epsilon)
            method = 'power'
            Ne, t = final_solver(N,args,Hamiltonian_args,dissipation_args,method)

            '''
            extract expectations 
            '''

            #t_conj = np.real(t) - np.imag(t) * 1j  # t conjugate
           # T = np.real(t_conj * t)  # Transmittivity T=|t|**2, take real part only

            '''
            save to value matrix
            '''

            Number[i] = Ne
            Trans[i] = t

            print('%d:\t \t number:\t\t %.8f \t\t T:\t\t%.8f' % (i, Number[i], Trans[i]))

    return Number,Trans


def epsilon_sweep(N,epsilon_list,omega_d,Omega_R0,dissipation_args):
    '''
    one order ,set N1,T1 matrix
    '''
    Number = np.zeros(len(epsilon_list))
    Trans = np.zeros(len(epsilon_list))
    i = 0

    for i in np.arange(len(epsilon_list)):
            args = get_param(Omega_R0, omega_d,epsilon_list[i])
            '''
            calculate evolution of red and blue detune
            '''
            Hamiltonian_args = get_Hamiltonian_param(Omega_R0, omega_d,epsilon_list[i])  # (Omega_R0,omega_d,epsilon)
            method = 'direct'
            Ne, t = final_solver(N,args,Hamiltonian_args,dissipation_args,method)

            '''
            extract expectations 
            '''

            #t_conj = np.real(t) - np.imag(t) * 1j  # t conjugate
           # T = np.real(t_conj * t)  # Transmittivity T=|t|**2, take real part only

            '''
            save to value matrix
            '''
            Number[i] = Ne
            Trans[i] = t

            print('%d:\t \t number:\t\t %.8f \t\t T:\t\t%.8f' % (i, Number[i], Trans[i]))

    return Number,Trans

def TwoD(N,epsilon,Omega_R0,omega_d_list,dissipation_args):
    '''
    one order ,set N1,T1 matrix
    '''
    Number = np.zeros(len(omega_d_list))
    Trans = np.zeros(len(omega_d_list))
    i = 0

    for i in np.arange(len(omega_d_list)):
            args = get_param(Omega_R0, omega_d_list[i], epsilon)
            '''
            calculate evolution of red and blue detune
            '''
            Hamiltonian_args = get_Hamiltonian_param(Omega_R0, omega_d_list[i], epsilon)  # (Omega_R0,omega_d,epsilon)
            method = 'direct'
            Number_expectation, t = final_solver(N,args,Hamiltonian_args,dissipation_args,method)

            '''
            extract expectations 
            '''

            t_conj = np.real(t) - np.imag(t) * 1j  # t conjugate
            T = np.real(t_conj * t)  # Transmittivity T=|t|**2, take real part only

            '''
            save to value matrix
            '''
            Number[i] = Number_expectation
            Trans[i] = T

            print('%d:\t \t number:\t\t %.8f \t\t T:\t\t%.8f' % (i, Number[i], Trans[i]))

    return Number,Trans




