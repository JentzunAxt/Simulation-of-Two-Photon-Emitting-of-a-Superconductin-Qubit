
Qubit Two Photon Laser



This program is intended to study the behavior of the system composed of a flux 
qubit under Rabi driving coupled to a reasonator.




1.parameters.py defines parameters of the system and functions to feed parameter to 
main_solver and other functions.






2.main_solver.py defines the function to solve the master equation, obtaining density
matrices at each time (store_states=True) and the expectation values of photon number
and transmittance.

intialization(N,omega_r,T) :  prepare a initial state with qubit in ground state and resonator in 
thermal equilibrium state.

evo_solver(N,timelist,args,Hamiltonian_args,dissipation_args,opts): solve the time-dependent 
master quation,  return density matrix and  expectations of N and T at each time.

final_solver(N,args,Hamiltonian_args,dissipation_args,method='direct'): use LU decomposition
to calculate the final state density matrix and  expectations of N and T.





3.plot_data.py defines functions to creat data for plots with parameters' sweep.





4.visualization1.py creates plots of FIG1, FIG2 in the article.(edit this)




5.visualization2.py creates plots of FIG3, FIG4 in the article.(edit this)




Tongliang and Yanjun

2021/12/15


