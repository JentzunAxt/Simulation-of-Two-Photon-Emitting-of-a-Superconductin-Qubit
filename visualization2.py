from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from parameters import  get_dissipation_param
from plot_data import  find12summit,TwoD
import matplotlib.cm as cm

T0=4.4e-9
N = 10
Omega_R0_list = np.linspace(0.001,30,100)
Omega_R0_list2 = np.linspace(0.1,5.1,100)
omega_d_list = np.linspace(1,15,100) 
omega_d_list2 = np.linspace(1,20,500)
epsilon_list = np.linspace(-0.25, 0.25, 100)  # ground energy bias, approximatley 0
dissipation_args = get_dissipation_param()

#Omega_R0-omega-d sweep

Nu,Tr = find12summit(N,0.06,omega_d_list,Omega_R0_list,dissipation_args)
Nu2,Tr2 = find12summit(N,0.06,omega_d_list,Omega_R0_list2,dissipation_args)
fig = plt.figure()
axiss = fig.add_subplot(121,projection='3d') 
X11, Y11 = np.meshgrid(omega_d_list,Omega_R0_list)
surf=axiss.plot_surface(X11,Y11, Nu,
                           cmap=cm.coolwarm, linewidth=0, antialiased=False)
axiss.set_xlabel('$\omega_d$')
axiss.set_ylabel('$\Omega_{R0}$')
axiss.set_zlabel('$N$')

axiss2 = fig.add_subplot(122,projection='3d')
X12, Y12 = np.meshgrid(omega_d_list,Omega_R0_list2)
surf=axiss2.plot_surface(X12, Y12, Tr2/T0,
                           cmap=cm.coolwarm, linewidth=0, antialiased=False)
axiss2.set_xlabel('$\omega_d$')
axiss2.set_ylabel('$\Omega_{R0}$')
axiss2.set_zlabel('$T$') 

#Section at Omega_R0 = 0.1
N2,T2 = TwoD(N,0.06,0.1,omega_d_list2,dissipation_args)
fig2d,ax = plt.subplots()
ax.plot(omega_d_list2,T2/T0)
ax.vlines(8.873, 0,3, colors = "r", linestyles = "dashed")
ax.vlines(11.355, 0,3, colors = "c", linestyles = "dashed")
plt.text(7, 2, '$\omega_{d} = 8.873(1 photon)$')
plt.text(10, 1, '$\omega_{d} = 11.355(2 photon)$')
ax.set_xlabel('$\omega_{d}$') 
ax.set_ylabel('T')

plt.show()






