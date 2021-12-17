from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from parameters import  get_dissipation_param,get_Hamiltonian_param,get_param
from plot_data import  epsilon_Omega_R0_sweep,find12summit,omega_d_sweep,Omega_R0_sweep,epsilon_sweep
import matplotlib.cm as cm
from main_solver import final_solver

T0=1e-4
N = 10
Omega_R0_list = np.linspace(0.5,12,200)
Omega_R0_list2 = np.linspace(2,15,400)
omega_d_list = np.linspace(6,7,50)
epsilon_list= np.linspace(-1, 1, 100)
epsilon_list2 = np.linspace(-0.6, 0.6, 50)  # ground energy bias, approximatley 0

dissipation_args = get_dissipation_param()

'''
search for maximum transmittivity bias (red detuning)
'''
epsilon_listr = np.linspace(-1, 1, 500)
Nr,Tr = epsilon_sweep(N,epsilon_listr,6.2,1.8,dissipation_args)
maxindexr=np.where(Tr==np.max(Tr))
print(-epsilon_listr[maxindexr],Tr[maxindexr]/T0)

fig2dr,ax = plt.subplots()
fig2dr.suptitle('$\omega_d =6.2GHz$, $\Omega_{R0}=1.8GHz $',fontsize=20)
ax.plot(epsilon_listr,Tr/T0)
ax.vlines(-epsilon_listr[maxindexr], -0.5, 5, colors = "r", linestyles = "dashed")
plt.text(0.5, 0.1, '$\epsilon = 0.3828$')

ax.set_xlabel('$\epsilon$')
ax.set_ylabel('T')
plt.show()
'''
cross section at maximum transmittvity
'''
epsilonm= -epsilon_listr[maxindexr].item()
Omega_R0_listr = np.linspace(1, 12, 400)

Nrm,Trm = Omega_R0_sweep(N,epsilonm,6.0,Omega_R0_listr,dissipation_args)
fig2drm,ax = plt.subplots()
fig2drm.suptitle(' $\epsilon=0.3828$, $\omega_d =6.2GHz$',fontsize=20)
ax.plot(Omega_R0_listr,Trm/T0)
ax.set_xlabel('$\Omega_{R0}$')
ax.set_ylabel('T')
plt.show()
'''
search for maximum transmittivity bias (blue detuning)
'''
epsilon_listb = np.linspace(-0.6, 0.6, 300)
Nb,Tb = epsilon_sweep(N,epsilon_listb,6.6,8,dissipation_args)
maxindexr=np.where(Tb==np.max(Tb))
print(-epsilon_listb[maxindexr],Tb[maxindexr]/T0)

fig2db,ax = plt.subplots()
fig2db.suptitle('$\omega_d =6.6GHz$, $\Omega_{R0}=8GHz$',fontsize=20)
ax.plot(epsilon_listb,Tb/T0)
ax.vlines(-epsilon_listb[maxindexr], -0.5, 2.0, colors = "r", linestyles = "dashed")
plt.text(0.4, 0.1, '$\epsilon = 0.2990$')

ax.set_xlabel('$\epsilon$')
ax.set_ylabel('T')
plt.show()

epsilonm= -epsilon_listr[maxindexr].item()
Omega_R0_listr = np.linspace(2, 16, 400)
Nrm,Trm = Omega_R0_sweep(N,epsilonm,6.6,Omega_R0_listr,dissipation_args)
fig2drm,ax = plt.subplots()
fig2drm.suptitle('$\epsilon = 0.2990$,$\omega_d =6.6GHz$',fontsize=20)
ax.plot(Omega_R0_listr,Trm/T0)
ax.set_xlabel('$\Omega_{R0}$')
ax.set_ylabel('T')
plt.show()

#plot  Omega_R0-epsilon sweep (red detune)

NumberR,TR=epsilon_Omega_R0_sweep(N,6.0,epsilon_list,Omega_R0_list, dissipation_args)

np.savetxt('NumberR.csv',NumberR,fmt='%.20f',delimiter=',')
np.savetxt('TR.csv',TR,fmt='%.20f',delimiter=',')
#NumberR=np.loadtxt('NumberR.csv',delimiter=',')
#TR=np.loadtxt('TR.csv',delimiter=',')

figR = plt.figure()
figR.suptitle('Red detuning',fontsize=28)

axisR1 = figR.add_subplot(121, projection='3d')
XR1, YR1 = np.meshgrid(Omega_R0_list,epsilon_list)
surfR1=axisR1.plot_surface(XR1,YR1, NumberR,
                           cmap=cm.coolwarm, linewidth=0, antialiased=False)
axisR1.set_xlabel('$\Omega_{R0}$')
axisR1.set_ylabel('$\epsilon$')
axisR1.set_zlabel('$N$')

axisR2 = figR.add_subplot(122, projection='3d')
XR2, YR2 = np.meshgrid(Omega_R0_list,epsilon_list)
surfR2=axisR2.plot_surface(XR2, YR2, TR/T0,
                           cmap=cm.coolwarm, linewidth=0, antialiased=False)
axisR2.set_xlabel('$\Omega_{R0}$')
axisR2.set_ylabel('$\epsilon$')
axisR2.set_zlabel('$T$')

#plot  Omega_R0-epsilon sweep (blue detune)

NumberB,TB=epsilon_Omega_R0_sweep(N,6.6,epsilon_list2,Omega_R0_list2,dissipation_args)


np.savetxt('NumberB.csv',NumberB,fmt='%.20f',delimiter=',')
np.savetxt('TB.csv',TB,fmt='%.20f',delimiter=',')
#NumberB=np.loadtxt('NumberB.csv',delimiter=',')
#TB=np.loadtxt('TB.csv',delimiter=',')
figB = plt.figure()
figB.suptitle('Blue detuning',fontsize=28)

axisB1 = figB.add_subplot(121, projection='3d')
XB1, YB1 = np.meshgrid(Omega_R0_list2,epsilon_list2)

surfB1=axisB1.plot_surface(XB1,YB1, NumberB,
                           cmap=cm.coolwarm, linewidth=0, antialiased=False)
axisB1.set_xlabel('$\Omega_{R0}$')
axisB1.set_ylabel('$\epsilon$')
axisB1.set_zlabel('$N$')

axisB2 = figB.add_subplot(122, projection='3d')
XB2, YB2 = np.meshgrid(Omega_R0_list2,epsilon_list2)
surfB2=axisB2.plot_surface(XB2, YB2, TB/T0,
                           cmap=cm.coolwarm, linewidth=0, antialiased=False)
axisB2.set_xlabel('$\Omega_{R0}$')
axisB2.set_ylabel('$\epsilon$')
axisB2.set_zlabel('$T$')

plt.show()








