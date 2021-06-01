import numpy as np
import Parameters as pm
import control as ctr
import control.matlab as sim
import matplotlib.pyplot as plt
import scipy as sp


concept = pm.ConceptParameters(0)
I = concept.MMOI()

#First system
#States: z, z', yaw (chi), yaw rate (chi')
#Inputs: T-mg (Fz), yaw rate (q)
#Outputs: z, yaw

n_input1 = 2
n_state1 = 4
n_output1 = 2

A1 = np.zeros((n_state1,n_state1))
A1[0,1] = 1
A1[2,3] = 1
print(A1)


B1 = np.zeros((n_state1,n_input1))
B1[1,0] = -1/concept.Mtot_concept
B1[2,1] = 1/I[0]


C1 = np.zeros((n_output1,n_state1))
C1[0,0] = 1
C1[1,2] = 1


D1 = np.zeros((n_output1,n_input1))

K = np.zeros((n_input1,n_output1))
K[1,1] = 100
K[0,0] = -0.5


ss1 = ctr.StateSpace(A1,B1,C1,D1)
ss1 = ss1.feedback(K)
ctr.pzmap(ss1)


T1 = np.arange(0,3000,0.1)
U1 = np.zeros((len(T1),n_input1))
U1[len(T1)//10:2*len(T1)//10,1] = 100
U1[100,0] = 100


y1,t1,x1 = sim.lsim(ss1,U1,T1)
#y1,t1 = sim.impulse(ss1,T1,input=0)



# Plotting functions Symmetric
fig,ax = plt.subplots(1,1,figsize=(5,10))
rx = ax.plot(T1,U1)
fig.suptitle('OEI Octorotor control input ')
ax.legend(iter(rx), ('Force in z','Yaw Torque'))

fig1, ax1 = plt.subplots(1, 2, figsize=(14, 10))
ax1 = ax1.ravel()
fig1.suptitle('Octorotor Control Response')
# ax1[0].plot(t, y[:, 3]*57.3)
# ax1[0].title.set_text(r"Roll angle $\phi$ over time")
# ax1[0].set_ylabel(r'$\phi$ [deg]')
# ax1[0].set_xlabel(r't [s]')
# ax1[1].plot(t, y[:, 4]*57.3)
# ax1[1].title.set_text(r"Pitch angle $\theta$ over time")
# ax1[1].set_ylabel(r'$\theta$ [deg]')
# ax1[1].set_xlabel(r't [s]')
ax1[1].plot(t1, y1[:, 1]*57.3)
ax1[1].title.set_text(r"Yaw angle $\chi$ over time")
ax1[1].set_ylabel(r'$\chi$ [deg]')
ax1[1].set_xlabel(r't [s]')
# ax1[3].plot(T, y[:,0])
# ax1[3].title.set_text(r"X position over time")
# ax1[3].set_ylabel(r'X [m]')
# ax1[3].set_xlabel(r't [s]')
# ax1[4].plot(T, y[:,1])
# ax1[4].title.set_text(r"Y position over time [m]")
# ax1[4].set_ylabel(r'Y [m]')
# ax1[4].set_xlabel(r't [s]')
ax1[0].plot(t1, y1[:,0])
ax1[0].title.set_text(r"Z position over time [m]")
ax1[0].set_ylabel(r'Z [m]')
ax1[0].set_xlabel(r't [s]')


# ax1[6].plot(T, U)
# ax1[6].title.set_text(r"Input Function")

plt.show()