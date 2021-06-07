import numpy as np
import Parameters as pm
import control as ctr
import control.matlab as sim
import matplotlib.pyplot as plt
import scipy as sp

concept = pm.ConceptParameters(0)
I = concept.MMOI()
nState = 12
nInput = 10
nOutput = 10

"""

State vector:
x, y, z, x', y', z', phi, theta, chi, phi', theta', chi'

Input Vector:
ft, fwx, fwy, fwz, Tx, Ty, Tz, Twx, Twy, Twz

Output Vector:
x, y, z, phi, theta, chi, x', y', phi', theta'

"""

#State Matrix
A = np.zeros((nState,nState))
A[0,3] = 1
A[1,4] = 1
A[2,5] = 1
A[3,7] = -concept.physics.g
A[4,6] = concept.physics.g
A[6,9] = 1
A[7,10] = 1
A[8,11] = 1

#Input Matrix
B = np.zeros((nState,nInput))
B[3,1] = 1/concept.Mtot_concept
B[4,2] = 1/concept.Mtot_concept
B[5,0] = -1/concept.Mtot_concept
B[5,3] = 1/concept.Mtot_concept
B[9,4] = 1/I[0]
B[9,7] = 1/I[0]
B[10,5] = 1/I[1]
B[10,8] = 1/I[1]
B[11,6] = 1/I[2]
B[11,9] = 1/I[2]

#Output matrix
C = np.zeros((nOutput,nState))
C[0,0] = 1
C[1,1] = 1
C[2,2] = 1
C[3,6] = 1
C[4,7] = 1
C[5,8] = 1
C[6,3] = 1
C[7,4] = 1
C[8,9] = 1
C[9,10] = 1

D = np.zeros((nOutput,nInput))

K = np.zeros((nInput,nOutput))
K_z_z = -10000

K_tx_y = 100
K_tx_ydot = 300
K_tx_phi = 6500
K_tx_phidot = 2500

K_ty_x = -10
K_ty_xdot = -100
K_ty_theta = 4000
K_ty_thetadot = 5000

K_tz_yaw  = 10000


K[0,2] = K_z_z

K[4,8] = K_tx_phidot
K[4,3] = K_tx_phi
K[4,7] = K_tx_ydot
K[4,1] = K_tx_y

K[5,9] = K_ty_thetadot
K[5,4] = K_ty_theta
K[5,6] = K_ty_xdot
K[5,0] = K_ty_x

K[6,5] = K_tz_yaw

ss = ctr.StateSpace(A,B,C,D)
ss = ss.feedback(K)

dt = 0.01
T = np.arange(0,350+dt,dt)
U = np.zeros((len(T), nInput))


az_OEI = np.loadtxt(r'Data/az.txt')
Mx_OEI = np.loadtxt(r'Data/Mx.txt')
My_OEI = np.loadtxt(r'Data/My.txt')

a_gust = np.loadtxt(r'Data/Gust/agust.txt')
F_gust = np.loadtxt(r'Data/Gust/Dgust.txt')
M_gust = np.loadtxt(r'Data/Gust/Mgust.txt')
alpha_gust = np.loadtxt(r'Data/Gust/alphagust.txt')
t_gust = np.loadtxt(r'Data/Gust/tgust.txt')
U_gust = np.loadtxt(r'Data/Gust/Ugust.txt')


# U[:len(az),0] = 0#az*concept.Mtot_concept
# U[:,4] = 300

U[len(t_gust):2*len(t_gust),1:4] = F_gust[:,:]
U[len(t_gust):2*len(t_gust),7:] = M_gust[:,:]



# U[:len(My),5] = -My
# U[:,8] = np.append([100*(np.sin(2*np.pi - t/250)) for t in range(5000)] , np.zeros(len(T)-5000))
# U[:,5] = np.append([2*np.pi - t**2 /1e5 for t in range(3000)] , np.zeros(len(T)-3000))
# U[len(T)//8:len(T)//7,8] = 100
# U[:,4] = [10*np.cos(np.pi+ t) for t in T]
y,t,x = sim.lsim(ss,U,T)






fig1, ax1 = plt.subplots(3,4,figsize=(15,7.5))
ax1 = ax1.ravel()
fig1.suptitle('OEI Octorotor Control Response')

ax1[0].plot(T, y[:,0])
ax1[0].title.set_text(r"X position over time")
ax1[0].set_ylabel(r'X [m]')
ax1[0].set_xlabel(r't [s]')

ax1[1].plot(T, y[:,6])
ax1[1].title.set_text(r"X velocity over time")
ax1[1].set_ylabel(r'u [m/s]')
ax1[1].set_xlabel(r't [s]')

ax1[2].plot(t, y[:, 4]*57.3)
ax1[2].title.set_text(r"Pitch angle $\theta$ over time")
ax1[2].set_ylabel(r'$\theta$ [deg]')
ax1[2].set_xlabel(r't [s]')

ax1[3].plot(t, y[:, 9]*57.3)
ax1[3].title.set_text(r"Pitch rate $\theta$' over time")
ax1[3].set_ylabel(r"$\theta$' [deg/s]")
ax1[3].set_xlabel(r't [s]')

ax1[4].plot(T, y[:,1])
ax1[4].title.set_text(r"Y position over time")
ax1[4].set_ylabel(r'Y [m]')
ax1[4].set_xlabel(r't [s]')

ax1[5].plot(T, y[:,7])
ax1[5].title.set_text(r"Y velocity over time")
ax1[5].set_ylabel(r'u [m/s]')
ax1[5].set_xlabel(r't [s]')

ax1[6].plot(t, y[:, 3]*57.3)
ax1[6].title.set_text(r"Roll angle $\phi$ over time")
ax1[6].set_ylabel(r'$\phi$ [deg]')
ax1[6].set_xlabel(r't [s]')

ax1[7].plot(t, y[:, 8]*57.3)
ax1[7].title.set_text(r"Roll rate $\phi$' over time")
ax1[7].set_ylabel(r"$\phi$' [deg/s]")
ax1[7].set_xlabel(r't [s]')

ax1[8].plot(T, y[:,2])
ax1[8].title.set_text(r"Z position over time [m]")
ax1[8].set_ylabel(r'Z [m]')
ax1[8].set_xlabel(r't [s]')

ax1[9].plot(t, y[:, 5]*57.3)
ax1[9].title.set_text(r"Yaw angle $\chi$ over time")
ax1[9].set_ylabel(r'$\chi$ [deg]')
ax1[9].set_xlabel(r't [s]')

ctr.pzmap(ss)

plt.subplots_adjust(0.065,0.075,0.985,0.91,0.375,0.52)




fig, ax = plt.subplots(1,1)
rx = ax.plot(t,U)
ax.legend(iter(rx), ('Force in z','Wind Force in x','Wind Force in y','Wind Force in z','Torque around x','Torque around y','Torque around z','Wind Torque around x','Wind Torque around y','Wind Torque around z'))
ax.title.set_text(r'Input function')



plt.show()






