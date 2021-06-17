import numpy as np
import Parameters as pm
import control as ctr
import control.matlab as sim
import matplotlib.pyplot as plt
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D

concept = pm.ConceptParameters(0)
I = concept.MMOI()
nState = 12
nInput = 10
nOutput = 10

"""

State vector:
x, y, z, x', y', z', phi, theta, chi, phi', theta', chi'
0  1  2  3   4   5   6    7      8    9     10      11

Input Vector:
ft, fwx, fwy, fwz, Tx, Ty, Tz, Twx, Twy, Twz
0    1    2    3   4    5   6  7    8    9
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
C = np.zeros((nOutput+4,nState))
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

D = np.zeros((nOutput+4,nInput))
D[10,4] = 1
D[11,5] = 1
D[12,6] = 1
D[13,0] = 1

K = np.zeros((nInput,nOutput+4))

K_z_z = -10000


# K_tx_y = 1500
# K_tx_ydot = 3000
# K_tx_phi = 15000
# K_tx_phidot = 13000


# K_ty_x = -75
# K_ty_xdot = -300
# K_ty_theta = 4000
# K_ty_thetadot = 4500

K_tx_y = 250_000
K_tx_ydot = 300_000
K_tx_phi = 750_000
K_tx_phidot = 70_000


K_ty_x = -250_000
K_ty_xdot = -300_000
K_ty_theta = 750_000
K_ty_thetadot = 70_000



K_tz_yaw  = 5000


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
T = np.arange(0,250+dt,dt)
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
T_cg = np.loadtxt(r'Data/Gust/cgtorque.txt')
t_cg = np.loadtxt(r'Data/Gust/cgtime.txt')

# M_gust[:,1] *= -1
# U[len(t_gust)//2:len(t_gust)+len(t_gust)//2,1:4] = F_gust[:,:]
# U[len(t_gust)//2:len(t_gust)+len(t_gust)//2,7:] = M_gust[:,:]
# U[len(t_cg)//2:len(t_cg)+len(t_cg)//2,8] = T_cg*2
# U[len(t_cg)//2:len(t_cg)+len(t_cg)//2,7] = T_cg*2
# U[len(T)//4:int(np.round(len(t_gust)*0.08))+len(T)//4,1] = 200
# U[len(T)//4:int(np.round(len(t_gust)*0.08))+len(T)//4,2] = 200
# U[len(T)//4:int(np.round(len(t_gust)*0.08))+len(T)//4,3] = 200
# U[len(T)//4:int(np.round(len(t_gust)*0.08))+len(T)//4,-3] = 60
# U[len(T)//4:int(np.round(len(t_gust)*0.08))+len(T)//4,-2] = -60
# U[len(T)//4:int(np.round(len(t_gust)*0.08))+len(T)//4,-1] = 60
U[len(T)//16:len(T)//16+len(T)//4,0] = np.sin(2*np.pi*T[:len(T)//4])*100
# U[len(T)//4:int(len(T)//4 + 60//dt), 7] = 2353
# U[len(T)//4:int(len(T)//4 + 60//dt), 8] = 1054

# U[len(az_OEI)//8:len(az_OEI)+len(az_OEI)//8,3] = az_OEI
# U[len(Mx_OEI)//8:len(Mx_OEI)+len(Mx_OEI)//8,-3] = Mx_OEI
# U[len(My_OEI)//8:len(My_OEI)+len(My_OEI)//8,-2] = My_OEI

# U[len(T)//4:,1] = 1
# U[len(T)//4:,2] = 1
# U[len(T)//4:,3] = 1
# U[len(T)//4:,-3] = 1
# U[len(T)//4:,-2] = -1
# U[len(T)//4:,-1] = 1

# U[:len(My),5] = -My
# U[:,8] = np.append([100*(np.sin(2*np.pi - t/250)) for t in range(5000)] , np.zeros(len(T)-5000))
# U[:,5] = np.append([2*np.pi - t**2 /1e5 for t in range(3000)] , np.zeros(len(T)-3000))
# U[len(T)//8:len(T)//7,8] = 100
# U[:,4] = [10*np.cos(np.pi+ t) for t in T]
# U[len(T)//4,:] = 1

y,t,x = sim.lsim(ss,U,T,[0,0,0,0,0,0,0,0,0,0,0,0])
#y,t,x = sim.lsim(ss,U,T,[1,1,1,1,1,1,1,1,1,1,1,1])


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
ax1[5].set_ylabel(r'v [m/s]')
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

# ax1[10].scatter(t,y[:,0],y[:,1])
# # ax1[10].plot(t, y[:, -1])
# ax1[10].title.set_text(r"Control torque")
# ax1[10].set_ylabel(r'$\tau$ [deg]')
# ax1[10].set_xlabel(r't [s]')


ctr.pzmap(ss)
fig1.subplots_adjust(top=0.92,
bottom=0.065,
left=0.065,
right=0.985,
hspace=0.405,
wspace=0.355)




fig, ax = plt.subplots(1,1)
rx = ax.plot(t,U)
ax.set_xlabel('Time t [s]')
ax.set_ylabel(r'Force [N], Torque [Nm]')
ax.legend(iter(rx), ('Force in z','Wind Force in x','Wind Force in y','Wind Force in z','Torque around x','Torque around y','Torque around z','Wind Torque around x','Wind Torque around y','Wind Torque around z'))
ax.title.set_text(r'Input function')

fig, ax = plt.subplots(1,1)
control = y[:,-4:]
print(np.amin(y[:,:]))
rx = ax.plot(t,y)
ax.set_xlabel('Time t [s]')
ax.set_ylabel(r'Force [N], Torque [Nm]')
ax.legend(iter(rx), ('Torque around x','Torque around y','Torque around z','Force in z'))
ax.title.set_text(r'Control inputs due to disturbance')

fig = plt.figure()
ax3d = Axes3D(fig)
ax3d.scatter(y[:,0],y[:,1],y[:,3])

plt.show()






