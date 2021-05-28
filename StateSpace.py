import numpy as np
import Parameters as pm
import control as ctr
import control.matlab as sim
import matplotlib.pyplot as plt

concept = pm.ConceptParameters(0)
I = concept.MMOI()
nState = 12
nInput = 4
nOutput = 6

d = ((concept.cabin.L_cabin/2)**2+(concept.cabin.W_cabin/2)**2)**0.5 + concept.propeller.D_prop/2
c = 1.5


A = np.zeros((nState,nState))
A[0,3] = 1
A[1,4] = 1
A[2,5] = 1
A[3,7] = -concept.physics.g
A[4,6] = concept.physics.g
A[6,9] = 1
A[7,10] = 1
A[8,11] = 1
print(A)
#A[0,3] = A[1,4] = A[2,5] = 1

B = np.zeros((nState,nInput))
B[5,0] = 1/concept.Mtot_concept
B[9,1] = 1/I[0]
B[10,2] = 1/I[1]
B[11,3] = 1/I[2]
print(B)
# B[3] = [0,d/I[0],0,-d/I[0]]
# B[4] = [d/I[1],0,-d/I[1],0]
# B[5] = [-c/I[2],c/I[2],-c/I[2],c/I[2]]


C = np.zeros((nOutput,nState))

C[0,0] = 1
C[1,1] = 1
C[2,2] = 1

C[3,6] = 1
C[4,7] = 1
C[5,8] = 1

D = np.zeros((nOutput,nInput))
# C[0,0] = 1
# C[1,1] = 1
# C[2,2] = 1

#U = [yaw torque,pitch torque, roll torque, T-mg]


ss = ctr.StateSpace(A,B,C,D)
T = np.arange(0,2.5,0.001)
U = np.zeros((len(T),4))
U[:,2] = -(1627+4438)/2
U[:,3] = -0.3*concept.physics.g*concept.Mtot_concept #T-mg

# U[len(U)//17:len(U)//16,2] = -10
# U[len(U)//10:len(U)//8,2] = 3
y,t,x = sim.lsim(ss,U,T)




# Plotting functions Symmetric
plt.plot(t,U)
plt.show()
fig1, ax1 = plt.subplots(3, 2, figsize=(15, 10))
ax1 = ax1.ravel()
ax1[0].plot(t, y[:, 0]*57.3)
ax1[0].title.set_text(r"Roll angle $\phi$")
ax1[1].plot(t, y[:, 1]*57.3)
ax1[1].title.set_text(r"Pitch angle $\theta$")
ax1[2].plot(t, y[:, 2]*57.3)
ax1[2].title.set_text(r"Yaw angle $\chi$")
ax1[3].plot(T, y[:,3])
ax1[3].title.set_text(r"X")
ax1[4].plot(T, y[:,4])
ax1[4].title.set_text(r"Y")
ax1[5].plot(T, y[:,5])
ax1[5].title.set_text(r"Z")

# ax1[6].plot(T, U)
# ax1[6].title.set_text(r"Input Function")

plt.show()

