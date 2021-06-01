import numpy as np
import Parameters as pm
import control as ctr
import control.matlab as sim
import matplotlib.pyplot as plt
import scipy as sp

concept = pm.ConceptParameters(0)
I = concept.MMOI()
nState = 12
nInput = 4
nOutput = 6

d = ((concept.cabin.L_cabin/2)**2+(concept.cabin.W_cabin/2)**2)**0.5 + concept.propeller.D_prop/2
c = 1.5


A = np.matrix(np.zeros((nState,nState)))
A[0,3] = 1
A[1,4] = 1
A[2,5] = 1
A[3,7] = -concept.physics.g
A[4,6] = concept.physics.g
A[6,9] = 1
A[7,10] = 1
A[8,11] = 1

#A[0,3] = A[1,4] = A[2,5] = 1

B = np.matrix(np.zeros((nState,nInput)))
B[5,0] = 1/concept.Mtot_concept
B[9,1] = 1/I[0]
B[10,2] = 1/I[1]
B[11,3] = 1/I[2]

print(B.shape)
# B[3] = [0,d/I[0],0,-d/I[0]]
# B[4] = [d/I[1],0,-d/I[1],0]
# B[5] = [-c/I[2],c/I[2],-c/I[2],c/I[2]]


C = np.matrix(np.zeros((nOutput,nState)))

C[0,0] = 1
C[1,1] = 1
C[2,2] = 1

C[3,6] = 1
C[4,7] = 1
C[5,8] = 1

D = np.matrix(np.zeros((nOutput,nInput)))
K = np.matrix(np.zeros((nInput,nOutput)))

az = np.loadtxt(r'Data/az.txt')
Mx = np.loadtxt(r'Data/Mx.txt')
My = np.loadtxt(r'Data/My.txt')

ss = ctr.StateSpace(A,B,C,D)
Ka = np.zeros((nInput,nOutput))
Ka[0,2] = -0.5
Ka[-1,-1] = 100
print(Ka)
notDone = True
i = 0
nBest = -1
bestK = K
imax = -1#9e9
nGood = 0

while nGood != 12 and i  < imax:
    i+=1

    Ka[:,:] = np.random.uniform(-1,1,(nInput,nOutput))*1000

    # kz = 0.0
    # K[0,2] = 0.
    # K[2,1] = 0.
    # K[0,2] = 0.
    # K[2,3] = 0. #impacts pitch
    # K[2,5] = 0.



    # K[1,3] = 5000
    # K[1,1] = 0.01


    # K[2,4] = 10000
    # K[2,1] = -0.01
    # K[2,5] = 0

    # K2 = ctr.place(A,B,[-0.1,-0.1,-0.2,-0.2,-0.3,-.3,-0.01,-0.01,-0.02,-0.02,-0.03,-0.03])


    # C[0,0] = 1
    # C[1,1] = 1
    # C[2,2] = 1

    #U = [T-mg,pitch torque, roll torque,yaw torque]


    # ss = ctr.minreal(ss)
    ssi = ss.feedback(Ka)
    poles = ssi.pole()
    nGood = 0

    for pole in poles:
        if pole.real <= 0:
            nGood += 1
        else:
            break

    if nGood > nBest:
        nBest = nGood
        bestK = Ka
        np.savetxt(r'Data/bestK.txt',bestK)
    if i%5000 ==0:
        print(nBest,i)



# ss= ss.append(K2.gain_matrix)
# ss = ctr.minreal(ss)


bestK = np.loadtxt(r'Data/bestK.txt')
ss = ss.feedback(bestK)
ctr.pzmap(ss)
tf = ctr.tf(ss)


# ss = ss.minreal()
T = np.arange(0,2500,0.01)
U = np.zeros((len(T),4))
U[10:50,1] = 1
U[10:50,2] = 1
# U[:len(My),2] = My
# U[:len(az),0] = az*concept.Mtot_concept #T-mg
#U[len(az):int(np.round(len(az)*39.75)),0] = az[-1]*concept.Mtot_concept


# U[len(U)//17:len(U)//16,2] = -10
# U[len(U)//10:len(U)//8,2] = 3
y,t,x = sim.lsim(ss,U,T)




# Plotting functions Symmetric
fig,ax = plt.subplots(1,1,figsize=(5,10))
rx = ax.plot(t,U)
fig.suptitle('OEI Octorotor control input ')
ax.legend(iter(rx), ('Force in z','Pitch Torque','Roll Torque','Yaw Torque'))

fig1, ax1 = plt.subplots(3, 2, figsize=(14, 10))
ax1 = ax1.ravel()
fig1.suptitle('OEI Octorotor Control Response')
ax1[0].plot(t, y[:, 3]*57.3)
ax1[0].title.set_text(r"Roll angle $\phi$ over time")
ax1[0].set_ylabel(r'$\phi$ [deg]')
ax1[0].set_xlabel(r't [s]')
ax1[1].plot(t, y[:, 4]*57.3)
ax1[1].title.set_text(r"Pitch angle $\theta$ over time")
ax1[1].set_ylabel(r'$\theta$ [deg]')
ax1[1].set_xlabel(r't [s]')
ax1[2].plot(t, y[:, 5]*57.3)
ax1[2].title.set_text(r"Yaw angle $\chi$ over time")
ax1[2].set_ylabel(r'$\chi$ [deg]')
ax1[2].set_xlabel(r't [s]')
ax1[3].plot(T, y[:,0])
ax1[3].title.set_text(r"X position over time")
ax1[3].set_ylabel(r'X [m]')
ax1[3].set_xlabel(r't [s]')
ax1[4].plot(T, y[:,1])
ax1[4].title.set_text(r"Y position over time [m]")
ax1[4].set_ylabel(r'Y [m]')
ax1[4].set_xlabel(r't [s]')
ax1[5].plot(T, y[:,2])
ax1[5].title.set_text(r"Z position over time [m]")
ax1[5].set_ylabel(r'Z [m]')
ax1[5].set_xlabel(r't [s]')


# ax1[6].plot(T, U)
# ax1[6].title.set_text(r"Input Function")

plt.show()
