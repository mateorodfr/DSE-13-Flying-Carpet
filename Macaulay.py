import numpy as np
import Parameters as pm
import matplotlib.pyplot as plt

def macaulay(z,z0,power):
    return ((z > z0) * (z>0))*(z-z0) ** power
def getConstant(concept,Ne,EI,L):
    Ze = 0
    b1 = (1/(2*EI))*Ne*macaulay(L,Ze,2)
    b2 = (1/(6*EI))*Ne*macaulay(L,Ze,3)
    b = np.array([[b1],[b2]])
    A = np.array([[1,0],[L,1]])
    return np.matmul(np.linalg.inv(A),b)
def beamDeflectY(z, Ne,concept,EI,L, macaulay,gc):
    Cs = gc(concept,Ne,EI,L)
    Vofx = Ne*macaulay(z,0,0)
    Mofx = Ne * macaulay(z,0,1)
    ThetaY = (-1/EI)*(0.5*Ne*macaulay(z,0,2))+Cs[0]
    DeflectionY = (-1/EI) * ((1/6)*Ne*macaulay(z,0,3)) + (Cs[0]*z + Cs[1])
    return Vofx, Mofx, ThetaY, DeflectionY
def getConstant2(concept,Teng,EI,L):
    b1 = (1/(EI))*Teng*macaulay(L,Ze,1)
    b2 = (1/(2*EI))*Teng*macaulay(L,Ze,2)
    b = np.array([[b1],[b2]])
    A = np.array([[1,0],[L,1]])
    return np.matmul(np.linalg.inv(A),b)






concept = pm.ConceptParameters(0)
EI = 1e6
L_arm = 1
L = L_arm + concept.propeller.D_prop/2
TW = 1
dz = 0.01
Ne = TW*concept.Mtot_concept*concept.physics.g/4
z = np.arange(0,L+dz,dz)

Vofy, Mofx, Thetay, DeflectionY = beamDeflectY(z,Ne,concept,EI,L,macaulay,getConstant)

fig1, ax1 = plt.subplots(2,2)
ax1 = ax1.ravel()
fig1.suptitle('Beam loading in Y direction')

ax1[0].plot(z,Vofy)
ax1[0].title.set_text(r"Shear force in Y")
ax1[0].set_ylabel(r'V [N]')
ax1[0].set_xlabel(r'z [m]')

ax1[1].plot(z,Mofx)
ax1[1].title.set_text(r"Moment around X")
ax1[1].set_ylabel(r'M [Nm]')
ax1[1].set_xlabel(r'z [m]')

ax1[2].plot(z,Thetay)
ax1[2].title.set_text(r"Deflection angle in Y")
ax1[2].set_ylabel(r'$\alpha$ [rad]')
ax1[2].set_xlabel(r'z [m]')

ax1[3].plot(z,DeflectionY)
ax1[3].title.set_text(r"Deflection in Y")
ax1[3].set_ylabel(r'v [m]')
ax1[3].set_xlabel(r'z [m]')

plt.show()



#def beamDeflectX(z, Ne,concept,EI,L_arm,TW, macaulay,gc):
 #   pass
  #  Cs = gc(concept,EI,L_arm,TW)
#    # Tofx = 0
#     Mofy = Teng
#     ThetaX = 
#     DeflectionX = 
    

#     return Vofy, Mofy, Thetay, DeflectionY