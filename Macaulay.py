import numpy as np
import Parameters as pm
import matplotlib.pyplot as plt



class Beam(object):

    " Material_Properties = [ sigma , tau , E-modulus, G-modulus] "
    " Crosssection_Properties = [A, Ixx, Iyy, xmax, ymax] "
    " ONLY SYMMETRIC CROSS SECTIONS "
    def __init__(self): 
        self.crossection_Ixx = 1
        self.crossection_Iyy = 1
        self.crossection_Area = 1
        self.crossection_h = 1
        self.crossection_w = 1
        self.crossection_E = 69 * 10 ** 6

#Function Library
def macaulay(z,z0,power):
    return ((z > z0) * (z>0))*(z-z0) ** power

def RotorbeamDeflect(z, Ne, concept, rotorarm,L, macaulay): # PLEASE NOTE that this analysis was conducted from the point of the PROPELLER, so the beam starts before the beam starts (if you will)
    Ze = 0 #Point of application of force, in this case at the base
    Te = 0 #Point of application of moment, in this case at the base
    EI = np.array([rotorarm.crossection_Ixx * rotorarm.crossection_E , rotorarm.crossection_Iyy * rotorarm.crossection_E]) #Making array of MMOI and E
    b1 = (1/(2*EI[0]))*Ne*macaulay(L,Ze,2)
    b2 = (1/(6*EI[0]))*Ne*macaulay(L,Ze,3)
    k1 = (1/(EI[1]))*concept.motor.Torque*macaulay(L,Ze,1)
    k2 = (1/(2*EI[1]))*concept.motor.Torque*macaulay(L,Ze,2)
    B = np.array([[b1],[b2]]) #Building B matrices
    K = np.array([[k1],[k2]])
    A = np.array([[1,0],[L,1]]) # A matrix is the same for both so we only make one
    Cs = np.matmul(np.linalg.inv(A),B) # Solving constants of integration
    Cs2 = np.matmul(np.linalg.inv(A),K)
    Vofy = Ne*macaulay(z,0,0) #Solving Shears & Moments
    Mofx = Ne * macaulay(z,0,1)
    Vofx = np.zeros(len(z))
    Mofy = concept.motor.Torque*np.ones(len(z))
    ThetaY = (-1/EI[0])*(0.5*Ne*macaulay(z,Ze,2))+Cs[0] #Solving angles of deflection and deflection
    DeflectionY = (-1/EI[0]) * ((1/6)*Ne*macaulay(z,Ze,3)) + (Cs[0]*z + Cs[1])
    ThetaX = (-1/EI[1])*(concept.motor.Torque*macaulay(z,Te,1)) + Cs2[0]
    DeflectionX = (-1/EI[1])*(0.5*concept.motor.Torque*macaulay(z,Te,2)) + (Cs2[0]*z + Cs2[1])
    return Vofy, Vofx, Mofy, Mofx, ThetaY, ThetaX, DeflectionY, DeflectionX

def RotorbeamStresses(Vofy, Vofx,  Mofy, Mofx, z, rotorarm):
    tau = (np.sqrt(Vofy**2+Vofx**2) / rotorarm.crossection_Area) * np.ones(len(z))
    sigma = (Mofx*rotorarm.crossection_h)/(2*rotorarm.crossection_Ixx)  + (Mofy*rotorarm.crossection_w)/(2*rotorarm.crossection_Iyy)
    return sigma , tau

def Rotorbeamfreq(rotorarm,Ne,L):
    return (1/(2*np.pi))* np.sqrt((3*rotorarm.crossection_E*rotorarm.crossection_Ixx)/(Ne*L**3))

### Parameters and values

concept = pm.ConceptParameters(0)
rotorarm = Beam()
L_arm = 1
L = L_arm + concept.propeller.D_prop/2
TW = 1
dz = 0.01
Ne = TW*concept.Mtot_concept*concept.physics.g/4
Teng = concept.motor.Torque
z = np.arange(0,L+dz,dz)

Vofy, Vofx, Mofy, Mofx, ThetaY, ThetaX, DeflectionY, DeflectionX = RotorbeamDeflect(z,Ne,concept,rotorarm,L,macaulay)
sigma, tau = RotorbeamStresses(Vofy, Vofx,  Mofy, Mofx, z, rotorarm)
Natfreq = Rotorbeamfreq(rotorarm,Ne,L)

print(Natfreq)
















#Plotting

fig0,ax0 = plt.subplots(2,1)
ax0 = ax0.ravel()
fig0.suptitle('Stresses in Beam')

ax0[0].plot(z,sigma)
ax0[0].title.set_text(r"Shear stress in X")
ax0[0].set_ylabel(r'stress [Pa]')
ax0[0].set_xlabel(r'z [m]')

ax0[1].plot(z,tau)
ax0[1].title.set_text(r"Shear force in Y")
ax0[1].set_ylabel(r'stress [Pa]')
ax0[1].set_xlabel(r'z [m]')


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

ax1[2].plot(z,ThetaY)
ax1[2].title.set_text(r"Deflection angle in Y")
ax1[2].set_ylabel(r'$\alpha$ [rad]')
ax1[2].set_xlabel(r'z [m]')

ax1[3].plot(z,DeflectionY)
ax1[3].title.set_text(r"Deflection in Y")
ax1[3].set_ylabel(r'v [m]')
ax1[3].set_xlabel(r'z [m]')


fig2, ax2 = plt.subplots(2,2)
ax2 = ax2.ravel()
fig2.suptitle('Beam loading in X direction')

ax2[0].plot(z,Vofx)
ax2[0].title.set_text(r"Shear force in X")
ax2[0].set_ylabel(r'V [N]')
ax2[0].set_xlabel(r'z [m]')

ax2[1].plot(z,Mofy)
ax2[1].title.set_text(r"Moment around Y")
ax2[1].set_ylabel(r'M [Nm]')
ax2[1].set_xlabel(r'z [m]')

ax2[2].plot(z,ThetaX)
ax2[2].title.set_text(r"Deflection angle in X")
ax2[2].set_ylabel(r'$\alpha$ [rad]')
ax2[2].set_xlabel(r'z [m]')

ax2[3].plot(z,DeflectionX)
ax2[3].title.set_text(r"Deflection in X")
ax2[3].set_ylabel(r'v [m]')
ax2[3].set_xlabel(r'z [m]')



plt.show()



