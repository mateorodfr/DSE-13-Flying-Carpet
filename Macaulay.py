import numpy as np
import Parameters as pm
import matplotlib.pyplot as plt



class Beam(object):
    #Only simmetric crossections
    def __init__(self): 
        self.crossection_Ixx = None
        self.crossection_Iyy = None
        self.crossection_Area = None
        self.crossection_h = None
        self.crossection_w = None
        self.crossection_E = None
        self.crossection_Density = None

#Function Library
def macaulay(z,z0,power):
    return ((z >= z0) * (z>=0))*(z-z0) ** power

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

def BridgeDeflect(z , BF_vec, BM_vec,  Weight_distrib_load,  People_distrib_load , bridge , L_bridge , macaulay):

    b1 = (-1/(bridge.crossection_E * bridge.crossection_Iyy)) * ( 0.5*BF_vec[0] * L_bridge **2 - BM_vec[1]*L_bridge) 
    b2 = (-1/(  bridge.crossection_E * bridge.crossection_Iyy)) * ( (1/6)*BF_vec[0] * L_bridge **3 - 0.5*BM_vec[1]*L_bridge**2)
    b3 = (-1/(bridge.crossection_E * bridge.crossection_Ixx)) * ( - BM_vec[0] * L_bridge - 0.5*BF_vec[1] * L_bridge**2 - np.sum(((Weight_distrib_load + People_distrib_load)*(1/6))* L_bridge**3) ) #Somebody please check this integration i'm not sure if its right but holy shit i cant figure out how to fix this bitch
    b4 = (-1/( bridge.crossection_E * bridge.crossection_Ixx)) * ( - 0.5* BM_vec[0] * L_bridge ** 2 - (1/3)*BF_vec[1] * L_bridge ** 3 - np.sum(((Weight_distrib_load + People_distrib_load)*(1/24))* L_bridge ** 4) ) 
    B = np.array([[b1],[b2],[b3],[b4]])
    A = np.array([[1,0,0,0],[L_bridge,1,0,0],[0,0,1,0],[0,0,L_bridge,1]])
    # print(B)
    Cs = np.matmul(np.linalg.inv(A),B)
    Weight_and_People_load = np.cumsum((People_distrib_load + Weight_distrib_load)*dz)
    Vofy = BF_vec[1]*macaulay(z,0,0) + Weight_and_People_load*macaulay(z,0,0)
    Moment_distrib = [] 
    for i in range(len(z)):
        Moment_distrib.append(np.sum( (People_distrib_load + Weight_distrib_load) * dz * (z[i]-z) * (z < z[i]) ))
    Mofx = - BF_vec[1]*macaulay(z,0,1) - BM_vec[0]*macaulay(z,0,0) - Moment_distrib # Not sure about this one
    Vofx = BF_vec[0]*macaulay(z,0,0) 
    Mofy = BF_vec[0]*macaulay(z,0,1) - BM_vec[1]*macaulay(z,0,0) 
    Nofz = BF_vec[2]*macaulay(z,0,0)
    Tofz = BM_vec[2]*macaulay(z,0,0)
    #### Le fail lies below




    Doubleint =  np.cumsum( dz*np.ones(len(z)) * np.cumsum( Moment_distrib * np.ones(len(z))*dz))
    ThetaX = (-1/(bridge.crossection_E * bridge.crossection_Iyy)) * ( 0.5*BF_vec[0] * z**2 - BM_vec[1] * z)  + Cs[0] 
    DeflectionX = (-1/( bridge.crossection_E * bridge.crossection_Iyy)) * ((1/6) * BF_vec[0] * z**3 - 0.5*BM_vec[1]*z**2)+ Cs[0]*z + Cs[1]
    ##### BEWARE ALL YE ENTER THE NEXT TWO LINES ARE DEBATABLE AT BEST AND DISASTROUS AT WORST
    ThetaY = (-1/(bridge.crossection_E * bridge.crossection_Ixx)) * ( - BM_vec[0] * z - 0.5* BF_vec[1] * z **2 - ((Weight_distrib_load + People_distrib_load)*(1/6))*macaulay(z,0,3)  )  + Cs[2] #FUCKFUCKFUCKFUCKFUCKFUCKFUCKFUCKFUCK
    DeflectionY = (-1/(bridge.crossection_E * bridge.crossection_Ixx)) * (- BM_vec[0] * 0.5 * z **2 - (1/6) * BF_vec[1] * z ** 3 - ((Weight_distrib_load + People_distrib_load)*(1/24))*macaulay(z,0,4)) + Cs[2] * z + Cs[3] #FUCKFUCKFUCKFUCKFUCKFUCKFUCKFUCKFUCK
    return Vofy, Vofx, Mofy, Mofx, Nofz,Tofz, ThetaY, ThetaX, DeflectionY, DeflectionX


def RotorbeamStresses(Vofy, Vofx,  Mofy, Mofx, z, rotorarm):
    tau = (np.sqrt(Vofy**2+Vofx**2) / rotorarm.crossection_Area) * np.ones(len(z))
    sigma = (Mofx*rotorarm.crossection_h)/(2*rotorarm.crossection_Ixx)  + (Mofy*rotorarm.crossection_w)/(2*rotorarm.crossection_Iyy)
    return sigma , tau

def Rotorbeamfreq(rotorarm,Ne,L):
    return (1/(2*np.pi))* np.sqrt((3*rotorarm.crossection_E*rotorarm.crossection_Ixx)/(Ne*L**3))











### CONTROL PANEL ###

initialiserotorarm = False
computerotorarm = False
plotrotorarm = False

initialisebridge = False
computebridge = False
plotbridge = False





### PARAMETERS
concept = pm.ConceptParameters(0)
Teng = concept.motor.Torque
TW = 1.1
Ne = TW*concept.Mtot_concept*concept.physics.g/4

if initialiserotorarm == True: # Putting rotorarm shit in beam class 
    rotorarm = Beam()
    rotorarm.crossection_Ixx = 1 
    rotorarm.crossection_Iyy = 1
    rotorarm.crossection_Area = 1
    rotorarm.crossection_h = 1
    rotorarm.crossection_w = 1
    rotorarm.crossection_E = 69 * 10 ** 6
    rotorarm.crossection_Density = 1
    L_arm = 1
    L = L_arm + concept.propeller.D_prop/2
    dz = 0.01
    z = np.arange(0,L+dz,dz)
else:
    pass

if initialisebridge == True:  # Putting bridge shit in beam class
    bridge = Beam()
    bridge.crossection_Ixx = 1 
    bridge.crossection_Iyy = 1
    bridge.crossection_Area = 1
    bridge.crossection_h = 1
    bridge.crossection_w = 1
    bridge.crossection_E = 69 * 10 ** 6
    bridge.crossection_Density = 1
    L_bridge = 1 #Length of this boinky bridge boi
    dz = 0.01
    z = np.arange(0,L_bridge+dz,dz)
    Weight_distrib_load = np.ones(len(z)) * bridge.crossection_Density * bridge.crossection_Area #Here we make array of unit weights for integration later
    People_distrib_load = np.ones(len(z)) * 10  # This heck of a boi is a placeholder for a distributed loading function we can alter to show the effect of passengers walking on the breedge ASSUMED CONSTANT BCS OTHERWISE ITS FUCKED
    BF_vec = [1,1,1]   # This vector contains the forces of the clamped condition, order - x y z
    BM_vec = [0,0,0]   # This Lad contains the moments of the clamped condition, order x y z ****IMPORTANT, TORQUES DEFINED POSITIVE IN AXIS OF INTERNAL LOADING


else:
    pass



### COMPUTATIONS

if computerotorarm == True:
    Vofy, Vofx, Mofy, Mofx, ThetaY, ThetaX, DeflectionY, DeflectionX = RotorbeamDeflect(z,Ne,concept,rotorarm,L,macaulay)
    sigma, tau = RotorbeamStresses(Vofy, Vofx,  Mofy, Mofx, z, rotorarm)
    Natfreq = Rotorbeamfreq(rotorarm,Ne,L)
else:
    pass


if computebridge == True:
    Vofy, Vofx, Mofy, Mofx, Nofz,Tofz, ThetaY, ThetaX, DeflectionY, DeflectionX = BridgeDeflect(z , BF_vec, BM_vec,  Weight_distrib_load,  People_distrib_load , bridge , L_bridge , macaulay)


else: 
    pass










#Plotting


if plotrotorarm == True:
    print("behold the natural frequency of our beam",(Natfreq))


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

else:
    pass



if plotbridge == True:
    #print("behold the natural frequency of our beam",(Natfreq))

    #fig0,ax0 = plt.subplots(2,1)
    #ax0 = ax0.ravel()
    #fig0.suptitle('Stresses in Beam')

    #ax0[0].plot(z,sigma)
    #ax0[0].title.set_text(r"Shear stress in X")
    #ax0[0].set_ylabel(r'stress [Pa]')
    #ax0[0].set_xlabel(r'z [m]')

    #ax0[1].plot(z,tau)
    #ax0[1].title.set_text(r"Shear force in Y")
    #ax0[1].set_ylabel(r'stress [Pa]')
    #ax0[1].set_xlabel(r'z [m]')


    fig1, ax1 = plt.subplots(2,2)
    ax1 = ax1.ravel()
    fig1.suptitle('Bridge Internal Loadings')

    ax1[0].plot(z,Vofy)
    ax1[0].title.set_text(r"Internal Shear force in Y")
    ax1[0].set_ylabel(r'Vy [N]')
    ax1[0].set_xlabel(r'z [m]')

    ax1[1].plot(z,Mofx)
    ax1[1].title.set_text(r"Internal Moment around X")
    ax1[1].set_ylabel(r'Mz [Nm]')
    ax1[1].set_xlabel(r'z [m]')

    ax1[2].plot(z,Vofx)
    ax1[2].title.set_text(r"Internal shear in X")
    ax1[2].set_ylabel(r'$Vx$ [N]')
    ax1[2].set_xlabel(r'z [m]')

    ax1[3].plot(z,Mofy)
    ax1[3].title.set_text(r"Internal Moment in Y")
    ax1[3].set_ylabel(r'My [Nm]')
    ax1[3].set_xlabel(r'z [m]')


    fig2, ax2 = plt.subplots(2,2)
    ax2 = ax2.ravel()
    fig2.suptitle('Bridge Deflections')

    ax2[0].plot(z,ThetaY)
    ax2[0].title.set_text(r"Deflection angle in Y")
    ax2[0].set_ylabel(r'V [N]')
    ax2[0].set_xlabel(r'z [m]')

    ax2[1].plot(z,ThetaX)
    ax2[1].title.set_text(r"Deflection angle in X")
    ax2[1].set_ylabel(r'M [Nm]')
    ax2[1].set_xlabel(r'z [m]')

    ax2[2].plot(z,DeflectionY)
    ax2[2].title.set_text(r"Deflection in Y")
    ax2[2].set_ylabel(r'$\alpha$ [rad]')
    ax2[2].set_xlabel(r'z [m]')

    ax2[3].plot(z,DeflectionX)
    ax2[3].title.set_text(r"Deflection in X")
    ax2[3].set_ylabel(r'v [m]')
    ax2[3].set_xlabel(r'z [m]')

    
    plt.show()

else:
    pass


"""
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
"""
