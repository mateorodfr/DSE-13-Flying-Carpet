import numpy as np
import matplotlib.pyplot as plt
import Parameters as pm



"General craft geometry library"

L_rotor = 2 # Rotor arm length
L_bridge = 2 # Bridge length
w_cabin = 1 # Cabin width
d_cabin = 1 # Cabin depth
h_cabin = 1 # Cabin height
alpha_1 = np.pi /3 #These are propeller rod angles, all defined clockwise from x axis on top of craft
alpha_2 = 2 * np.pi /3
alpha_3 = 4 * np.pi /3
alpha_4 = 5 *  np.pi /3
Beta = np.pi /6 # This is bridge angle going up. This is relative to the negative x axis

Ixx_init = 1
Iyy_init = 1
Izz_init = 1
mass = 1

"Desired accellerations"

angular_acc_x = 1
angular_acc_y = 1
angular_acc_z = 1
acc_z = 1

"Array Library"

#Unit vectors of forces
L_unit = np.array([0,0,-1])
N_unit = np.array([0,0,-1])
B_unit = np.array([0,0,-1])
W_unit = np.array([0.37,0.61,0.71])
#Position vectors of forces
L_1_pos = np.array([d_cabin/2 + L_rotor*np.cos(alpha_1), w_cabin/2 + L_rotor*np.sin(alpha_1),0])
L_2_pos = np.array([-d_cabin/2 + L_rotor*np.cos(alpha_2), w_cabin/2 + L_rotor*np.sin(alpha_2),0])
L_3_pos = np.array([-d_cabin/2 + L_rotor*np.cos(alpha_3), -w_cabin/2 + L_rotor*np.sin(alpha_3),0])
L_4_pos = np.array([d_cabin/2 + L_rotor*np.cos(alpha_4), -w_cabin/2 + L_rotor*np.sin(alpha_4),0])
N_1_pos = np.array([d_cabin/2,w_cabin/2,-h_cabin])
N_2_pos = np.array([-d_cabin/2,w_cabin/2,-h_cabin])
N_3_pos = np.array([-d_cabin/2,-w_cabin/2,-h_cabin])
N_4_pos = np.array([d_cabin/2,-w_cabin/2,-h_cabin])
c_g_pos = np.array([0,0,1])
B_pos = np.array([d_cabin/2 + L_bridge*np.cos(Beta),0,-h_cabin-L_bridge*np.sin(Beta)])

"Scalar Library"

B = 1 # Bridge force
W = 1 # Weight
N = np.array([1,1,1,1]) #Normal force for each leg, defined clockwise starting at front right (looking at x axis, it's the same as for the rotors)

"Here we want to solve for thrust and accellerations"

#def DynamicAnalysis(        type):

#A_1 = np.cross(N_unit,N_1_pos).append(1)
#A_2 = np.cross(N_unit,N_2_pos).append(1)
#A_3 = np.cross(N_unit,N_3_pos).append(1)
#A_4 = np.cross(N_unit,N_4_pos).append(1)

Ttot = acc_z*mass

A = np.vstack((np.append(np.cross(L_unit,L_1_pos),1),np.append(np.cross(L_unit,L_2_pos), 1),np.append(np.cross(L_unit,L_3_pos),1),np.append(np.cross(L_unit,L_4_pos),1)))

A = np.matrix.T

print((np.append(np.cross(L_unit,L_1_pos), 1))

# B = np.matrix([])

   # L = np.linalg.inv(A).dot(B)

    #return 

# Fz : -m*alpha_z = mg*cos(phi)*cos(theta) - (L1+L2+L3+L4) + (1 if boarding else 0) * [(N21+N22+N23+N24) + (Bx)]  # This one i'd like to get total T from
#"If possible with the equation above, I would use it to get Ttot from"

# Fy :  m*alpha_y = w *cos(phi)*sin(theta) - (1 if boarding else 0) * (By + sum_of(Ny))  #The problem with this one and the one below is that once both bridge and 
# Fx :  m*alpha_x = w *sin(phi)*sin(theta) - (1 if boarding else 0) * (Bx + sum_of(Ny))

# Mz : I_zz*alpha_z = sum_of(T_eng) + (1 if boarding else 0) * [sum_of(+/- w/2 * Nx) + sum_of(+/- d/2 * Ny) - By*(d/2 + L_B * cos(beta))]


# W






#"Here all of the 5 main assemblies are defined, for each some properties are given to use in the analysis"

class Rotorarm(object):

    " Material_Properties = [ sigma , tau , E-modulus, G-modulus] "
    " Crosssection_Properties = [A, Ixx, Iyy, xmax, ymax] "
    " ONLY SYMMETRIC CROSS SECTIONS "


    Material_Properties0 = []
    Crosssection_Properties0 = []
    Material_Properties1 = []
    Crosssection_Properties1 = []


    def __init__(self,keyC=0,keyM=0):
        if keyC == 0:
            self.CrosssectionArea = self.Crosssection_Properties0[0]
            self.Ixx = self.Crosssection_Properties0[1]
            self.Iyy = self.Crosssection_Properties0[2]
            self.ymax = self.Crosssection_Properties0[3]
            self.zmax = self.Crosssection_Properties0[4]

        if keyM == 0:
            self.sigma = self.Material_Properties0[0]
            self.tau = self.Material_Properties0[1]
            self.Emodulus = self.Material_Properties0[2]
            self.Gmodulus = self.Material_Properties0[3] 


        elif keyC == 1:
            self.CrosssectionArea = self.Crosssection_Properties1[0]
            self.Ixx = self.Crosssection_Properties1[1]
            self.Iyy = self.Crosssection_Properties1[2]
            self.ymax = self.Crosssection_Properties1[3]
            self.zmax = self.Crosssection_Properties1[4]

        elif keyM == 1:
            self.sigma = self.Material_Properties1[0]
            self.tau = self.Material_Properties1[1]
            self.Emodulus = self.Material_Properties1[2]
            self.Gmodulus = self.Material_Properties1[3] 

class Bridge(object):

    " Material_Properties = [ sigma , tau , E-modulus, G-modulus] "
    " Crosssection_Properties = [A, Ixx, Iyy, ymax, zmax] "


    Material_Properties0 = []
    Sectionbridge_Properties0 = []
    Material_Properties1 = []
    Sectionbridge_Properties1 = []


    def __init__(self,keyC=0,keyM=0):
        if keyS == 0:
            self.CrosssectionArea = self.Crosssection_Properties0[0]
            self.Ixx = self.Sectionbridge_Properties0[1]
            self.Iyy = self.Sectionbridge_Properties0[2]
            self.ymax = self.Sectionbridge_Properties0[3]
            self.zmax = self.Sectionbridge_Properties0[4]

        if keyM == 0:
            self.sigma = self.Material_Properties0[0]
            self.tau = self.Material_Properties0[1]
            self.Emodulus = self.Material_Properties0[2]
            self.Gmodulus = self.Material_Properties0[3] 


        elif keyS == 1:
            self.CrosssectionArea = self.Crosssection_Properties1[0]
            self.Ixx = self.Sectionbridge_Properties1[1]
            self.Iyy = self.Sectionbridge_Properties1[2]
            self.ymax = self.Sectionbridge_Properties1[3]
            self.zmax = self.Sectionbridge_Properties1[4]

        elif keyM == 1:
            self.sigma = self.Material_Properties1[0]
            self.tau = self.Material_Properties1[1]
            self.Emodulus = self.Material_Properties1[2]
            self.Gmodulus = self.Material_Properties1[3] 


class Cabin(object):

    " Material_Properties = [ sigma , tau , E-modulus, G-modulus] "
    " Crosssection_Properties = [A, Ixx, Iyy, ymax, zmax] "


    Material_Properties0 = []
    Cabin_Properties0 = []
    Material_Properties1 = []
    Cabin_Properties1 = []


    def __init__(self,keyC=0,keyM=0):
        if keyL == 0:
            self.CrosssectionArea = self.Crosssection_Properties0[0]
            self.Ixx = self.Cabin_Properties0[1]
            self.Iyy = self.Cabin_Properties0[2]
            self.ymax = self.Cabin_Properties0[3]
            self.zmax = self.Cabin_Properties0[4]

        if keyM == 0:
            self.sigma = self.Material_Properties0[0]
            self.tau = self.Material_Properties0[1]
            self.Emodulus = self.Material_Properties0[2]
            self.Gmodulus = self.Material_Properties0[3] 


        elif keyL == 1:
            self.CrosssectionArea = self.Crosssection_Properties1[0]
            self.Ixx = self.Cabin_Properties1[1]
            self.Iyy = self.Cabin_Properties1[2]
            self.ymax = self.Cabin_Properties1[3]
            self.zmax = self.Cabin_Properties1[4]

        elif keyM == 1:
            self.sigma = self.Material_Properties1[0]
            self.tau = self.Material_Properties1[1]
            self.Emodulus = self.Material_Properties1[2]
            self.Gmodulus = self.Material_Properties1[3] 

class LandingGear(object):

    " Material_Properties = [ sigma , tau , E-modulus, G-modulus] "
    " Crosssection_Properties = [A, Ixx, Iyy, ymax, zmax] "


    Material_Properties0 = []
    LandingGear_Properties0 = []
    Material_Properties1 = []
    LandingGear_Properties1 = []


    def __init__(self,keyC=0,keyM=0):
        if keyL == 0:
            self.CrosssectionArea = self.Crosssection_Properties0[0]
            self.Ixx = self.LandingGear_Properties0[1]
            self.Iyy = self.LandingGear_Properties0[2]
            self.ymax = self.LandingGear_Properties0[3]
            self.zmax = self.LandingGear_Properties0[4]

        if keyM == 0:
            self.sigma = self.Material_Properties0[0]
            self.tau = self.Material_Properties0[1]
            self.Emodulus = self.Material_Properties0[2]
            self.Gmodulus = self.Material_Properties0[3] 


        elif keyL == 1:
            self.CrosssectionArea = self.LandingGear_Properties1[0]
            self.Ixx = self.LandinGgear_Properties1[1]
            self.Iyy = self.LandingGear_Properties1[2]
            self.ymax = self.LandingGear_Properties1[3]
            self.zmax = self.LandingGear_Properties1[4]

        elif keyM == 1:
            self.sigma = self.Material_Properties1[0]
            self.tau = self.Material_Properties1[1]
            self.Emodulus = self.Material_Properties1[2]
            self.Gmodulus = self.Material_Properties1[3] 

class TopPart(object):

    " Material_Properties = [ sigma , tau , E-modulus, G-modulus] "
    " Crosssection_Properties = [A, Ixx, Iyy, ymax, zmax] "


    Material_Properties0 = []
    TopPart_Properties0 = []
    Material_Properties1 = []
    TopPart_Properties1 = []


    def __init__(self,keyC=0,keyM=0):
        if keyL == 0:
            self.CrosssectionArea = self.Crosssection_Properties0[0]
            self.Ixx = self.TopPart_Properties0[1]
            self.Iyy = self.TopPart_Properties0[2]
            self.ymax = self.TopPart_Properties0[3]
            self.zmax = self.TopPart_Properties0[4]

        if keyM == 0:
            self.sigma = self.Material_Properties0[0]
            self.tau = self.Material_Properties0[1]
            self.Emodulus = self.Material_Properties0[2]
            self.Gmodulus = self.Material_Properties0[3] 


        elif keyL == 1:
            self.CrosssectionArea = self.LandingGear_Properties1[0]
            self.Ixx = self.TopPart_Properties1[1]
            self.Iyy = self.TopPart_Properties1[2]
            self.ymax = self.TopPart_Properties1[3]
            self.zmax = self.TopPart_Properties1[4]

        elif keyM == 1:
            self.sigma = self.Material_Properties1[0]
            self.tau = self.Material_Properties1[1]
            self.Emodulus = self.Material_Properties1[2]
            self.Gmodulus = self.Material_Properties1[3] 







def internal_loading_beam(FSX,FSY,FSZ,MSX,MSY,MSZ,zz,plot=False):
    "ONLY WORKS FOR NON THIN WALLED BEAMS CLAMPED AT ONE SIDE"
    "Force-> F (S/E) (X/Y/Z), Moment -> M (S/E) (X/Y/Z), this means forces & moments are defined by S or E for start or end of the beam, and XYZ for the direction)"
    "This analysis and loading diagram go from unclamped edge to clamped edge"
    print(FSZ)
    NormalZ = FSZ
    ShearY = FSY
    ShearX = FSX
    MomentX = MSX + FSY * zz
    MomentY = MSY - FSX * zz
    TorqueZ = MSZ

    if plot == True:
        plt.subplot(231)
        plt.ylim(-np.max(NormalZ),np.max(NormalZ))
        #plt.xlim(-zz[len(zz)],zz[len(zz)])
        plt.plot(zz,NormalZ)
        plt.title("NormalZ")
        plt.grid()
        plt.subplot(232)
        plt.ylim(-np.max(ShearY),np.max(ShearY))
        #plt.xlim(-zz[len(zz)],zz[len(zz)])
        plt.plot(zz,ShearY)
        plt.title("ShearY")
        plt.grid()
        plt.subplot(233)
        plt.ylim(-np.max(ShearX),np.max(ShearX))
        #plt.xlim(-zz[len(zz)],zz[len(zz)])
        plt.plot(zz,ShearX)
        plt.title("ShearX")
        plt.grid()
        plt.subplot(234)
        plt.ylim(-np.max(MomentX),np.max(MomentX))
        #plt.xlim(-zz[len(zz)],zz[len(zz)])
        plt.plot(zz,MomentX)
        plt.title("MomentX")
        plt.grid()
        plt.subplot(235)
        plt.ylim(-np.max(MomentY),np.max(MomentY))
        #plt.xlim(-zz[len(zz)],zz[len(zz)])
        plt.plot(zz,MomentY)
        plt.title("MomentY")
        plt.grid()
        plt.subplot(236)
        plt.ylim(-np.max(TorqueZ),np.max(TorqueZ))
        #plt.xlim(-zz[len(zz)],zz[len(zz)])
        plt.plot(zz,TorqueZ)
        plt.title("TorqueZ")
        plt.grid()
        plt.show()
        
        return zz, NormalZ, ShearX, ShearY, MomentX, MomentY, TorqueZ
    elif plot == False:
        return zz, NormalZ, ShearX, ShearY, MomentX, MomentY, TorqueZ



def internal_stresses_deflection_beam(zz,NormalZ,ShearY,ShearX,MomentX,MomentY,TorqueZ,A,Ixx,Iyy,E,G,L,w,h):
    
    "Here we check every point of the internal loading for stresses & deflections"
    
    sigma_Z = ((MomentX)* h) / (2 * Ixx) + (MomentY * w) / (2 * Iyy) + NormalZ / A
    tau_arr_X = ShearX / A
    tau_arr_Y = ShearY / A
    deflection_X = (-1/(E * I)) * MomentY*zz*zz
    deflection_Y = (-1/(E * I)) * MomentX*zz*zz


    if plot == True:
            plt.subplot(231)
            plt.ylim(-np.max(NormalZ),np.max(NormalZ))
            #plt.xlim(-zz[len(zz-1)],zz[len(zz-1)])
            plt.plot(zz,NormalZ)
            plt.title("NormalZ")
            plt.grid()
            plt.subplot(232)
            plt.ylim(-np.max(ShearY),np.max(ShearY))
            #plt.xlim(-zz[len(zz-1)],zz[len(zz-1)])
            plt.plot(zz,ShearY)
            plt.title("ShearY")
            plt.grid()
            plt.subplot(233)
            plt.ylim(-np.max(ShearX),np.max(ShearX))
            #plt.xlim(-zz[len(zz-1)],zz[len(zz-1)])
            plt.plot(zz,ShearX)
            plt.title("ShearX")
            plt.grid()
            plt.subplot(234)
            plt.ylim(-np.max(MomentX),np.max(MomentX))
            #plt.xlim(-zz[len(zz-1)],zz[len(zz-1)])
            plt.plot(zz,MomentX)
            plt.title("MomentX")
            plt.grid()
            plt.subplot(235)
            plt.ylim(-np.max(MomentY),np.max(MomentY))
            #plt.xlim(-zz[len(zz-1)],zz[len(zz-1)])
            plt.plot(zz,MomentY)
            plt.title("MomentY")
            plt.grid()
            plt.subplot(236)
            plt.ylim(-np.max(TorqueZ),np.max(TorqueZ))
            #plt.xlim(-zz[len(zz-1)],zz[len(zz-1)])
            plt.plot(zz,TorqueZ)
            plt.title("TorqueZ")
            plt.grid()
            plt.show()
    
    "Not sure if the next one's right, note to self"
    #deflection_Y = (1/6) * (ShearY / (E * Iyy)) * zz ** 3 + tau_arr_Y / G

    return(maxsigma,maxtau,maxVX,maxVY)



#### INPUT PARAMETERS ####


L = 1
P1 = 1
Tengine = 1



##### INDIVIDUAL LOAD FUNCTION ARRAYS #########
zz = np.linspace(0,L,1000,True)

def Flatload(F,zz):
    return F + 0 * zz
def Linearload(F,zz):
    return F * zz
def Zeroload(zz):
    return np.zeros(len(zz))



#### LOADING CASES ####

#Top rotor loading cases

E1 = [Zeroload(zz),-Flatload(P1,zz),Zeroload(zz),Zeroload(zz),Zeroload(zz),-Flatload(Tengine,zz)] #One engine inoperative
E2 = [Zeroload(zz),-Flatload(P1,zz),Zeroload(zz),Zeroload(zz),Zeroload(zz),Zeroload(zz)] #Normal operation
E3 = [100,-P1,0,0,0,-Tengine] # Rotorarm stuck to something, induces normal force
E4 = [0,-P1+ 100 ,0,0,0,-Tengine] # Something is hanging on the arm (debris)


#### INPUTTING AND TESTING ####

Ein = E1
[] , [] , [] , [] , [] , [], [] = internal_loading_beam(Ein[0],Ein[1],Ein[2],Ein[3],Ein[4],Ein[5],zz,True)[0] , internal_loading_beam(Ein[0],Ein[1],Ein[2],Ein[3],Ein[4],Ein[5],zz,True)[0] , internal_loading_beam(Ein[0],Ein[1],Ein[2],Ein[3],Ein[4],Ein[5],zz,True)[0]
[] = internal_stresses_deflection_beam() 


