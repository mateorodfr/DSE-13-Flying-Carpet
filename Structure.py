import numpy as np
import matplotlib.pyplot as plt






class Rotorarm(object):

    " Material_Properties = [ sigma , tau , E-modulus, G-modulus] "
    " Crosssection_Properties = [A, Ixx, Iyy, ymax, zmax] "


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


def internal_loading_beam(L,FSX,FSY,FSZ,MSX,MSY,MSZ,plot=False):
    "ONLY WORKS FOR NON THIN WALLED BEAMS CLAMPED AT ONE SIDE"
    "Force-> F (S/E) (X/Y/Z), Moment -> M (S/E) (X/Y/Z), this means forces & moments are defined by S or E for start or end of the beam, and XYZ for the direction)"
    "This analysis and loading diagram go from unclamped edge to clamped edge"
    zz = np.linspace(0,L,1000,True)
    NormalZ = - FSZ(zz)
    ShearY = FSY(zz)
    ShearX = FSX(zz)
    MomentX = MSX(zz)
    MomentY = MSY(zz)
    TorqueZ = MSZ(zz)

    if plot:
        plt.subplot(231)
        plt.plot(zz,NormalZ)
        plt.title("NormalZ")
        plt.subplot(232)
        plt.plot(zz,ShearY)
        plt.title("ShearY")
        plt.subplot(233)
        plt.plot(zz,ShearX)
        plt.title("ShearX")
        plt.subplot(234)
        plt.plot(zz,MomentX)
        plt.title("MomentX")
        plt.subplot(235)
        plt.plot(zz,MomentY)
        plt.title("MomentY")
        plt.subplot(236)
        plt.plot(zz,TorqueZ)
        plt.title("TorqueZ")
        plt.show()
        
    return zz, NormalZ, ShearX, ShearY, MomentX, MomentY, TorqueZ




def internal_stresses_&_deflection_beam(zz,NormalZ,ShearY,ShearX,MomentX,MomentY,TorqueZ,A,Ixx,Iyy,E,G):
    "Here we check every point of the internal loading for stresses & deflections"
    
    





    return(maxsigma,maxtau,maxVX,maxVY)






