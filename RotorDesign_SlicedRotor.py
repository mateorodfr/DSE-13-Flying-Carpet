import numpy as np

AoA = np.load(r"airfoils\naca0012aero\AoA.npy")
Cd  = np.load(r"airfoils\naca0012aero\cd.npy")
Cl  = np.load(r"airfoils\naca0012aero\cl.npy")
Cm  = np.load(r"airfoils\naca0012aero\cm.npy")
Cp  = np.load(r"airfoils\naca0012aero\cp.npy")

ClA = np.vstack((Cl, AoA)).T
print(Cl[0])

class Blade:
    def __init__(self, length: float, Nslice: int = 1):
        self.length = length
        self.Nslice = Nslice
    
    def setSliceNum(self, Nslice: int):
        self.Nslice = Nslice

class BladeSlice(object):
    def __init__(self, local_r: float, local_dr: float):
        self.r  = local_r
        self.dr = local_dr

    def getLocalVel(self, ang_vel: float) -> float:
        '''
        v = omega * r
        '''
        return ang_vel*self.r
    
    def getLocalCl(self, AoA) -> float:
        '''
        Cl depends on the angle of attack
        '''
        return ...
    def getLocalArea(self) -> float:
        return self.dr * ...

    def getLocalLift(self) -> float:
        '''
        L = Cl * 1/2 * rho * V * V * S
        '''
        Cl  = self.getLocalCl(...)
        rho = 1.225 # HARD CODED, has to be addapted for altitude
        V   = self.getLocalVel(ang_vel)
        S   = self.getLocalArea()
        
        return Cl * .5 * rho * V * V * S
