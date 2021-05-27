import numpy as np

AoA = np.load(r"C:\Users\marvd\Documents\GitHub\DSE-13-Flying-Carpet\airfoils\naca2412aero\AoA.npy")
Cd  = np.load(r"C:\Users\marvd\Documents\GitHub\DSE-13-Flying-Carpet\airfoils\naca2412aero\cd.npy")
Cl  = np.load(r"C:\Users\marvd\Documents\GitHub\DSE-13-Flying-Carpet\airfoils\naca2412aero\cl.npy")
Cm  = np.load(r"C:\Users\marvd\Documents\GitHub\DSE-13-Flying-Carpet\airfoils\naca2412aero\cm.npy")
Cp  = np.load(r"C:\Users\marvd\Documents\GitHub\DSE-13-Flying-Carpet\airfoils\naca2412aero\cp.npy")

class Blade:
    def __init__(self, length: float):
        self.length = length
        self.chord  = self.length * .1
        self.t      = self.chord  * .1

class BladeSlice:
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
