import numpy as np
import RotorDesign

AoA = np.load(r"airfoils\naca0012aero\AoA.npy")
Cd  = np.load(r"airfoils\naca0012aero\cd.npy")
Cl  = np.load(r"airfoils\naca0012aero\cl.npy")
Cm  = np.load(r"airfoils\naca0012aero\cm.npy")
Cp  = np.load(r"airfoils\naca0012aero\cp.npy")

ClA = np.vstack((Cl, AoA)).T
CdA = np.vstack((Cd, AoA)).T
CmA = np.vstack((Cm, AoA)).T
CpA = np.vstack((Cp, AoA)).T

class Blade:
    def __init__(self, length: float, Nslice: int = 1):
        self.b = length
        self.Nslice = self.setSliceNum(Nslice)
    
    def setSliceNum(self, Nslice: int) -> None:
        self.Nslice = Nslice
        self.genSlices()
    
    def genSlices(self):
        self.Slices = [BladeSlice(self.b/N, (self.b/N)*(N-.5), ) for N in range(self.Nslice)]


class BladeSlice(object):
    def __init__(self, local_dr: float, local_r: float, local_pitch: float = 0):
        self.r     = local_r
        self.dr    = local_dr
        self.chord = 1
        self.pitch = local_pitch

    def getLocalVel(self, ang_vel: float) -> float:
        '''
        v = omega * r
        '''
        return ang_vel*self.r

    def getLocalAoA(self) -> float:
        return ...
    
    def getLocalCl(self) -> float:
        '''
        Cl depends on the angle of attack
        '''
        closestAoAindex = np.abs(AoA - self.getLocalAoA()).argmin()
        return ClA[closestAoAindex]
    
    def getLocalArea(self) -> float:
        return self.dr * self.chord

    def getLocalLift(self) -> float:
        '''
        L = Cl * 1/2 * rho * V * V * S
        '''
        Cl  = self.getLocalCl(...)
        rho = 1.225 # HARD CODED, has to be addapted for altitude
        V   = self.getLocalVel(ang_vel)
        S   = self.getLocalArea()
        
        return Cl * .5 * rho * V * V * S
