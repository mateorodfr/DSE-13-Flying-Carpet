import numpy as np
# import RotorDesign


def main():
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
            self.b       = length
            self.Nslice  = self.setSliceNum(Nslice)
        
        def setSliceNum(self, Nslice: int) -> None:
            self.Nslice = Nslice
            self.genSlices()
        
        def genSlices(self):
            self.Slices = [BladeSlice(self.b/N, (self.b/N)*(N-.5), 20) for N in range(1, self.Nslice+1)]


    class BladeSlice(Blade):
        def __init__(self, local_dr: float, local_r: float, local_pitch: float = 0):
            self.r     = local_r
            self.dr    = local_dr
            self.chord = .25
            self.pitch = local_pitch

            self.ang_vel = 31.666666666666668 # HARD CODED, has to be adapted function

            self.dL    = self.getLocalCl()

        def getLocalVel(self) -> float:
            """
            v = omega * r
            """
            return self.ang_vel*self.r

        def getLocalAoA(self) -> float:
            return self.pitch # Assumption for now. To-Do: calc/estimate actual AoA
        
        def getLocalCl(self) -> float:
            """
            Cl depends on the angle of attack
            """
            closestAoAindex = np.abs(AoA - self.getLocalAoA()).argmin()
            return ClA[closestAoAindex][0]
        
        def getLocalArea(self) -> float:
            return self.dr * self.chord

        def getLocalLift(self) -> float:
            """
            L = Cl * 1/2 * rho * V * V * S
            """
            Cl  = self.getLocalCl()
            rho = 1.225 # HARD CODED, has to be adapted for altitude
            V   = self.getLocalVel()
            S   = self.getLocalArea()
            
            return Cl * .5 * rho * V * V * S
    
    ROTOR_BLADE = Blade(2.5, 5)
    NN = 1
    print(f"Local Velocity of slice {NN} equals \t {ROTOR_BLADE.Slices[NN].getLocalVel()} [m/s]"
         +f"\nLocal radius of slice {NN} equals \t {ROTOR_BLADE.Slices[NN].r} [m]"
         +f"\nLocal Thrust of slice {NN} equals \t {ROTOR_BLADE.Slices[NN].getLocalLift()} [N]\n\n")












if __name__ == "__main__":
    main()