import numpy as np
# import RotorDesign


def main():
    AoA = np.load(r"airfoils\naca633418aero\AoA.npy")
    Cd  = np.load(r"airfoils\naca633418aero\cd.npy")
    Cl  = np.load(r"airfoils\naca633418aero\cl.npy")
    Cm  = np.load(r"airfoils\naca633418aero\cm.npy")
    Cp  = np.load(r"airfoils\naca633418aero\cp.npy")

    ClA = np.vstack((Cl, AoA)).T
    CdA = np.vstack((Cd, AoA)).T
    CmA = np.vstack((Cm, AoA)).T
    CpA = np.vstack((Cp, AoA)).T

    class Blade:
        def __init__(self, span: float, Nslice: int = 1):
            self.span    = span
            self.Nslice  = self.setSliceNum(Nslice)
            self.Lift    = self.getLift()
        
        def setSliceNum(self, Nslice: int) -> None:
            self.Nslice = Nslice
            self.genSlices()
        
        def genSlices(self):
            self.Slices = [BladeSlice(self.span/self.Nslice, (self.span/self.Nslice)*(N-.5), 13) for N in range(1, self.Nslice+1)]

        def getLift(self) -> float:
            return np.sum([sl.getLocalLift() for sl in self.Slices])


    class BladeSlice(Blade):
        def __init__(self, local_dr: float, local_r: float, local_pitch: float = 0):
            self.r     = local_r
            self.dr    = local_dr
            self.chord = .25*.5
            self.pitch = local_pitch

            self.ang_vel = 87.5 # HARD CODED, has to be adapted function

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
    
    ROTOR_BLADE = Blade(2.5/2, 10)
    NN = 0
    print(f"\n\n\nLocal v of slice {NN} equals \t {ROTOR_BLADE.Slices[NN].getLocalVel()} [m/s]"
         +f"\nLocal Cl of slice {NN} equals \t {ROTOR_BLADE.Slices[NN].getLocalCl()} [-]"
         +f"\nLocal r of slice {NN} equals \t {ROTOR_BLADE.Slices[NN].r} [m]"
         +f"\nLocal T of slice {NN} equals \t {ROTOR_BLADE.Slices[NN].getLocalLift()} [N]")

             
    print(f"\nMax Weight equals {ROTOR_BLADE.Lift*8*4/9.80665} [kg]\n\n\n")





if __name__ == "__main__":
    main()