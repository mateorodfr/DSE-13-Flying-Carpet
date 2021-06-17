"""Small script to determine spacing of stringers on the flanges to prevent shear buckling"""
import numpy as np
import matplotlib.pyplot as plt

class ShearBuckling:

    def __init__(self, stress, length, Ks, Emod, thickness):
        self.stress = stress
        self.length = length
        self.Ks = Ks
        self.Emod = Emod
        self.thickness = thickness
        self.minimal_pitch = self.find_minimal_pitch()
        self.n_stringers = self.find_required_stringers()

    def find_minimal_pitch(self):
        return self.thickness / np.sqrt(self.stress/(self.Ks*self.Emod))

    def find_required_stringers(self):
        return np.ceil(self.length / self.minimal_pitch)

    def plot_findings(self):
        fig, ax = plt.subplots(figsize=(8, 6))

        ax.scatter(self.stress / 1e6, self.minimal_pitch, marker="x", color="red")
        ax.grid()
        ax.set_xlabel("Shear stress [MPa]")
        ax.set_ylabel("Required stringer spacing [m]")

        plt.show()

def main():
    stresses = np.linspace(15e6, 25e6, 11)
    buckling_analysis = ShearBuckling(stresses, 2, 6.4, 70e9, 1e-3)
    print(f"Current stress value: \t\t  {stresses/1e6}")
    print(f"Required number of stringers: {buckling_analysis.n_stringers}")
    print(f"Minimum stringer spacing: \t {buckling_analysis.minimal_pitch}")
    buckling_analysis.plot_findings()


if __name__ == "__main__":
    main()
