from xfoil import XFoil
from xfoil.model import Airfoil
import numpy as np
import matplotlib.pyplot as plt
import os

xf = XFoil()

def read_airfoil(filename):
    """Opens a file in the format of two columns and creates Airfoil object"""
    data = np.genfromtxt(filename, skip_header=1).T
    print(data)
    bad = np.isnan(data).any()
    if not bad:
        x_arr, y_arr = data[0], data[1]
        return Airfoil(x_arr, y_arr)
    else:
        raise ValueError("airfoil file contains NaN")


def compare_airfoils():
    """Soon outdated"""
    path_str = r"airfoils/all_airfoils/coord_seligFmt"
    all_files = os.listdir(path_str)
    start_idx = all_files.index('naca64206.dat')
    files = all_files[start_idx:start_idx + 5]
    # reading airfoil geometry from txt file
    for filename in files:
        foil = read_airfoil(os.path.join(path_str, filename))
        xf.airfoil = foil
        xf.Re = 1e6
        xf.max_iter = 40
        a, cl, cd, cm, cp = xf.aseq(-5, 20, 2)
        plt.plot(a, cl, label=filename.split('.')[0])

    plt.legend()
    plt.grid()
    plt.show()


def analyse_airfoil(airfoilname, astart, afinal, astep):
    path_str = r"airfoils/all_airfoils/coord_seligFmt"
    foil = read_airfoil(os.path.join(path_str, airfoilname + '.dat'))
    xf.airfoil = foil
    xf.Re = 1e6
    xf.max_iter = 40
    a, cl, cd, cm, cp = xf.aseq(astart, afinal, astep)
    return a, cl, cd, cm, cp


class Airfoildata:
    """Contains airfoil data and methods to export"""

    def __init__(self, foilname):
        path_str = r"airfoils/all_airfoils/coord_seligFmt"
        self.foilname = foilname
        self.filepath = os.path.join(path_str, self.foilname + '.dat')

    def ReadGeometry(self, plot=False):
        """Opens a file in the format of two columns and creates Airfoil object"""
        data = np.genfromtxt(self.filepath, skip_header=1).T
        bad = np.isnan(data).any()
        if not bad:
            self.x_geom, self.y_geom = data[0], data[1]
        else:
            raise ValueError("airfoil file contains NaN")
        if plot:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.plot(self.x_geom, self.y_geom, color='red')
            ax.set_xlabel("x/c")
            ax.set_ylabel("t/c")
            ax.axis("equal")
            ax.grid()
            fig.tight_layout()
            plt.show()

    def getThicknessToCordRatio(self):
        """compute t/c from the geometry"""
        ncoord = len(self.y_geom)
        if ncoord % 2 == 0:
            top, bottom = self.y_geom[:ncoord//2], self.y_geom[ncoord//2::]
        else:
            top, bottom = self.y_geom[:ncoord//2 + 1], self.y_geom[ncoord//2::]
        thickness = top - bottom
        self.ttocratio, self.maxthickloc = np.max(thickness), self.x_geom[np.argmax(thickness)]

    def Analyse(self, astart, afinal, astep, Re=1e6, max_iter=40):
        """Run the XFoil script for this airfoil"""
        xf.airfoil = Airfoil(self.x_geom, self.y_geom)
        xf.Re, xf.max_iter = Re, max_iter
        self.a, self.cl, self.cd, self.cm, self.cp = xf.aseq(astart, afinal, astep)

    def getFromFile(self):
        """construct airfoil data from binary files"""
        data_dir = f"airfoils/{self.foilname}aero"
        self.a = np.load(os.path.join(data_dir, "AoA.npy"))
        self.cl = np.load(os.path.join(data_dir, "cl.npy"))
        self.cd = np.load(os.path.join(data_dir, "cd.npy"))
        self.cm = np.load(os.path.join(data_dir, "cm.npy"))
        self.cp = np.load(os.path.join(data_dir, "cp.npy"))

    def getLiftSlope(self):
        """Lift slope based on linear regime AoA [-6, 6]"""
        idx_arr = np.where(np.abs(self.a) < 6)
        self.cla = np.average(np.gradient(self.cl[idx_arr], self.a[idx_arr]))

    def PlotAeroData(self):
        """Plot the damn thing"""
        fig, axes = plt.subplots(2, 2, figsize=(8, 6))
        axes[0, 0].plot(self.a, self.cl, color='red')
        axes[0, 0].set_xlabel(r"$\alpha$ [deg]")
        axes[0, 0].set_ylabel("CL [-]")
        axes[0, 0].grid()

        axes[0, 1].plot(self.cl, self.cd, color='blue')
        axes[0, 1].set_xlabel("CL [-]")
        axes[0, 1].set_ylabel("CD [-]")
        axes[0, 1].grid()

        axes[1, 0].plot(self.a, self.cm, color='cyan')
        axes[1, 0].set_xlabel(r"$\alpha$ [deg]")
        axes[1, 0].set_ylabel("Cm [-]")
        axes[1, 0].grid()

        axes[1, 1].plot(self.a, self.cl/self.cd, color='deeppink')
        axes[1, 1].set_xlabel(r"$\alpha$ [deg]")
        axes[1, 1].set_ylabel("L/D [-]")
        axes[1, 1].grid()

        fig.tight_layout()
        plt.show()

    def ExportData(self):
        target_dir = f"airfoils/{self.foilname}aero"
        if not os.path.isdir(target_dir):
            os.mkdir(target_dir)
        np.save(os.path.join(target_dir, "AoA.npy"), self.a)
        np.save(os.path.join(target_dir, "cl.npy"), self.cl)
        np.save(os.path.join(target_dir, "cd.npy"), self.cd)
        np.save(os.path.join(target_dir, "cm.npy"), self.cm)
        np.save(os.path.join(target_dir, "cp.npy"), self.cp)


if __name__ == "__main__":
    airfoil = Airfoildata("naca642215")
    airfoil.ReadGeometry(plot=True)
    airfoil.getThicknessToCordRatio()
    airfoil.Analyse(-20, 20, 0.5)
    airfoil.getLiftSlope()
    print(airfoil.ttocratio, airfoil.maxthickloc, airfoil.cla)
    # airfoil.PlotAeroData()
