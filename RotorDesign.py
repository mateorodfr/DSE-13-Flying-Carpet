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
    bad = np.isnan(data).any() # TODO fix check to find NaN values
    if not bad:
        x_arr, y_arr = data[0], data[1]
        return Airfoil(x_arr, y_arr)
    else:
        raise ValueError("airfoil file contains NaN")

def compare_airfoils():
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


# TODO build software to analyse airfoils in rotor design

if __name__ == "__main__":
    # select airfoils from database
    foil = read_airfoil(r"airfoils/naca23018.csv")
    xf.airfoil = foil
    xf.Re = 1e6
    xf.max_iter = 40
    a, cl, cd, cm, cp = xf.aseq(-20, 20, 0.5)
    # plt.plot(a, cl)
    # plt.show()

    # TODO make this a function
    np.save("airfoils/naca23018aero/AoA.npy", a)
    np.save("airfoils/naca23018aero/cl.npy", cl)
    np.save("airfoils/naca23018aero/cd.npy", cd)
    np.save("airfoils/naca23018aero/cm.npy", cm)
    np.save("airfoils/naca23018aero/cp.npy", cp)
