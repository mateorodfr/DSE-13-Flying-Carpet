import numpy as np
from pybemt.solver import Solver
import os
import matplotlib.pyplot as plt
import numpy as np
import configparser
from scipy.interpolate import interp1d


def setINI(n, *args) -> None:
        if not args:
            coaxial = input("Set coaxial:\t")
            rpm = input("Set RPM:\t")
            v_inf = input("Set v_inf:\t")
            dz = input("Set dz:\t\t")
            
            nblades = input("Set nblades:\t")
            diameter = input("Set diameter:\t")
            radius_hub = input("Set radius_hub:\t")
            section = input("Set airfoil:\t")
            chord = input("Set chord:\t")
            pitch = input("Set pitch:\t")

            rho = input("Set rho:\t")
            mu  = input("Set mu:\t\t")

            dr = float(diameter)/(2*n)

        lines  = f"[case]\n"
        lines += f"coaxial = {coaxial}\n"
        lines += f"rpm = {rpm}\n"
        lines += f"rpm2 = {rpm}\n"
        lines += f"v_inf = {v_inf}\n"
        lines += f"dz = {dz}\n"
        lines += f"\n"
        lines += f"[rotor]\n"
        lines += f"nblades = {nblades}\n"
        lines += f"diameter = {diameter}\n"
        lines += f"radius_hub = {radius_hub}\n"
        lines += f"section = " + " ".join([str(section) for _ in range(n)]) + "\n"
        lines += f"radius = " + " ".join([(str(np.round(dr * (_+.5), 3))) for _ in range(n)]) + "\n"
        lines += f"chord = " + " ".join([chord for _ in range(n)]) + "\n"
        lines += f"pitch = " + " ".join([pitch for _ in range(n)]) + "\n"
        lines += f"\n"
        if coaxial:
            lines += f"[rotor2]\n"
            lines += f"nblades = {nblades}\n"
            lines += f"diameter = {diameter}\n"
            lines += f"radius_hub = {radius_hub}\n"
            lines += f"section = " + " ".join([str(section) for _ in range(n)]) + "\n"
            lines += f"radius = " + " ".join([(str(np.round(dr * (_+.5), 3))) for _ in range(n)]) + "\n"
            lines += f"chord = " + " ".join([chord for _ in range(n)]) + "\n"
            lines += f"pitch = " + " ".join([pitch for _ in range(n)]) + "\n"
            lines += f"\n"
        lines += f"[fluid]\n"
        lines += f"rho = {rho}\n"
        lines += f"mu = {mu}"

        fpath = os.path.join(os.getcwd(), "08062021rot_v2.ini")
        with open(fpath, "w+") as f:
            f.write(lines)

def main() -> None:
    def CT_CP(Thrust, Power, rho, D, f):
        return (Thrust/(rho * D*D*D*D * f * f), Power/(rho * D*D*D*D*D * f * f * f))
    
    
    solver = Solver("08062021rot_v2.ini")

    
    T,Q,P,dfU,T2,Q2,P2,dfL = solver.run()
    try:
        TSR, CT, CP = solver.turbine_coeffs(T, Q, P)
        TSR2, CT2, CP2 = solver.turbine_coeffs(T2, Q2, P2)
    except ZeroDivisionError:
        TSR, CT, CP = np.nan, *CT_CP(T, P, float(config["fluid"]["rho"]), float(config["rotor"]["diameter"]), float(config["case"]["rpm"])/60)
        TSR2, CT2, CP2 = np.nan, *CT_CP(T2, P2, float(config["fluid"]["rho"]), float(config["rotor2"]["diameter"]), float(config["case"]["rpm2"])/60)

    init_M = 1667 # [kg]
    T_W = 1.1
    M_needed = init_M * (T_W * 4) / 8

    P_avg = (P+P2)/2

    print(f"\n\n{M_needed = } [kg]")
    print(f"\nMass carry-able = {(T+T2)*4 / G} [kg]\n")
    print(f"{P_avg = }\n")
    print(f"Rotor1: {CT = }, {CP = }, Ct/Cp = {CT / CP}\n")
    print(f"Rotor2: {CT2 = }, {CP2 = }, Ct/Cp = {CT2 / CP2}\n")

    


if __name__ == '__main__':
    '''
    For v2 following settings apply for sealevel:
    1.0g -> rpm = 958
    1.1g -> rpm = 1003
    1.8g -> rpm = 1282
    '''

    G = 9.80665

    ini_file_path = os.path.join(os.getcwd(), "08062021rot_v2.ini")
    config = configparser.ConfigParser()
    config.read(ini_file_path)
    
    # setINI(15)

    main()

    # print(np.round(np.linspace(0.195, 0.175, 15), 4))

    # interpo = interp1d([0.125, 0.375, 0.625, 0.875, 1.125, 1.375], [10.0, 15.0, 17.0, 15.0, 13.0, 11.0], kind = "cubic", fill_value="extrapolate")
    # rng = np.asarray([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45])
    # print(np.round(interpo(rng), 3))