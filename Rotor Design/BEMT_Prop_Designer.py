import numpy as np
from pybemt.solver import Solver
import os
import matplotlib.pyplot as plt
import numpy as np
import pprint
import configparser

def main():
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
        lines += f"radius = " + " ".join([(str(np.round(dr * (_+.5), 2))) for _ in range(n)]) + "\n"
        lines += f"chord = " + " ".join([chord for _ in range(n)]) + "\n"
        lines += f"pitch = " + " ".join([pitch for _ in range(n)]) + "\n"
        lines += f"\n"
        if coaxial:
            lines += f"[rotor2]\n"
            lines += f"nblades = {nblades}\n"
            lines += f"diameter = {diameter}\n"
            lines += f"radius_hub = {radius_hub}\n"
            lines += f"section = " + " ".join([str(section) for _ in range(n)]) + "\n"
            lines += f"radius = " + " ".join([(str(np.round(dr * (_+.5), 2))) for _ in range(n)]) + "\n"
            lines += f"chord = " + " ".join([chord for _ in range(n)]) + "\n"
            lines += f"pitch = " + " ".join([pitch for _ in range(n)]) + "\n"
            lines += f"\n"
        lines += f"[fluid]\n"
        lines += f"rho = {rho}\n"
        lines += f"mu = {mu}"

        with open("rot_sl_OEI.ini.ini", "w+") as f:
            f.write(lines)
    def CT_CP(Thrust, Power, rho, D, f):
        return Thrust/(rho * D*D*D*D * f * f), Power/(rho * D*D*D*D*D * f * f * f)
    
    
    # setINI(6)

    # solver = Solver(r"C:\Users\marvd\Documents\GitHub\DSE-13-Flying-Carpet\Rotor Design\rot.ini")
    solver = Solver(os.path.join(os.getcwd(), r"rot_sl_OEI.ini"))
    # df, df_section = solver.run_sweep("rpm", 40, 900, 1300)
    # pp.pprint(df)
    
    T,Q,P,dfU,T2,Q2,P2,dfL = solver.run()
    try:
        a1 = solver.turbine_coeffs(T, Q, P)
        a2 = solver.turbine_coeffs(T2, Q2, P2)
    except ZeroDivisionError:
        a1 = (np.nan, *CT_CP(T, P, float(config["fluid"]["rho"]), float(config["rotor"]["diameter"]), float(config["case"]["rpm"])/60))
        a2 = (np.nan, *CT_CP(T2, P2, float(config["fluid"]["rho"]), float(config["rotor2"]["diameter"]), float(config["case"]["rpm2"])/60))

    print(f"\n{(T+T2)*4 / 9.80665} [kg]")
    print(f"\n rotor1: Ct = {a1[1]}, Cp = {a1[2]}, Ct/Cp = {a1[1] / a1[2]}")
    print(f"\n rotor2: Ct = {a2[1]}, Cp = {a2[2]}, Ct/Cp = {a2[1] / a2[2]}")



if __name__ == '__main__':
    pp = pprint.PrettyPrinter(indent=4)
    config = configparser.ConfigParser()
    config.read("rot_sl_OEI.ini")
    main()
