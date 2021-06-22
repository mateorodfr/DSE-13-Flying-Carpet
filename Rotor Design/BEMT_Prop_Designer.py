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
        
    
    solver = Solver(r"C:\Users\marvd\Documents\GitHub\DSE-13-Flying-Carpet\Rotor Design\08062021rot_v2.ini")

    Ts=[]
    T2s=[]
    Qs=[]
    Q2s=[]
    Ps=[]
    P2s=[]
    def run_sweep_coaxial(n, vinf):
        rpms, rpms2 = np.linspace(300, 1400, n), np.linspace(300, 1400, n)
        for rpm_now, rpm2_now in zip(rpms, rpms2):
            solver.rpm = rpm_now
            solver.rpm2 = rpm2_now
            T,Q,P,dfU,T2,Q2,P2,dfL = solver.run()

            Ts.append(T)
            T2s.append(T2)
            Qs.append(Q)
            Q2s.append(Q2)
            Ps.append(P)
            P2s.append(P2)

        fig, ax = plt.subplots(2, 2, figsize=(11, 11), dpi = 135, sharex=True)

        ax[0][0].set_title(r"Thrust versus RPMs @ $V_{\infty} = -1~[ms^{-1}]$")
        ax[0][0].set_xlim(0, 1500)
        ax[0][0].set_xlabel("RPMs [-]")
        ax[0][0].set_xticks(np.arange(0, 1500, 200))
        ax[0][0].set_ylim(0, 7500)
        ax[0][0].set_ylabel("Thrust [N]")
        ax[0][0].grid(ls = "-.")
        ax[0][0].spines["right"].set_visible(False)
        ax[0][0].spines["top"].set_visible(False)
        ax[0][0].plot(rpms, np.asarray(Ts), c="black", label="Top Propeller")
        ax[0][0].plot(rpms, np.asarray(T2s), c="black", ls="--", label="Bottom Propeller")
        ax[0][0].fill_between([1300, 1500], 0, 8000, hatch = "//", fc="#FF000055", edgecolor="white", linewidth=0.0)
        ax[0][0].legend(loc = "upper left", ncol=2, fancybox=True, framealpha=1, shadow=True, borderpad=1)

        ax[0][1].set_title(r"Power versus RPMs @ $V_{\infty} = -1~[ms^{-1}]$")
        ax[0][1].set_xlim(0, 1500)
        ax[0][1].set_xlabel("RPMs [-]")
        ax[0][1].set_xticks(np.arange(0, 1500, 200))
        ax[0][1].set_ylim(0, 160)
        ax[0][1].set_ylabel("Power [kW]")
        ax[0][1].fill_between([0, 1500], 120, 160, hatch = "//", fc="#FF000055", edgecolor="white", linewidth=0.0)
        ax[0][1].fill_between([1300, 1500], 0, 120, hatch = "//", fc="#FF000055", edgecolor="white", linewidth=0.0)
        ax[0][1].grid(ls = "-.")
        ax[0][1].spines["right"].set_visible(False)
        ax[0][1].spines["top"].set_visible(False)
        ax[0][1].plot(rpms, np.asarray(Ps)*1e-3, c="black", label="Top Propeller")
        ax[0][1].plot(rpms, np.asarray(P2s)*1e-3, c="black", ls="--", label="Bottom Propeller")

        ax[1][0].set_title(r"Torque versus RPMs @ $V_{\infty} = -1~[ms^{-1}]$")
        ax[1][0].set_xlim(0, 1500)
        ax[1][0].set_xlabel("RPMs [-]")
        ax[1][0].set_xticks(np.arange(0, 1500, 200))
        ax[1][0].set_ylim(0, 1100)
        ax[1][0].set_ylabel("Torque [Nm]")
        ax[1][0].fill_between([0, 1500], 930, 1100, hatch = "//", fc="#FF000055", edgecolor="white", linewidth=0.0)
        ax[1][0].fill_between([1300, 1500], 0, 930, hatch = "//", fc="#FF000055", edgecolor="white", linewidth=0.0)
        ax[1][0].grid(ls = "-.")
        ax[1][0].spines["right"].set_visible(False)
        ax[1][0].spines["top"].set_visible(False)
        ax[1][0].plot(rpms, Qs, c="black", label="Top Propeller")
        ax[1][0].plot(rpms, Q2s, c="black", ls="--", label="Bottom Propeller")

        ax[1][1].set_title(r"T/W versus RPMs @ $V_{\infty} = -1~[ms^{-1}]$")
        ax[1][1].set_xlim(0, 1500)
        ax[1][1].set_xlabel("RPMs [-]")
        ax[1][1].set_xticks(np.arange(0, 1500, 200))
        ax[1][1].set_ylim(0, 3.25)
        ax[1][1].set_ylabel("Thrust-to-Weight [-]")
        ax[1][1].fill_between([1300, 1500], 0, 3.25, hatch = "//", fc="#FF000055", edgecolor="white", linewidth=0.0)
        ax[1][1].grid(ls = "-.")
        ax[1][1].spines["right"].set_visible(False)
        ax[1][1].spines["top"].set_visible(False)
        ax[1][1].plot(rpms, 4*(np.asarray(Ts)+np.asarray(T2s))/(1667*G), c="black", ls="-", label="All Propellers", zorder=10)
        ax[1][1].axhline(y = 1.00, xmin = 0, xmax = 230/1500, c="#00A6D6", ls="-.")
        ax[1][1].axhline(y = 1.00, xmin = 370/1500, xmax = 1300/1500, c="#00A6D6", ls="-.")
        ax[1][1].text(300, 1.02, "Hover", fontsize = 7, horizontalalignment="center", verticalalignment="center", c="#00A6D6", weight="bold")
        ax[1][1].axhline(y = 1.10, xmin = 0, xmax = 1300/1500, c="#0066A2", ls="--")
        ax[1][1].text(300, 1.2, "Nominal Acceleration", fontsize = 7, horizontalalignment="center", verticalalignment="center", c="#0066A2", weight="bold")
        ax[1][1].axhline(y = 0.90, xmin = 0, xmax = 1300/1500, c="#0066A2", ls="--")
        ax[1][1].text(300, 0.80, "Nominal Deceleration", fontsize = 7, horizontalalignment="center", verticalalignment="center", c="#0066A2", weight="bold")
        ax[1][1].axhline(y = 1.95, xmin = 0, xmax = 1300/1500, c="#C3312F", ls="-.")
        ax[1][1].text(300, 1.85, "Critical Acceleration", fontsize = 7, horizontalalignment="center", verticalalignment="center", c="#C3312F", weight="bold")
        ax[1][1].legend(loc = "upper left", fancybox=True, framealpha=1, shadow=True, borderpad=1)

        fig.tight_layout()

        try:
            fig.savefig(rf"Rotor Design/T-RPM_P-RPM_at(Vinf {vinf}).png")
        except Exception:
            pass

        plt.show()

    # run_sweep_coaxial(20, -1)
            
    T,Q,P,dfU,T2,Q2,P2,dfL = solver.run()
    try:
        TSR, CT, CP = solver.turbine_coeffs(T, Q, P)
        TSR2, CT2, CP2 = solver.turbine_coeffs(T2, Q2, P2)
    except ZeroDivisionError:
        TSR, CT, CP = np.nan, *CT_CP(T, P, float(config["fluid"]["rho"]), float(config["rotor"]["diameter"]), float(config["case"]["rpm"])/60)
        TSR2, CT2, CP2 = np.nan, *CT_CP(T2, P2, float(config["fluid"]["rho"]), float(config["rotor2"]["diameter"]), float(config["case"]["rpm2"])/60)

    init_M = 1331 # [kg] -> 1497 -> 1488 -> 1385 -> 1344 -> 1334 -> 1331
    # T_W = 0.9
    # M_needed = init_M * (T_W * 4) / 8

    omega = float(config["case"]["rpm"])*2*np.pi / 60
    omega2 = float(config["case"]["rpm2"])*2*np.pi / 60

    P_avg = (P+P2)/2

    # print(f"\n\n{M_needed = } [kg]")
    print(f"\nMass carry-able = {(T+T2)*4 / G} [kg]\n")
    print(f"{P_avg = }\n")
    print(f"Rotor1: {CT = }, CQ = {CP/omega}, {CP = }, Ct/Cp = {CT / CP}\n") 
    print(f"Rotor2: {CT2 = }, CQ2 = {CP2/omega2}, {CP2 = }, Ct2/Cp2 = {CT2 / CP2}\n")

    


if __name__ == '__main__':
    '''
    For v2 following settings apply for sealevel:
    1.0g -> rpm = 958
    1.1g -> rpm = 1003
    1.8g -> rpm = 1282
    '''

    #1685

    G = 9.80665

    ini_file_path = os.path.join(os.getcwd(), "08062021rot_v2.ini")
    config = configparser.ConfigParser()
    config.read(r"C:\Users\marvd\Documents\GitHub\DSE-13-Flying-Carpet\Rotor Design\08062021rot_v2.ini")
    
    # setINI(15)

    main()

    # print(np.round(np.linspace(0.195, 0.175, 15), 4))

    # interpo = interp1d([0.125, 0.375, 0.625, 0.875, 1.125, 1.375], [10.0, 15.0, 17.0, 15.0, 13.0, 11.0], kind = "cubic", fill_value="extrapolate")
    # rng = np.asarray([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45])
    # print(np.round(interpo(rng), 3))