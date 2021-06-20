import pprint
import numpy as np
import matplotlib.pyplot as plt

def noise_empi_piston(P_br: float) -> float:
    return 57 + 14*np.log10(P_br)

def noise_empi_prop(P_br: float, D: float, M_t: float, B: int, N_p: int) -> float:
    return 83.4 + 15.3*np.log10(P_br) - 20*np.log10(D) + 38.5*M_t - 3*(B-2) + 10*np.log10(N_p)

def M_t(D: float, a: float, RPM: float) -> float:
    return (np.pi*D*RPM)/(a*60)

def dB2I(noisedB: float) -> float:
    I_min = 10**(-12)
    return (10**(noisedB/10))*I_min

def I2dB(I: float) -> float:
    I_min = 10**(-12)
    return 10*np.log10(I/I_min)

def inv_sq_law(r: float) -> float:
    try:
        return 1/(r*r)
    except ZeroDivisionError:
        return np.nan

def main() -> None:

    pp = pprint.PrettyPrinter(indent = 4)

    P_br1 = 29_329.2  # W
    P_br2 = 31_114.1  # W
    dist  = np.arange(0, 1000, 1)

    # I1 = dB2I(noise_empi_piston(P_br1)) * inv_sq_law(dist)
    # I2 = dB2I(noise_empi_piston(P_br2)) * inv_sq_law(dist)

    M_t1 = M_t(3.0, 347, 810)
    M_t2 = M_t(3.0, 347, 846)
    I1 = dB2I(noise_empi_prop(P_br1*1e-3, 3.0, M_t1, 3, 4)) * inv_sq_law(dist)
    I2 = dB2I(noise_empi_prop(P_br2*1e-3, 3.0, M_t2, 3, 4)) * inv_sq_law(dist)

    SPL1 = I2dB(I1)
    SPL2 = I2dB(I2)

    fig, ax = plt.subplots(1, figsize=(11, 11), dpi = 135)

    # ax.set_title(r"Noise versus Distance")
    ax.set_xlabel("Distance from source [m]", fontsize=16)
    ax.set_ylabel("Sound Pressure Level [dB]", fontsize=16)
    ax.grid(ls = "-.")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xticks(np.arange(0, 1100, 100))
    ax.set_yticks(np.arange(0, 120, 5))
    ax.set_xlim(0, 1000)
    ax.set_ylim(52.5, 115)
    ax.fill_between([0, 1000], 85, 115, hatch = "//", fc="#FF000055", edgecolor="white", linewidth=4.0)
    ax.plot(dist, SPL1, c="black", label="Top Propeller")
    ax.plot(dist, SPL2, c="black", ls = "--", label="Bottom Propeller")
    ax.tick_params(labelsize=16)
    ax.legend(loc = "upper right", ncol=2, fancybox=True, framealpha=1, shadow=True, borderpad=1, fontsize=13)

    fig.tight_layout()
    plt.show()

    try:
        fig.savefig(rf"Rotor Design/Noise_vs_Dist.pdf")
    except Exception:
        pass

if __name__ == "__main__":
    main()