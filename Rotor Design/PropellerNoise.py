import pprint
import numpy as np
import matplotlib.pyplot as plt

def noise_empi(P_br: float) -> float:
    return 57 + 14*np.log10(P_br)

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

    I1 = dB2I(noise_empi(P_br1)) * inv_sq_law(dist)
    I2 = dB2I(noise_empi(P_br2)) * inv_sq_law(dist)

    SPL1 = I2dB(I1)
    SPL2 = I2dB(I2)

    fig, ax = plt.subplots(1, figsize=(11, 11), dpi = 135)

    ax.set_title(r"Noise versus Distance")
    ax.set_xlabel("Distance from source [m]")
    ax.set_ylabel("Sound Pressure Level [dB]")
    ax.grid(ls = "-.")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.plot(dist, SPL1, c="black", label="Top Propeller")
    ax.plot(dist, SPL2, c="black", ls = "--", label="Bottom Propeller")
    ax.legend(loc = "upper right", ncol=2, fancybox=True, framealpha=1, shadow=True, borderpad=1)

    fig.tight_layout()
    plt.show()

    # try:
    #     fig.savefig(rf"Rotor Design/Noise_vs_Dist.pdf")
    # except Exception:
    #     pass

if __name__ == "__main__":
    main()