"""A calculation tool for disk loading of the rotors,
including comparison to other rotorcraft"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar

rholst = [1.225, 0.923247]  # kg/m^3
alst = [340.294, 329.174]  # m/s
MTOW = 1667  # kg
DiskLoad = np.linspace(100, 1000, 91)  # N/m^2  !!!!
g = 9.80665
Cp = 0.05  # normal range between 0.01 and 0.4
maxrpm_motor = 1300
design_D = 3    # m

colors = ["red", "green"]

fig, ax1 = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

ThrustReq = 1.9 * MTOW * g
TperEngine = ThrustReq / 8
Areq = TperEngine / DiskLoad
Dreq = np.sqrt(4*Areq/np.pi)

Design_diskload = DiskLoad[np.argmin(np.abs(Dreq - design_D))] / g

for rho, a, col in zip(rholst, alst, colors):
    Preq = TperEngine * np.sqrt(DiskLoad / (2 * rho))
    RPMreq = 60*np.cbrt(Preq / (rho * Cp * Dreq**5))  # rev/min


    Design_Ct = TperEngine / (rho * (maxrpm_motor / 60) ** 2 * design_D ** 4)
    maxCttoCpratio = 4.25
    maxCp = Design_Ct / maxCttoCpratio
    Design_maxCp = 120e3 / (rho * (1300 / 60) ** 3 * design_D ** 5)

    print(f"Ct = {Design_Ct} \n"
          f"Cp = {Design_maxCp} \n"
          f"maxCp = {maxCp} \n")


    # RPMreq_interp = interp1d(RPMreq, DiskLoad, fill_value="extrapolate")
    # maxdiskload = root_scalar(lambda x: RPMreq_interp(x) - maxrpm_motor, x0=800, x1=1100).root
    # print(maxdiskload)
    ax1[0].set_ylabel('Power required [kW]')
    ax1[0].plot(DiskLoad/10, Preq/1e3, label=f"rho = {rho}", color=col)
    ax1[0].tick_params(axis='y')
    # ax1[0].axvline(maxdiskload / g, 0.05, 0.95, linestyle='--', color="black")

    ax1[1].set_xlabel('Disk Loading (kg/m^2)')
    ax1[1].set_ylabel('RPM required [rev/min]')
    ax1[1].plot(DiskLoad / g, RPMreq, color=col)
    ax1[1].plot(DiskLoad / g, (60 * 0.6 * a / np.pi) / Dreq, linestyle='--', label=f"maxrpm {rho}", color=col) # tip speed reaches mach 0.7
    ax1[1].axhline(maxrpm_motor, label="engine limit", linestyle="--", color="cyan")
    # ax1[1].scatter(maxdiskload / g, maxrpm_motor, marker='x', color='black')

ax2 = ax1[0].twinx()  # instantiate a second axes that shares the same x-axis

color = "blue"
ax2.set_ylabel('Diameter required [m]', color=color)  # we already handled the x-label with ax1
ax2.plot(DiskLoad / g, Dreq, color=color)  # , linestyle="--")
ax2.axvline(Design_diskload, color="black", linestyle="--")
ax2.tick_params(axis='y', labelcolor=color)

ax1[1].axvline(Design_diskload, color="black", linestyle="--")
ax1[0].legend()

fig.tight_layout()
plt.show()
