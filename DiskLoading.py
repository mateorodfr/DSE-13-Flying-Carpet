"""A calculation tool for disk loading of the rotors,
including comparison to other rotorcraft"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar

rholst = [1.225, 1.0, 0.88]  # kg/m^3
alst = [343, 336, 326]  # m/s
MTOW = 2200  # kg
DiskLoad = np.linspace(100, 1000, 91)  # N/m^2  !!!!
g = 9.80665
Cp = 0.05  # normal range between 0.01 and 0.4
maxrpm_motor = 1300

colors = ["red", "blue", "green"]

fig, ax1 = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

for rho, a, col in zip(rholst, alst, colors):
    ThrustReq = 1.8 * MTOW * g
    TperEngine = ThrustReq / 8
    Preq = TperEngine * np.sqrt(DiskLoad / (2 * rho))
    Areq = TperEngine / DiskLoad
    Dreq = np.sqrt(4*Areq/np.pi)
    RPMreq = 60*np.cbrt(Preq / (rho * Cp * Dreq**5))  # rev/min

    # RPMreq_interp = interp1d(RPMreq, DiskLoad, fill_value="extrapolate")
    # maxdiskload = root_scalar(lambda x: RPMreq_interp(x) - maxrpm_motor, x0=800, x1=1100).root
    # print(maxdiskload)
    ax1[0].set_xlabel('Disk Loading (kg/m^2)')
    ax1[0].set_ylabel('Power required [kW]')
    ax1[0].plot(DiskLoad/10, Preq/1e3, label=f"rho = {rho}", color=col)
    ax1[0].tick_params(axis='y')
    # ax1[0].axvline(maxdiskload / g, 0.05, 0.95, linestyle='--', color="black")

    ax1[1].set_ylabel('RPM required [rev/min]')
    ax1[1].plot(DiskLoad / g, RPMreq, color=col)
    ax1[1].plot(DiskLoad / g, (36 * a / np.pi) / Dreq, linestyle='--', label=f"maxrpm {rho}", color=col)
    ax1[1].axhline(maxrpm_motor, label="engine limit", linestyle="--", color="cyan")
    # ax1[1].scatter(maxdiskload / g, maxrpm_motor, marker='x', color='black')

ax2 = ax1[0].twinx()  # instantiate a second axes that shares the same x-axis

color = "red"
ax2.set_ylabel('Diameter required [m]', color=color)  # we already handled the x-label with ax1
ax2.plot(DiskLoad / g, Dreq, color=color)  # , linestyle="--")
ax2.tick_params(axis='y', labelcolor=color)

ax1[0].legend()

fig.tight_layout()
plt.show()
