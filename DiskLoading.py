"""A calculation tool for disk loading of the rotors,
including comparison to other rotorcraft"""

import numpy as np
import matplotlib.pyplot as plt

MTOW = 2200 # kg
DiskLoad = np.linspace(10, 100, 100) # kg/m^2
# RPM = 1300 # rev/min
rho = 1.225 # kg/m^3

ThrustReq = 1.2 * MTOW * 9.81
TperEngine = ThrustReq / 8
Preq = TperEngine * np.sqrt(DiskLoad * (2 * rho))
Areq = TperEngine / DiskLoad
Dreq = np.sqrt(4*Areq/np.pi)

fig, ax1 = plt.subplots(figsize=(8, 6))

color = 'tab:red'
ax1.set_xlabel('Disk Loading (kg/m^2)')
ax1.set_ylabel('Power required [W]', color=color)
ax1.plot(DiskLoad, Preq, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Diameter required [m]', color=color)  # we already handled the x-label with ax1
ax2.plot(DiskLoad, Dreq, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()
plt.show()
