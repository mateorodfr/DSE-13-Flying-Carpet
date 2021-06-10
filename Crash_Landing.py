"""File to analyse crash landing Drop"""
import numpy as np
import matplotlib.pyplot as plt


class CrashLanding:

    def __init__(self, thrusttoweight, time, height, kmperh=False):
        self.cd = 0.8
        self.area = 12.04
        self.rho = 1.225
        self.mass = 1995
        self.g = 9.80665
        self.height = height
        self.thrusttoweight = thrusttoweight
        self.time = time
        self.vterm = self.terminal_velocity()
        self.dt = np.average(np.diff(self.time))
        self.velocity, self.altitude = self.solve_ODE()
        self.ground_contact = self.get_ground_contact()
        if kmperh:
            self.velocity *= 3.6
            self.vterm *= 3.6

    def terminal_velocity(self):
        return np.sqrt((2 * self.mass / self.area * self.cd * self.rho) * (1 - self.thrusttoweight) * self.g)

    def solve_ODE(self):
        v = 0
        z = self.height
        vlst = []
        zlst = []
        for _ in self.time:
            a = self.g * (1 - self.thrusttoweight) - (0.5 * self.area * self.cd * self.rho) * v*v / self.mass
            v += a * self.dt
            z -= v * self.dt + 0.5 * a * self.dt**2
            vlst.append(v)
            zlst.append(z)

        return np.array(vlst), np.array(zlst)

    def get_ground_contact(self):
        return self.time[np.argmin(np.abs(self.altitude))]

    def set_altitude_velocity_nonegative(self):
        self.velocity[np.where(self.altitude < 0)], self.altitude[np.where(self.altitude < 0)] = np.nan, np.nan

    def compute_impact_momentum(self):
        pass    # TODO

def main():
    thrustratios = [0, 0.8, 0.9, 0.95, 1]
    t_arr = np.linspace(0, 50, 1000)
    alt = 400

    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    for Tratio in thrustratios:
        crash = CrashLanding(Tratio, t_arr, alt, kmperh=True)
        crash.set_altitude_velocity_nonegative()
        tt = crash.time
        vv, zz = crash.velocity, crash.altitude
        ax[0].plot(tt, vv, label=f"Terminal = {np.round(crash.vterm, 2)}")
        ax[0].set_ylabel("velocity [km/h]")
        ax[0].axvline(crash.ground_contact, color="black", linestyle='--')
        ax[0].legend()

        ax[1].plot(tt, zz)
        ax[1].axhline(color="black")
        ax[1].axvline(crash.ground_contact, color="black")
        ax[1].set_ylabel("altituede [m]")
        ax[1].set_xlabel("time [s]")

    fig.tight_layout()
    plt.show()



if __name__ == "__main__":
    main()
