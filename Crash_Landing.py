"""File to analyse crash landing Drop"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz


class CrashLanding:

    def __init__(self, thrusttoweight, time, height, kmperh=False, triangle=True):
        self.cd = 0.8
        self.area = 12.04
        self.rho = 0.92   # 1.225
        self.mass = 1667
        self.g = 9.80665
        self.compression = 0.65
        self.height = height
        self.thrusttoweight = thrusttoweight
        self.time = time
        self.vterm = self.terminal_velocity()
        self.dt = np.average(np.diff(self.time))
        self.velocity, self.altitude = self.solve_ODE()
        self.ground_contact = self.get_ground_contact()
        self.contact_speed = self.velocity[np.where(self.time == self.ground_contact)][0]
        self.impact_duration = self.compute_impact_duration()
        self.impact_time = np.linspace(0, self.impact_duration, 1000)
        self.peakacc = self.compute_peak_acceleration()
        self.avgacc = 0.5 * self.peakacc
        if triangle:
            self.acceleration_arr = self.triangle_function()
        else:
            self.acceleration_arr = self.haversine_function()
        self.velocity_arr = self.crash_velocity()
        self.distance_arr = self.crash_distance()
        self.jerk_arr = self.crash_jerk()
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

    def compute_impact_duration(self):
        return self.compression / (0.5 * self.contact_speed)

    def compute_peak_acceleration(self):
        return 2 * self.contact_speed / self.impact_duration

    def triangle_function(self):
        return -1 * ((2 * self.peakacc / self.impact_duration) * self.impact_time
                     * (self.impact_time <= self.impact_duration / 2) * (self.impact_time >= 0)
                     + (2 * self.peakacc - (2 * self.peakacc / self.impact_duration) * self.impact_time)
                     * (self.impact_time >= self.impact_duration / 2) * (self.impact_time <= self.impact_duration))

    def haversine_function(self):
        return -1 * (self.peakacc * (np.sin(np.pi * self.impact_time / self.impact_duration))**2)

    def crash_velocity(self):
        return cumtrapz(self.acceleration_arr, self.impact_time, initial=0) + self.contact_speed

    def crash_distance(self):
        return cumtrapz(self.velocity_arr, self.impact_time, initial=0)

    def crash_jerk(self):
        return np.gradient(self.acceleration_arr, self.impact_time)


def main():
    thrustratios = [0.8, 0.9, 0.95]
    t_arr = np.linspace(0, 50, 1000)
    alt = 350

    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    fig2, ax2 = plt.subplots(3, 1, figsize=(8, 6), sharex=True)

    styles = [("darkgray", "dashdot"), ("dimgrey", "dashed"), ("black", "solid")]

    for Tratio, style in zip(thrustratios, styles):
        color, line = style
        crash = CrashLanding(Tratio, t_arr, alt, kmperh=False, triangle=True)
        crash.set_altitude_velocity_nonegative()
        print(crash.impact_duration, crash.peakacc / crash.g)

        tt = crash.time
        vv, zz = crash.velocity, crash.altitude

        tt_crash = crash.impact_time
        tt_crash -= np.average(tt_crash)
        jj_crash, aa_crash = crash.jerk_arr, crash.acceleration_arr
        vv_crash, ss_crash = crash.velocity_arr, crash.distance_arr

        ax[0].plot(tt, vv, color=color, linestyle=line, label=f"T/W = {Tratio}")
        ax[0].set_ylabel("velocity [m/s]")
        # ax[0].axvline(crash.ground_contact, color=color, linestyle=line)
        ax[0].legend()

        ax[1].plot(tt, zz, color=color, linestyle=line)
        # ax[1].axhline(color="black")
        # ax[1].axvline(crash.ground_contact, color="black")
        ax[1].set_ylabel("altituede [m]")
        ax[1].set_xlabel("time [s]")

        major_xticks = np.linspace(0, np.max(crash.time), 6)
        minor_xticks = np.linspace(0, np.max(crash.time), 51)

        major_vticks = np.linspace(0, 25, 6)
        minor_vticks = np.linspace(0, 25, 11)

        major_altticks = np.linspace(0, 350, 8)
        minor_altticks = np.linspace(0, 350, 15)

        ax[0].set_xticks(major_xticks)
        ax[0].set_xticks(minor_xticks, minor=True)
        ax[0].set_yticks(major_vticks)
        ax[0].set_yticks(minor_vticks, minor=True)
        ax[0].grid(which='minor', alpha=0.2)
        ax[0].grid(which='major', alpha=0.5)

        ax[1].set_xticks(major_xticks)
        ax[1].set_xticks(minor_xticks, minor=True)
        ax[1].set_yticks(major_altticks)
        ax[1].set_yticks(minor_altticks, minor=True)
        ax[1].grid(which='minor', alpha=0.2)
        ax[1].grid(which='major', alpha=0.5)

        # ax2[0].plot(tt_crash, ss_crash, label=f"T/W = {Tratio}")
        # # ax2[0].axhline(crash.compression, color="black", linestyle="--")
        # ax2[0].set_ylabel("distance [m]")
        # ax2[0].legend()

        ax2[0].plot(tt_crash, vv_crash, color=color, linestyle=line, label=f"T/W = {Tratio}")
        ax2[0].set_ylabel("velocity [m/s]")
        ax2[0].legend()

        ax2[1].plot(tt_crash, aa_crash / crash.g, color=color, linestyle=line)
        ax2[1].set_ylabel("Acceleration [G]")

        ax2[2].plot(tt_crash, jj_crash / crash.g, color=color, linestyle=line)
        ax2[2].set_ylabel("Jerk [G/s]")
        ax2[2].set_xlabel("time [s]")

    major_tticks = np.linspace(-0.06, 0.06, 7)
    minor_tticks = np.linspace(-0.06, 0.06, 31)

    major_v_crash_ticks = np.linspace(0, 25, 6)
    minor_v_crash_ticks = np.linspace(0, 25, 11)

    major_accticks = np.linspace(-100, 0, 6)
    minor_accticks = np.linspace(-100, 0, 11)

    major_jerkticks = np.linspace(-4000, 4000, 6)
    minor_jerkticks = np.linspace(-4000, 4000, 11)

    ax2[0].set_xticks(major_tticks)
    ax2[0].set_xticks(minor_tticks, minor=True)
    ax2[0].set_yticks(major_v_crash_ticks)
    ax2[0].set_yticks(minor_v_crash_ticks, minor=True)
    ax2[0].grid(which='minor', alpha=0.2)
    ax2[0].grid(which='major', alpha=0.5)

    ax2[1].set_xticks(major_tticks)
    ax2[1].set_xticks(minor_tticks, minor=True)
    ax2[1].set_yticks(major_accticks)
    ax2[1].set_yticks(minor_accticks, minor=True)
    ax2[1].grid(which='minor', alpha=0.2)
    ax2[1].grid(which='major', alpha=0.5)

    ax2[2].set_xticks(major_tticks)
    ax2[2].set_xticks(minor_tticks, minor=True)
    ax2[2].set_yticks(major_jerkticks)
    ax2[2].set_yticks(minor_jerkticks, minor=True)
    ax2[2].grid(which='minor', alpha=0.2)
    ax2[2].grid(which='major', alpha=0.5)

    # fig.tight_layout()
    fig2.tight_layout()
    plt.show()
    # plt.savefig(r"figures/Crash_OEI.pdf")


if __name__ == "__main__":
    main()
