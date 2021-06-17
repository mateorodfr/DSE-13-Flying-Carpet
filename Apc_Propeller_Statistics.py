import os

import matplotlib.pyplot as plt
import numpy as np
import pprint
# import matplotlib.pyplot as plt

# def is_all_n(s):
#     """Checks if a string is a float"""
#     try:
#         for n in s:
#             float(n)
#         return True
#     except ValueError:
#         return False

# def contains_NaN(s):
#     """Takes a list and checks if it contains a NaN concatenated with a float like -0.81-NaN"""
#     newa = []
#     for it in s:
#         try:
#             newa.append(float(it))
#         except Exception:
#             newa.append(float(it[0:len(it)-4]))
#             newa.append(np.NaN)
#     return newa

# pp = pprint.PrettyPrinter(indent=4)
# path = r"APC_Propellers\PERFILES_WEB\PERFILES2\\"
# files = sorted(os.listdir(path)) # get all files we want to read

# with open(path+f"{files[3]}") as f:
# # with open(os.path.join(path, "PER3_5x4R-RH.dat")) as f:
#     alllines = f.readlines()
#     fileout = [_.strip().split() for _ in alllines][1:] # first line is always part of the header

# if fileout[0]: # deal with fact that there are two types of headers
#     RPM_start = float(fileout[17][-1])
#     fileout = fileout[20:]
# else:
#     RPM_start = float(fileout[10][-1])
#     fileout = fileout[13:]

# a = [line for line in fileout if line] # drop empty lines
# final = []
# for i, line in enumerate(a): # filter out all the in between headers
#     if line[0][0].isalpha() or line[0].startswith("("):
#         pass  # lines starting with letters are skipped
#     else:
#         final.append(contains_NaN(line))

# final = np.asarray(final, dtype=float) # convert to numpy array

# # data = {}
# # for i in range(final.shape[0]//30):
# #     data[int(RPM_start * (i + 1))] = final[i:i + 30, :]
# # pp.pprint(data)
# final3D = np.asarray([final[i*30:(i+1)*30, :] for i in range(final.shape[0]//30)], dtype=float)

# pp.pprint(final3D[1].shape)

pp = pprint.PrettyPrinter(indent=4, depth=1)
path = r"APC_Propellers\PERFILES_WEB\PERFILES2\\"
allfiles = ["".join((path, fi)) for fi in sorted(os.listdir(path))]

def contains_NaN(s):
    """Takes a list and checks if it contains a NaN concatenated with a float like -0.81-NaN"""
    newa = []
    for it in s:
        try:
            newa.append(float(it))
        except Exception:
            newa.append(float(it[0:len(it)-4]))
            newa.append(np.NaN)
    return newa

def read_individual_prop():
    for file_0 in allfiles:
        with open(r"APC_Propellers\PERFILES_WEB\PERFILES2\PER3_5x3E.dat", "r") as f:
            lines = f.readlines()[1:-1]
            content = [c.strip() for c in lines]

            if content[17].startswith("PROP"):
                # pp.pprint("17")
                N = 16
            elif content[12].startswith("PROP"):
                # pp.pprint("12")
                N = 11
            else:
                pp.pprint("Different Format")
                break
            STEP = 37
            content = content[N:]
            chunks = len(content[N-1:])//STEP + 1
            # pp.pprint(content)
            final3D = np.asarray([content[index*STEP:(index+1)*STEP] for index in range(chunks)])[:, 5:-2]
            try:
                final3D = np.array(np.char.split(final3D).tolist(), dtype=float)
            except ValueError:
                pass # TODO implement fix for NaN values
        break

    pp.pprint(final3D[0])

def pitch_to_angle(Diameter, pitch):
    """Pitch at 75% of radius in inch to blade angle"""
    pg = pitch * 2.54 / 100
    D = Diameter * 2.54 / 100
    beta_75 = np.degrees(np.arctan((4 * pg)/(3 * np.pi * D)))
    return beta_75

def decipher_tag(tag):
    """Get the diameter and pitch from a propeller tag"""
    # The tags are given in the format DxP, some contain letters or dashes at the end
    diameter, pitch = tag.split('x')
    if "-" in pitch:
        pitch = pitch.split("-")[0]
    # strip the last letter until the remaining string is a number or empty string
    while (not pitch.isdigit()) or (len(pitch) == 0):
        pitch = pitch[:-1:]
    return float(diameter), float(pitch)


class ApcPropellerData:

    def __init__(self):
        self.filedir = r"APC_Propellers\PERFILES_WEB"
        self.propfiles = r"APC_Propellers\PERFILES_WEB\PERFILES2"
        self.thrustfile = "PER2_STATIC-2.DAT"
        self.rpmfile = "PER2_RPMRANGE.DAT"
        self.rpmranges = self.get_rpmranges()
        self.data = self.read_static_thrust()
        self.pitch_data, self.final_data, self.below_cpreq, self.above_ctcp = self.get_data_per_propeller()

    def get_rpmranges(self):
        """Read the rpm file and get the range in a dictionary"""
        with open(os.path.join(self.filedir, self.rpmfile), "r") as f:
            lines = f.readlines()
            lines = [line.split() for line in lines]
            rpm_ranges = {}
            is_tenfold = lambda x: x % 10 == 0
            for tag, min, max in lines:
                tag = self.tag2proper(tag)
                min, max = int(min), int(max)
                if not is_tenfold(int(max)):    # bodged fix for some rpms ending in 999
                    max += 1
                rpm_ranges[tag] = (min, max)
            return rpm_ranges

    def tag2proper(self, tag):
        properfile = os.path.join(self.propfiles, f"PER3_{tag}.dat")
        with open(properfile, "r") as pf:
            lines = pf.readlines()
            pline = " ".join([line.split() for line in lines][0])
            proper = pline.split()[0]
        return proper
    
    def read_static_thrust(self):
        """Read thrust file and slap everything into one large array"""
        with open(os.path.join(self.filedir, self.thrustfile), "r") as f:
            lines = f.readlines()[2:]

        datalines = [line.split() for line in lines if line.strip()[:4].isdigit() and not line.strip().endswith(".dat")]
        return np.asarray(datalines, dtype=float)    # Why won't you work

    def get_data_per_propeller(self):
        """Propeller data for each propeller for 1000 RPM"""
        cp, ct = self.data[:, 4], self.data[:, 5]

        ctovercp = ct/cp

        above_ctcp = []
        pitch_angles = []
        nested_data = []
        i = 0
        # Loop over the dictionary of rpm ranges
        for tag, vals in self.rpmranges.items():
            minrpm, maxrpm = vals
            # extract pitch angle from tag
            diameter, pitch = decipher_tag(tag)
            pitch_angle = pitch_to_angle(diameter, pitch)
            pitch_angles.append(pitch_angle)

            # Extract rotor data from vals and ct/cp
            nvals = (maxrpm - minrpm) // 1000
            cpi, cti, ctovercpi = cp[i], ct[i], ctovercp[i]
            if minrpm != 1000:
                cpi, cti, ctovercpi = np.nan, np.nan, np.nan
            nested_data.append([cpi, cti, ctovercpi])
            above_ctcp.append((cti/cpi > 3.202))
            if cti/cpi > 3.202 and cti > 0.09789676:
                print(tag)
            
            i += nvals
        nested_data = np.asarray(nested_data)
        below_cpreq = [(ct_val < 0.09789676) for ct_val in nested_data[:, 1]]

        return np.array(pitch_angles), nested_data, below_cpreq, above_ctcp


    def plot_propdata(self):
        fig, ax = plt.subplots(figsize=(16, 12), dpi=135)
        ax.grid(ls="-.", zorder=1, alpha=.5)
        for _ in range(len(self.pitch_data)):
            if self.below_cpreq[_]:
                ax.scatter(self.pitch_data[_], self.final_data[_, -1], color="purple", marker="x", zorder=10, alpha = .5)
            else:
                ax.scatter(self.pitch_data[_], self.final_data[_, -1], color="red", marker=".", zorder=10)
        
        ax.scatter([12.823], [3.202], s=60, color="blue", marker="$\circ$", zorder=10)
        ax.set_yticks(np.arange(0, 5.0, 0.5))
        ax.set_xticks(np.arange(0, 35, 2.5))
        ax.set_xlabel(r"Blade angle ($\beta_{0.75}$) [deg]", fontsize=18)
        ax.set_ylabel(r"$C_{T}/C_{P}$ [-]", fontsize=18)
        ax.tick_params(labelsize=18)
        ax.set_ylim(0, 4.25)
        ax.set_xlim(0, 32.5)

        plt.show()
        fig.savefig("CT-CPvsBeta075.pdf")



if __name__ == "__main__":
    # print(pitch_to_angle(*decipher_tag("15x55MR")))
    # print(pitch_to_angle(*decipher_tag("20x15(WCAR-T6)")))
    data = ApcPropellerData()
    data.plot_propdata()
