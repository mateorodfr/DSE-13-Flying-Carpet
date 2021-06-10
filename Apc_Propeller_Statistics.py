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
        self.thrustfile = "PER2_STATIC-2.DAT"
        self.rpmfile = "PER2_RPMRANGE.DAT"
        self.rpmranges = self.get_rpmranges()
        self.data = self.read_static_thrust()
        self.pitch_data, self.final_data = self.get_data_per_propeller()

    def get_rpmranges(self):
        """Read the rpm file and get the range in a dictionary"""
        with open(os.path.join(self.filedir, self.rpmfile), "r") as f:
            lines = f.readlines()
            lines = [line.split() for line in lines]
            rpm_ranges = {}
            is_tenfold = lambda x: x % 10 == 0
            for tag, min, max in lines:
                min, max = int(min), int(max)
                if not is_tenfold(int(max)):    # bodged fix for some rpms ending in 999
                    max += 1
                rpm_ranges[tag] = (min, max)
            return rpm_ranges

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
            i += nvals

            # TODO Fix decimals in pitch reading (1225x48 = 12.25x4.8)

            # TODO check minimal Ct of 0.098

        return np.array(pitch_angles), np.array(nested_data)


    def plot_propdata(self):
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(self.pitch_data, self.final_data[:, -1], color="red", marker=".")
        ax.set_yticks(np.linspace(0, 5, 6))
        ax.set_xlabel("Blade angle [deg]")
        ax.set_ylabel("Ct/Cp [-]")

        plt.show()



if __name__ == "__main__":
    # print(pitch_to_angle(*decipher_tag("15x55MR")))
    # print(pitch_to_angle(*decipher_tag("20x15(WCAR-T6)")))

    data = ApcPropellerData()
    data.plot_propdata()
