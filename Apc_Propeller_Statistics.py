import os
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

def get_rpmranges():
    with open(r"APC_Propellers\PERFILES_WEB\PER2_RPMRANGE.DAT", "r") as f:
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

def read_static_thrust():
    with open(r"APC_Propellers\PERFILES_WEB\PER2_STATIC-2.DAT", "r") as f:
        lines = f.readlines()[2:]

    datalines = [line.split() for line in lines if line.strip()[:4].isdigit() and not line.strip().endswith(".dat")]
    return np.asarray(datalines, dtype=float)    # Why won't you work

def get_data_per_propeller(datalines, rpm_ranges):
    for tag, vals in rpm_ranges.values():
        diameter, pitch = tag.split('x')
        if "-" in pitch:
            pitch = pitch.split("-")[0]
        while (not pitch.isdigit()) or (len(pitch) == 0):
            pitch = pitch[:-1:]
        pitch = float(pitch)




print(get_rpmranges()["9x4"])
print(read_static_thrust())
