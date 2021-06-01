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
