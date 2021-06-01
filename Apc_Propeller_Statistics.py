import os
import numpy as np
import pprint

def is_all_n(s):
    """Checks if a string is a float"""
    try:
        for n in s:
            float(n)
        return True
    except ValueError:
        return False

def contains_NaN(s):
    """Takes a list and checks if it contains a NaN concatenated with a float like -0.81-NaN"""
    newa = []
    for it in s:
        try:
            newa.append(float(it))
        except Exception:
            # print(it)
            newa.append(float(it[0:len(it)-4]))
            newa.append(np.NaN)
    return newa

pp = pprint.PrettyPrinter(indent=4)
path = r"APC_Propellers\PERFILES_WEB\PERFILES2\\"
files = sorted(os.listdir(path)) # get all files we want to read

# with open(path+f"{files[3]}") as f:
with open(os.path.join(path, "PER3_5x4R-RH.dat")) as f:
    alllines = f.readlines()
    fileout = [_.strip().split() for _ in alllines][1:] # first line is always part of the header

if fileout[0]: # deal with fact that there are two types of headers
    RPM_start = float(fileout[17][-1])
    fileout = fileout[20:]
else:
    RPM_start = float(fileout[10][-1])
    fileout = fileout[13:]

a = [line for line in fileout if line] # drop empty lines
final = []
for i, line in enumerate(a): # filter out all the in between headers
    if line[0][0].isalpha() or line[0].startswith("("):
        pass  # lines starting with letters are skipped
    else:
        final.append(contains_NaN(line))

final = np.asarray(final, dtype=float) # convert to numpy array
pp.pprint(final)
print(final.shape)
data = {}
for i in range(final.shape[0]//30):
    data[RPM_start * (i + 1)] = final[i:i + 30, :]

pp.pprint(data)
