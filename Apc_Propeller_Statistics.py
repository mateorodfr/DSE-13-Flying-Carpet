import os
import numpy as np
import pprint

def is_all_n(s):
    try:
        for n in s:
            float(n)
        return True
    except ValueError:
        return False

def contains_NaN(s):
    newa = []
    for it in s:
        try:
            newa.append(float(it))
        except:
            print(it)
            newa.append(float(it[0:len(it)-4]))
            newa.append(np.NaN)
    return newa

pp = pprint.PrettyPrinter(indent=4)
path = r"APC_Propellers\PERFILES_WEB\PERFILES2\\"
files = sorted(os.listdir(path))

with open(path+f"{files[0]}") as f:
    alllines = f.readlines()
    fileout = [_.strip().split() for _ in alllines][17:]


a = [line for line in fileout if line]
a = [line for i, line in enumerate(a) if (i+1)%31 != 0 and (i+1)%32 != 0 and (i+1)%33 != 0]
pp.pprint(a)