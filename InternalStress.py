import Parameters as pm
import numpy as np
from BridgeDeflection import section,Nz_int,Mx_int,My_int

Ixx = section.Ix
Iyy = section.Iy

ys = [-section.h/2,section.h/2,section.h/2,-section.h/2]
xs = [-section.w/2,section.w/2,-section.w/2,section.w/2]





sigmas = []
maxsigma = []
for i in range(len(ys)):
    sigma_z = Nz_int/section.A + Mx_int/Ixx * ys[i] + (My_int/Iyy) * xs[i]
    sigmas.append(sigma_z)

    if np.sum(np.abs(sigma_z)) > np.sum(np.abs(maxsigma)):
        maxsigma =  np.sum(np.abs(sigma_z))
        sigma_final = sigma_z

f = lambda x: x**2
print(f(2))