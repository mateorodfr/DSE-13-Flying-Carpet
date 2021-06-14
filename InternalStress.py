import Parameters as pm
import numpy as np
from BridgeDeflection import Vy_int,Vx_int, Tz_int
import matplotlib.pyplot as plt

taus = []
maxtau = []
prop_circle = [0.05]
prop_square = [0.05,0.05]
t_circ = [0.001]
t_square = [0.001,0.0015]
section = pm.CrossSectionParameters('square',prop_square,t_square)
Vymax, Vxmax , Tmax = np.amax(np.abs(Vy_int)) , np.amax(np.abs(Vx_int)), np.amax(np.abs(Tz_int)) ######CHANGE MEEEEE########


tau, taumax, tauxyz, xyarr = internalShear(section,Vymax,Vxmax,Tmax)

