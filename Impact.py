import numpy as np
import Parameters as pm
import matplotlib.pyplot as plt


concept = pm.ConceptParameters(0)

n = 10
nh = 500


h = np.linspace(0,400,nh)
V = np.linspace(0,20,nh)

ds = 0.1
s = np.linspace(0.2,concept.cabin.L_cabin,n)


print(s)
Fxys,Fzs = [],[]
for i in range(len(s)):

    Fxy = (concept.Mtot_concept*0.5*V*V)/s[i]
    Fz = (concept.Mtot_concept*concept.physics.g*h)/s[i]
    Fxys.append(Fxy)
    Fzs.append(Fz)
Fxys = np.transpose(np.array(Fxys))
Fzs = np.transpose(np.array(Fzs))



fig, ax = plt.subplots(1,1)
rx = plt.plot(V, Fxys)
ax.legend(iter(rx), [str(np.round(j,2)) +' [m]' for j in s])
ax.title.set_text(r'Impact force as a fucntion of velocity (Sideways)')


fig1, ax1 = plt.subplots(1,1)
rx = plt.plot(V, Fzs)
ax1.legend(iter(rx), [str(np.round(j,2)) +' [m]' for j in s])
ax1.title.set_text(r'Impact force as a fucntion of velocity (gravity)')
plt.show()

