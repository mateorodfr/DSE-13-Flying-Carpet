import numpy as np
import matplotlib.pyplot as plt


def linearMotion(initial, t_low,t_high, dt):

    initial = np.array(initial)
    t_out = np.arange(t_low,t_high+dt,dt)
    output = []
    for i in range(len(t_out)):
        output.append(initial)
        initial = initial + np.array([0,initial[0],initial[1]+2*initial[0]*t])*dt
    return output,t_out
def getTrajectory(init,As,ts):
    output, ts_out = [],[]
    for i in range(0,len(As)):

        init,t_out = linearMotion([As[i],init[-1][1],init[-1][2]],ts[i],ts[i+1],dt)
        output.extend(init)
        ts_out = np.hstack((ts_out,t_out))

    output = np.array(output)
    ts_out = np.array(ts_out)

    return output, ts_out



dt = 1
t = 250
t2 = 100

ts = [0,100,300,450,700,750,850,1000]
As = [-0.5,1,-0.8,-0.4,0.6,-0.5,0.6]
init = [[0,0,0]]

output,ts_out = getTrajectory(init,As,ts)


# initial = [-0.5,10,100]
# output,t_out = linearMotion(initial,0,t,dt)

# initial1 = [0.6,output[-1][1],output[-1][2]]
# output1, t_out1 = linearMotion(initial1,t,t2+t,dt)


# y = np.vstack((output,output1))
# t = np.concatenate((t_out,t_out1))



plt.plot(ts_out,output[:,2])
plt.show()