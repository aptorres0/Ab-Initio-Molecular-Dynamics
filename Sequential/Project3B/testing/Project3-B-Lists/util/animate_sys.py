#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

data = np.loadtxt('../build/pos.txt') # 2000 x 3*108

# want to convert the data above into 3 matrices where each column is the value at a timepoint
x_columns = [i*3 for i in range(108)]
y_columns = [i*3+1 for i in range(108)]
z_columns = [i*3+2 for i in range(108)]
X = data[:, x_columns]
Y = data[:, y_columns]
Z = data[:, z_columns]

# create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = X[0]
y = Y[0]
z = Z[0]

points, = ax.plot(x, y, z,'o')
L = np.cbrt(108.0/0.8)
ax.set(
    xlim=(0,L), 
    ylim=(0,L), 
    zlim=(0,L)
    )
txt = fig.suptitle('Time Step: 0')

def update_points(num,points):
    txt.set_text(f'Time Step: {num}') # for debug purposes
    
    # Update the points
    points.set_data(X[num], Y[num])
    points.set_3d_properties(Z[num])
    
    # return the updated points
    return points, txt

ani=animation.FuncAnimation(fig, update_points, frames=2000, fargs=(points,), interval=1, blit=False)    

#plt.show()
# enable keyboard inturrpts to close figures
plt.show(block=False)
plt.pause(1)
input()
plt.close()
