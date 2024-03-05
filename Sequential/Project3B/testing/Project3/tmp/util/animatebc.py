#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import mpl_toolkits.mplot3d.axes3d as p3

data = np.loadtxt("../build/steps.txt", float)
print(data.shape)

r = np.array([data[:, i : (i + 2)] for i in range(108 * 3)])
print(r.shape)
r1 = data[:, 0:2]
r2 = data[:, 3:5]
r3 = data[:, 6:8]
r4 = data[:, 9:11]

print(data.shape)
x = data[:, [[0], [3], [6], [9]]]
y = data[:, [[1], [4], [7], [10]]]
z = data[:, [[2], [5], [8], [11]]]

print(np.squeeze(x).shape)
print(np.squeeze(x[0].T))
print(np.squeeze(x[1].T))
print(np.squeeze(y[0].T))
print(np.squeeze(z[0].T))

# q = [
#    [-4.32, -2.17, -2.25, 4.72, 2.97, 1.74],
#    [2.45, 9.73, 7.45, 4.01, 3.42, 1.80],
#    [-1.40, -1.76, -3.08, -9.94, -3.13, -1.13],
# ]
# v = [
#    [0.0068, 0.024, -0.014, -0.013, -0.0068, -0.04],
#    [0.012, 0.056, -0.022, 0.016, 0.0045, 0.039],
#    [-0.0045, 0.031, 0.077, 0.0016, -0.015, -0.00012],
# ]

# xt = np.array(q[0])
# print(xt)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
points = ax.plot(np.squeeze(x[0].T), np.squeeze(y[0].T), np.squeeze(z[0].T), "*")[0]
ax.set(
    xlim=([0, 5.2]),
    ylim=([0, 5.2]),
    zlim=([0, 5.2]),
)
txt = fig.suptitle("")


def update_points(i, x, y, z, points):
    txt.set_text("num={:d}".format(i))  # for debug purposes

    points.set_data(np.squeeze(x[i].T), np.squeeze(y[i].T))
    points.set_3d_properties(np.squeeze(z[i].T), "z")


anim = FuncAnimation(fig, update_points, frames=data.shape[0], fargs=(x, y, z, points))

# plt.draw()
plt.show()
