from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

np.set_printoptions(threshold=np.inf)
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')


def init():
    global points
    points = np.loadtxt('../build/pos.txt')


def animate(i):
    ax1.clear()
    ax1.plot(points[:i, 0], points[:i, 1], points[:i, 2])

ani = animation.FuncAnimation(fig, animate, init_func=init, interval=1000)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')
plt.show()