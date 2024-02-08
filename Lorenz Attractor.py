"Lorenz Attractor"

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

X, Y, Z = [], [], []

def generator():
    sigma, beta = 10.0, 8 / 3
    rou = 28

    delta = 0.01
    x, y, z = 0.1, 0, 0
    for i in range(10000):
        _x = x + delta * sigma * (y - x)
        _y = y + delta * (x * (rou - z) - y)
        _z = z + delta * (x * y - beta * z)

        x, y, z = _x, _y, _z
        X.append(x)
        Y.append(y)
        Z.append(z)

mpl.rcParams["legend.fontsize"] = 18
fig = plt.figure()
axes = Axes3D(fig)

generator()
axes.plot(X, Y, Z, label = "Lorenz Attractor")
axes.legend()

plt.show()
