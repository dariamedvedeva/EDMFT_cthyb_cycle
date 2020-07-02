import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

kpt = 64

#filename = "X_loc_v0.dat"
filename = "X_loc_v0.dat"
D = np.loadtxt(filename)

x           = D[:,0]
y           = D[:,1]
dispersion  = D[:,2]

area    = dispersion
colors  = dispersion

plt.figure()
plt.title("NbS2 (t = [0.0614, 0.0975, 0.0151]; V = [0.6368, 0.4733, 0.4386])")
plt.ylabel("k_y")
plt.xlabel("k_x")
plt.scatter(x, y, c=colors, cmap='Spectral', marker = 's')
plt.colorbar()
plt.show()
