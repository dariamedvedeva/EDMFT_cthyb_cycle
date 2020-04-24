import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

kpt = 

filename = "t_dispersion.dat"
D = np.loadtxt(filename)

dispersion = D[:, 2]
matrix = [ [0.0] * np.int(kpt) for i in range(np.int(kpt))]

for k_y in range(kpt):
    for k_x in range(kpt):
        matrix[k_y][k_x] = dispersion[k_y * kpt + k_x]

plt.figure()
plt.imshow(matrix, origin="lower", cmap=cm.rainbow)
plt.colorbar()
plt.show()
