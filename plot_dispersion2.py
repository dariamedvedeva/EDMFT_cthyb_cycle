import matplotlib.pyplot as plt
import numpy as np
#from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.transforms as mtransforms

kpt = 64

filename = "X_loc_v0.dat"
D = np.loadtxt(filename)

x           = D[:,0]
y           = D[:,1]
dispersion  = D[:,2]

size = []
for i in range(len(dispersion)):
    size.append(2.0)
    
#area    = size
colors  = dispersion
colormap = 'Spectral' # 'viridis' #'CMRmap'  # 'Spectral'
plt.figure()
#plt.title("NbS2 (t = [0.0614, 0.0975, 0.0151]; V = [0.6368, 0.4733, 0.4386])")
plt.title("NbS2")
#plt.ylabel("k_y")
#plt.xlabel("k_x")
plt.axis('off')
m = '.'
plt.scatter(x, -4.*np.pi /np.sqrt(3) + y, c=colors, cmap=colormap, marker = m) # IV
plt.scatter(-2.*np.pi + x, -6.*np.pi /np.sqrt(3) + y, c=colors, cmap=colormap, marker = m) # III
plt.scatter(x, y, c=colors, cmap=colormap, marker = m) # I
plt.scatter(-2.*np.pi + x, -2.*np.pi /np.sqrt(3) + y, c=colors, cmap=colormap, marker = m) # II
# points
plt.plot(0, 0, color='black', marker='o')
plt.plot(4. * np.pi / 3, 0, color='black', marker='o')
plt.plot(0, 2 * np.pi / np.sqrt(3), color='black', marker='o')
plt.plot(2 * np.pi / 3, 2 * np.pi / np.sqrt(3), color='black', marker='o')
# lines
# 1
x = [2. * np.pi / 3, 4. * np.pi / 3]
y = [2 * np.pi / np.sqrt(3), 0]
plt.plot(x, y, color='black')
# 2
x = [-2. * np.pi / 3, 2. * np.pi / 3]
y = [2 * np.pi / np.sqrt(3), 2 * np.pi / np.sqrt(3)]
plt.plot(x, y, color='black')
# 3
x = [-2. * np.pi / 3, -4. * np.pi / 3]
y = [2 * np.pi / np.sqrt(3), 0]
plt.plot(x, y, color='black')
# 4
x = [-2. * np.pi / 3, -4. * np.pi / 3]
y = [-2 * np.pi / np.sqrt(3), 0]
plt.plot(x, y, color='black')
# 5
x = [-2. * np.pi / 3, 2. * np.pi / 3]
y = [-2 * np.pi / np.sqrt(3), -2 * np.pi / np.sqrt(3)]
plt.plot(x, y, color='black')
# 6
x = [2. * np.pi / 3, 4. * np.pi / 3]
y = [-2 * np.pi / np.sqrt(3), 0]
plt.plot(x, y, color='black')

plt.colorbar()
plt.show()

def get_image(x_min, x_max, y_min, y_max):
    D = np.loadtxt(filename)
    x           = D[:,0]
    y           = D[:,1]
    dispersion  = D[:,2]
    
    kpt = 64
    delta_x = (2.*np.pi) / (kpt)
    delta_y = (4.*np.pi /np.sqrt(3)) / (kpt)
    x = np.arange(x_min, x_max, delta_x)
    y = np.arange(y_min, y_max, delta_y)
    print(str(x.shape))
    X, Y = np.meshgrid(x, y)
    i=0
    j=0
    Z = [ [0.0] * np.int(kpt) for i in range(np.int(kpt)) ]
    while (i < kpt):
        while (j < kpt):
            Z[j][i] = float(dispersion[kpt * i + j])
            j += 1
        else:
            j  = 0
            i += 1
    return Z

def do_plot(ax, Z, transform, extention):
    im = ax.imshow(Z, interpolation='none', extent = extention, origin='lower', clip_on=True)

    trans_data = transform + ax.transData
    im.set_transform(trans_data)

    # display intended extent of the image
    ax.set_xlim(-2.*np.pi, 2.*np.pi)
    ax.set_ylim(-13, 13)

# prepare image and figure
fig, (ax2) = plt.subplots(1, 1)

Z1 = get_image(0, 2.*np.pi, 0, 4.*np.pi /np.sqrt(3))
Z2 = get_image(-2.*np.pi, 0, -2.*np.pi /np.sqrt(3), 0)
Z3 = get_image(-2.*np.pi, 0, -6.*np.pi /np.sqrt(3), 0)
Z4 = get_image(0, 2.*np.pi, -4.*np.pi /np.sqrt(3), 0)

# image rotation
#do_plot(ax1, Z1, mtransforms.Affine2D().rotate_deg(0),  [0, 2.*np.pi, 0, 2.*np.pi /np.sqrt(3)])

# image skew
a = 0
b = 15
do_plot(ax2, Z1, mtransforms.Affine2D().skew_deg(a, b), [0, 2.*np.pi, 0, 6.*np.pi /np.sqrt(3)])
do_plot(ax2, Z1, mtransforms.Affine2D().skew_deg(a, b), [-2.*np.pi, 0, 0, 6.*np.pi /np.sqrt(3)])
do_plot(ax2, Z3, mtransforms.Affine2D().skew_deg(a, b), [-2.*np.pi, 0, -6.*np.pi /np.sqrt(3), 0])
do_plot(ax2, Z4, mtransforms.Affine2D().skew_deg(a, b), [0, 2.*np.pi, -6.*np.pi /np.sqrt(3), 0])

#plt.show()
