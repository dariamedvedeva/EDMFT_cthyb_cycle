# -*- coding: utf-8 -*-
# 3 columns: x, y and intencity
import sys
import os
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
global kpt
if len (sys.argv) > 1:
    print ("Ploting...".format (sys.argv[1]) )
    input_filename = sys.argv[1]
    kpt = int(sys.argv[2])
    print(kpt)
else:
    print ("Enter filename")

#save picture
def save(name = '', fmt = 'png'):
    pwd = os.getcwd()
    plt.savefig('%s.%s' % (name, fmt), fmt = 'png')
    os.chdir(pwd)


def picture():
    
    ##############
    # read file
    ##############
    f = open(input_filename)
    data = f.readlines()

    global kpt
    intensity = [ [0.0] * np.int(kpt) for i in range(np.int(kpt))]

    i=0 
    j=0
    for str in data:
        value = str.split()
        value = [float(x) for x in value]  
        intensity[i][j] = value[2]
    i +=1
        if(i==kpt): 
            j += 1
            i  = 0
    f.close()



    #############################
    # PLOT
    ############################
    plt.figure()
    # plt.title('Static magnetic susceptibility X(q)\n', color='k', size=20)

    size_mesh = kpt * kpt
    X = [[0.0] for i in range(np.int(kpt*kpt))]
    Y = [[0.0] for i in range(np.int(kpt*kpt))]
    num = 0
    x_0 = 0.0
    y_0 = 0.0
    # Mesh = [ [0] * 2 for i in range(size_mesh)]
    #  for x in range(0, kpt, 1):  # along x
    #    for y in range(0, kpt, 1):  # along y
    #
    #       X[num] = x_0 - 0.5 * np.float(y)
    #       Y[num] = round(y_0 + np.float(y) * (np.sqrt(3.0) / 2.0), 6)
    #
    #   num += 1
    #   x_0 += 1.0
    #   y_0 = 0.0

    # (X, Y) = np.meshgrid(X, Y)
    # vmin=1.17, vmax=1.67

    plt.imshow(intensity, origin="lower", cmap=cm.rainbow)

    # plt.xticks((0, 49), ('$\Gamma$', 'K'), color='k', size=20)
    # plt.yticks((0, 49), ('$\Gamma$', 'K'), color='k', size=20)

    plt.grid(False)
    plt.axis('off')
    plt.colorbar()
    plt.draw()


    filename = 'out'
 #   save(name=filename, fmt='pdf')
 #   save(name=filename, fmt='png')
    print ("Finish.")
    plt.show()
    exit()

picture()
