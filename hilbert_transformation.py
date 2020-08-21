# Medvedeva 04.05.2020
# the Hylbert transformation from rho(w) to G(iw_n) and Delta(iw_n)
# and the same on real frequencies

import numpy as np
from scipy import integrate
import os, sys
import matplotlib.pyplot as plt

# Parameters
E_f = 0.0

#print("python script.py dos.dat E_f number_of_mats_points beta")

if len (sys.argv) > 1:
    fileDOS = sys.argv[1]
    E_f = float(sys.argv[2])
    number_of_mats_freq = int(sys.argv[3])
    beta = float(sys.argv[4])
    print("Filename is {}. Number od mats freqs is {}. Beta is {}".format(fileDOS, number_of_mats_freq, beta))
else:
    print ("\n Parameters:\n 1 - this_script_name.py \n 2 - file_with_rho_from_Lda.dat \n 3 - E_f \n 4 - n \n 5 - b. \n\n Here E_f - Fermi energy, n is number of matsubara points,  b is beta.\n")
    sys.exit()

ro          = []
omega_real  = []

# Read bare DOS (from LDA)
f = open(fileDOS, "r")
lines_in_file = f.readlines()
ind = 0
for line in lines_in_file:
    line_data = line.split()
    omega_real.append(float(line_data[0]))
    ro.append(float(line_data[1]))
f.close()
len_ro = ro.__len__()

# -- scipy.integrate.trapz(y, x=None, dx=1.0, axis=-1)[source]
# -- Integrate along the given axis using the composite trapezoidal rule.
# -- Integrate y (x) along given axis.

integral = np.round(integrate.trapz(ro, omega_real), 4)
print("Integral ro = {}".format(integral))
if (integral == np.round(1.0, 4)):
    print("DOS is normalized")
else:
    ro /= integral
    # after normalization
    integral2 = np.round(integrate.trapz(ro, omega_real), 4)
    print("Integral ro = {} after normalization".format(integral2))

# Calculate the energy level of the orbital(s), E_d
# How to calculate:
# \int { E * ro(w) dw } = E_0 \int { ro(w) dw }
# == >
# \int{(E - E_0) * ro(w) * dw} = 0
# in other words: E_d = \int { x * f(x) dx } for normalized DOS

x_fx = np.zeros(len_ro, np.float)
for i in range(len_ro):
    x_fx[i] = ro[i] * omega_real[i]
E_d = np.round(integrate.trapz(x_fx, omega_real), 4)
print("E_d = {}".format(E_d))

# Calculate G(iw_n)

## create an array of mat. freq-s
mats_freq = np.zeros(number_of_mats_freq, np.complex)
for n in range(number_of_mats_freq):
    mats_freq[n] = 0.0 + 1j * (2 * n + 1) * np.pi / beta
    
## transformation

print("E_F = {}".format(E_f))
G_iw = np.zeros(number_of_mats_freq, np.complex)
for iwn in range(number_of_mats_freq):
    for omega in range(len_ro):
        G_iw[iwn] += ro[omega] / (mats_freq[iwn] - omega_real[omega] + E_f)
G_iw *= (omega_real[1] - omega_real[0])

G_real  = np.zeros(len_ro, np.complex)
for omega1 in range(len_ro):
    for omega2 in range(len_ro):
        G_real[omega1] += ro[omega2] / (omega_real[omega1] - omega_real[omega2] + 0.1 * 1j + E_f)
G_real *= (omega_real[1] - omega_real[0])

## Delta
Delta_iw = np.zeros(number_of_mats_freq, np.complex)
for iwn in range(number_of_mats_freq):
    Delta_iw[iwn] = mats_freq[iwn] - E_d - ( 1.0 / G_iw[iwn])

Delta_real = np.zeros(len_ro, np.complex)
Delta_real = -(1.0 / G_real).imag


# # # # # # # # # # # # # #
#   S A V E  F I L E S    #
# # # # # # # # # # # # # #
filename1 = 'Green_mats.out'
filename2 = 'Delta.dat'
f1  = open(filename1, "w")
f2  = open(filename2, "w")
for i in range(number_of_mats_freq):
    f1.write("{}\t {}\t {}\n".format(mats_freq.imag[i], G_iw.real[i], G_iw.imag[i]))        # G_iw
    f2.write("{}\t {}\t {}\n".format(mats_freq.imag[i], Delta_iw.real[i], Delta_iw.imag[i]))    # D_iw
f1.close()
f2.close()

filename1 = 'Green_re.out'
filename2 = 'Delta_re.out'
f1 = open(filename1, "w")
f2 = open(filename2, "w")
for i in range(len_ro):
    f1.write("{}\t {}\t {}\n".format(omega_real[i], Delta_real.real[i],  Delta_real.imag[i]))    # D_w
    f2.write("{}\t {}\t {}\n".format(omega_real[i], Delta_real.real[i],  Delta_real.imag[i]))    # D_w
f1.close()
f2.close()
