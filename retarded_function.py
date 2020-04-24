import numpy as np
from math import cos, sin, pi, floor, ceil
import sys
import scipy.optimize as optimization
import iteration_cycle
import h5py as h5
#############################
#
#     RETARDED FUNCTION
#
#############################
"""
From CT-HYB SEGMENT help's notation:
N   - number of imaginary time discretization points
N_W - number of bosonic Matsubara frequencies
"""

def compute_function(nmats, N, beta, lambda_filename):

    print ("Fourier transform Lambda and creation of K")
    
    # data & parameters
    data_charge_, data_spin_ = iteration_cycle.read_2_functions_in_file(lambda_filename)
    
    c_k_w    = np.zeros(nmats, np.complex)
    c_kp_w   = np.zeros(nmats, np.complex)
    sz_k_w   = np.zeros(nmats, np.complex)
    sz_kp_w  = np.zeros(nmats, np.complex)

    for n_w in range(nmats):
        if n_w != 0:
            # Avoids divide by zero. Constant term will be shifted, in k_tau -= k_tau[0]
            # The w=0 term gives ~tau**2/beta*data_[0] to k_tau and tau*2/beta*data_[0] to kp_tau, which is added later
            W = 2 * np.pi / beta * n_w
            c_k_w[n_w]   = data_charge_[n_w] * -1. / (W * W) #(1/(iW)^2)=-1/W^2
            c_kp_w[n_w]  = data_charge_[n_w] * 1.j / W
            sz_k_w[n_w]  = data_spin_[n_w] * -1.0/(W*W) # (1/(iW)^2)=-1/W^2
            sz_kp_w[n_w] = data_spin_[n_w] * 1j/W

    number_of_tau_points = N + 1
    c_k_tau     = np.zeros(number_of_tau_points, np.float)
    c_kp_tau    = np.zeros(number_of_tau_points, np.float)
    sz_k_tau    = np.zeros(number_of_tau_points)
    sz_kp_tau   = np.zeros(number_of_tau_points)

    for t in range(number_of_tau_points):
        tau     = t * beta / N
        tmp1_c  = 0.0
        tmp2_c  = 0.0
        tmp1_sz = 0.0
        tmp2_sz = 0.0
        for n in range(nmats):
            tmp1_c  += c_k_w[n].real  * cos( tau*(2*n) * pi / beta )  + c_k_w[n].imag  * sin( tau*(2*n) * pi / beta )
            tmp2_c  += c_kp_w[n].real * cos( tau*(2*n) * pi / beta )  + c_kp_w[n].imag * sin( tau*(2*n) * pi / beta )
            tmp1_sz += sz_k_w[n].real*cos( tau*(2*n)*pi/beta ) + sz_k_w[n].imag*sin( tau*(2*n)*pi/beta )
            tmp2_sz += sz_kp_w[n].real*cos( tau*(2*n)*pi/beta ) + sz_kp_w[n].imag*sin( tau*(2*n)*pi/beta )

        tmp1_c  += data_charge_[0].real*( (tau - beta/2.)**2 - beta**2/4. )/4.
        tmp2_c  += data_charge_[0].real*(tau - beta/2.)/2.
        tmp1_sz += data_spin_[0].real*( (tau - beta/2.)**2 - beta**2/4. )/4.
        tmp2_sz += data_spin_[0].real*(tau - beta/2.)/2.

        c_k_tau[t]  = tmp1_c * 2. / beta
        c_kp_tau[t] = tmp2_c * 2. / beta
        sz_k_tau[t] = tmp1_sz * 2./beta
        sz_kp_tau[t]= tmp2_sz * 2./beta

    print ("Symmetry check of k_tau (charge): ", c_k_tau[0],  c_k_tau[-1]) # First and last element
    print ("Symmetry check of k_tau (spin): ",   sz_k_tau[0], sz_k_tau[-1]) # First and last element
    c_k_tau -= c_k_tau[0]   # Ensures K(0) = 0 by fixing the integration constant
    sz_k_tau -= sz_k_tau[0] # Ensures K(0) = 0 by fixing the integration constant
    f = open("K_tau.dat","w")
    for n in range(len(c_k_tau)):
        #f.write("%e %e %e  %e %e  \n"%(n * beta / N, c_k_tau[n], c_kp_tau[n], sz_k_tau[n], sz_kp_tau[n]))
        f.write("%e %e %e  \n"%(n * beta / N, c_k_tau[n], c_kp_tau[n]))
    f.close()
