# # # # # # # # # # # # # # # # # # # # # # # #
#         Medvedeva D.S. 2018 - 2020          #
# For questions: medvedeva.ds@gmail.com       #
# Skype: daryacooper                          #
# # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import subprocess
import os.path
import sys
import os
import time
from scipy import integrate
import set_parameters
import h5py
import discrete_fourier
import iteration_cycle as it_cyc
import subprocess

def read_real_function(filename):
    # read file
    D = np.loadtxt(filename)
    argument = D[:,0]
    function = D[:,1]
    return argument, function
    
def create_input_file_maxent(beta, data_file_name_for_maxent, num_of_data_points, PARTICLE_HOLE_SYMMETRY, min_w, max_w, max_iters_for_fitting):
    file_with_parameters = open("in.param", "w")
    # inverse temperature
    file_with_parameters.write("BETA=")
    file_with_parameters.write(str(beta))
    file_with_parameters.write("\n")

    # 0 || 1
    file_with_parameters.write("PARTICLE_HOLE_SYMMETRY=")
    file_with_parameters.write(str(bool(PARTICLE_HOLE_SYMMETRY)))
    file_with_parameters.write("\n")
    
    # num of data points
    file_with_parameters.write("NDAT=")
    file_with_parameters.write(str(num_of_data_points))
    file_with_parameters.write("\n")

    # num of output frequencies
    file_with_parameters.write("NFREQ=")
    file_with_parameters.write(str(1000))
    file_with_parameters.write("\n")

    # G(iw)
    file_with_parameters.write("DATASPACE=frequency")
    file_with_parameters.write("\n")

    #fermionic|bosonic values
    file_with_parameters.write("KERNEL=bosonic")
    file_with_parameters.write("\n")

    # location of data file
    file_with_parameters.write("DATA=")
    file_with_parameters.write(data_file_name_for_maxent)
    file_with_parameters.write("\n")

    # Minimum frequency
#    file_with_parameters.write("OMEGA_MIN=")
#    file_with_parameters.write(str(min_w))
#    file_with_parameters.write("\n")

    # Maximum frequency
#    file_with_parameters.write("OMEGA_MAX=")
#    file_with_parameters.write(str(max_w))
#    file_with_parameters.write("\n")
    
    # Type of frequency grid (default value: Lorentzian) Quadratic
    file_with_parameters.write("FREQUENCY_GRID=Lorentzian")
    file_with_parameters.write("\n")

    # log_min for log grid (default value: 0.0001)
#    file_with_parameters.write("LOG_MIN=")
#    file_with_parameters.write(str(0.0001))
#    file_with_parameters.write("\n")

    # Default model for entropy (default value: flat) "Gaussian"
#    file_with_parameters.write("DEFAULT_MODEL=\"double Gaussian\"")     # <== For Susceptibility
#    file_with_parameters.write("DEFAULT_MODEL=\"Gaussian\"")           # <== For DOS
#    file_with_parameters.write("\n")

    # stddev - For Gaussian models
#    file_with_parameters.write("SIGMA=0.5")
#    file_with_parameters.write("\n")
    

    # shift of a model (default value: 0)
#    file_with_parameters.write("SHIFT=5.5")
#    file_with_parameters.write("SHIFT=1.0")
#    file_with_parameters.write("\n")

    # Maximum Iterations for the fitting routine (default value: 1000)
    file_with_parameters.write("MAX_IT=")
    file_with_parameters.write(str(max_iters_for_fitting))
    file_with_parameters.write("\n")

    # Number of alpha samples (default value: 60)
#    file_with_parameters.write("N_ALPHA=100")
#    file_with_parameters.write("\n")
    
#    file_with_parameters.write("ALPHA_MIN=0.005")
#    file_with_parameters.write("\n")
#
#    file_with_parameters.write("ALPHA_MAX=5")
#    file_with_parameters.write("\n")

    # true to print verbose output (default value: false)
    file_with_parameters.write("VERBOSE=0")
    file_with_parameters.write("\n")
    
#    file_with_parameters.write("NORM=0.374355")
#    if(NORM > 0.0):
    file_with_parameters.write("NORM=")
    file_with_parameters.write(str(NORM))
    file_with_parameters.write("\n")
    
    file_with_parameters.close()

def construct_Xw_file_for_maxent(beta, filename):
    freq, func = it_cyc.read_freq_function("Xw.dat")
    static_error_re = 0.0001
    static_error_im = 10**(-16)
    error_array_re = np.zeros(func.shape, np.float)
    error_array_im = np.zeros(func.shape, np.float)
    for i in range(len(error_array_re)):
        error_array_re[i] = static_error_re
        error_array_im[i] = static_error_im
        
    np.savetxt(filename, np.column_stack((freq.imag, np.abs(func.real), error_array_re, np.abs(func.imag), error_array_im)))

def run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting):
    
    # in.param
    create_input_file_maxent(beta, filename_for_maxent, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting)
    
    construct_Xw_file_for_maxent(beta, filename_for_maxent)
    
    subprocess.call([path_to_maxent, "in.param"])


#############################################
#                                           #
#            TUNE OF EXECUTION              #
#                                           #
#############################################

server          = False
path_to_maxent  = '/Users/witcher/workspace/CT_HYB_SEGMENT/Maxent/build2/maxent'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
#    C T - H Y B    S E G M E N T   C A L C U L A T I O N     #
#                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

print (" ")
print ("\t++++++++++++++++++++")
print ("\t    MAXENT  ")
print ("\t++++++++++++++++++++")
print (" ")

                # - - - - - - - - - - - - - - - - - #
                #  (1) S E T  P A R A M E T E R S   #
                # - - - - - - - - - - - - - - - - - #
lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm = \
set_parameters.set_model_parameters()
# run CT-HYB SEGMENT solver
# off and make a Fouriet transform into the \tau - space by the first frequencies with smooth data.
number_of_bosonic_frequencies           = 52
start_time = time.time()

filename_for_maxent = 'Xw_for_maxent.dat'
# local = True  --> use Gloc.dat (local    Green's function)
# local = False --> use Gw.dat   (impurity Green's function)
local = False
min_w = -10.0
max_w = 10.0
max_iterations_for_fitting = 4000

                            # - - - - - -  #
                            #  (2) R U N   #
                            # - - - - - -  #

NORM = 0.14509 # Look at output and put it here. I am still not sure how it works
run(path_to_maxent, beta, filename_for_maxent, local, number_of_bosonic_frequencies, particle_hole_symm, min_w, max_w, max_iterations_for_fitting)
os.system("mv in.out.maxspec.dat in.out.maxspec_Xw.dat")

                # - - - - - - - - - - - - - - - - - -  #
                #  (3) C O N S T R U C T X_LOC(w, q)   #
                # - - - - - - - - - - - - - - - - - -  #

print ("************ Calculation is finished. ************")
print("Time = {} min".format(np.round((time.time() - start_time)/60),2))
print ("************         *(~.~)*          ************")
