#!/opt/intel/intelpython2/bin/python
import numpy as np
import subprocess
import maxent
import os.path
import sys
import os
import time
from scipy import integrate

def read_real_function(filename):
    # read file
    D = np.loadtxt(filename)
    argument = D[:,0]
    function = D[:,1]
    return argument, function

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

# run CT-HYB SEGMENT solver
number_of_fermionic_freqs               = 1024
number_of_fermionic_freqs_for_fourier   = 512   # Because of noise we cut the tail of Delta (fermionic hybr. function)
# off and make a Fouriet transform into the \tau - space by the first frequencies with smooth data.
number_of_bosonic_frequencies           = 1024
number_of_discrete_tau_points           = 4096  # Friedrich - 4096


number_of_iterations = 1
start_time = time.time()

for iteration in range(0, number_of_iterations, 1):
    print (" ")
    print ("++++++++++++++++++++++++++")
    print ("ITERATION NUMBER ", str(iteration))
    print ("++++++++++++++++++++++++++")
    print (" ")
  
    # Maxent for Gw.
    filename_for_maxent = 'Gw_for_maxent.dat'
    local = False
    min_w = -5.0
    max_w = 5.0
    max_iterations_for_fitting = 1000000
     # +++++++++++++++++++++ #
    NORM = 1.0 # look at output
    maxent.run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting, NORM)

    maxent_filename = 'in.out.maxspec.dat'
    if (os.path.exists(maxent_filename)):
        w, dos = iteration_cycle.read_real_function(maxent_filename)
    else:
        print("File >> {} << doesn't exist.".format(maxent_filename))
        sys.exit()

    integral = np.round(integrate.trapz(dos, w), 4)
    print("Integral dos = {}".format(integral))
    if (integral == np.round(1.0, 4)):
        print("DOS is normalized")
    else:
        dos /= integral
        # after normalization
        integral2 = np.round(integrate.trapz(dos, w), 4)
        print("Integral dos = {} after normalization".format(integral2))
    os.system("mv in.out.maxspec.dat in.out.maxspec_Gw.dat ")

    print("Time for one iteration {} min".format(np.round((time.time() - start_time)/60),2))
        
print ("************ Calculation is finished. ************")
print("Time = {} min".format(np.round((time.time() - start_time)/60),2))
print ("************         *(~.~)*          ************")
