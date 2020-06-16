# # # # # # # # # # # # # # # # # # # # # # # #
#         Medvedeva D.S. 2018 - 2020          #
# For questions: medvedeva.ds@gmail.com       #
# Skype: daryacooper                          #
# # # # # # # # # # # # # # # # # # # # # # # #

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

start_time = time.time()

print (" ")
print ("\t++++++++++++++++++++")
print ("\t    MAXENT  ")
print ("\t++++++++++++++++++++")
print (" ")

                # - - - - - - - - - - - - - - - - - #
                #  (1) S E T  P A R A M E T E R S   #
                # - - - - - - - - - - - - - - - - - #

# Maxent for Gw. The script will create a file 'Gw_for_maxent.dat' and after will use it
# to run maxent.
filename_for_maxent = 'Gw_for_maxent.dat'
# local = True  --> use Gloc.dat (local    Green's function)
# local = False --> use Gw.dat   (impurity Green's function)
local = False
min_w = -5.0
max_w = 5.0
max_iterations_for_fitting = 1000000


                            # - - - - - -  #
                            #  (2) R U N   #
                            # - - - - - -  #

NORM = 1.0 # Look at output and put it here. I am still not sure how it works
maxent.run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting, NORM)


                # - - - - - - - - - - - - - - - - - -  #
                #  (3) D O S  &&  O C C U P A T I O N  #
                # - - - - - - - - - - - - - - - - - -  #

# 'in.out.maxspec.dat' - the name of the file is maxent's output and reserved.
maxent_filename = 'in.out.maxspec.dat'
if (os.path.exists(maxent_filename)):
    w, dos = iteration_cycle.read_real_function(maxent_filename)
else:
    print("File >> {} << doesn't exist.".format(maxent_filename))
    sys.exit()

# I check if DOS is normalized on 1. (For one orbital)
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

        
print ("************ Calculation is finished. ************")
print("Time = {} min".format(np.round((time.time() - start_time)/60),2))
print ("************         *(~.~)*          ************")
