#!/opt/intel/intelpython2/bin/python
# edmft.py
import iteration_cycle
import discrete_fourier
import retarded_function
import ct_hyb
import numpy as np
import subprocess
import maxent
import initial_hybr_functions
import os.path
import sys
import os
import time
import temptale as tmp
import pathes

#############################################
#                                           #
#            TUNE OF EXECUTION              #
#                                           #
#############################################

server          = False
num_omp_threads = 1

type_of_calc    = "dmft"

main_cycle          = True      # Only solver execution and self-consistency
maxent_run          = False     # MaxEnt execution to get DoS of the real frequencies

#############################################
#                                           #
#               PATHES IF .exe              #
#                                           #
#############################################

if server:
    path_to_exec_file, num_mpi_threads, path_to_maxent = pathes.get_server_run()
else:
    path_to_exec_file, num_mpi_threads, path_to_maxent = pathes.get_local_run()

#############################################
#                                           #
#        SPECIFICATION OF A MODEL           #
#                                           #
#############################################
lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm = \
    set_parameter.set_model_parameters()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
#    C T - H Y B    S E G M E N T   C A L C U L A T I O N     #
#                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# run CT-HYB SEGMENT solver
number_of_fermionic_freqs               = 1024
number_of_fermionic_freqs_for_fourier   = 512   # Because of noise we cut the tail of Delta (fermionic hybr. function)
# off and make a Fouriet transform into the \tau - space by the first frequencies with smooth data.
number_of_bosonic_frequencies           = 64
number_of_discrete_tau_points           = 4096  # Friedrich - 4096


start_time = time.time()
       
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
#            M A X E N T   C A L C U L A T I O N              #
#                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if (maxent_run):
    # doen't work properly because of parameters
    local = False
    
    if (local):
        filename_for_maxent = 'Gloc_for_maxent.dat'
    else:
        filename_for_maxent = 'Gw_for_maxent.dat'

    min_w = -10.0
    max_w = 10.0
    max_iterations_for_fitting = 100000
    # +++++++++++++++++++++ #
    maxent.run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting)

print ("************ Calculation is finished. ************")
print("Time = {} min".format(np.round((time.time() - start_time)/60),2))
print ("************         *(~.~)*          ************")
