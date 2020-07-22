#!/opt/intel/intelpython2/bin/python
# edmft.py
import iteration_cycle
import fourier
import retarded_function
import ct_hyb
import numpy as np
import subprocess
import maxent
import construct_discrete_hybr_functions
import os.path
import sys
import os
import time
import files_and_folders as tmp
import pathes
import parameters
from scipy import integrate

#############################################
#                                           #
#            TUNE OF EXECUTION              #
#                                           #
#############################################

server          = False
num_omp_threads = 1
type_of_calc    = "edmft"

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
    parameters.set_model_parameters()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
#    C T - H Y B    S E G M E N T   C A L C U L A T I O N     #
#                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# run CT-HYB SEGMENT solver
number_of_fermionic_freqs               = 1024
number_of_fermionic_freqs_for_fourier   = 1024   # Because of noise we cut the tail of Delta (fermionic hybr. function)
# off and make a Fouriet transform into the \tau - space by the first frequencies with smooth data.
number_of_bosonic_frequencies           = 1024
number_of_discrete_tau_points           = 4096  # Friedrich - 4096


number_of_iterations = 1
start_time = time.time()

for iteration in range(0, number_of_iterations, 1):
    print (" ")
    print ("++++++++++++++++++++++++++")
    print ("+  EDMFT - consistency   + ")
    print ("++++++++++++++++++++++++++")
    print (" ")
    
#    #-------------------------------------------------------#
#    #                 4. New Delta function                 #
#    #-------------------------------------------------------#
#    iteration_cycle.Gloc(mu, Nk, t, lattice_type, U)
#    mixing_parameter = 0.25
#    iteration_cycle.new_delta(mixing_parameter)
#
#    #-------------------------------------------------------#
#    #               5. New Lambda function                  #
#    #-------------------------------------------------------#
#    if (type_of_calc == "edmft"):
#        interaction = Coulomb
#        iteration_cycle.X_loc(beta, interaction, Nk, lattice_type)
#        mixing_parameter = 0.1
#        iteration_cycle.new_lambda(mixing_parameter)
#    else:
#        print("New Lambda function is not calculated.")

    #-------------------------------------------------------#
    #                    6. try fourier                     #
    #-------------------------------------------------------#
#    file = 'Delta_new_minimized.dat'
#    tmp.check_delta_file_exist(file)
#    discrete_fourier.compute(beta, number_of_discrete_tau_points, number_of_fermionic_freqs, number_of_fermionic_freqs_for_fourier, 'Delta_new_minimized')
#
#    lambda_file_name = 'Lambda_new_smooth.dat'
#    retarded_function.compute_function(number_of_bosonic_frequencies, number_of_discrete_tau_points, beta, lambda_file_name)


    #-------------------------------------------------------#
    #             6. rename files for new_iteration         #
    #-------------------------------------------------------#
#    iteration_cycle.rename_files(iteration, type_of_calc)

    #-------------------------------------------------------#
    #                    7. plot                            #
    #-------------------------------------------------------#
    os.system("./plot_second_part.sh")
    os.system("mv plot.pdf plot_{}.pdf".format(str('sec_part')))

    #-------------------------------------------------------#
    #             6. Copy the results into folder           #
    #-------------------------------------------------------#
#    tmp.create_dir_with_files(type_of_calc, iteration)


    print("Time for one iteration {} min".format(np.round((time.time() - start_time)/60),2))
        
print ("************ Calculation is finished. ************")
print("Time = {} min".format(np.round((time.time() - start_time)/60),2))
print ("************         *(~.~)*          ************")
