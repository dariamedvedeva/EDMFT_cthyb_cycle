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
import set_parameters

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
    set_parameters.set_model_parameters()


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

    #-------------------------------------------------------#
    #               1. Delta(w) -> Delta (\tau)             #
    #-------------------------------------------------------#

    file = 'Delta.dat'
    tmp.check_delta_file_exist(file)
    discrete_fourier.compute(beta, number_of_discrete_tau_points, number_of_fermionic_freqs, number_of_fermionic_freqs_for_fourier, 'Delta')

    #------------------------------------------------------#
    #               2. Lambda(w) -> Lambda(\tau)           #
    #------------------------------------------------------#
    if(type_of_calc != "dmft"):
        lambda_file_name = 'Phi.dat'
        retarded_function.compute_function(number_of_bosonic_frequencies, number_of_discrete_tau_points, beta, lambda_file_name)

        file = 'K_tau.dat'
        tmp.check_ktau_file_exists('K_tau.dat')
    else:
        print ("DMFT Calculation --> we don't constract Phi.dat and K_tau.dat files.")

    #-------------------------------------------------------#
    #                       3. Solver                       #
    #-------------------------------------------------------#
    if (os.path.exists(path_to_exec_file)):
        if (hartree_shift != 0.0):
            print("\n > > >  W A R N I N G  < < <")
            print("mu = U/2 is implemented in this code and is calculated automatically. ")
            print("You are trying to calculate somethig away from half-filling, that means the sign problem "
                  "can occure.")
            sys.exit()
        else:
             # run solver if executable file exists
            print("CT-QMC calculation began.")
        ct_hyb.run_ct_hyb(path_to_exec_file, num_omp_threads, num_mpi_threads, beta, U, hartree_shift, number_of_bosonic_frequencies, number_of_discrete_tau_points, number_of_fermionic_freqs)
    else:
        print("Error! Check the path to the solver execution file. \(^o^)/ ")
        sys.exit()

    #-------------------------------------------------------#
    #                 4. New Delta function                 #
    #-------------------------------------------------------#
    # 1. Gloc
    iteration_cycle.Gloc(mu, Nk, t, lattice_type, U)
    # 2. Delta
    mixing_parameter = 0.25
    iteration_cycle.new_delta(mixing_parameter)

    #-------------------------------------------------------#
    #               5. New Lambda function                  #
    #-------------------------------------------------------#
    if (type_of_calc == "edmft"):
        # 1. Xloc
        interaction = Coulomb
        iteration_cycle.X_loc(beta, interaction, Nk, lattice_type)
        # 2. Lambda
        mixing_parameter = 0.25
        iteration_cycle.new_lambda(mixing_parameter)
    else:
        print("New Lambda function is not calculated.")

    #-------------------------------------------------------#
    #                    6. plot                            #
    #-------------------------------------------------------#
    os.system("./plot.sh")
    os.system("mv plot.pdf plot_{}.pdf".format(str(iteration)))

    #-------------------------------------------------------#
    #             7. Copy the results into folder           #
    #-------------------------------------------------------#
    tmp.create_dir_with_files(type_of_calc,iteration)

    #-------------------------------------------------------#
    #             8. rename files for new_iteration         #
    #-------------------------------------------------------#
#    iteration_cycle.rename_files(str(iteration), type_of_calc)

    #-------------------------------------------------------#
    #             9. rename files for new_iteration         #
    #-------------------------------------------------------#
    tmp.prepare_files_for_new_it(type_of_calc, create_dir_with_files)


    print("Time for one iteration {} min".format(np.round((time.time() - start_time)/60),2))
        
print ("************ Calculation is finished. ************")
print("Time = {} min".format(np.round((time.time() - start_time)/60),2))
print ("************         *(~.~)*          ************")
