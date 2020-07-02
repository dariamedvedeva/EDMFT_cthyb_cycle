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
    set_parameters.set_model_parameters()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
#    C T - H Y B    S E G M E N T   C A L C U L A T I O N     #
#                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# run CT-HYB SEGMENT solver
number_of_fermionic_freqs               = 1024
number_of_fermionic_freqs_for_fourier   = 1000   # Because of noise we cut the tail of Delta (fermionic hybr. function)
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
    
#    # 1. Gloc
#    iteration_cycle.Gloc(mu, Nk, t, lattice_type, beta, U)
#
#    # 2. maxent Gloc(iw) -> Gloc(w)
#    filename_for_maxent = 'Gloc_for_maxent.dat'
#    local = True
#    min_w = -5.0
#    max_w = 5.0
#    max_iterations_for_fitting = 1000000
#    # +++++++++++++++++++++ #
#    NORM = 1.0 # look at output
#    maxent.run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting, NORM)
#
#    # 3. normalization of Gloc(w)
#    maxent_filename = 'in.out.maxspec.dat'
#    if (os.path.exists(maxent_filename)):
#        w, dos = iteration_cycle.read_real_function(maxent_filename)
#    else:
#        print("File >> {} << doesn't exist.".format(maxent_filename))
#        sys.exit()
#
#    integral = np.round(integrate.trapz(dos, w), 4)
#    print("Integral dos = {}".format(integral))
#    if (integral == np.round(1.0, 4)):
#        print("DOS is normalized")
#    else:
#        dos /= integral
#        # after normalization
#        integral2 = np.round(integrate.trapz(dos, w), 4)
#        print("Integral dos = {} after normalization".format(integral2))
#
#    # 4. mu half filling from Gloc
#    mu_half_filling = set_parameters.get_shift_half_filling(dos, 2*w[-1], abs(w[0]-w[1]))
#    os.system("mv in.out.maxspec.dat in.out.maxspec_G_loc_before.dat ")
#
##    mu_half_filling = 0.351889
##     5. G_loc with half filling
#    iteration_cycle.Gloc(mu_half_filling + mu, Nk, t, lattice_type, beta, U)
#    maxent.run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting, NORM)
#
#    # 6. normalization of Gloc(w)
#    maxent_filename = 'in.out.maxspec.dat'
#    if (os.path.exists(maxent_filename)):
#        w, dos = iteration_cycle.read_real_function(maxent_filename)
#    else:
#        print("File >> {} << doesn't exist.".format(maxent_filename))
#        sys.exit()
#
#    integral = np.round(integrate.trapz(dos, w), 4)
#    print("Integral dos = {}".format(integral))
#    if (integral == np.round(1.0, 4)):
#        print("DOS is normalized")
#    else:
#        dos /= integral
#        # after normalization
#        integral2 = np.round(integrate.trapz(dos, w), 4)
#        print("Integral dos = {} after normalization".format(integral2))
#
#    # 7. mu half filling from Gloc
##    mu_half_filling = set_parameters.get_shift_half_filling(dos, 2*w[-1], abs(w[0]-w[1]))
#    os.system("mv in.out.maxspec.dat in.out.maxspec_G_loc_before.dat ")
  
#    # 8. Maxent for Gw.
#    filename_for_maxent = 'Gw_for_maxent.dat'
#    local = False
#    min_w = -5.0
#    max_w = 5.0
#    max_iterations_for_fitting = 1000000
#     # +++++++++++++++++++++ #
#    NORM = 1.0 # look at output
#    maxent.run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting, NORM)
#
#    maxent_filename = 'in.out.maxspec.dat'
#    if (os.path.exists(maxent_filename)):
#        w, dos = iteration_cycle.read_real_function(maxent_filename)
#    else:
#        print("File >> {} << doesn't exist.".format(maxent_filename))
#        sys.exit()

#    integral = np.round(integrate.trapz(dos, w), 4)
#    print("Integral dos = {}".format(integral))
#    if (integral == np.round(1.0, 4)):
#        print("DOS is normalized")
#    else:
#        dos /= integral
#        # after normalization
#        integral2 = np.round(integrate.trapz(dos, w), 4)
#        print("Integral dos = {} after normalization".format(integral2))
#    os.system("mv in.out.maxspec.dat in.out.maxspec_Gw.dat ")

    #-------------------------------------------------------#
    #                 4. New Delta function                 #
    #-------------------------------------------------------#
    iteration_cycle.Gloc(mu, Nk, t, lattice_type, U)
    mixing_parameter = 0.05
#    iteration_cycle.new_delta(mixing_parameter)
    
    #-------------------------------------------------------#
    #               5. New Lambda function                  #
    #-------------------------------------------------------#
    if (type_of_calc == "edmft"):
        interaction = Coulomb
        iteration_cycle.X_loc(beta, interaction, Nk, lattice_type)
        mixing_parameter = 0.45
        iteration_cycle.new_lambda(mixing_parameter)
    else:
        print("New Lambda function is not calculated.")

    #-------------------------------------------------------#
    #                    6. try fourier                     #
    #-------------------------------------------------------#
    file = 'Delta_new_minimized.dat'
    tmp.check_delta_file_exist(file)
    discrete_fourier.compute(beta, number_of_discrete_tau_points, number_of_fermionic_freqs, number_of_fermionic_freqs_for_fourier, 'Delta_new_minimized')
     
    lambda_file_name = 'Lambda_new_smooth.dat'
    retarded_function.compute_function(number_of_bosonic_frequencies, number_of_discrete_tau_points, beta, lambda_file_name)


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
