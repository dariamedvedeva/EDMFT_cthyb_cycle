# # # # # # # # # # # # # # # # # # # # # # # #
#         Medvedeva D.S. 2018 - 2020          #
# For questions: medvedeva.ds@gmail.com       #
# Skype: daryacooper                          #
# # # # # # # # # # # # # # # # # # # # # # # #
print("\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #")
print("#    SELF-CONSISTENCY CYCLE FOR EDMFT, BASED ON CT-HYB SEGMENT    # ")
print("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n")
print("Developed by Medvedeva D. 2018 - 2020")
print("For any questions: medvedeva.ds@gmail.com  ")
print("Skype: daryacooper")
print("Information about calculation will be saved in log.file(a).\n")
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n\n")

import iteration_cycle
import fourier
import retarded_function
import ct_hyb
import numpy as np
import sys
import os
import time
import files_and_folders as tmp
import parameters

def get_date_time():
    t = time.localtime()
    current_time =  time.strftime("%Y-%m-%d %H:%M:%S", t)
    return current_time

#############################################
#            TUNE OF EXECUTION              #
#############################################
# You can preset 2 ways to execute: for local execution and on server.
server          = False
num_omp_threads = 1
type_of_calc    = "dmft"

#############################################
#                SET PATHES                 #
#            >  See pathes.py  <            #
#############################################
if server:
    path_to_exec_file, num_mpi_threads, path_to_maxent = parameters.get_server_run()
else:
    path_to_exec_file, num_mpi_threads, path_to_maxent = parameters.get_local_run()

#################################################
#        SPECIFICATION OF THE MODEL             #
# (1) Copy the file set_parameters_change.py    #
# (2) Change the name as parameters.py      #
# (3) Set parameters                            #
#################################################
parameters.set_model_parameters()
lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm, \
sweeps, time_limit, delta_mix, lambda_mix, number_of_max_it, start_it =  parameters.get_model_parameters()

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


start_time = time.time()
log_file = open("log.file", "a")
log_file.write("\n\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
log_file.write("\nThe calculations was started at {}\n".format(get_date_time()))
log_file.write("from {} till {} iterations.\n".format(start_it, number_of_max_it))
log_file.write("\nPARAMETERS\n")
log_file.write("lattice_type:\t" + str(lattice_type) + "\n")
log_file.write("beta =\t" + str(beta)+ "\n")
log_file.write("U =\t" + str(U)+ "\n")
log_file.write("mu Anderson =\t" + str(hartree_shift)+ "\n")
log_file.write("Nk =\t" + str(Nk)+ "\n")
log_file.write("NN =\t" + str(num_of_neighbours)+ "\n")
log_file.write("Hopping =\t{}\n". format(t))
log_file.write("Coulomb =\t{}\n".format(Coulomb))
log_file.write("mu lattice =\t{}\n".format(Coulomb))
if (particle_hole_symm == 1):
    log_file.write("Particle - hole symmetry - yes.\n")
else:
    log_file.write("Particle - hole symmetry - no.\n")
log_file.write("sweeps (N_MEAS = 2000) =\t{}\n".format(sweeps))
log_file.write("time_limit  =\t{} hours\n".format(time_limit/3600))
log_file.write("Mixing\n")
log_file.write("Delta mixing parameter  =\t{}\n".format(delta_mix))
log_file.write("Lambda mixing parameter  =\t{}\n\n".format(lambda_mix))
log_file.close()
# - - - - - - - - - - - - -
# start_it, number_of_iterations you can tune in parameters.py
for iteration in range(start_it, number_of_max_it, 1):
    print (" ")
    print ("+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +")
    print ("\t\t       ITERATION NUMBER ", str(iteration), "")
    print ("+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +")
    print (" ")

    log_file = open("log.file", "a")
    log_file.write("Iteration # {} is started at {}\n".format(iteration, get_date_time()))
    log_file.close()

    #-------------------------------------------------------#
    #               1. Delta(w) -> Delta (\tau)             #
    #-------------------------------------------------------#

    file = 'Delta.dat'
    tmp.check_delta_file_exist(file)
    fourier.compute(beta, number_of_discrete_tau_points,
                    number_of_fermionic_freqs, number_of_fermionic_freqs_for_fourier, 'Delta')

    #------------------------------------------------------#
    #               2. Lambda(w) -> Lambda(\tau)           #
    #------------------------------------------------------#
    if(type_of_calc != "dmft"):
        lambda_file_name = 'Phi.dat'
        retarded_function.compute_function(number_of_bosonic_frequencies,
                                           number_of_discrete_tau_points, beta, lambda_file_name)

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
        else:
             # run solver if executable file exists
            print ("\n+ - - - - - - - - - - - - - - - +")
            print ("!   Run CT-HYB SEGMENT solver   !")
            print ("+ - - - - - - - - - - - - - - - +\n")
        log_file = open("log.file", "a")
        log_file.write("CT-HYB was started on {} threads at {}.\n".format(num_mpi_threads, get_date_time()))
        log_file.close()
        ct_hyb.run_ct_hyb(path_to_exec_file, num_omp_threads, num_mpi_threads, beta, U, hartree_shift,
                          number_of_bosonic_frequencies, number_of_discrete_tau_points,
                          number_of_fermionic_freqs, type_of_calc, sweeps, time_limit)
    else:
        print("Error! Check the path to the solver execution file. \(^o^)/ ")
        sys.exit()

    #-------------------------------------------------------#
    #                 4. New Delta function                 #
    #-------------------------------------------------------#
    # 1. Gloc
    iteration_cycle.Gloc(mu, Nk, t, lattice_type, U)
    # 2. Delta
    mixing_parameter = delta_mix
    iteration_cycle.new_delta(mixing_parameter)

    #-------------------------------------------------------#
    #               5. New Lambda function                  #
    #-------------------------------------------------------#
    if (type_of_calc == "edmft"):
        # 1. Xloc
        interaction = Coulomb
        iteration_cycle.X_loc(beta, interaction, Nk, lattice_type)
        # 2. Lambda
        mixing_parameter = lambda_mix
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
    tmp.prepare_files_for_new_it(type_of_calc, iteration)

    print("Time for one iteration = {} min".format(np.round((time.time() - start_time)/60),2))
    
    parameters.save_param_file(lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours,
                               t, Coulomb, mu, particle_hole_symm, sweeps, time_limit, delta_mix, lambda_mix)

    log_file = open("log.file", "a")
    log_file.write("Time for the iteration = {} min".format(np.round((time.time() - start_time)/60),2))
    log_file.write("Iteration was finished at {}\n".format( get_date_time()))
    log_file.write("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    log_file.close()

print ("************ Calculation is finished. ************")
print("Time = {} min".format(np.round((time.time() - start_time)/60),2))
print ("************         *(~.~)*          ************")

log_file = open("log.file", "a")
log_file.write("Time for the calculation = {} min".format(np.round((time.time() - start_time)/60),2))
log_file.write("Calculation was finished at {}\n".format(get_date_time()))
log_file.close()
