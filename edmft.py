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

#############################################
#                                           #
#            TUNE OF EXECUTION              #
#                                           #
#############################################

server          = False
num_omp_threads = 1

type_of_calc    = "edmft"

main_cycle      = True     # Only solver execution and self-consistency
maxent_run      = False     # MaxEnt execution to get DoS of the real frequencies

#############################################
#                                           #
#               PATHES IF .exe              #
#                                           #
#############################################

path_to_maxent  = '/Users/witcher/workspace/CT_HYB_SEGMENT/Maxent/build2/maxent'

if server:
    # server
    path_to_exec_file = '/storage/praha1/home/medvedeva/workspace/other_codes/CT_HYB_SEGMENT/CT-HYB-SEGMENT/build/alps_cthyb'
    num_mpi_threads = 42
else:
    # local
    path_to_exec_file = '/Users/witcher/workspace/CT_HYB_SEGMENT/CT-HYB-SEGMENT/build/alps_cthyb'
    num_mpi_threads = 4

#############################################
#                                           #
#        SPECIFICATION OF A MODEL           #
#                                           #
#############################################

lattice_type        = 'square' # write square || triangular
beta                = 100.     # inversive temperature as \beta = 1./T
U                   = 3.0      # local (inter-site) Hubbard repulsion
mu                  = U / 2.0  # for a half filling U / 2. In the tutorial it is written
# that mu = U/2 isn't implemented, but it is (!!!). Automatically mu = U/2, for half-filling.
# The sign problem can occure away from half-filling. Don't touch.
hartree_shift       = 0.0      # Hartree shift (\mu in ct-hyb). Should be 0.0, because it is
# the chemical potential, but it id
Nk                  = 64       # num. of kpoints in each direction, 64 is better for qmc (Friedrich K.)
num_of_neighbours   = 3
particle_hole_symm  = 1

#############################################
#                                           #
#      STATIC PARAMETERS OF A MODEL         #
#                                           #
#############################################
# t         - value of a hopping integral
# Coulomb   - value of non-local (intra-site) Coulomb interaction.
# In papers it figurates as V.

if lattice_type == 'square':
    t       = -0.25
    Coulomb = 0.5
elif lattice_type == 'triangular':
    t    = np.empty(num_of_neighbours, dtype=np.float)

    t[0] =  -0.233     # C2F
    t[1] =  0.00
    t[2] =  0.00

    Coulomb     = np.empty(num_of_neighbours, dtype=np.float)
    Coulomb[0]  = 2.46
    Coulomb[1]  = 0.0
    Coulomb[2]  = 0.0

print ("Lattice type is ", lattice_type)
print ("5*beta*U/(2pi) ~ ", str(int(5*beta*U/(2.*np.pi))))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
#    C T - H Y B    S E G M E N T   C A L C U L A T I O N     #
#                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# run CT-HYB SEGMENT solver
number_of_fermionic_freqs               = 128
number_of_fermionic_freqs_for_fourier   = 120   # Because of noise we cut the tail of Delta (fermionic hybr. function) off and make a Fouriet transform into the \tau - space by the first frequencies with smooth data.
number_of_bosonic_frequencies           = 64
number_of_discrete_tau_points           = 4096  # Friedrich - 4096


number_of_iterations = 1
if (main_cycle):
    for iteration in range(0, number_of_iterations, 1):
        print (" ")
        print ("++++++++++++++++++++++++++")
        print ("ITERATION NUMBER ", str(iteration))
        print ("++++++++++++++++++++++++++")
        print (" ")
#-------------------------------------------------------#
#         0. Initial  hybridization functions.          #
#-------------------------------------------------------#
        if iteration == 0:
            
        # square
            Ek  = [-0.5, 0.0, 0.5]
            Vk  = [0.32, 0.32, 0.32]
            w0p = [1.1, 0.5]
            Wp  = [0.75, 0.35]
            
            initial_hybr_functions.construct_hybr_functions(beta, number_of_fermionic_freqs, number_of_bosonic_frequencies, U, Ek, Vk, w0p, Wp)
        
#-------------------------------------------------------#
#               1. Delta(w) -> Delta (\tau)             #
#-------------------------------------------------------#
        discrete_fourier.compute(beta, number_of_discrete_tau_points, number_of_fermionic_freqs, number_of_fermionic_freqs_for_fourier, 'Delta')
        
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
        file = 'Delta.dat'                                                                              #
        if (os.path.exists(file)):                                                                      #
            print("File " + str(file) + " is constructed.\n")                                           #
        else:                                                                                           #
            print("There is no file " + str(file) + ". The CT-HYB needs Delta_tau.dat to run DMFT.")    #
            sys.exit()                                                                                  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
        
#------------------------------------------------------#
#               2. Lambda(w) -> Lambda(\tau)           #
#------------------------------------------------------#
        if(type_of_calc != "dmft"):
            if iteration == 0:
                # Firstly constructed file will have no noise.
                lambda_file_name = 'Phi.dat'
            else:
                # Use file with reduced noise after other iterations.
                lambda_file_name = "Lambda_new_smooth.dat"
                if (os.path.exists(lambda_file_name)):
                    print("Calculation K_tau.dat from " + str(lambda_file_name) + " with removed noise  ")
                else:
                    print("There is no file ", str(lambda_file_name))
                    sys.exit()
                    
            retarded_function.compute_function(number_of_bosonic_frequencies, number_of_discrete_tau_points, beta, lambda_file_name)
            
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
            file = 'K_tau.dat'                                                                          #
            if (os.path.exists(file)):                                                                  #
                print("File " + str(file) + " is constructed.\n")                                       #
            else:                                                                                       #
                print("There is no file " + str(file) +                                                 #
                ". The CT-HYB needs Delta_tau.dat and K_tau.dat to run EDMFT.")                         #
                sys.exit()                                                                              #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
        else:
            print ("DMFT Calculation --> we don't constract Phi.dat and K_tau.dat files.")
        
#-------------------------------------------------------#
#                       3. Solver                       #
#-------------------------------------------------------#
        if (os.path.exists(path_to_exec_file)):
            if (hartree_shift != 0.0):
                print("\n > > >  W A R N I N G  < < <")
                print("mu = U/2 is implemented in this code and is calculated automatically. ")
                print("You are trying to calculate somethig away from half-filling, that means the sign problem can occure.")
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
        iteration_cycle.new_delta(mu, Nk, t, lattice_type, beta, U)

#-------------------------------------------------------#
#               5. New Lambda function                  #
#-------------------------------------------------------#
        interaction = Coulomb
        iteration_cycle.new_lambda(beta, interaction, Nk, lattice_type)

#-------------------------------------------------------#
#             6. rename files for new_iteration         #
#-------------------------------------------------------#
        iteration_cycle.rename_files(iteration)

#-------------------------------------------------------#
#                    7. plot                            #
#-------------------------------------------------------#
        os.system("./plot.sh")
        os.system("mv plot.pdf plot_{}.pdf".format(str(iteration)))
        
#-------------------------------------------------------#
#             6. Copy the results into folder           #
#-------------------------------------------------------#
        folder = "results_iteration_" + str(iteration)
        os.system("mkdir " + folder)
        os.system("cp ./*.dat " + folder )
        os.system("rm -fr *.dat")
        os.system("cp " + folder + "/Phi.dat")
        os.system("cp " + folder + "/Delta.dat")
#        os.system("cp " + folder + "/Lambda.dat")
        os.system("cp " + folder + "/Lambda_new_smooth.dat")

        
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

    min_w = -7.0
    max_w = 7.0
    max_iterations_for_fitting = 100000
    # +++++++++++++++++++++ #
    maxent.run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting)

print ("************ Calculation is finished. ************")
print ("************         *(~.~)*          ************")
