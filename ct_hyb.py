# # # # # # # # # # # # # # # # # # # # # # # #
#         Medvedeva D.S. 2018 - 2020          #
# For questions: medvedeva.ds@gmail.com       #
# Skype: daryacooper                          #
# # # # # # # # # # # # # # # # # # # # # # # #
import subprocess
import os


def run_ct_hyb(path, num_omp_threads, num_mpi_threads, beta, U, mu, N_W, N,
               number_of_fermionic_freqs, type_of_calculation, sweeps, time_limit):
    
# Simple check of parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Common parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    os.putenv('OMP_NUM_THREADS', str(int(num_omp_threads)))
    path_for_delta_function = 'Delta_tau_ct_hyb.dat'
    arg_mpi = ['mpirun', '-np', str(int(num_mpi_threads))]
    arg = []
    arg.append('--BETA')		   	    # inverse temperature
    arg.append(str(beta))
    arg.append('--L') 			        # Number of Brillouin Zone points. It was 512
    arg.append(str(64))
    arg.append('--U')   			    # interaction value. Only specify if you are not reading an U matrix
    arg.append(str(U))
    arg.append('--MU')			        # chemical potential / orbital energy values
    arg.append(str(mu))
    arg.append('--N')                   # number of imaginary time discretization points
    arg.append(str(N))
    arg.append('--cthyb.N_MEAS')	    # number of updates per measurement
    arg.append('2000')
    arg.append('--cthyb.N_SHIFT')       # number of shift updates
    arg.append('10')
    arg.append('--cthyb.N_SWAP')        # number of shift swap updates
    arg.append('10')
    arg.append('--cthyb.N_LEGENDRE')    # number of legendre coefficients
    arg.append('25')


# SOLVER control parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    arg.append('--FLAVORS')               # number of spin-orbitals (sometimes called flavors)
    arg.append(str(2))
    arg.append('--cthyb.SWEEPS')          # total number of Monte Carlo sweeps to be done 10^9
    arg.append(str(sweeps))
    arg.append('--cthyb.THERMALIZATION')  # thermalization sweeps
    arg.append('100')
    arg.append('--SEED')                  # PRNG seed
    arg.append(str(0))
    arg.append('--MAX_TIME')              # maximum solver runtime. 10800 - 3 hours
    arg.append(str(time_limit))
    arg.append('--cthyb.MEASURE_freq')    # measure in frequency domain
    arg.append(str(1))
    arg.append('--cthyb.MEASURE_nn')
    arg.append(str(1))
    if (type_of_calculation == "edmft"):
        arg.append('--cthyb.MEASURE_nnw')
        arg.append(str(1))
    
    
# DELTA parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    arg.append(str('--NMATSUBARA'))   	    # number of matsubara coefficients
    arg.append(str(number_of_fermionic_freqs))
    arg.append('--Tmax')              	    # maximum time to check if simulation has finished
    arg.append(str(600))
    arg.append('--cthyb.J')            	    # Hund
    arg.append(str(0.0))
    arg.append('--cthyb.DELTA_IN_HDF5') 	# true if hybridization function file is in hdf5 format
    arg.append('0')
    arg.append('--cthyb.DELTA')         	#  path for hybridization function file
    arg.append(path_for_delta_function)


# RETARDED INTERACTION parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (type_of_calculation == "edmft"):
        arg.append(str('--cthyb.N_W'))          # number of bosonic Matsubara frequencies
        arg.append(str(N_W))
        arg.append('--cthyb.K_IN_HDF5')         # set to true if retarded interaction K is stored in hdf5
        arg.append(str(0))
        arg.append(str('--cthyb.RET_INT_K'))    # file with the retarted interaction information. See doc for format.
        arg.append("K_tau.dat")
    
# parameters of DATA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    arg.append('--cthyb.MEASURE_freq') 	    # measure in frequency domain
    arg.append(str(1))
    if (type_of_calculation == "edmft"):
        arg.append('--cthyb.MEASURE_nnw')  	    # meaure density-density correlation functions in frequency domain (susceptibility NN)
        arg.append(str(1))
    arg.append('--cthyb.TEXT_OUTPUT')
    arg.append('1')

# print arg
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    with open("params", "w") as fff:
        iii = 1
        for ar in arg:
            ar = ar.replace('--', '')
            fff.write(ar)
            if iii % 2 == 0 and iii > 0:
                fff.write("\n")
            else:
                fff.write("=")
            iii += 1

# run
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    arg = arg_mpi + [ path, "params"]
    subprocess.call(arg)
