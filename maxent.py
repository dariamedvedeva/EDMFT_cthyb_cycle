#!/opt/intel/intelpython2/bin/python
import numpy as np
import h5py
import discrete_fourier
import iteration_cycle as it_cyc
import subprocess

def create_input_file_maxent(beta, data_file_name_for_maxent, num_of_data_points, PARTICLE_HOLE_SYMMETRY, min_w, max_w, max_iters_for_fitting, NORM):
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
    if (PARTICLE_HOLE_SYMMETRY == 0):
        file_with_parameters.write(str(num_of_data_points*2))
    if (PARTICLE_HOLE_SYMMETRY == 1):
        file_with_parameters.write(str(num_of_data_points))
    file_with_parameters.write("\n")

    # num of output frequencies
    file_with_parameters.write("NFREQ=")
    file_with_parameters.write(str(2001))
    file_with_parameters.write("\n")

    # G(iw)
    file_with_parameters.write("DATASPACE=frequency")
    file_with_parameters.write("\n")

    #fermionic|bosonic values
    file_with_parameters.write("KERNEL=fermionic")
    file_with_parameters.write("\n")

    # location of data file
    file_with_parameters.write("DATA=")
    file_with_parameters.write(data_file_name_for_maxent)
    file_with_parameters.write("\n")

    # Minimum frequency
    file_with_parameters.write("OMEGA_MIN=")
    file_with_parameters.write(str(min_w))
    file_with_parameters.write("\n")

    # Maximum frequency
    file_with_parameters.write("OMEGA_MAX=")
    file_with_parameters.write(str(max_w))
    file_with_parameters.write("\n")
    
    # Type of frequency grid (default value: Lorentzian)
    file_with_parameters.write("FREQUENCY_GRID=Quadratic")
    file_with_parameters.write("\n")

    # log_min for log grid (default value: 0.0001)
    file_with_parameters.write("LOG_MIN=")
    file_with_parameters.write(str(0.0001))
    file_with_parameters.write("\n")

    # Default model for entropy (default value: flat) "Gaussian"
    file_with_parameters.write("DEFAULT_MODEL=\"double Gaussian\"")
    file_with_parameters.write("\n")

    # stddev - For Gaussian models
    file_with_parameters.write("SIGMA=0.5")
    file_with_parameters.write("\n")
    
#    file_with_parameters.write("GAMMA=0.5")
#    file_with_parameters.write("\n")

    # shift of a model (default value: 0)
#    file_with_parameters.write("SHIFT=5.5")
    file_with_parameters.write("SHIFT=1.0")
    file_with_parameters.write("\n")

    # Maximum Iterations for the fitting routine (default value: 1000)
    file_with_parameters.write("MAX_IT=")
    file_with_parameters.write(str(max_iters_for_fitting))
    file_with_parameters.write("\n")

    # Number of alpha samples (default value: 60)
    file_with_parameters.write("N_ALPHA=100")
    file_with_parameters.write("\n")
    
    file_with_parameters.write("ALPHA_MIN=0.005")
    file_with_parameters.write("\n")
    
    file_with_parameters.write("ALPHA_MAX=5")
    file_with_parameters.write("\n")

    # true to print verbose output (default value: false)
    file_with_parameters.write("VERBOSE=0")
    file_with_parameters.write("\n")
    
#    file_with_parameters.write("NORM=0.374355")
    if(NORM > 0.0):
        file_with_parameters.write("NORM=")
        file_with_parameters.write(str(NORM))
        file_with_parameters.write("\n")
    
    file_with_parameters.close()

def construct_Gw_file_for_maxent(beta, filename):
    # D = np.loadtxt(filename)
    #w       = D[:,0]
    #Re_Gw   = (D[:,1] + D[:,3]) / 2.
    #Im_Gw   = (D[:,2] + D[:,4]) / 2.
    
    f = h5py.File("sim.h5",'r')
    
    Re_Gw_0 = f['simulation/results/gw_re_0/mean/value'][()]
    Im_Gw_0 = f['simulation/results/gw_im_0/mean/value'][()]
    
    Re_Gw_1 = f['simulation/results/gw_re_1/mean/value'][()]
    Im_Gw_1 = f['simulation/results/gw_im_1/mean/value'][()]
    
    #Re_Error_0 = f['simulation/results/gw_re_0/mean/error'].value
    #Im_Error_0 = f['simulation/results/gw_im_0/mean/error'].value
    
    #Re_Error_1 = f['simulation/results/gw_re_1/mean/error'].value
    #Im_Error_1 = f['simulation/results/gw_im_1/mean/error'].value
   
    Re_Error_0 = np.zeros(Re_Gw_0.shape, np.float)
    Im_Error_0 = np.zeros(Re_Gw_0.shape, np.float)
    Re_Error_1 = np.zeros(Re_Gw_0.shape, np.float)
    Im_Error_1 = np.zeros(Re_Gw_0.shape, np.float)
    
    static_error = 0.01
    for i in range(len(Re_Error_0)):
        Re_Error_0[i] = static_error
        Im_Error_0[i] = static_error
        Re_Error_1[i] = static_error
        Im_Error_1[i] = static_error
    
    
    
    Re_Gw    = (Re_Gw_0 + Re_Gw_1) / 2.
    Re_Error = (Re_Error_0 + Re_Error_1) / 2.
    Im_Gw    = (Im_Gw_0 + Im_Gw_1) / 2.
    Im_Error = (Im_Error_0 + Im_Error_1) / 2.

    w = []
    for i in range(len(Re_Gw)):
        w.append(discrete_fourier.fermionic_mats(i, beta))
    
    np.savetxt(filename, np.column_stack((w, Re_Gw, Re_Error, Im_Gw, Im_Error)))

def construct_Gloc_file_for_maxent(input_file, output_file, beta):
    Gloc = it_cyc.read_function(input_file)
    w = []
    for i in range(len(Gloc)):
        w.append(discrete_fourier.fermionic_mats(i, beta))

    Re_Error = np.zeros(Gloc.shape, np.float)
    Im_Error = np.zeros(Gloc.shape, np.float)

    static_error = 0.01
    for i in range(len(Re_Error)):
        Re_Error[i] = static_error
        Im_Error[i] = static_error

    np.savetxt(output_file, np.column_stack((w, Gloc.real, Re_Error, Gloc.imag, Im_Error)))

def run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting, NORM):
    
    ###
    # 1. Create input file "in.param" for MaxEnt
    ###
    create_input_file_maxent(beta, filename_for_maxent, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting, NORM)
    
    if (local):
        ###
        # 2.2 Construct Gloc for MaxEnt
        ###
        construct_Gloc_file_for_maxent("G_loc.dat", filename_for_maxent, beta)
    else:
        ###
        # 2.1 Construct Gw for MaxEnt
        ###
        construct_Gw_file_for_maxent(beta, filename_for_maxent)
    
    ###
    # 3. Construct Gw for MaxEnt
    ###
    subprocess.call([path_to_maxent, "in.param"])
