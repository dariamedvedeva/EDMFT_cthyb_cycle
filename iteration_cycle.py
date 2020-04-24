import subprocess
import numpy as np
import h5py
import os
import discrete_fourier
import shutil
import smooth_function

global sqrt32, sqrt3
sqrt32 = np.sqrt(3.)/2.
sqrt3 = np.sqrt(3.)

def read_G_imp_frequencies(filename):
    # read file
    D = np.loadtxt(filename)
    frequencies = 1.j*D[:,0]
    function = D[:,1] + 1j * D[:,2]
    return frequencies, function

def read_function(filename):
    # read file
    D = np.loadtxt(filename)
    function = D[:,1] + 1j * D[:,2]
    return function

def read_2_functions_in_file(filename):
    # read file
    D = np.loadtxt(filename)
    function1 = D[:,1] + 1j * D[:,2]
    function2 = D[:,3] + 1j * D[:,4]
    return function1, function2

def read_real_function(filename):
    # read file
    D = np.loadtxt(filename)
    argument = D[:,0]
    function = D[:,1]
    return argument, function

def self_energy(freq, G_0, G):
    Sigma = 1./G_0 - 1./G
    np.savetxt("Sigma.dat", np.column_stack((freq.imag, Sigma.real, Sigma.imag)))
    return Sigma 

def g_0(freq, chem_pot, delta):
    g_0 = 1./(freq + chem_pot - delta)
    np.savetxt("G_0.dat", np.column_stack((freq.imag, g_0.real, g_0.imag)))
    return g_0    

def read_correlator_file(filename):
    D = np.loadtxt(filename)
    freq = D[:,0] * 1j
    nnw_00 = D[:,1]
    nnw_10 = D[:,2]
    nnw_11 = D[:,3]
    return freq, nnw_00, nnw_10, nnw_11

def ro():
    f = h5py.File("sim.h5",'r+')
    # h5ls sim.h5/simulation/results/density_0
    # List all groups
    #print("Keys: %s" % f.keys())
    #a_group_key = list(f.keys())[0]

    # Get the data
    #list
    #[u'count', u'mean', u'tau', u'timeseries', u'error_bins']
    density_0 = f['simulation/results/density_0/mean/value'].value
    density_1 = f['simulation/results/density_1/mean/value'].value
    print (density_0, density_1)
    charge_density  = density_0 + density_1
    spin_density    = density_0 - density_1
    return charge_density, spin_density

def bosonic_frequencies(number_of_freqs, beta):
    freqs_array = []
    for i in range(number_of_freqs):
        freqs_array.append(0.0 + (2. * i * np.pi / beta) * 1j)
    return freqs_array

def susceptibility(correlator_filename, beta):
    freq, nn_00, nn_10, nn_11 = read_correlator_file(correlator_filename)
    #chiw_charge = (nn_00 + nn_10) * 2.0
    #chiw_spin  = (nn_00 - nn_10) * 2.0
    chiw_charge = nn_00 + 2. * nn_10 + nn_11
    chiw_spin 	= nn_00 - 2. * nn_10 + nn_11
    np.savetxt("X_not_0.dat", np.column_stack((freq.imag, chiw_charge, chiw_spin)))
    charge_density, spin_density = ro()
    chiw_charge[0]  -= beta * charge_density * charge_density
    chiw_spin[0]    -= beta * spin_density * spin_density
    # bosonic_freq = bosonic_frequencies(len(chiw_charge), beta)
    chiw_charge_ = -chiw_charge + 0.0 * 1j
    chiw_spin_   = -chiw_spin   + 0.0 * 1j
    np.savetxt("Xw.dat", np.column_stack((freq.imag, chiw_charge_.real, chiw_charge_.imag, chiw_spin_.real, chiw_spin_.imag)))
    return freq, chiw_charge_, chiw_spin_

#def delta_from_sigma(mu, Nk, t, lattice_type):
#    global sqrt32, sqrt3
#    frequencies, G_imp_omega = read_G_imp_frequencies("Gw.dat")
#    Delta = read_function("Delta.dat")
#    Sigma = read_function("Sw.dat")
#    G_loc       = np.zeros(G_imp_omega.shape, dtype=np.complex64)
#    Delta_new   = np.zeros(G_imp_omega.shape, dtype=np.complex64)
#
#    print ("**Compute local 1PGF")
#
#
#    G_loc += 1.0/ (frequencies + mu - t_k - Sigma)
#    G_loc /= Nk ** 2
#    np.savetxt("G_loc.dat", np.column_stack((frequencies.imag, G_loc.real, G_loc.imag)))
#
#    G_0 = g_0(frequencies, mu, Delta)
#
#    #####################
#    # New Delta
#    #####################
#    print ("**Compute new delta function")
#    Delta_2 = frequencies + mu - Sigma - 1./G_loc
#    np.savetxt("Delta_from_Sigma.dat", np.column_stack((frequencies.imag, Delta_2.real, Delta_2.imag)))

def retarded_function_ct_hyb(filename):
    time, Lt = read_real_function(filename + ".dat")
    dLdt = numerical_derivative(time, Lt)
    time.pop(len(time) - 1)
    Lt.pop(len(Lt) - 1)
    subprocess.call("rm", filename + "_ct_hyb.dat")
    f = open(filename + "_ct_hyb.dat", "w")
    for i in range(len(time) - 1):
        f.write(str(i))
        f.write(' ')
        f.write(str(Lt[i]))
        f.write(' ')
        f.write(str(dLdt))
        if(i < len(tau) - 1):
            f.write('\n')
    f.close()

def numerical_derivative(x, y):
    dydx = []
    if (len(x) == len(y)):
        for i in range(len(y) - 1):
            value = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            dydx.append(value)
    return dydx

def compute_numerical_derivative(x, y):
    size = len(y)
    res = np.zeros(size, 'd') # 'd' for double
    # centered differences
    for idx in range(1, size-1):
       res[idx] = (y[idx+1] - y[idx-1]) / (x[idx+1] - x[idx-1])
    # one-sided differences
    res[0] = (y[1] - y[0]) / (x[1] - x[0])
    res[-1] = (y[size-1] - y[size-2]) / (x[size-1] - x[size-2])
    return res

def calculate_b_vectors(a1, a2, a3):
    
    vol = a1[0]*(a2[1]*a3[2]-a2[2]*a3[1])- a1[1]*(a2[0]*a3[2]-a3[0]*a2[2]) + a1[2]*(a2[0]*a3[1]-a2[1]*a3[0])
    print("vol = ", vol)
    
#     %%%%%%%%%%%%%%%%%%%%%%%
#     reciprocal space cell
#     %%%%%%%%%%%%%%%%%%%%%%%
    tpi = 2.0 * np.pi
    b1 = [0.0, 0.0, 0.0]
    b2 = [0.0, 0.0, 0.0]
    b3 = [0.0, 0.0, 0.0]
    
    b1[0]=tpi*(a2[1]*a3[2]-a3[1]*a2[2])/vol
    b1[1]=tpi*(a2[2]*a3[0]-a2[0]*a3[2])/vol
    b1[2]=tpi*(a2[0]*a3[1]-a2[1]*a3[0])/vol
    b2[0]=tpi*(a3[1]*a1[2]-a3[2]*a1[1])/vol
    b2[1]=tpi*(a3[2]*a1[0]-a3[0]*a1[2])/vol
    b2[2]=tpi*(a3[0]*a1[1]-a1[0]*a3[1])/vol
    b3[0]=tpi*(a1[1]*a2[2]-a2[1]*a1[2])/vol
    b3[1]=tpi*(a1[2]*a2[0]-a1[0]*a2[2])/vol
    b3[2]=tpi*(a1[0]*a2[1]-a2[0]*a1[1])/vol
    
    return b1, b2, b3

def interaction_dispersion(filename, Nk, lattice_type, param):
    
    kpoints = np.array([[[x, y] for x in range(Nk)] for y in range(Nk)]).reshape((Nk * Nk, 2))
    
    kpoints_coordinates = np.zeros((Nk*Nk, 2), dtype = np.float)
    
    interaction = np.zeros(Nk * Nk, dtype = np.float)
    
    b1, b2, b3 = get_b_vectors(lattice_type)

    print('b1 = ', b1 )
    print('b2 = ', b2 )
    print('b3 = ', b3 )

    kstep_x, kstep_y = get_kstep(lattice_type, Nk )

    for k in range(Nk * Nk):
        vec_ = [0.0, 0.0]
        for i in range(2):
            kpoints_coordinates[k][i] = kpoints[k][0] * b1[i]/ (Nk - 1) + kpoints[k][1] * b2[i] / (Nk - 1)
    
    # print(kpoints_coordinates)
    
    dispersion = open(filename + ".dat", "w")
    
    i = 0
    for k in kpoints:
#        print(k)
        k_x = k[0]
        k_y = k[1]
        
        # SQUARE
        if lattice_type == 'square':
            interaction[i] = -2. * param * np.sum(np.cos(-np.pi + k * kstep_x))
#            interaction[i] = -2. * param * (np.cos(-np.pi + k_x * kstep_x) + np.cos(-np.pi + k_y * kstep_y))
            
        # TRIANG
        if lattice_type == 'triangular':
            interaction_01 = interaction[0]
            interaction_02 = interaction[1]
            interaction_03 = interaction[2]
                # the first order neighbour
            interaction[i] = -2. * param[0] * (np.cos(-np.pi + (k_x) * kstep_x) + 2. * np.cos(-np.pi + (0.5 * k_x) * kstep_x) * np.cos(-np.pi + (k_y * sqrt32) * kstep_y))
                # the second order neighbours
            interaction[i] += -2. * param[1] * (np.cos(-np.pi + (sqrt3 * k_y) * kstep_y) + 2. * np.cos(-np.pi + (1.5 * k_x) * kstep_x) * np.cos(-np.pi + (sqrt32 * k_y) * kstep_y))
                # the third order neighbours
            interaction[i] += -2. * param[2] * (np.cos(-np.pi + (2. * k_x) * kstep_x) + 2. * np.cos(-np.pi + (k_x) * kstep_x) * np.cos(-np.pi + (sqrt3 * k_y) * kstep_y))
            
        dispersion.write(str(np.round(k_x, 6)))
        dispersion.write('\t')
        dispersion.write(str(np.round(k_y,6)))
        dispersion.write('\t')
        dispersion.write(str(np.round(interaction[i],6)))
        dispersion.write('\n')
        i += 1
        
    dispersion.close()
    
    return interaction, kpoints_coordinates
     
def new_delta(mu, Nk, t, lattice_type, beta, U):
    print ("Construct new delta function")
    
    # read
    global sqrt32, sqrt3
    print ("**Reading...")
    frequencies, G_imp_omega = read_G_imp_frequencies("Gw.dat")
    Delta = read_function("Delta.dat")
    Sigma = read_function("Sw.dat")

    # init arrays for local quantities
    G_loc       = np.zeros(G_imp_omega.shape, dtype=np.complex64)
    Delta_new   = np.zeros(G_imp_omega.shape, dtype=np.complex64)
   
    # if (Delta.Shape() != G_imp_omega.Shape()):
    #	print "Error: Array sizes of Delta and Gw(impurity) functions are not equal!"
    
    print ("**Compute local 1PGF...")
    print ("\nFor dimension 2!!\n")

    # Corrections
    corrections(frequencies, mu, Delta, Sigma, beta, U)
    
    t_k, kpoints_coordinates = interaction_dispersion("t_k", Nk, lattice_type, t)

#   We don't use + mu in that formula, because mu is implemented in the solver
#   and you get it as a real part of Sw.dat (Self-energy).
    
    for i in range(t_k.__len__()):
        G_loc += 1.0 / (frequencies + mu - t_k[i] - Sigma)
    G_loc /= Nk ** 2
    
    np.savetxt("G_loc.dat", np.column_stack((frequencies.imag, G_loc.real, G_loc.imag)))
    # G_loc_smooth = smooth_function.smooth_data(frequencies.imag, G_loc, 400, "G_loc_smooth")
    
    ################################################################################
    #                              New Delta
    ################################################################################
    print ("**Compute new delta function")
    mix_par = 0.75
    print("Mixing:")
    print('{} of old Delta and {} of new Delta'.format(int((1.0 - mix_par) * 100),  int(mix_par * 100)))
    intermixed_part = (1. / G_loc - 1. / G_imp_omega)
    Delta_1 = (1.0 - mix_par) * Delta + mix_par * intermixed_part
    np.savetxt("Delta_new.dat", np.column_stack((frequencies.imag, Delta_1.real, Delta_1.imag)))
    print("File Delta_new.dat is created")

def corrections(iw, mu, Delta_omega, Sigma_omega, beta, U):
    # 1. Compute g_omega
    g_omega = 1.0 / (iw + mu - Delta_omega - Sigma_omega)
#    np.savetxt("G_OMEGA.dat", np.column_stack((iw.imag, g_omega.real, g_omega.imag)))
    # 2. High frequency corrections
    X = discrete_fourier.chi(iw, Delta_omega)
    analytical_part = np.zeros((iw.shape),dtype=np.complex)
#    analytical_part = X / iw
    g_omega_minus_analytical_part = g_omega - analytical_part
    # 3. Compute g_tau(tau = 0)
    g_tau_0 = 0.0
    for i in range(len(g_omega_minus_analytical_part)):
        g_tau_0 += g_omega_minus_analytical_part[i] / beta
    print ('g (tau = 0) =', g_tau_0)

def get_b_vectors(lattice_type):
    if lattice_type == 'square':
       # real. vectors a1 & a2 & a3
        a1 = [ 1.0, 0.0 , 0.0]
        a2 = [ 0.0, 1.0 , 0.0]
        a3 = [ 0.0, 0.0 , 10.0]
        # recipr. vectors b1 & b2 & b3
        b1, b2, b3 = calculate_b_vectors(a1, a2, a3)
           
    if lattice_type == 'triangular':
        # real. vectors a1 & a2 & a3
        a1 = [0.866025, -0.5, 0.0]
        a2 = [0.866025, 0.5,  0.0]
        a3 = [0.0, 0.0, 10.0]
        # recipr. vectors b1 & b2 & b3
        b1, b2, b3 = calculate_b_vectors(a1, a2, a3)
    return b1, b2, b3

def get_kstep(lattice_type, Nk):
    # Nk - number of points
    # Nk - 1 - number of "slices" of the Mesh
    b1, b2, b3 = get_b_vectors(lattice_type)
    if lattice_type == 'square':
        # compute step in the grid
        kstep_x = 2.0 * (2. * np.pi) / (Nk - 1)
        kstep_y = 2.0 * (2. * np.pi) / (Nk - 1)
        # kstep_x = b1[0] / (Nk - 1)
        # kstep_y = b2[1] / (Nk - 1)
        # kstep_z = b3[2] / (Nk - 1)
        print("kstep_x = " , str(kstep_x), "\nkstep_y = " , str(kstep_y))
              
    if lattice_type == 'triangular':
        # compute step in the grid
        kstep_y = ((2. * np.pi / sqrt3) * 2.) / (Nk - 1)
        kstep_x = (2. * (2. * np.pi)) / (Nk - 1)
        print("kstep_x = " , str(kstep_x), "\nkstep_y = " , str(kstep_y))
    return kstep_x, kstep_y

def new_lambda(beta, interaction, Nk, lattice_type):
    
    print ("Construct new lambda function")
    
    global sqrt32, sqrt3
    mix_par = 0.75
    print("mix_par = ", str(mix_par))
    
    # read
    Sigma = read_function("Sw.dat")
    Lambda = read_function("Phi.dat")
    frequencies, X_charge, X_spin = susceptibility("nnw.dat", beta)
    lambda_charge, lambda_spin = read_2_functions_in_file("Phi.dat")
    
    # init arrays for local quantities
    X_charge_loc        = np.zeros(Lambda.shape, dtype=np.complex64)
    Lambda_new          = np.zeros(Lambda.shape, dtype=np.complex64)
    Lambda_spin_new     = np.zeros(Lambda.shape, dtype=np.complex64)
    
    print ("**Compute local 2PGF")
    print(interaction)
    V_k, kpoints_coordinates = interaction_dispersion("V_k", Nk, lattice_type, interaction)
   
    for i in range(V_k.__len__()):
        X_charge_loc += 1.0/ (1.0 / X_charge + lambda_charge - V_k[i])
    X_charge_loc /= Nk ** 2
    
    np.savetxt("X_loc.dat", np.column_stack((frequencies.imag, X_charge_loc.real, X_charge_loc.imag)))
    print("File X_loc.dat is constructed.\n")
    
    #####################
    # Check the phase
    #####################
    
    if ((Sigma[0] - Sigma[1]) > 0.0):
        metal = True
        print ("Impurity is in the metal phase")
    elif ((Sigma[0] - Sigma[1]) < 0.0):
        metal = False
        print ("Impurity is in the insulator phase")
    else:
        print ("WTF with your Sigma? (New Lambda)")

    #########################################################################################################
    #                                               New Lambda
    #########################################################################################################
    
    if (metal):
        # In metallic regimes a fast convergence at low energies is more important
        intermixed_part = (X_charge - X_charge_loc)
    else:
        # Mainly convergency is at high frequencies
        intermixed_part = (1. / X_charge - 1. / X_charge_loc)

    Lambda_new = Lambda + mix_par * intermixed_part

    np.savetxt("Lambda_new.dat", np.column_stack((frequencies.imag, Lambda_new.real, Lambda_new.imag, Lambda_spin_new.real, Lambda_spin_new.imag)))

    Lambda_new_smooth = np.zeros(Lambda_new.shape, np.complex)
    Lambda_new_smooth = smooth_function.smooth_data(frequencies.imag, Lambda_new, 10000, "Lambda_new_smooth")

    print (Lambda_new_smooth[0])
    np.savetxt("Lambda_new_smooth.dat", np.column_stack(
        (frequencies.imag, Lambda_new_smooth.real, Lambda_new_smooth.imag, Lambda_spin_new.real, Lambda_spin_new.imag)))


def rename_files(number_of_iteration):
    
    print ("Rename files.")
    
    # Delta.dat     -> Delta_iteration.dat
    shutil.copy("Delta.dat", "Delta_" + str(number_of_iteration) + ".dat")
    os.remove("Delta.dat")
    
    # Delta_new.dat -> Delta.dat
    shutil.copy("Delta_new.dat", "Delta.dat")

    # Phi.dat   -> Phi_iteration.dat
    shutil.copy("Phi.dat", "Phi_" + str(number_of_iteration) + ".dat")
    os.remove("Phi.dat")
    
    # Lambda_new.dat + (spin_lambda = 0) -> Phi.dat
    shutil.copy("Lambda_new.dat", "Phi.dat")
