import subprocess
import numpy as np
import h5py
import os
import sys
import discrete_fourier
import shutil
import smooth_function
from scipy.interpolate import UnivariateSpline
from pylab import *
from scipy.optimize import curve_fit
import delta_min

global sqrt32, sqrt3
sqrt32 = np.sqrt(3.)/2.
sqrt3 = np.sqrt(3.)
global triang_points

def set_num_of_tr_points(number):
    global triang_points
    triang_points = number

def get_num_of_tr_points():
    global triang_points
    return triang_points

def read_G_imp_frequencies(filename):
    # read file
    D = np.loadtxt(filename)
    frequencies = 1.j*D[:,0]
    function = D[:,1] + 1j * D[:,2]
    return frequencies, function
    
def read_freq_function(filename):
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
#    print("vol = ", vol)
    
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
#        a1 = [0.866025, -0.5, 0.0]
#        a2 = [0.866025, 0.5,  0.0]
        a1 = [ 1.0, 0.0, 0.0]
        a2 = [ -0.5 , 0.866025 , 0.0]
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
#        print("kstep_x = " , str(kstep_x), "\nkstep_y = " , str(kstep_y))
              
    if lattice_type == 'triangular':
        # compute step in the grid
        kstep_x = (2. * (4.0 * np.pi / 3. ))  / (Nk - 1)
        kstep_y = (2. * (2. * np.pi / sqrt3)) / (Nk - 1)
#        print("kstep_x = " , str(kstep_x), "\nkstep_y = " , str(kstep_y))
    return kstep_x, kstep_y

def interaction_dispersion(filename, Nk, lattice_type, param):
    
    kpoints = np.array([[[x, y] for x in range(Nk)] for y in range(Nk)]).reshape((Nk * Nk, 2))
    kpoints_coordinates = np.zeros((Nk*Nk, 2), dtype = np.float)
    interaction = np.zeros(Nk * Nk, dtype = np.float)
    
    b1, b2, b3 = get_b_vectors(lattice_type)
    
#    print('b1 = ', b1 )
#    print('b2 = ', b2 )
#    print('b3 = ', b3 )

    kstep_x, kstep_y = get_kstep(lattice_type, Nk)
    
    for k in range(Nk * Nk):
        vec_ = [0.0, 0.0]
        for i in range(2):
            if (lattice_type == "square"):
                kpoints_coordinates[k][i] = kpoints[k][0] * b1[i]/ (Nk - 1) + kpoints[k][1] * b2[i] / (Nk - 1)
            elif (lattice_type == "triangular"):
                kpoints_coordinates[k][i] = kpoints[k][0] * b1[i]/ (Nk - 1) + kpoints[k][1] * b2[i] / (Nk - 1)
    
    dispersion = open(filename + ".dat", "w")
    
    if lattice_type == 'triangular':
        abs_k = sqrt3
        abs_b = 4. * np.pi * sqrt3 / 3.
        
#        x_boundary = -4. * np.pi / 3.
#        y_boundary = -2. * np.pi / sqrt3
# OR
        x_boundary = 0.0
        y_boundary = 0.0
#        print("shift for G point s = ( {}, {} )".format(x_boundary, y_boundary))
    
    m = 0
    i = 0
    for k in kpoints_coordinates:
        
        # SQUARE
        if lattice_type == 'square':
            interaction[i] = -2. * param * np.sum(np.cos(-np.pi + k))
        # TRIANG
        if lattice_type == 'triangular':
            k_x = k[0]
            k_y = k[1]
            x_coord = x_boundary + k_x
            y_coord = y_boundary + k_y
        
            # first
            interaction[i]  = 2. * param[0] * ( np.cos(x_coord) + np.cos(0.5 * x_coord + sqrt32 * y_coord) + np.cos(0.5 * x_coord - sqrt32 * y_coord) )
            # second
            interaction[i] += 2. * param[1] * ( np.cos(sqrt3 * y_coord) + np.cos(3./2. * x_coord + sqrt32 * y_coord) + np.cos(3./2. * x_coord - sqrt32 * y_coord) )
            # third
            interaction[i] += 2. * param[2] * ( np.cos(2. * x_coord) + np.cos(x_coord + sqrt3 * y_coord) + np.cos(x_coord - sqrt3 * y_coord) )
            
#            interaction[i]  = -2. * param[0] * (np.cos(x_coord) + 2. * np.cos(0.5 * x_coord) * np.cos(y_coord * sqrt32))
#            interaction[i] += -2. * param[1] * (np.cos(sqrt3 * y_coord) + 2. * np.cos(3./2. * x_coord) * np.cos(y_coord * sqrt32))
#            interaction[i] += -2. * param[2] * (np.cos(2. * x_coord) + 2. * np.cos(x_coord) * np.cos(sqrt3 * y_coord ))
            m +=1
            
        dispersion.write(str(np.round(k_x, 6)))
        dispersion.write('\t')
        dispersion.write(str(np.round(k_y,6)))
        dispersion.write('\t')
        dispersion.write(str(np.round(interaction[i],6)))
        dispersion.write('\n')
        i += 1
    dispersion.close()
    
    set_num_of_tr_points(m)
    return interaction, kpoints_coordinates
     
def Gloc(mu, Nk, t, lattice_type, U):
   
    global sqrt32, sqrt3
    frequencies, G_imp_omega = read_G_imp_frequencies("Gw.dat")
    Delta = read_function("Delta.dat")
    Sigma = read_function("Sw.dat")
#    print(Sigma)
    
#    if (Sigma[1].imag - Sigma[0].imag > 0.0):
#        metal = True
#    else:
#        metal = False
        
    G_loc       = np.zeros(G_imp_omega.shape, dtype=np.complex64)
    G_exact     = np.zeros(G_imp_omega.shape, dtype=np.complex64)
    Delta_new   = np.zeros(G_imp_omega.shape, dtype=np.complex64)
   
#    if (Delta.Shape() != G_imp_omega.Shape()):
#        print "Error: Array sizes of Delta and Gw(impurity) functions are not equal!"

    # Corrections
#    corrections(frequencies, mu, Delta, Sigma, beta, U)
    
    ################################################################################
    #                              New G_loc
    ################################################################################
    
    t_k, kpoints_coordinates = interaction_dispersion("t_k", Nk, lattice_type, t)
    
#   We don't use + mu in that formula, because mu is implemented in the solver
#   and you get it as a real part of Sw.dat (Self-energy).?
    if (lattice_type == "square"):
        for i in range(t_k.__len__()):
            G_loc += 1.0 / (frequencies + mu - t_k[i] - Sigma)
        G_loc /= Nk ** 2
    elif (lattice_type == "triangular"):
#        Sigma = smooth_function.smooth_data(frequencies.imag, Sigma, 400, "Sigma_smooth")
        for i in range(get_num_of_tr_points()):
#            G_loc += 1.0 / (frequencies + mu - t_k[i] - Sigma)
            # because of constant shift
            G_loc += 1.0 / (frequencies + mu - t_k[i] - (Sigma - U/2.))
        G_loc /= get_num_of_tr_points()
    np.savetxt("G_loc.dat", np.column_stack((frequencies.imag, G_loc.real, G_loc.imag)))
    # G_loc_smooth = smooth_function.smooth_data(frequencies.imag, G_loc, 400, "G_loc_smooth")
    return 0
    

#    /*
# Non interacting function
#    print("m = {}".format(get_num_of_tr_points()))
#    for i in range(get_num_of_tr_points()):
#        G_exact += 1.0 / (frequencies + mu - t_k[i])
#    G_exact /= get_num_of_tr_points()
#    np.savetxt("G_0.dat", np.column_stack((frequencies.imag, G_exact.real, G_exact.imag)))
#    */
    
def new_delta(mix_par):
    filename = 'Delta_new.dat'
    print ("**Compute new delta function (write in file {}).".format(filename))

    frequencies, G_imp_omega = read_G_imp_frequencies('Gw.dat')
# Tryings to smooth. All of them works
#    G_imp_omega = smooth_function.smooth_data(frequencies.imag, G_imp_omega, 400, "G_imp_omega")
#    G_imp_omega = smooth_function.reduce_noise(frequencies, G_imp_omega, 15, "Gw_smooth")
    Delta = read_function('Delta.dat')
    G_loc = read_function('G_loc.dat')
    
    print("Mixing:")
    print("{} of old Delta and {} of new Delta".format(int((1.0 - mix_par) * 100),  int(mix_par * 100)))
    intermixed_part = (1. / G_imp_omega - 1. / G_loc)
    Delta_new = (1.0 - mix_par) * Delta + mix_par * intermixed_part
    np.savetxt(filename, np.column_stack((frequencies.imag, Delta_new.real, Delta_new.imag)))

    print("\n+ + + + + + + + + + + + + + + + + + ")
    print("+  Reduce noise in Delta function +")
    print("+ + + + + + + + + + + + + + + + + + \n")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   (1) crutches && bicycles --> works only fro high frequencies
#    print(" >> crutches && bicycles var.")
#    # first i frequencies
#    i = 50
#    chi = (Delta_new[i]*frequencies[i])
#
#    def analytical_approx_re(arg, c0, c2, c4, c6):
#        return c0 - c2/(arg**2) + c4/(arg**4) + c6/(arg**6)
#
#    def analytical_approx_im(arg, c0, c1, c3, c5):
##        return c0 - c1/arg - c3/(arg**3) + c5/(arg**5)
#        return c0 - c1/arg + c3/(arg**3) + c5/(arg**5)
#
#    # real part
#    c_start = [Delta_new.real[0], chi.real, chi.real, chi.real]
#    c, cov = curve_fit(analytical_approx_re, frequencies.imag[:i], Delta_new.real[:i], c_start, maxfev = 500000)
#    print(c)
#    Delta_new.real = analytical_approx_re(frequencies.imag, c[0], c[1], c[2], c[3])
#
#    # imag part
#    g_start = [Delta_new.imag[0], chi.imag, chi.imag, chi.imag]
#    g, cov = curve_fit(analytical_approx_im, frequencies.imag[:i], Delta_new.imag[:i], g_start, maxfev = 500000)
#    print(g)
#    Delta_new.imag = analytical_approx_im(frequencies.imag, g[0], g[1], g[2], g[3])
#
#    filename = 'Delta_new_extrapolation.dat'
#    np.savetxt(filename, np.column_stack((frequencies.imag, Delta_new.real, Delta_new.imag)))
#    print("Reconstructed function was saved in file " + filename + ".\n" )
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
#    # high frequency decomposition --> some strange behaviour of the function can occur
#    print(" >> high frequency decomposition")
#    # first i frequencies
#    num_of_used_freqs = 50
#    bath_size = 5
#    filename_t2g = 'Delta_new.dat'
#    D = delta_min.DeltaMin(filename_t2g, bath_size, num_of_used_freqs)
#    coeffs = D.minimize3("delta") # coeffs[E ... , V ...]
#    coeffs[bath_size : bath_size *2] = abs(coeffs[bath_size : bath_size *2])
#    print("Result values:")
#    print(*np.around(coeffs, decimals=3), sep = ",")
#
#    def high_freq_decomp(E, V, bath_size, freq):
#        # for 5 poles
#        func = np.zeros(len(freq), np.complex64)
#        D = []
#        for i in range(bath_size):
#            D.append(V[i]**2)
#        j = 0
#        D_I = D_II = D_III = D_IV = 0.0
#        for i in range(bath_size):
#            D_I   += D[i]
#            D_II  -= D[i] * E[i]
#            D_III += D[i] * E[i]**2
#            D_IV  += D[i] * E[i]**3
#        for w in freq:
#            func[j] -= 1j * (D_I) * 1/w
#            func[j] += D_II * (1/w)**2
#            func[j] += 1j * D_III * (1/w)**3
#            func[j] += D_IV * (1/w)**4
#            j += 1
#        return func
#
#    decomposition = high_freq_decomp(coeffs[:bath_size], coeffs[bath_size:bath_size*2], bath_size, frequencies.imag)
#    filename = 'Delta_new_decomposition.dat'
#    np.savetxt(filename, np.column_stack((frequencies.imag, decomposition.real, decomposition.imag, )))
#    print("Decomposed function was saved in file " + filename + ".\n" )
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    # minimization --> perfect.
    # Should be done twice.
    print(" >> minimization")
    num_of_used_freqs = 25
    bath_size = 9
    filename = 'Delta_new.dat'
    params = [-5.0,-2.5,-1.5,-0.5,0.0,0.5,1.5,2.5,5.0,0.5,0.1,0.5,0.3,0.5,0.4,0.6,0.72,0.5]
#    params = [-44.569,-50.354,-10.446,-4.818,-0.211,0.0,0.261,6.253,71.805,0.0,0.0,0.0,0.205,0.067,0.072,0.183,0.0,0.81]
    D2 = delta_min.DeltaMin(filename, bath_size, num_of_used_freqs, params)
    coeffs2 = D2.minimize("delta") # coeffs[E ... , V ...]
    coeffs2[bath_size : bath_size *2] = abs(coeffs2[bath_size : bath_size *2])
    print("Result values:")
    print(*np.around(coeffs2, decimals=3), sep = ",")
    
    minimized_function  = D2.delta_model(frequencies.imag, coeffs2)
    filename = 'Delta_new_minimized.dat'
    np.savetxt(filename, np.column_stack((frequencies.imag, minimized_function.real, minimized_function.imag)))
    print("Minimized Delta function was saved into the file " + filename + ".\n" )
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    return 0

def corrections(iw, mu, Delta_omega, Sigma_omega, beta, U):
    # 1. Compute g_omega
    g_omega = 1.0 / (iw + mu - Delta_omega - Sigma_omega)
#    np.savetxt("G_OMEGA.dat", np.column_stack((iw.imag, g_omega.real, g_omega.imag)))
    # 2. High frequency corrections
    X = discrete_fourier.chi(iw, Delta_omega)
    analytical_part = np.zeros((iw.shape),dtype=np.complex)
    analytical_part = X / iw
    g_omega_minus_analytical_part = g_omega - analytical_part
    # 3. Compute g_tau(tau = 0)
    g_tau_0 = 0.0
    for i in range(len(g_omega_minus_analytical_part)):
        g_tau_0 += g_omega_minus_analytical_part[i] / beta
    print ('g (tau = 0) =', g_tau_0)

def ro():
    f = h5py.File("sim.h5",'r+')
    # h5ls sim.h5/simulation/results/density_0
    # List all groups
    #print("Keys: %s" % f.keys())
    #a_group_key = list(f.keys())[0]

    # Get the data
    #list
    #[u'count', u'mean', u'tau', u'timeseries', u'error_bins']
    density_0 = f['simulation/results/density_0/mean/value'][()]
    density_1 = f['simulation/results/density_1/mean/value'][()]
#    print (density_0, density_1)
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
    chiw_spin   = nn_00 - 2. * nn_10 + nn_11
 #   np.savetxt("X_not_0.dat", np.column_stack((freq.imag, chiw_charge, chiw_spin)))
    charge_density, spin_density = ro()
#    print(type(beta), type(charge_density))
    chiw_charge[0] -= beta * charge_density * charge_density
    chiw_spin[0]   -= beta * spin_density * spin_density
    # bosonic_freq = bosonic_frequencies(len(chiw_charge), beta)
    chiw_charge_ = -chiw_charge + 0.0 * 1j
    chiw_spin_   = -chiw_spin   + 0.0 * 1j
    np.savetxt("Xw.dat", np.column_stack((freq.imag, chiw_charge_.real, chiw_charge_.imag, chiw_spin_.real, chiw_spin_.imag)))
    return freq, chiw_charge_, chiw_spin_


def X_loc(beta, interaction, Nk, lattice_type):
    global sqrt32, sqrt3
    filename_chi_loc_v0 = "X_loc_v0.dat"
    # read
    frequencies, X_charge, X_spin = susceptibility("nnw.dat", beta)
    
    V_k, kpoints_coordinates = interaction_dispersion("V_k", Nk, lattice_type, interaction)
    
    # init arrays for local quantities
    lambda_charge, lambda_spin = read_2_functions_in_file("Phi.dat")
    X_charge_loc    = np.zeros(X_charge.shape, dtype=np.complex64)
    X_charge_loc_v0 = np.zeros(V_k.shape, dtype=np.complex64)
    
    chi_loc_v0 = open(filename_chi_loc_v0, "w")
    
    for j in range(X_charge_loc.__len__()):
        for i in range(get_num_of_tr_points()):
            X_charge_loc[j] += 1.0/ (1.0 / X_charge[j] + lambda_charge[j] - V_k[i])
            if (j == 0):
                X_charge_loc_v0[i] = 1.0/ (1.0 / X_charge[j] + lambda_charge[j] - V_k[i])
                chi_loc_v0.write(str(np.round(kpoints_coordinates[i][0], 6)))
                chi_loc_v0.write('\t')
                chi_loc_v0.write(str(np.round(kpoints_coordinates[i][1], 6)))
                chi_loc_v0.write('\t')
                chi_loc_v0.write(str(np.round(X_charge_loc_v0.real[i],6)))
#                chi_loc_v0.write('\t')
#                chi_loc_v0.write(str(np.round(X_charge_loc_v0.imag[i],6)))
                chi_loc_v0.write('\n')
    
    X_charge_loc /= get_num_of_tr_points()
    
    chi_loc_v0.close()
#    np.savetxt("X_loc_v0.dat", np.column_stack((kpoints_coordinates, X_charge_loc.real, X_charge_loc.imag)))
    np.savetxt("X_loc.dat", np.column_stack((frequencies.imag, X_charge_loc.real, X_charge_loc.imag)))
    print("File X_loc.dat is constructed.\n")
    
    
def new_lambda(mix_par):
    Sigma        = read_function("Sw.dat")
    Lambda       = read_function("Phi.dat")
    frequencies, X_charge     = read_freq_function("Xw.dat")
    X_charge_loc = read_function("X_loc.dat")
        
    Lambda_new      = np.zeros(X_charge_loc.shape, dtype=np.complex64)
    Lambda_spin_new = np.zeros(X_charge_loc.shape, dtype=np.complex64)
    
    print(X_charge_loc.__len__(), X_charge.__len__(), Lambda.__len__())
    
    # Check the phase
    if ((Sigma[0] - Sigma[1]) > 0.0):
        metal = True
        print ("Impurity is in the metal phase")
    elif ((Sigma[0] - Sigma[1]) < 0.0):
        metal = False
        print ("Impurity is in the insulator phase")
    elif (Sigma[0] == Sigma[1]):
        print( "No interaction. Sigma == 0.0" )
    else:
        print ("WTF with your Sigma? (New Lambda)")
    
    # New Lambda
#    if (metal):
#        # In metallic regimes a fast convergence at low energies is more important
#        intermixed_part = (X_charge - X_charge_loc)
#    else:
#        # Mainly convergency is at high frequencies
#        intermixed_part = (1. / X_charge - 1. / X_charge_loc)

    intermixed_part = (1. / X_charge - 1. / X_charge_loc)
    
    for j in range(X_charge.__len__()):
        Lambda_new[j] = Lambda[j] + mix_par * intermixed_part[j]

    np.savetxt('Lambda_new.dat', np.column_stack((frequencies.imag, Lambda_new.real, Lambda_new.imag, Lambda_spin_new.real, Lambda_spin_new.imag)))

    Lambda_new_smooth = np.zeros(Lambda_new.shape, np.complex)
    num_of_points = 100
    # the last parameter - type of smooth: 1 - filter, 2 - minimization
    Lambda_new_smooth = smooth_function.smooth_data(frequencies.imag, Lambda_new, num_of_points, 'Lambda_new_smooth', 2)

    print (Lambda_new_smooth[0])
    np.savetxt("Lambda_new_smooth.dat", np.column_stack(
        (frequencies.imag, Lambda_new_smooth.real, Lambda_new_smooth.imag, Lambda_spin_new.real, Lambda_spin_new.imag)))


def rename_files(number_of_iteration, type_of_calc):
    print ("Rename files.")
    
#    print(" Delta.dat \t -> \t Delta_iteration.dat " )
    shutil.copy("Delta.dat", "Delta_" + str(number_of_iteration) + ".dat")
    os.remove("Delta.dat")
    
    # Delta_new.dat -> Delta.dat
#    print("Delta_new_extrapolation.dat \t -> \t Delta.dat ")
    # shutil.copy("Delta_new.dat", "Delta.dat")
    shutil.copy("Delta_new_extrapolation.dat", "Delta.dat")
    

    if (type_of_calc == "edmft"):
#        print ("Phi.dat \t -> \t Phi_iteration.dat ")
        shutil.copy("Phi.dat", "Phi_" + str(number_of_iteration) + ".dat")
        os.remove("Phi.dat")
    
        # Lambda_new.dat + (spin_lambda = 0) -> Phi.dat
        shutil.copy("Lambda_new.dat", "Phi.dat")
