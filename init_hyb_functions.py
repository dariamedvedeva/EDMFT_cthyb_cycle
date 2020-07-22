def init_delta(Ek, Vk, number_of_freqs, beta):

    Delta   = np.zeros(number_of_freqs, np.complex)
    w       = np.zeros(number_of_freqs, np.complex)

    if(len(Ek) == len(Vk)):
        K = len(Ek)
    else:
        print ("The lenghts of Ek and Vk are different.")
        os.abort()

    for ind in range(number_of_freqs):
        w[ind] = fourier.fermionic_mats(ind, beta) * 1j
        for orbital_index in range(K):
            Delta[ind] += (np.abs(Vk[orbital_index])**2) / (w[ind] - Ek[orbital_index])

    return w, Delta

def init_lambda(w0p, Wp, number_of_freqs, beta):

    Lambda   = np.zeros(number_of_freqs, np.complex)
    v        = np.zeros(number_of_freqs, np.complex)

    if(len(w0p) == len(Wp)):
        P = len(Wp)
    else:
        print ("The lenghts of Ek and Vk are different.")
        os.abort()

    for ind in range(0, number_of_freqs, 1):
        v[ind] = fourier.bosonic_mats(ind, beta) * 1j
        for orbital_index in range(P):
            Lambda[ind] += 2. * (np.abs(Wp[orbital_index])**2 * w0p[orbital_index]) / ((v[ind]**2 - w0p[orbital_index]**2))
    return v, Lambda


def rewrite_Phi_file_for_cthyb_solver(freq, lambda_charge):
    n = len(lambda_charge)
    lambda_spin = np.zeros(n, np.complex)
    np.savetxt("Phi.dat", np.column_stack((freq.imag, lambda_charge.real, lambda_charge.imag, lambda_spin.real, lambda_spin.imag)))
    print("File Phi.dat saved.")


def construct_hybr_functions(beta, ferm, bos, U, Ek, Vk, w0p, Wp):

    frequencies, initial_delta =  init_delta(Ek, Vk, ferm, beta)
    np.savetxt("Delta.dat", np.column_stack((frequencies.imag, initial_delta.real, initial_delta.imag)))
    print("File Delta.dat saved.")

    frequencies2, initial_lambda = init_lambda(w0p, Wp, bos, beta)
    rewrite_Phi_file_for_cthyb_solver(frequencies2, initial_lambda)

import parameters
import numpy as np
import fourier

parameters.set_model_parameters()
lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm, \
sweeps, time_limit, delta_mix, lambda_mix, number_of_iterations, start_from_it \
    = parameters.get_model_parameters()

beta = 100.
# Because of noise we cut the tail of Delta (fermionic hybr. function) off and make a Fouriet transform into the \tau - space by the first frequencies with smooth data.
number_of_fermionic_freqs               = 1024
number_of_fermionic_freqs_for_fourier   = 1024   # Because of noise we cut the tail of Delta (fermionic hybr. function)
# off and make a Fouriet transform into the \tau - space by the first frequencies with smooth data.
number_of_bosonic_frequencies           = 1024
number_of_discrete_tau_points           = 4096  # Friedrich - 4096

Ek  = [-0.5, 0.0, 0.5]
Vk  = [0.32, 0.32, 0.32]
w0p = [1.1, 0.5]
Wp  = [0.75, 0.35]

construct_hybr_functions(beta, number_of_fermionic_freqs, number_of_bosonic_frequencies, U, Ek, Vk, w0p, Wp)