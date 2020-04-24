import numpy as np
import discrete_fourier
import os

def init_delta(Ek, Vk, number_of_freqs, beta):

    Delta   = np.zeros(number_of_freqs, np.complex)
    w       = np.zeros(number_of_freqs, np.complex)
    
    if(len(Ek) == len(Vk)):
        K = len(Ek)
    else:
        print ("The lenghts of Ek and Vk are different.")
        os.abort()

    for ind in range(number_of_freqs):
        w[ind] = discrete_fourier.fermionic_mats(ind, beta) * 1j
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

    for ind in range(number_of_freqs):
        v[ind] = discrete_fourier.bosonic_mats(ind, beta) * 1j
        for orbital_index in range(P):
            Lambda[ind] += 2. * (np.abs(Wp[orbital_index])**2 * w0p[orbital_index]) / ((v[ind]**2 - w0p[orbital_index]**2))

    return v, Lambda


def rewrite_Phi_file_for_cthyb_solver(freq, lambda_charge):
    n = len(lambda_charge)
    lambda_spin = np.zeros(n, np.complex)
    
    np.savetxt("Phi.dat", np.column_stack((freq.imag, lambda_charge.real, lambda_charge.imag, lambda_spin.real, lambda_spin.imag)))
    np.savetxt("Phi_0.dat", np.column_stack((freq.imag, lambda_charge.real, lambda_charge.imag, lambda_spin.real, lambda_spin.imag)))


def construct_hybr_functions(beta, ferm, bos, U, Ek, Vk, w0p, Wp):

    frequencies, initial_delta =  init_delta(Ek, Vk, ferm, beta)
    np.savetxt("Delta.dat", np.column_stack((frequencies.imag, initial_delta.real, initial_delta.imag)))
    np.savetxt("Delta_0.dat", np.column_stack((frequencies.imag, initial_delta.real, initial_delta.imag)))

    frequencies, initial_lambda = init_lambda(w0p, Wp, bos, beta)
    rewrite_Phi_file_for_cthyb_solver(frequencies, initial_lambda)
