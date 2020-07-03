import initial_hybr_functions
import set_parameters
lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm = \
set_parameters.set_model_parameters()

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
    
initial_hybr_functions.construct_hybr_functions(beta, number_of_fermionic_freqs, number_of_bosonic_frequencies, U, Ek, Vk, w0p, Wp)

