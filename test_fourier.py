import discrete_fourier

# discrete_fourier.compute(100.0, "Delta")
number_of_discrete_tau_points           = 4096
discrete_fourier.compute(100.0, number_of_discrete_tau_points, 1024, 250, "Lambda_new_smooth")
