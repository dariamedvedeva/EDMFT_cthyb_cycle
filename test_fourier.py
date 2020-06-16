# # # # # # # # # # # # # # # # # # # # # # # #
#         Medvedeva D.S. 2018 - 2020          #
# For questions: medvedeva.ds@gmail.com       #
# Skype: daryacooper                          #
# # # # # # # # # # # # # # # # # # # # # # # #

#
# The easiest way to do the fourier transform
#

import discrete_fourier

# discrete_fourier.compute(100.0, "Delta")
number_of_discrete_tau_points           = 4096
#discrete_fourier.compute(40.0, number_of_discrete_tau_points, 1024, 250, "Lambda_new_smooth")
#discrete_fourier.compute(40.0, number_of_discrete_tau_points, 1024, 250, "Delta")
# number of tau points, number of points for fourier?
discrete_fourier.compute(40.0, number_of_discrete_tau_points, 1024, 120, "Delta")
