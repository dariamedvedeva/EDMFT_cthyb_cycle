# # # # # # # # # # # # # # # # # # # # # # # #
#         Medvedeva D.S. 2018 - 2020          #
#                    ~ ~ ~                    #
#      Smooth / Minimized Lambda function     #
# # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
from scipy.signal import lfilter
import lambda_min

# - - - - - - - - - - - - - - # - - - - - - - - - - - - - - # - - - - - - - - - - - - - - - - #
#  Smooth Lambda_new function by filtfilt procedure.                                          #
#  Will not work in case, if the system is close to the Charge ordering phase.                #
#  Close to CO it is better to use minimization by large number of peaks.                     #
# - - - - - - - - - - - - - - # - - - - - - - - - - - - - - # - - - - - - - - - - - - - - - - #

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filtfilt(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def reduce_noise(x, y, n, filename_to_save):
#    documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html
#    n = 15  # the larger n is, the smoother curve will be
    a = 1
#    n = 15
    b = [1.0 / n] * n
    y_filtred_real = filtfilt(b,a,y.real)
    y_filtred_imag = filtfilt(b,a,y.imag)
    np.savetxt(filename_to_save + ".dat", np.column_stack((x.imag, y_filtred_real.real, y_filtred_imag.real)))
    func = y_filtred_real + 1j * y_filtred_imag
    return func

def smooth_data(x, function, cutoff, filename_to_save, type_of_smooth):
    if (type_of_smooth == 1):
        # filter
        real = function.real
        imag = function.imag
        fs = 5000
        smooth_real = butter_lowpass_filtfilt(real, cutoff, fs)
        smooth_imag = butter_lowpass_filtfilt(imag, cutoff, fs)
        np.savetxt(filename_to_save + ".dat", np.column_stack((x, smooth_real, smooth_imag)))
        result_function = smooth_real + 1j * smooth_imag
    elif (type_of_smooth == 2):
        # minimization
        print(" >> minimization")
        num_of_freqs = 25
        bath_size = 2
        params = [1.0, 2.0, 0.5, 0.5]
        L = lambda_min.LambdaMin('Lambda_new.dat', bath_size, num_of_freqs, params, verb=False)
        Bose_bath = L.minimize()
        print("Result values:")
        print(*np.around(Bose_bath, decimals=3), sep = ",")
        minimized_function  = L.lambda_model(x, Bose_bath)
        np.savetxt(filename_to_save + ".dat", np.column_stack((x, minimized_function.real, minimized_function.imag)))
        print("Minimized Delta function was saved into the file " + filename_to_save + ".\n" )
        result_function = minimized_function
    return result_function
