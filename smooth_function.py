# # # # # # # # # # # # # # # # # # # # # # # #
#         Medvedeva D.S. 2018 - 2020          #
#                                             #
#               Smooth function               #
# # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
from scipy.signal import lfilter

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filtfilt(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def smooth_data(x, function, cutoff, filename_to_save):
    real = function.real
    imag = function.imag
    fs = 5000
    smooth_real = butter_lowpass_filtfilt(real, cutoff, fs)
    smooth_imag = butter_lowpass_filtfilt(imag, cutoff, fs)
#    print (smooth_real[0])
    np.savetxt(filename_to_save + ".dat", np.column_stack((x, smooth_real, smooth_imag)))
    result_function = smooth_real + 1j * smooth_imag
    return result_function

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
