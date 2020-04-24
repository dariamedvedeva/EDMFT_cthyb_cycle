import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

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
    fs = 500000
    smooth_real = butter_lowpass_filtfilt(real, cutoff, fs)
    smooth_imag = butter_lowpass_filtfilt(imag, cutoff, fs)
    print (smooth_real[0])
    #np.savetxt(filename_to_save + ".dat", np.column_stack((x, smooth_real, smooth_imag)))
    result_function = smooth_real + 1j * smooth_imag
    return result_function
