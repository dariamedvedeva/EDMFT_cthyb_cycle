# Discrete Fourier Transform (DFT)
# Medvedeva, Iskakov

import numpy as np
import sys
import os


Delta = False
Lambda = False

def DFT(function, frequencies, tau):
    # from time to frequency
    N = len(function)
    dt = tau[1] - tau[0]
    FmList = np.zeros(frequencies.shape, dtype=np.complex)
    for iw, w in enumerate(frequencies):
        for it, t in enumerate(tau):
            FmList[iw] += dt * function[it] * np.exp( w * t)
    return FmList

def InverseDFT(function, frequencies, tau, beta):
    # from frequency to time
    N = len(function)
    fnList = np.zeros((tau.shape),dtype=np.complex)
    for it, t in enumerate(tau):
        for iw, w in enumerate(frequencies):
            fnList[it] += function[iw] * np.exp(- w * t)
    return fnList/beta

def fermionic_mats(n, beta):
    return (2 * n + 1) * np.pi / beta

def bosonic_mats(n, beta):
    return 2 * n * np.pi / beta

def read_delta_frequencies(filename, beta, total_number_of_freqs, freqs_no_noise):
    # read file
    global Delta, Lambda
    D = np.loadtxt(filename)
    if D[0,0] == 0.0:
        Lambda = True
    else:
        Delta  = True
     
    # frequencies from file
    freq_right  = 1.j*D[:freqs_no_noise,0]
    freq_left   = -freq_right[::-1]

    # function from file
    func_right  = D[:freqs_no_noise,1] + 1j * D[:freqs_no_noise,2]
    func_left   = func_right[::-1].real - 1j * func_right[::-1].imag

    # const for tail
    const_for_tail = freq_right[-1].imag * func_right[-1].imag

    if (len(freq_right) == len(func_right)):
        number_of_start_freqs = len(freq_right)
    else:
        print ("Lengths of frequencies and a function are different")
        os.exit()
    frequencies_tail    = np.zeros(total_number_of_freqs - freqs_no_noise, np.complex)
    function_tail       = np.zeros(total_number_of_freqs - freqs_no_noise, np.complex)

    i = 0
    for n in range(number_of_start_freqs, total_number_of_freqs, 1):
        frequencies_tail[i]  = fermionic_mats(n, beta) * 1j
        function_tail[i]     = 1j * const_for_tail / (fermionic_mats(n, beta))
        i = i + 1

    function_left_add_part = function_tail[::-1].real - 1j * function_tail[::-1].imag
    frequencies = np.concatenate((freq_left[:], freq_right[:]))
    function    = np.concatenate((func_left[:], func_right[:]))
    return frequencies, function

def chi(frequencies, function):
    # Only for Delta and 1PGF
    chi = function[-1]  * frequencies[-1]
    return chi

def compute(beta, num_of_time_points, number_of_freqs, number_of_freqs_for_fourier, filename):
    global Delta, Lambda
    num_of_time_points += 1
    
    print ("Fourier transform of Delta")
    print ("__________________________")
    
    # +++++++++++++++++ #
    # 1. read file Delta(w) and construct negative frequency part
    # +++++++++++++++++ #
    frequencies, delta_freq = read_delta_frequencies(filename + ".dat", beta, number_of_freqs, number_of_freqs_for_fourier)
    
    print ("The number of w points \t= ", len(frequencies))
    
    # +++++++++++++++++ #
    # 2. tau array
    # +++++++++++++++++ #
    # We need much more tau-points to get proper transform back to matsubara
    print ("The number of tau points = ", num_of_time_points)
    tau = np.array(range(num_of_time_points))*beta/(num_of_time_points-1)
    print ("Time step = ", np.round(tau[1] - tau[0],6))

    # +++++++++++++++++ #
    # 3. FT w -> t
    # +++++++++++++++++ #
    
    # 3.1 FT(Delta - tail)
    X = chi(frequencies, delta_freq)
    print ("Tail const = ", X)
    analytical_part = np.zeros((frequencies.shape),dtype=np.complex)
    if (Delta):
        analytical_part = X / frequencies
#        analytical_part = X / frequencies + X / (frequencies**2)
    if (Lambda):
        analytical_part = 0.0 +1j * 0.0
    delta_minus_chi_part = delta_freq - analytical_part
    print ("FT omega -> tau")
    delta_time_first = InverseDFT(delta_minus_chi_part, frequencies, tau, beta)

    # 3.2 FT
    FT_analytical_part = - X * 0.5
#    FT_analytical_part =  X * (2*tau - beta)/4.0
    # 3.3 Delta(tau) = Delta_minus_anylytical_part(tau) - analytical_part(tau)
    delta_time = delta_time_first + FT_analytical_part
    
    # 3.4 save function
#    np.savetxt(filename + "_tau.dat", np.column_stack((tau, delta_time.real, delta_time.real)))
#    np.savetxt(filename + "_tau2.dat", np.column_stack((tau, delta_time.imag, delta_time.imag)))

    # 3.5 Check that Delta_tau < 0
    error_text = 'Numerical problem: Delta(tau) > 0.0. \n' \
                 'Recommendations: use smaller mixing parameter ||  increase the number of sweeps || increase number of points for minimization || decrease number of frequencies for the fourier transform. \n'
    def positive_check(val, error):
        if (val >= 0.0):
            print(error_text)
            log_file = open("log.file", "a")
            log_file.write(error_text)
            log_file.close()
            os.system("./plot_fourier.sh")
            sys.exit()

    f = open(filename + "_tau_ct_hyb.dat", "w")
    for i in range(len(tau)):
        f.write(str(i))
        f.write(' ')
        f.write(str(delta_time[i].real))
        f.write(' ')
        f.write(str(delta_time[i].real)) #????? 26.06.2020
 #       f.write(str(delta_time[i].imag))
        if(i < len(tau) - 1):
            f.write('\n')
    f.close()
    
    for value in delta_time.real:
        positive_check(value, error_text)

    # +++++++++++++++++ #
    # 4. FT t -> w
    # +++++++++++++++++ #
#    print ("FT tau -> omega")
#    delta_freq2 = DFT(delta_time.real, frequencies, tau)
##    np.savetxt(filename + "_check_FT.dat", np.column_stack((frequencies.imag, delta_freq2.real, delta_freq2.imag)))

def gt_tail(beta, m, tau):
    """returns the Fourier transform G(tau)=(1/BETA)*\sum_{n=-\infty}^\infty e^{-i w_n*tau } dtau
       of G(i*w_n) = 1/(i*w_n)^m;
    //note: the Fourier of 1=1/(i*w_n)^0 is delta(tau)"""
    if(m>9): raise NameError("Fourier transform of the tail not implemented for orders larger than 9")
    elif m==0: raise NameError("cannot handle Fourier transform of 1")
    elif m==1: return -0.5
    elif m==2: return (2*tau - beta)/4.0
    elif m==3: return (beta*tau - tau**2)/4.0
    elif m==4: return (beta**3 - 6.0*beta*tau**2 + 4.0*tau**3)/48.0
    elif m==5: return -tau*(beta**3 - 2.0*beta*tau**2 + tau**3)/48.0
    elif m==6: return -1.0*(beta - 2.*tau)*(beta**2 + beta*tau - tau**2)**2/480.0
    elif m==7: return (3.0*(beta**5)*tau - 5.0*(beta*tau)**3 + 3.0*tau**5 - tau**6)/1440.0
    elif m==8: return (17.0*beta**7 - 84.0*(beta**5)*tau**2 + 70.0*beta**3*tau**4 - 28.0*beta*tau**6 + 8.0*tau**7)/80640.0
    elif m==9: return -tau*(17.0*beta**7 - 28.0*(beta**5)*tau**2 + 14.0*beta**3*tau**4 - 4.0*beta*tau**6 + tau**7)/80640.0
    else: return 0.0
