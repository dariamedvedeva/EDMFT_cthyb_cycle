import discrete_fourier
import numpy as np

filename = 'chi_charge.dat'
data = np.loadtxt(filename)
beta = data[-1, 0]
num_of_tau_points  = len(data[:])
num_of_freq_points = int((num_of_tau_points - 1) / 2.)
tau_points = []
func_tau = []

for i in range(num_of_tau_points):
    tau_points.append(float(data[i, 0]))
    func_tau.append(float(data[i, 1]))

freqs = np.zeros(num_of_freq_points, np.complex)
for n in range(num_of_freq_points):
    freqs[n] = 1j*discrete_fourier.bosonic_mats(n, beta)
    
#func_omega = np.fft(func_tau)
func_omega = discrete_fourier.DFT(func_tau, freqs, tau_points)
np.savetxt("chi_charge_omega.dat", np.column_stack((freqs.imag, func_omega.real, func_omega.imag)))
