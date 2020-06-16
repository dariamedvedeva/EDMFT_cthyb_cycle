
import numpy as np
from pylab import *
from scipy.optimize import leastsq, minimize


#from scipy.optimize import least_squares
import math
import sys
left_lim = -5.
right_lim = 5.

class DeltaMin:
  Lambda = []
  omega = []
  X0 = []
  X = []
  func = []
  odd = False
  Verbose = False
  random_start_params = False
  chit = False
  params_fr = []

  def __init__(self, filename, bath_size, num_of_freqs, params, verb=False):
    iwmax = num_of_freqs
    self.bath = bath_size
    self.Verbose = verb
#    self.iwmax = 100
    # continuous function file
    f = open(filename, "r")
    file_lines = f.readlines()

    self.omega.clear()
    self.func.clear()

    for i in range(min(len(file_lines), iwmax)):
      l = file_lines[i]
      line_data = l.split()
      self.omega.append(float(line_data[0]))
      self.func.append(float(line_data[1]) + 1.j * float(line_data[2]))

    # create arrays
    self.X0.append([])
    self.X.append([])
    
    # start parameters
#    Epsk = []
#    Vk  = []
#    for i in range(bath_size):
#        # idea: from -5 to 5 with step 10/bath_size
#        step = (right_lim - left_lim) / (bath_size - 1)
#        Epsk.append(left_lim + i* step)
#        Vk.append(0.3)
#    params  = Epsk + Vk


    for i in range(len(params)):
      self.X0[0].append(float(params[i]))

  def delta_model(self, t, coeffs):
    Delta = 0.0
    for i in range(self.bath):
      Delta = Delta + coeffs[i + self.bath] * coeffs[i + self.bath] / (t * 1.j - coeffs[i])
    return Delta

  def delta_model_serias(self, coeffs):
    t = self.omega
    Delta = np.zeros(len(t), np.complex64)
    for omega in range(len(t)):
        Delta[omega] = coeffs[0] + coeffs[1]/(1j*t[omega]) + coeffs[2]/((1j*t[omega])**2) + coeffs[3]/((1j*t[omega])**3) + coeffs[4]/((1j*t[omega])**4) + coeffs[5]/((1j*t[omega])**5) + coeffs[6]/((1j*t[omega])**6)
    return Delta
    
  def delta_model_sprint(self, t, coeffs):
    Delta = coeffs[0] + coeffs[1]/(1j*t) + coeffs[2]/((1j*t)**2) + coeffs[3]/((1j*t)**3) + coeffs[4]/((1j*t)**4) + coeffs[5]/((1j*t)**5) + coeffs[6]/((1j*t)**6)
    return Delta

  def delta_model2(self, coeffs):
    t = self.omega
    Delta = np.zeros(len(t), np.complex64)
    for i in range(self.bath):
      for omega in range(len(t)):
        Delta[omega] = Delta[omega] + coeffs[i + self.bath] * coeffs[i + self.bath] / (t[omega] * 1.j - coeffs[i])
    return Delta
  
  def delta_model_frozen_orbitals(self, coeffs):
    t = self.omega
    num_of_frozen_orbs = int(len(self.params_fr) / 2)
    Delta = np.zeros(len(t), np.complex64)
    for omega in range(len(t)):
      for j in range(num_of_frozen_orbs):
        # j - frozen orbitals
        Delta[omega] = Delta[omega] + self.params_fr[j + num_of_frozen_orbs] * self.params_fr[j + num_of_frozen_orbs] / (t[omega] * 1.j - self.params_fr[j])
      for i in range(self.bath):
          # i - changed orbitals
        Delta[omega] = Delta[omega] + coeffs[i + self.bath] * coeffs[i + self.bath] / (t[omega] * 1.j - coeffs[i])
    return Delta
    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #----------------------------------#
  #    R e s i d u a l   v a l u e   #
  #----------------------------------#
   
  def  residuals(self, coeffs, y, t):
    a = (y - self.delta_model(t, coeffs))
    return a.imag ** 2 + a.real ** 2

  def residuals2(self, coeffs, y, t):
    a = (y - self.delta_model(t, coeffs))
    b = (a.real ** 2 + a.imag **2) / (len(t) + 1)
    return b

  def residuals3(self, coeffs, y, t):
    a = (y - self.delta_model2(coeffs))
    b = (a.real ** 2 + a.imag **2) / (len(t))
    return b

  def residuals3_1(self, coeffs):
    t = self.omega
    a = (self.func - self.delta_model2(coeffs))
    b = (a.real ** 2 + a.imag **2) / (np.sqrt(t))
    c = sum(b)
    return c

  def residuals3_2(self, coeffs):
    t = self.omega
    a = (self.func - self.delta_model_frozen_orbitals(coeffs))
    b = (a.real ** 2 + a.imag **2) / (np.sqrt(t))
    c = sum(b)
    return c

  def residuals4(self, coeffs):
    t = self.omega
    a = (self.func - self.delta_model2(coeffs))
    b = (a.real ** 2 + a.imag **2) / t
    c = sum(b)
    return c

  def residuals5(self, coeffs):
    t = self.omega
    a = (self.func - self.delta_model2(coeffs))
    b = (a.real ** 2 + a.imag **2) / (np.power(t, 0.25))
    c = sum(b)
    return c

  def residuals7(self, coeffs):
    t = self.omega
    a = (self.func - self.delta_model_serias(coeffs))
    b = (a.real ** 2 + a.imag **2)
    c = sum(b)
    return c

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  
  #----------------------------------#
  #      M i n i m i z a t i o n     #
  #----------------------------------#
  def minimize(self, postfix):
    # minimize only | a.imag ** 2 + a.real ** 2 |
    t = np.array(self.omega, dtype=float)
    y = np.array(self.func, dtype=np.complex)
    x0 = np.array(self.X0[0][:2 * self.bath], dtype=float)
    print ("Start values:\t", x0)
    x, flags = leastsq(self.residuals, x0, args=(y, t), maxfev=700000)
    self.X[0] = x
    name = "Delta_" + postfix + ".dat"
    self.Delta_output(name, t, x)
    return self.X[0]

  def minimize2(self, postfix):
    # N = number of w => use (N + 1)
    t = np.array(self.omega, dtype=float)
    y = np.array(self.func, dtype=np.complex)
    x0 = np.array(self.X0[0][:2 * self.bath], dtype=float)
    print ("Start values:\t", x0)
    x, flags = leastsq(self.residuals2, x0, args=(y, t), maxfev=700000)
    self.X[0] = x
    name = "Delta_" + postfix + ".dat"
    print(name)
    self.Delta_output(name, t, x)
    return self.X[0]

  def minimize3(self, postfix):
    # sqrt( w )
    t = np.array(self.omega, dtype=float)
    y = np.array(self.func, dtype=np.complex)
    x0 = np.array(self.X0[0][:2 * self.bath], dtype=float)
    print ("Start values:\t", x0)
    full_info = minimize(self.residuals3_1, x0)
    result = full_info.x
    self.X[0] = result
    name = "Delta_" + postfix + ".dat"
    print(name)
    self.Delta_output(name, t, result)
    return self.X[0]

  def minimize4(self, postfix):
    # w ^ (1)
    t = np.array(self.omega, dtype=float)
    y = np.array(self.func, dtype=np.complex)
    x0 = np.array(self.X0[0][:2 * self.bath], dtype=float)

    print ("Start values:\t", x0)

    full_info = minimize(self.residuals4, x0)
    result = full_info.x
    self.X[0] = result

    name = "Delta_" + postfix + ".dat"
    print(name)
    self.Delta_output(name, t, result)

    return self.X[0]

  def minimize5(self, postfix):
    # w ^ (1/4)
    t = np.array(self.omega, dtype=float)
    y = np.array(self.func, dtype=np.complex)
    x0 = np.array(self.X0[0][:2 * self.bath], dtype=float)

    print ("Start values:\t", x0)

    full_info = minimize(self.residuals5, x0)
    result = full_info.x
    self.X[0] = result

    name = "Delta_" + postfix + ".dat"
    print(name)
    self.Delta_output(name, t, result)

    return self.X[0]

  def minimize6(self, postfix):
      #frozen peaks
        t = np.array(self.omega, dtype=float)
        y = np.array(self.func, dtype=np.complex)
        x0 = np.array(self.X0[0][:2 * self.bath], dtype=float)
        
        print ("Start values:\t", x0)
        
        full_info = minimize(self.residuals3_2, x0)
        result = full_info.x
        self.X[0] = result
        
        name = "Delta_" + postfix + "_with_frozen_orbitals.dat"
        print(name)
        self.Delta_output(name, t, result)
        
        return self.X[0]

  def minimize7(self, postfix):
      #frozen peaks
        t = np.array(self.omega, dtype=float)
        y = np.array(self.func, dtype=np.complex)
        # 7 init values
        X = y[-1] * t[-1]
        x0 = [y[0], X, X, X, X, X, X]
#        x0 = np.array(self.X0[0][:2 * self.bath], dtype=float)
        
        print ("Start values:\t", x0)
        
        full_info = minimize(self.residuals7, x0)
        result = full_info.x
        self.X[0] = result
        
        name = "Delta_" + postfix + "_serias.dat"
        print(name)
     #   self.Delta_output_s(name, t, result)
        
        return self.X[0]
        
 # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
  def Delta_output(self, filename, omega, ar_of_params):
    f = open(str(filename), "w")
    print(len(omega))
    for t_ in omega:
      f.write("{0:.10f}".format(t_))
      f.write("\t")
      f.write("{0:.10f} {1:.10f}".format(self.delta_model(t_, ar_of_params).real,self.delta_model(t_, ar_of_params).imag))
      f.write("\n")
    f.close()
    
  def Delta_output_s(self, filename, omega, ar_of_params):
      f = open(str(filename), "w")
      print(len(omega))
      for t_ in omega:
        f.write("{0:.10f}".format(t_))
        f.write("\t")
        f.write("{0:.10f} {1:.10f}".format(self.delta_model_sprint(t_, ar_of_params).real,self.delta_model_serias(t_, ar_of_params).imag))
        f.write("\n")
      f.close()
