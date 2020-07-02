
import sys
import numpy as np
from pylab import *
from scipy.optimize import least_squares
import math
import sys


class LambdaMin:
  omega     = []
  X0        = []
  X         = []
  function  = []
  bath      = 0
  Beta      = 0

  def __init__(self, data_file, bath, iwmax, params, verb=False):
    # open Lambda file
    f = open(data_file, "r")
    file_lines = f.readlines()
    self.function.append([])
    # read data
    for i in range(min(len(file_lines), iwmax)):
      l = file_lines[i]
      line_data = l.split()
      self.omega.append(float(line_data[0]))
#      self.Lambda.append(float(line_data[1]) + 1j * float(line_data[2]))
      self.function[0].append(float(line_data[1]))
      
    self.bath = bath
    # initialization of bath states
    # first we have bath number of w0 then bath number of W
    self.X0.append([])
    self.X.append([])
#    tmp = bath_line.split()
    for i in range(len(params)):
       self.X0[0].append(float(params[i]))

  def lambda_model(self, t, coeffs):
    Lambda = 0.0
    for i in range(self.bath):
        Lambda = Lambda - 2.0*coeffs[i]*coeffs[self.bath+i]*coeffs[self.bath+i]/(t*t + coeffs[i]*coeffs[i])
    return Lambda

  def residuals(self, coeffs, y, t):
    return y - self.lambda_model(t, coeffs)

  def minimize(self):
    t = np.array(self.omega, dtype=float)
    y = np.array(self.function[0], dtype=float)
    x0 = np.array(self.X0[0][:2*self.bath], dtype=float)
    print(x0)
    x = least_squares(self.residuals, x0, args=(y, t), max_nfev=100000, bounds=(0, 40)).x
    self.X[0] = x
    return self.X[0]
