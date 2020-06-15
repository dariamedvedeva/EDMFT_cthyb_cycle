#!/usr/bin/env python

#import h5py
import numpy as np

import sys
import subprocess

global sqrt32, sqrt3
sqrt32 = np.sqrt(3.)/2.
sqrt3 = np.sqrt(3.)

#
# Impurity solver wrapper
#
class ImpuritySolver(object):
  __param_file__ = "holstein.param"

  #
  # Impurity solver wrapper constructor
  #
  def __init__(self, path, avg, U, xmu, V, t, beta, BosBit, Nev, Nk=64):
    self.path = path
    self.avg = avg
    self.U = U
    self.xmu = xmu
    self.V = V
    self.t = t
    self.beta = beta
    self.BosBit = BosBit
    self.Nev = Nev
    # initialize reciprocal space mesh
    # number of k points in each direction
    self.Nk = Nk
    # compute step in the grid
    kpoints = np.zeros((Nk*Nk, 2), dtype = np.float)
    kpoints_coordinates = np.array([[[x, y] for x in range(Nk)] for y in range(Nk)]).reshape((Nk * Nk, 2))
    # recipr. vectors b1 & b2
    b1 = np.array([-2. * np.pi, 2. * np.pi / sqrt3])
    b2 = np.array([-2. * np.pi, -2. * np.pi / sqrt3])
    
    for i, k in enumerate(kpoints_coordinates):
        kpoints[i] = b1 * k[0] / Nk + b2 * k[1] / Nk
    self.kpoints = kpoints
    #print kpoints
    
    # compute step in the grid
    #kstep_y = (2. * np.pi / sqrt3 * 2.) / Nk
    #kstep_x = (4. * np.pi) / Nk
    #self.kstep_y = (2. * np.pi / sqrt3 * 2.) / Nk
    #self.kstep_x = (4. * np.pi) / Nk
    #self.kstep = 2 * np.pi / Nk
    # fill k point grid
    #self.kpoints = np.array([[[x, y] for x in range(Nk)] for y in range(Nk)]).reshape((Nk * Nk, 2))

  def init_params(self, bath_len, bbath_len):
    f = open(self.__param_file__, "w+")
    f.write("NSITES=" + str(bath_len + 1) + "\n")
    f.write("NBBITS=" + str(self.BosBit) + "\n")
    f.write("NBLEVEL=" + str(bbath_len) + "\n")
    f.write("NSPINS=2" + "\n")
    f.write("INPUT_FILE=input.h5" + "\n")

    f.write("\n[storage]" + "\n")
    f.write("MAX_DIM=2643325" + "\n")
    f.write("MAX_SIZE=80000000" + "\n")
    f.write("EIGENVALUES_ONLY=0" + "\n")

    f.write("\n[arpack]" + "\n")
    f.write("NEV=" + str(self.Nev) + "\n")
    f.write("NCV=35" + "\n")
    f.write("SECTOR=FALSE" + "\n")

    f.write("\n[lanc]" + "\n")
    f.write("BOLTZMANN_CUTOFF=1e-9" + "\n")
    f.write("NOMEGA=8024" + "\n")
    f.write("BETA=" + str(self.beta) + "\n")
    f.write("EMIN=-17.0" + "\n")
    f.write("EMAX=17.0"  + "\n")
    f.close()

  def prepare_input(self, Vk, Epsk, W, w0):
    Eps0 = np.array([0., 0.])
    Nk = len(Vk)
    Ns = len(Eps0) + len(Epsk)
    sectors = np.array([[0, 0], ])
    data = h5py.File("input.h5", "w");
    data.create_dataset("BETA", shape=(), dtype='f', data=self.beta)

    hop_g = data.create_group("sectors")
    hop_g.create_dataset("values", data=sectors)

    bath = data.create_group("Bath")

    if (Epsk.shape != Vk.shape):
      raise "Incorrect shape for Hybridisation and Epsk"
    Epsk_g = bath.create_group("Epsk")
    Epsk_g.create_dataset("values", shape=(len(Epsk), 2,), data=Epsk, dtype=np.float)
    Vk_g = bath.create_group("Vk")
    Vk_g.create_dataset("values", data=np.array(Vk), dtype=np.float)

    w0_g = bath.create_group("w0")
    w0_g.create_dataset("values", shape=(len(w0),), data=w0, dtype=np.float)
    W_g = bath.create_group("W")
    W_g.create_dataset("values", data=np.array(W), dtype=np.float)

    hop_g = data.create_group("Eps0")
    hop_g.create_dataset("values", data=Eps0)
    data.create_dataset("U", shape=(), data=self.U)
    data.create_dataset("AVG", shape=(), data=self.avg)
    data.create_dataset("mu", shape=(), data=self.xmu)
    data.close()
    self.init_params(len(Vk), len(w0))

  def extract_delta(self, Vk, Epsk, W, w0):
    global sqrt32, sqrt3
    # read impurity solver results
    sim = h5py.File("sim.h5", "r");
    # extract occupation number
    avg = sim["results/static_observables/N"][()]
    self.avg = avg
    # extract 1PGF
    G_omega = np.array(sim["results/G_omega/data"][()]).astype(np.float32).view(np.complex64)
    G_omega_r = np.array(sim["results/G_omega_r/data"][()]).astype(np.float32).view(np.complex64)
    f_omega = np.array(sim["results/G_omega/mesh/1/points"][()])
    f_omega_r = np.array(sim["results/G_omega_r/mesh/1/points"][()])
    # extract 2PGF
    Chi_omega = np.array(sim["results/ChiN_omega/data"][()]).astype(np.float32).view(np.complex64)
    b_omega = np.array(sim["results/ChiN_omega/mesh/1/points"][()])
    # init arrays for local quantities
    G_loc = np.zeros(G_omega.shape, dtype=np.complex64)
    G_loc_r = np.zeros(G_omega_r.shape, dtype=np.complex64)
    Delta = np.zeros(G_omega.shape, dtype=np.complex64)
    Delta_r = np.zeros(G_omega_r.shape, dtype=np.complex64)
    Chi_loc = np.zeros(Chi_omega.shape, dtype=np.complex64)
    Lambda = np.zeros(Chi_omega.shape, dtype=np.complex64)
    # Compute fermionic hybridization function // imag frequencies
    for i in range(len(G_omega)):
      for j in range(len(Vk)):
        for s in range(2):
          Delta[i, :, s, :] += Vk[j, s] * Vk[j, s] / (f_omega[i] * 1.j - Epsk[j, s])
    # Compute fermionic hybridization function // real frequencies
    for i in range(len(G_omega_r)):
      for j in range(len(Vk)):
        for s in range(2):
          Delta_r[i, :, s, :] += Vk[j, s] * Vk[j, s] / (f_omega_r[i] * 1.j - Epsk[j, s])      
    # compute local 1PGF
    dispersion = open("t_dispersion.dat", "w")
    for k in self.kpoints:
      k_x = k[0]
      k_y = k[1]
      
      # compute dispersion
      #t_k = -2. * self.t * np.sum(np.cos(-np.pi + k * self.kstep))
      t_01 = self.t[0]
      t_02 = self.t[1]
      t_03 = self.t[2]
      
      # 1st neighbour
      t_k = 2. * t_01 * (np.cos((0.5 * k_x + sqrt32 * k_y)) + np.cos((0.5 * k_x - sqrt32 * k_y)) + np.cos(k_x))
      # 2nd neighbours
      t_k += 2. * t_02 * (np.cos(1.5 * k_x + sqrt32 * k_y) + np.cos(sqrt3 * k_y) + np.cos(1.5 * k_x - sqrt32 * k_y))
      # 3rd neighbours
      t_k += 2. * t_03 * (np.cos(2. * k_x) + np.cos(k_x + sqrt3 * k_y) + np.cos(k_x - sqrt3 * k_y))

      dispersion.write(str(k_x))
      dispersion.write('\t')
      dispersion.write(str(k_y))
      dispersion.write('\t')
      dispersion.write(str(t_k))
      dispersion.write('\n')

      G_loc += 1.0 / (1.0 / G_omega + Delta - t_k)
      G_loc_r += 1.0 / (1.0 / G_omega_r + Delta_r - t_k)
    G_loc /= self.Nk ** 2
    G_loc_r /= self.Nk ** 2
    
    dispersion.close()
    
    # save impurity 1PGF
    np.savetxt("G_imp.dat", np.column_stack([f_omega.real, G_omega[:, 0, 0].real, G_omega[:, 0, 0].imag]))
    # save local 1PGF
    np.savetxt("G_loc.dat", np.column_stack([f_omega.real, G_loc[:, 0, 0].real, G_loc[:, 0, 0].imag]))
    np.savetxt("G_loc_r.dat", np.column_stack([f_omega_r.real, G_loc_r[:, 0, 0].real, G_loc_r[:, 0, 0].imag]))
    # compute and save new fermionic hybridization function
    Delta_new = Delta + 0.5 * (1 / G_omega - 1 / G_loc)
    np.savetxt("delta.dat", np.column_stack(
      [f_omega.real, Delta_new[:, 0, 0, 0].real, Delta_new[:, 0, 0, 0].imag, Delta_new[:, 0, 1, 0].real,
       Delta_new[:, 0, 1, 0].imag]))
    # compute bosonic hybridization function
    for i in range(len(Chi_omega)):
      for j in range(len(W)):
        Lambda[i, :, :] += 2 * w0[j] * (W[j] ** 2) / (-b_omega[i] ** 2 - w0[j] ** 2)
    # compute local 2PGF
    for k in self.kpoints:
        k_x = k[0]
        k_y = k[1]
        # V_q = -2 * self.V * np.sum(np.cos(-np.pi + q * self.kstep))
        interaction_01 = self.V[0]
        interaction_02 = self.V[1]
        interaction_03 = self.V[2]
        # the first order neighbour
        interaction_k = 2. * interaction_01 * (np.cos((0.5 * k_x + sqrt32 * k_y)) + np.cos((0.5 * k_x - sqrt32 * k_y)) + np.cos(k_x))
        # the second order neighbours
        interaction_k += 2. * interaction_02 * (np.cos(1.5 * k_x + sqrt32 * k_y) + np.cos(sqrt3 * k_y) + np.cos(1.5 * k_x - sqrt32 * k_y))
        # the third order neighbours
        interaction_k += 2. * interaction_03 * (np.cos(2. * k_x) + np.cos(k_x + sqrt3 * k_y) + np.cos(k_x - sqrt3 * k_y))
        Chi_loc += 1.0 / (1.0 / Chi_omega + Lambda - interaction_k)
    Chi_loc /= self.Nk ** 2
    # save impurity 2PGF
    np.savetxt("Chi_imp.dat", np.column_stack([b_omega.real, Chi_omega[:, 0, 0].real, Chi_omega[:, 0, 0].imag])) 
    # save local 2PGF
    np.savetxt("Chi_loc.dat", np.column_stack([b_omega.real, Chi_loc[:, 0, 0].real, Chi_loc[:, 0, 0].imag]))
    # compute and save new bosonic hybridization function
    Lambda_new = Lambda + ((Chi_loc - Chi_omega) / (Chi_loc * Chi_omega))
    np.savetxt("lambda.dat", np.column_stack([b_omega.real, Lambda_new[:, 0, 0].real, Lambda_new[:, 0, 0].imag]))

  def solve(self, Vk, Epsk, W, w0, realGF=False):
    self.prepare_input(Vk, Epsk, W, w0)
    subprocess.call(self.path + " " + self.__param_file__ + " --REALGF={}".format(realGF).upper(), shell=True)
    self.extract_delta(Vk, Epsk, W, w0)
