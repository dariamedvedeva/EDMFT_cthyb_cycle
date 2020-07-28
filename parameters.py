import numpy as np
#############################################
#                                           #
#           MODEL PARAMETERS                #
#                                           #
#############################################

#  Default parameters. To set parameters - change values in set_model_parameters()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lattice_type        = 'triangular'  # write square || triangular
beta                = 1.0           # inversive temperature as \beta = 1./T
U                   = 1.0           # local (inter-site) Hubbard repulsion
hartree_shift       = 0.0           # Hartree shift (\mu in ct-hyb). for a half filling U / 2. In the tutorial it is written
                                    # that mu = U/2 isn't implemented, but it is (!!!).
Nk                  = 64            # num. of kpoints in each direction, 64 is better for qmc (Friedrich K.)
num_of_neighbours   = 3

#  Hubbard model parameters
t    = np.empty(num_of_neighbours, dtype=np.float)
t[0] = t[1] = t[2] = 0.0
Coulomb     = np.empty(num_of_neighbours, dtype=np.float)
Coulomb[0]  = Coulomb[1] = Coulomb[2] = 0.0

#  ct_hub parameters
sweeps     = 10**10          # 10**8 .. 10**10 is enough
hours_max  = 3               # max hours
time_limit = hours_max*60*60 # in seconds

#  mixing_parameters
delta_mix  = 0.1
lambda_mix = 0.6

#  SET PARAMETERS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def set_model_parameters():
    global lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours
    global t, Coulomb, mu, particle_hole_symm, sweeps
    global time_limit, delta_mix, lambda_mix
    global max_it_num, start_from_it
    lattice_type        = 'triangular' # write square || triangular
    beta                = 100.     # inversive temperature as \beta = 1./T
    U                   = 5.16      # local (inter-site) Hubbard repulsion
    #mu                  = U/2.   # for a half filling U / 2. In case of square lattice it should be mu = U/2. !!!!
    #mu                  = 0.8 * t
    hartree_shift       = 0.0      # Hartree shift (\mu in ct-hyb). for a half filling U / 2. In the tutorial it is written
    # that mu = U/2 isn't implemented, but it is (!!!). Automatically mu = U/2, for half-filling.
    # The sign problem can occure away from half-filling. Don't touch.
    Nk                  = 64       # num. of kpoints in each direction, 64 is better for qmc (Friedrich K.)
    num_of_neighbours   = 3

    #ct_hub parameters
    sweeps     = 10**10          # 10**8 .. 10**11
    hours_max  = 3               # max hours
    time_limit = hours_max*60*60 # in seconds

    #mixing_parameters
    delta_mix  = 0.00
    lambda_mix = 0.70

    # iterations
    max_it_num = 4
    start_from_it = 1
    
    if (start_from_it > max_it_num):
        print ("The number of start iteration is less than the number of the last one.")
    #############################################
    #                                           #
    #      STATIC PARAMETERS OF A MODEL         #
    #                                           #
    #############################################
    # t         - value of a hopping integral
    # Coulomb   - value of non-local (intra-site) Coulomb interaction.
    # In papers it figurates as V.

    t    = np.empty(num_of_neighbours, dtype=np.float)
    t[0] = -0.233
    t[1] = 0.0
    t[2] = 0.0

    Coulomb     = np.empty(num_of_neighbours, dtype=np.float)
    Coulomb[0]  = 1.0 #2.46
    Coulomb[1]  = 0.0
    Coulomb[2]  = 0.0

    mu = 0.0
    particle_hole_symm  = 0
    save_param_file(lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours,
                    t, Coulomb, mu, particle_hole_symm, sweeps, time_limit, delta_mix, lambda_mix)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def get_model_parameters():
    global lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours
    global t, Coulomb, mu, particle_hole_symm, sweeps
    global time_limit, delta_mix, lambda_mix
    global max_it_num, start_from_it

    print ("Lattice type is ", lattice_type)
    print ("5*beta*U/(2pi) ~ ", str(int(5*beta*U/(2.*np.pi))))
    print ("Required number of N_tau >= {}".format(5. * U * beta))
    print ("mu = {}".format(mu))
    if (particle_hole_symm == 1):
        print("Particle - hole symmetry - yes.")
    else:
        print("Particle - hole symmetry - no.")

    print("Hopping is {}".format(t))
    print("Coulomb is {}".format(Coulomb))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    return lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm, \
           sweeps, time_limit, delta_mix, lambda_mix, max_it_num, start_from_it

def save_param_file(lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm,
                    sweeps, time_limit,  delta_mix, lambda_mix):
    # -- save parameters in file for the iteration --
    f = open("model_params.dat", "w")
    f.write("lattice_type:\t" + str(lattice_type) + "\n")
    f.write("beta =\t" + str(beta)+ "\n")
    f.write("U =\t" + str(U)+ "\n")
    f.write("mu Anderson =\t" + str(hartree_shift)+ "\n")
    f.write("Nk =\t" + str(Nk)+ "\n")
    f.write("NN =\t" + str(num_of_neighbours)+ "\n")
    f.write("Hopping =\t{}\n". format(t))
    f.write("Coulomb =\t{}\n".format(Coulomb))
    f.write("mu lattice =\t{}\n".format(Coulomb))
    if (particle_hole_symm == 1):
        f.write("Particle - hole symmetry - yes.\n")
    else:
        f.write("Particle - hole symmetry - no.\n")
    f.write("sweeps (N_MEAS = 2000) =\t{}\n".format(sweeps))
    f.write("time_limit  =\t{} hours\n".format(time_limit/3600))
    f.write("Mixing\n")
    f.write("Delta mixing parameter  =\t{}\n".format(delta_mix))
    f.write("Lambda mixing parameter  =\t{}\n".format(lambda_mix))
    f.close()

def get_shift_half_filling(dos, Erange, dE):
  if (dE < 0.0):
    dE = abs(dE)
  weight = 0.0
  Eshift = 0.0
  NE = len(dos)
  for n in range(NE):
    weight += dos[n]*dE
#    print (-Erange+n*dE,weight)
    if weight >= 0.5:
      if weight==0.5: Eshift=-Erange+n*dE
      else: Eshift=-Erange+(n-0.5)*dE
      break
  print ("Energy shift to obtain half-filling: E_shift = %f"%(Eshift))
  return Eshift


def get_van_Hove_filling(dos, Erange, dE):
  #find maximum
  NE=len(dos)
  Evl = [-Erange+x*dE for x in range(0,NE+1)]
  max=dos[0]
  maxn=0
  for n in range (0,NE):
    if dos[n] >= max:
      max=dos[n]
      maxn=n
  print ("maximum of dos @ E=%f"%Evl[maxn])
  vhf=0.0
  N=len(dos)
  for n in range (0,maxn):
    vhf += 2*dos[n]*dE #van Hove doping; factor 2 is for spin
  print ("optimal hole doping: delta ~ %f"%(1.0-vhf))
  return vhf

def get_server_run():
    path_to_exec_file   = '/storage/praha1/home/medvedeva/workspace/other_codes/CT_HYB_SEGMENT/CT-HYB-SEGMENT/build/alps_cthyb'
    num_mpi_threads     = 16
    path_to_maxent      = ''
    print_sources(path_to_exec_file, num_mpi_threads, path_to_maxent)
    return path_to_exec_file, num_mpi_threads, path_to_maxent

def get_local_run():
    path_to_exec_file   = '/Users/witcher/workspace/CT_HYB_SEGMENT/CT-HYB-SEGMENT/build/alps_cthyb'
    num_mpi_threads     = 3
    path_to_maxent      = '/Users/witcher/workspace/CT_HYB_SEGMENT/Maxent/build2/maxent'
    print_sources(path_to_exec_file, num_mpi_threads, path_to_maxent)
    return path_to_exec_file, num_mpi_threads, path_to_maxent

def print_sources(path_to_exec_file, num_mpi_threads, path_to_maxent):
   print("Path to the solver \t: ", path_to_exec_file)
   print("Num of MPI threads \t: ", num_mpi_threads)
   print("Path to the MaxEnt \t: ", path_to_maxent)
   print("\n")

