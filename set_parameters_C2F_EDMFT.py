import numpy as np
#############################################
#                                           #
#        SPECIFICATION OF A MODEL           #
#                                           #
#############################################
def set_model_parameters():
    lattice_type        = 'triangular' # write square || triangular
    beta                = 100.     # inversive temperature as \beta = 1./T
    U                   = 5.16      # local (inter-site) Hubbard repulsion
    #mu                  = U/2.   # for a half filling U / 2. In case of square lattice it should be mu = U/2. !!!!
    #mu                  = 0.8 * t
    hartree_shift       = -0.0      # Hartree shift (\mu in ct-hyb). for a half filling U / 2. In the tutorial it is written
    # that mu = U/2 isn't implemented, but it is (!!!). Automatically mu = U/2, for half-filling.
    # The sign problem can occure away from half-filling. Don't touch.
    Nk                  = 64       # num. of kpoints in each direction, 64 is better for qmc (Friedrich K.)
    num_of_neighbours   = 3

    #ct_hub parameters
    sweeps = 10**10
       
    #############################################
    #                                           #
    #      STATIC PARAMETERS OF A MODEL         #
    #                                           #
    #############################################
    # t         - value of a hopping integral
    # Coulomb   - value of non-local (intra-site) Coulomb interaction.
    # In papers it figurates as V.

    if lattice_type == 'square':
        t       = -0.25
        Coulomb = 0.5
        mu = U / 2.
        particle_hole_symm  = 1
        
    elif lattice_type == 'triangular':
        t    = np.empty(num_of_neighbours, dtype=np.float)
        t[0] = -0.233
        t[1] =  0.006
        t[2] = -0.021

        Coulomb     = np.empty(num_of_neighbours, dtype=np.float)
        Coulomb[0]  = 2.46
        Coulomb[1]  = 1.66
        Coulomb[2]  = 1.46
        
#        mu =  0.8 * t[0]
        mu = 0.0
        
        particle_hole_symm  = 0
        
    save_param_file(lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm, sweeps)
    
    return lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm, sweeps

def save_param_file(lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm, sweeps):
    # -- save parameters in file for the iteration --
    print ("Lattice type is ", lattice_type)
    print ("5*beta*U/(2pi) ~ ", str(int(5*beta*U/(2.*np.pi))))
    print ("mu = {}".format(mu))
    if (particle_hole_symm == 1):
        print("Particle - hole symmetry - yes.")
    else:
        print("Particle - hole symmetry - no.")

    print("Hopping is {}".format(t))
    print("Coulomb is {}".format(Coulomb))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    f = open("save_model_params.dat", "w")
    f.write("lattice_type:\t" + str(lattice_type))
    f.write("beta =\t" + str(beta))
    f.write("U =\t" + str(U))
    f.write("mu Anderson =\t" + str(hartree_shift))
    f.write("Nk =\t" + str(Nk))
    f.write("NN =\t" + str(num_of_neighbours))
    f.write("Hopping =\t{}". format(t))
    f.write("Coulomb =\t{}".format(Coulomb))
    f.write("mu lattice =\t{}".format(Coulomb))
    if (particle_hole_symm == 1):
        f.write("Particle - hole symmetry - yes.")
    else:
        f.write("Particle - hole symmetry - no.")
    f.write("sweeps (N_MEAS = 2000) =\t{}".format(sweeps))
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
