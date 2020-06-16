# # # # # # # # # # # # # # # # # # # # # # # #
#         Medvedeva D.S. 2018 - 2020          #
# For questions: medvedeva.ds@gmail.com       #
# Skype: daryacooper                          #
# # # # # # # # # # # # # # # # # # # # # # # #

#################################################
#                                               #
#        SPECIFICATION OF A MODEL               #
#                                               #
# (1) Copy the file set_parameters_change.py    #
# (2) Change the name as set_parameters.py      #
# (3) Set parameters                            #
#                                               #
#################################################

import numpy as np
#############################################
#                                           #
#        SPECIFICATION OF A MODEL           #
#                                           #
#############################################
def set_model_parameters():
    lattice_type        = 'triangular' # write square || triangular
    beta                = 50.       # inversive temperature as \beta = 1./T. Kind of estimated temperatures for qmc are in
                                    # range [0+ .. 100] stable, beta > 100 - unstable behavior of the solver can occur.
                                    
    U                   = 1.0       # local (inter-site) Hubbard repulsion
    
    hartree_shift       = -0.0      # Hartree shift (dctype = double counting = ... = \mu in ct-hyb)
                                    # for a half filling U / 2. In the tutorial it is written
                                    # that mu = U/2 isn't implemented, but it is (!!!). Automatically mu = U/2 for
                                    # half-filling. The sign problem can occure away from half-filling. Don't touch.
                                    
    Nk                  = 64        # number of kpoints in each direction, 64 is better for qmc (Friedrich K.)
    num_of_neighbours   = 3         # It is better to save it like this and use 0.0 values for the neighbours of the
                                    # larger order.

    #############################################
    #                                           #
    #      STATIC PARAMETERS OF A MODEL         #
    #                                           #
    #############################################
    # t         - value of the hopping integrals.
    
    # Coulomb   - value of non-local (intra-site) Coulomb interaction.
    # In papers it figurates as V.

    if lattice_type == 'square':
        # for square lattice everything is easier
        t       = -0.25
        Coulomb = 0.5
        mu = U / 2.
        particle_hole_symm  = 1     # Because the lattice is symmetrical
        
    elif lattice_type == 'triangular':
        t    = np.empty(num_of_neighbours, dtype=np.float)

        t[0] =  0.0614
        t[1] =  0.0975
        t[2] =  0.0151

        Coulomb     = np.empty(num_of_neighbours, dtype=np.float)
        Coulomb[0]  = 0.0
        Coulomb[1]  = 0.0
        Coulomb[2]  = 0.0
            
        mu = 0.0                    # this is the chemical potential (NOT IN THE IMPURITY MODEL)
        
        particle_hole_symm  = 0     # Because the lattice is frustrated
        
    print ("Lattice type is ", lattice_type)
    print ("5*beta*U/(2pi) ~ ", str(int(5*beta*U/(2.*np.pi))))
    print ("mu = {}".format(mu))
    if (particle_hole_symm == 1):
        print("Particle - hole symmetry - yes.")
    else:
        print("Particle - hole symmetry - no.")
    
    print("Hopping is {}".format(t))
    print("Coulomb is {}".format(Coulomb))
    return lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, Coulomb, mu, particle_hole_symm 
    
def get_shift_half_filling(dos, Erange, dE):
  if (dE < 0.0):
    dE = abs(dE)
  weight = 0.0
  Eshift = 0.0
  NE = len(dos)
  for n in range(NE):
    weight += dos[n]*dE
    if weight >= 0.5:
      if weight==0.5: Eshift=-Erange+n*dE
      else: Eshift=-Erange+(n-0.5)*dE
      break
  print ("Energy shift to obtain half-filling: E_shift = {}".format(Eshift))
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
