import numpy as np
import parameters
import iteration_cycle

global sqrt32, sqrt3
sqrt32 = np.sqrt(3.)/2.
sqrt3 = np.sqrt(3.)

def interaction_GMKG(filename, Nk, lattice_type, param):
    kpoints = np.array([[[x, y] for x in range(Nk)] for y in range(Nk)]).reshape((Nk * Nk, 2))
    kpoints_coordinates = np.zeros((Nk*Nk, 2), dtype = np.float)
    interaction = np.zeros(Nk * Nk, dtype = np.float)
    
    b1, b2, b3 = iteration_cycle.get_b_vectors(lattice_type)
    
#    print('b1 = ', b1 )
#    print('b2 = ', b2 )
#    print('b3 = ', b3 )
    
    kstep_x, kstep_y = iteration_cycle.get_kstep(lattice_type, Nk)
    
    for k in range(Nk * Nk):
        vec_ = [0.0, 0.0]
        for i in range(2):
            if (lattice_type == "square"):
                kpoints_coordinates[k][i] = kpoints[k][0] * b1[i]/ (Nk - 1) + kpoints[k][1] * b2[i] / (Nk - 1)
            elif (lattice_type == "triangular"):
                kpoints_coordinates[k][i] = kpoints[k][0] * b1[i]/ (Nk - 1) + kpoints[k][1] * b2[i] / (Nk - 1)
    
    GM = open(filename + "GM.dat", "w")
    MK = open(filename + "MK.dat", "w")
    KG = open(filename + "KG.dat", "w")
    x_boundary = 0.0
    y_boundary = 0.0
    i = 0
    for k in kpoints_coordinates:
        # TRIANG
        if lattice_type == 'triangular':
            k_x = k[0]
            k_y = k[1]
            
            x_coord = x_boundary + k_x
            y_coord = y_boundary + k_y
            
#             G -> M
            if (x_coord == 0.0) and (y_coord >= 0.0) and (y_coord <= 2.*np.pi/sqrt3):
                interaction[i]  = 2. * param[0] * ( np.cos(x_coord) + np.cos(0.5 * x_coord + sqrt32 * y_coord) + np.cos(0.5 * x_coord - sqrt32 * y_coord) )
                interaction[i] += 2. * param[1] * ( np.cos(sqrt3 * y_coord) + np.cos(3./2. * x_coord + sqrt32 * y_coord) + np.cos(3./2. * x_coord - sqrt32 * y_coord) )
                interaction[i] += 2. * param[2] * ( np.cos(2. * x_coord) + np.cos(x_coord + sqrt3 * y_coord) + np.cos(x_coord - sqrt3 * y_coord) )
            
                GM.write(str(np.round(k_x, 6)))
                GM.write('\t')
                GM.write(str(np.round(k_y,6)))
                GM.write('\t')
                GM.write(str(np.round(interaction[i],6)))
                GM.write('\n')
                i += 1
#             M -> K and (x_coord >= 0.0) and (x_coord <= 2.*np.pi/3.)
            elif (abs(y_coord - 2.* np.pi/sqrt3) < 0.01) and (x_coord >= 0.0) and (x_coord <= 2.*np.pi/3.):
                interaction[i]  = 2. * param[0] * ( np.cos(x_coord) + np.cos(0.5 * x_coord + sqrt32 * y_coord) + np.cos(0.5 * x_coord - sqrt32 * y_coord) )
                interaction[i] += 2. * param[1] * ( np.cos(sqrt3 * y_coord) + np.cos(3./2. * x_coord + sqrt32 * y_coord) + np.cos(3./2. * x_coord - sqrt32 * y_coord) )
                interaction[i] += 2. * param[2] * ( np.cos(2. * x_coord) + np.cos(x_coord + sqrt3 * y_coord) + np.cos(x_coord - sqrt3 * y_coord) )
            
                MK.write(str(np.round(k_x, 6)))
                MK.write('\t')
                MK.write(str(np.round(k_y,6)))
                MK.write('\t')
                MK.write(str(np.round(interaction[i],6)))
                MK.write('\n')
                i += 1
#           K -> G  (y_coord >= 0.0) and (y_coord <= 2.*np.pi/sqrt3) and (x_coord >= 0.0) and (x_coord <= 2.*np.pi / 3.) and
            elif (y_coord >= 0.0) and (y_coord <= 2.*np.pi/sqrt3) and (x_coord >= 0.0) and (x_coord <= 2.*np.pi / 3.) and ( x_coord * sqrt3 - y_coord < 0.000001):
                interaction[i]  = 2. * param[0] * ( np.cos(x_coord) + np.cos(0.5 * x_coord + sqrt32 * y_coord) + np.cos(0.5 * x_coord - sqrt32 * y_coord) )
                interaction[i] += 2. * param[1] * ( np.cos(sqrt3 * y_coord) + np.cos(3./2. * x_coord + sqrt32 * y_coord) + np.cos(3./2. * x_coord - sqrt32 * y_coord) )
                interaction[i] += 2. * param[2] * ( np.cos(2. * x_coord) + np.cos(x_coord + sqrt3 * y_coord) + np.cos(x_coord - sqrt3 * y_coord) )
           
                KG.write(str(np.round(k_x, 6)))
                KG.write('\t')
                KG.write(str(np.round(k_y,6)))
                KG.write('\t')
                KG.write(str(np.round(interaction[i],6)))
                KG.write('\n')
                i += 1
    GM.close()
    MK.close()
    KG.close()
    
    return interaction, kpoints_coordinates

"""
    main part
"""

filename = 'Xw.dat'
freq, X_w = iteration_cycle.read_real_function(filename)
number_of_real_points = len(X_w)

parameters.set_model_parameters()
lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, V, mu, particle_hole_symm, sweeps, time_limit,  delta_mix, lambda_mix, number_of_iterations, start_from_it = \
parameters.get_model_parameters()
# 1. Selection G - M - K - G
test, kpoints_coordinates = interaction_GMKG('test', Nk, lattice_type, V)


    
