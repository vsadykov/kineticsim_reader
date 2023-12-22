# This python file is used to read the fields files and extract the high-cadence properties
# (temperatures, magnetic field energies, etc)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import kineticsim_reader as kr
import pickle
import os

# location of simulation files (assumed to be all in one folder)
simfiles_folder = '../'
# location to save the results of readings
savefiles_folder = './processing_results/'
# simulation filename, number of species, proton specie indexes, he specie indexes
simfiles_spec = [['fields.d10_A0.5Hepp_beta0.5eps1e-4_256',2,[0],[1]],\
                 ['fields.d10_A0.75Hepp_beta1_256',2,[0],[1]],\
                 ['fields.d10_E11Ap3.3Aa2.0Vd0.42',2,[0],[1]],\
                 ['fields.d10_E11Ap4.3Aa1.6',2,[0],[1]],\
                 ['fields.d10_E11Ap4.3Aa1.6Vd0.32',2,[0],[1]],\
                 ['fields.d10_He++v2_256_iden0eps1e-4t600',2,[0],[1]],\
                 ['fields.d10_pv1.5_128_64_iden0eps1e-4_dx0.75_long',2,[0],[1]],\
                 ['fields.d10_pv2a_128x3_iden0eps1e-4_dx0.75',3,[0,1],[2]]]


# running throught the models
for simfile_spec in simfiles_spec:
    # reading specifications for the filename
    simfile = simfiles_folder + simfile_spec[0]
    simfile_sh = simfile_spec[0]
    kspi = simfile_spec[1]
    kspi_pr = simfile_spec[2]
    kspi_he = simfile_spec[3]
    # header information and understanding simulations
    print(" ")
    print("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
    print("-> SIMULATION: "+simfile)
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    dx, dy, tim, bx, by, bz, ex, ey, ez, rvxh, rvyh, rvzh, rdnh, tpal, tper, me_perp, me_tot = kr.read_fieldsfile(simfile, kspi=kspi)    
    # saving properties
    np.save(savefiles_folder + simfile_sh + '.timing.npy', tim)
    np.save(savefiles_folder + simfile_sh + '.scales.npy', [dx,dy])
    np.save(savefiles_folder + simfile_sh + '.tpal.npy', tpal)
    np.save(savefiles_folder + simfile_sh + '.tper.npy', tper)
    np.save(savefiles_folder + simfile_sh + '.me_perp.npy', me_perp)
    np.save(savefiles_folder + simfile_sh + '.me_tot.npy', me_tot)
