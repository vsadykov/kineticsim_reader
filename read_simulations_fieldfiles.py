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
simfiles_spec = [\
['fields.d10_A0.5Hepp_beta0.5eps1e-4_256',2,[0],[1]],\
['fields.d10_A0.75Hepp_beta1_256',2,[0],[1]],\
['fields.d10_E11Ap3.3Aa2.0Vd0.42',2,[0],[1]],\
['fields.d10_E11Ap4.3Aa1.6',2,[0],[1]],\
['fields.d10_E11Ap4.3Aa1.6Vd0.32',2,[0],[1]],\
['fields.d10_E12Ap1.86Aa1.0Vd0.32_256_256x256',2,[0],[1]],\
['fields.d10_E12Ap1.86Aa1.0Vd0.32_512_256x256',2,[0],[1]],\
['fields.d10_He++A10_256_iden0eps0',2,[0],[1]],\
['fields.d10_He++v2_256_iden0eps1e-4t600',2,[0],[1]],\
['fields.d10_He++vd1.5_256_iden0eps1e-4',2,[0],[1]],\
['fields.d10_pv1.5_128_64_iden0eps1e-4_dx0.75_long',2,[0,1],[]],\
['fields.d10_pv1Ap2Apb2betac0.214betab0.858_128_128x2_dx0.75_t3000',2,[0,1],[]],\
['fields.d10_pv2a_128x3_iden0eps1e-4_dx0.75',3,[0,1],[2]],\
['fields.d10_pv2Ap1Ab1betac0.429betab0.858_128_128x2_dx0.75_t3000',2,[0,1],[]],\
['fields.d10_pv2Ap1Ab2betac0.429betab0.858_128_128x2_dx0.75_t3000',2,[0,1],[]],\
['fields.d10_pv2Ap2Apb2betac0.214betab0.858_128_128x2_dx0.75_t3000',2,[0,1],[]],\
['fields.d10_pv2av2.3_128x3_iden0eps1e-4_dx0.75',4,[0,1],[2,3]],\
['fields.d10_pv2av2Ap1Aa1beta0.429_128_128x2_dx0.75_t3000',4,[0,1],[2,3]],\
['fields.d10_pv2av2_rdna0.03375_128x3_iden0eps1e-4_dx0.75_t6000',4,[0,1],[2,3]],\
['fields.d10_vap1.2Ap1Aa0.75_rdna_0.05',2,[0],[1]],\
['fields.d10_vap1.2Ap3.35Aa2.05rdna_0.007',2,[0],[1]],\
['fields.d10_vap1.5Ap1.5Aa1rdna_0.007',2,[0],[1]],\
['fields.d10_e260945ap1.30.5_1',2,[0,1],[]],\
['fields.d10_e260955ap2.20.4_2',2,[0,1],[]],\
['fields.d10_e261013ap1.50.6_3',2,[0,1],[]],\
['fields.d10_e261016ap1.70.6_4',2,[0,1],[]],\
['fields.d10_e261019ap1.50.4_5',2,[0,1],[]],\
['fields.d10_e261022ap1.40.4_6',2,[0,1],[]],\
['fields.d10_e261040ap1.40.4_7',2,[0,1],[]],\
['fields.d10_pv1.4av2Ap2Apb2betac0.214betab0.858_128_128x4_dx0.75',4,[0,1],[2,3]],\
['fields.d10_pv2Ap2Apb2beta_pb0.429_128_128x2_dx0.75',2,[0,1],[]],\
['fields.d10_pv2av1.4Ap1Aa1betac0.214betab0.858_128_128x4_dx0.75',4,[0,1],[2,3]],\
['fields.d10_pv2av1.4Ap2Aa2betac0.214betab0.858_128_128x4_dx0.75',4,[0,1],[2,3]],\
['fields.d10_pv2av1.4Ap2Ab2beta0.429_128_128x4_dx0.75_t3000',4,[0,1],[2,3]],\
['fields.d10_pv2av2Ap1Aa1beta0.429_128_128x4_dx0.75_t3000',4,[0,1],[2,3]],\
['fields.d10_pv2av2Ap2Aa2beta0.429_128_128x4_dx0.75_t3000',4,[0,1],[2,3]]\
]

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
