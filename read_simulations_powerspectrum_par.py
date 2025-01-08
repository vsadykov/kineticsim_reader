# This python file is used to read the fields files and extract the high-cadence properties
# (temperatures, magnetic field energies, etc)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import kineticsim_reader as kr
import pickle
import os
import multiprocessing as mp


def process_simulation_ps(simfile_spec):
    # reading specifications for the filename
    simfile = '../' + simfile_spec[0]
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
    # there are several aspects to do after we have read the entire fields file:
    # 1. Resampling the electric field array to the timing of 2.5 gyroperiods
    # 2. Extending the solution for 30 gyroperiods back and forth from the time interval
    # 3. Computing the power spectrum at each timing moment averaged over 50 gyroperiods and centered at the moment of measurements.
    timing_nonint = tim[:,1]
    timing_int = np.arange(timing_nonint[0]-30.0, timing_nonint[-1]+30.0, 2.5)
    ez_int = np.zeros([timing_int.shape[0],ez.shape[1],ez.shape[2]], dtype=float)
    for i in range (0, ez.shape[1], 1):
        for j in range (0, ez.shape[2], 1):
            ez_nonint = ez[:,i,j]
            ez_int[:,i,j] = np.interp(timing_int, timing_nonint, ez[:,i,j])
    # power spectrum after interpolation and extension
    ps_total = np.zeros([21,timing_int.shape[0]-24], dtype=float)
    for t in range (2, timing_int.shape[0]-22, 1):
        print("Current step " + str(t-2) + " out of " + str(timing_int.shape[0]-24) + " for " + simfile_sh)
        for i in range (0, ez.shape[1], 1):
            for j in range (0, ez.shape[2], 1):
                ez_ts = ez_int[t:t+21,i,j]
                ps = np.abs(np.fft.fft(ez_ts))**2
                freqs = np.fft.fftfreq(ps.size, 2.5)
                idx = np.argsort(freqs)
                ps_total[:,t-2] += ps[idx]
    ps_total /= ez.shape[1]*ez.shape[2]    
    # saving simulations
    np.save('./processing_results/' + simfile_sh + '.ps_ps.npy', ps_total)
    np.save('./processing_results/' + simfile_sh + '.ps_freqs.npy', freqs[idx])
    np.save('./processing_results/' + simfile_sh + '.ps_inttiming.npy', timing_int[12:-12])


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
['fields.d10_e261040ap1.40.4_7',2,[0,1],[]]\
]

with mp.Pool(28) as p:
    p.map(process_simulation_ps, simfiles_spec)
