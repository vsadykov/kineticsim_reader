# This python file is used to read the simulations and to preprocess them into ML-compatible format

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import kineticsim_reader as kr
import pickle
import os

simfile = '../particles.d11_E11Ap4.3Aa1.6Vd0.32'
foldername = './outputs/outputs.d11_E11Ap4.3Aa1.6Vd0.32/'

# print("-------------- HEADER AND FILE INFORMATION -----------------")
# header = kr.show_fileinfo(simfile)

print("-------------- SAVING FILES INTO THE FOLDER -----------------")
kr.save_simulationframes(simfile, 'distributions.', foldername, framestart = 0, frameend = 47)

print("-------------- READING PROTON VDFS -----------------")
fnames = []
for fname in os.listdir(foldername):
    if (fname[0:5] != 'distr'): continue
    fnames.append(fname)

fnames.sort()
anisotropies_all = []
moments_all = []
times_all = []
for fname in fnames:
    timfrm, timep, xp, yp, vxp, vyp, vzp = kr.retrieve_simulationframe(foldername + fname)
    # generating histograms for proton distributions
    hist, vx_edges, vy_edges = kr.generate_histogram(vxp, vyp, 0.02, 0, 2)
    kr.visualize_histogram(hist, vx_edges, vy_edges, 'Proton VDF', to_image = True, imfile = foldername + 'Proton_VDF.' + fname[0:-4] + '.png')
    data = {'timeframe':timfrm, 'timep':timep, 'hist':hist, 'vx_edges':vx_edges, 'vy_edges':vy_edges}
    histfile = foldername + 'Proton_VDF.' + fname[0:-4] + '.pkl'
    with open(histfile, "wb") as outfile:
        pickle.dump(data, outfile)
    anisotropy, moments = kr.calculate_anisotropy_moments(vxp, vyp, vzp, 0, 2)
    print("Proton Anisotropy for t = " + str(timep[0]) + " is " + str(anisotropy))
    anisotropies_all.append(anisotropy)
    moments_all.append(moments)
    times_all.append(timep)
outputs_file = foldername + 'Proton_parameters.pkl'
data = {'moments_all':moments_all, 'anisotropies_all':anisotropies_all, 'times_all':times_all}
with open(outputs_file, "wb") as outfile:
    pickle.dump(data, outfile)

print("-------------- READING HELIUM VDFS -----------------")
anisotropies_all = []
moments_all = []
times_all = []
for fname in fnames:
    timfrm, timep, xp, yp, vxp, vyp, vzp = kr.retrieve_simulationframe(foldername + fname)
    # generating histograms for proton distributions
    hist, vx_edges, vy_edges = kr.generate_histogram(vxp, vyp, 0.02, 2, 4)
    kr.visualize_histogram(hist, vx_edges, vy_edges, 'He VDF', to_image = True, imfile = foldername + 'He_VDF.' + fname[0:-4] + '.png')
    data = {'timeframe':timfrm, 'timep':timep, 'hist':hist, 'vx_edges':vx_edges, 'vy_edges':vy_edges}
    histfile = foldername + 'He_VDF.' + fname[0:-4] + '.pkl'
    with open(histfile, "wb") as outfile:
        pickle.dump(data, outfile)
    anisotropy, moments = kr.calculate_anisotropy_moments(vxp, vyp, vzp, 2, 4)
    print("Proton Anisotropy for t = " + str(timep[0]) + " is " + str(anisotropy))
    anisotropies_all.append(anisotropy)
    moments_all.append(moments)
    times_all.append(timep)
outputs_file = foldername + 'He_parameters.pkl'
data = {'moments_all':moments_all, 'anisotropies_all':anisotropies_all, 'times_all':times_all}
with open(outputs_file, "wb") as outfile:
    pickle.dump(data, outfile)

#os.system('rm -rf ' + foldername + 'distributions*')
