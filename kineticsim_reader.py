"""
__author__ = "Viacheslav M Sadykov"
__version__ = 0.0.1
__status__ = "Testing"
"""

# The current library presents a set of functions to read and pre-process the
# results of hybrid kinetic simulations. Check the description of indivisual
# routines for more details.


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import pickle
from scipy.stats import moment


def show_fileinfo(filename, kspi=4, nsim_out = False):
    """
    The routine presents the header information and informs about the number
    of snapshots available for the analysis. The routine also allows one to
    read the header of the file.
    Input:
    filename - the relative path to the file containing the result of hybrid
    kinetic simulations
    kspi - number of species; default number is 4
    Output:
    header - the header array of the simulations
    """
    # reading the file header
    with open(filename, 'rb') as file:
        # setting up some parameters
        ihparm = 19     # number of parameters stored in the header
        ihdr = 19 + 8*kspi
        # reading the header information
        header = np.empty([ihdr], dtype=np.float32)
        header = np.fromfile(file, dtype=np.float32, count=len(header))
        # presenting header information
        print("Header information:")
        print("HEADER[0]:          dt:", header[0])
        print("HEADER[1]:          dx:", header[1])
        print("HEADER[2]:          dy:", header[2])
        print("HEADER[3]:       itmax:", header[3])
        print("HEADER[4]:        ifld:", header[4])
        print("HEADER[5]:      betaen:", header[5])
        print("HEADER[6]:         nis:", header[6])
        print("HEADER[7]:          nx:", header[7])
        print("HEADER[8]:          ny:", header[8])
        print("HEADER[9]:         npt:", header[9])
        print("HEADER[10]:    ifields:", header[10])
        print("HEADER[11]: iparticles:", header[11])
        print("HEADER[12]:    ienergy:", header[12])
        print("HEADER[13]:    ifilter:", header[13])
        print("HEADER[14]:        bx0:", header[14])
        print("HEADER[15]:        by0:", header[15])
        print("HEADER[16]:        bz0:", header[16])
        print("HEADER[17]:      icont:", header[17])
        print("HEADER[18]:      xlres:", header[18])
        # understanding how many frames are in the file
        nis = int(header[6])
        npts = np.empty([nis], dtype=np.int32)
        for _is in range (0, nis, 1):
            ihs = ihparm + 8*_is
            npts[_is] = np.int32(header[ihs + 7])
        nptsm=np.int32(np.amax(npts))
        xp = np.empty([nptsm,nis], dtype=np.float32)
        yp = np.empty([nptsm,nis], dtype=np.float32)
        vxp = np.empty([nptsm,nis], dtype=np.float32)
        vyp = np.empty([nptsm,nis], dtype=np.float32)
        vzp = np.empty([nptsm,nis], dtype=np.float32)
        framecur = 1
        while (True):
            timfrm = np.fromfile(file, dtype=np.float32, count=1)
            timep = np.fromfile(file, dtype=np.float32, count=1)
            for i in range (0, nis, 1):
                var = np.empty([5*npts[i]], dtype=np.float32)
                var = np.fromfile(file, dtype=np.float32, count=len(var))
                xp[0:npts[i],i] = np.copy(var[0::5])
                yp[0:npts[i],i] = np.copy(var[1::5])
                vxp[0:npts[i],i] = np.copy(var[2::5])
                vyp[0:npts[i],i] = np.copy(var[3::5])
                vzp[0:npts[i],i] = np.copy(var[4::5])
            print('Frame record in file:', framecur, ', Frame number:', timfrm[0], ', Frame timing:', timep[0])
            framecur += 1
            if (timfrm[0] + header[11] >= header[3]): break
    print('Total number of frames:', framecur-1)
    # returning the header information
    if (nsim_out == True):
        return framecur-1, header
    else:
        return header


def return_selectedframe(filename, framenumber, kspi=4):
    """
    The routine retuns the velocity and position distributions for the selected
    frame record within the file.
    The routine is updated and utilizes 'seek' command to skip to the appropriate
    position in the file.
    Input:
    filename - the relative path to the file containing the result of hybrid
    kinetic simulations
    framenumber - the frame record number within the file
    kspi - number of species; default number is 4
    Output:
    timfrm - frame number
    timep - timing of the frame
    xp - X-positions of particles
    yp - Y-positions of particles
    vxp - X-velocities of particles
    vyp - Y-velocities of particles
    vzp - Z-velocities of particles
    """
    # reading the file header
    with open(filename, 'rb') as file:
        # setting up some parameters
        ihparm = 19     # number of parameters stored in the header
        ihdr = 19 + 8*kspi
        # reading the header information
        header = np.empty([ihdr], dtype=np.float32)
        header = np.fromfile(file, dtype=np.float32, count=len(header))
        # understanding how many frames are in the file
        nis = int(header[6])
        npts = np.empty([nis], dtype=np.int32)
        for _is in range (0, nis, 1):
            ihs = ihparm + 8*_is
            npts[_is] = np.int32(header[ihs + 7])
        nptsm=np.int32(np.amax(npts))
        xp = np.empty([nptsm,nis], dtype=np.float32)
        yp = np.empty([nptsm,nis], dtype=np.float32)
        vxp = np.empty([nptsm,nis], dtype=np.float32)
        vyp = np.empty([nptsm,nis], dtype=np.float32)
        vzp = np.empty([nptsm,nis], dtype=np.float32)
        # skipping the appropriate number of bytes before the frame of interest
        skipbytes = 2*4
        for i in range (0, nis, 1): skipbytes += 5*npts[i]*4
        file.seek(skipbytes*(framenumber-1), 1)
        # reading the timeframe of interest
        timfrm = np.fromfile(file, dtype=np.float32, count=1)
        timep = np.fromfile(file, dtype=np.float32, count=1)
        for i in range (0, nis, 1):
            var = np.empty([5*npts[i]], dtype=np.float32)
            var = np.fromfile(file, dtype=np.float32, count=len(var))
            xp[0:npts[i],i] = np.copy(var[0::5])
            yp[0:npts[i],i] = np.copy(var[1::5])
            vxp[0:npts[i],i] = np.copy(var[2::5])
            vyp[0:npts[i],i] = np.copy(var[3::5])
            vzp[0:npts[i],i] = np.copy(var[4::5])
    return timfrm, timep, xp, yp, vxp, vyp, vzp


def save_simulationframes(filename, outfile_stamp, foldername, framestart = 0, frameend = -1, kspi=4):
    """
    The routine saves the particle coordinate distributions and velocity distributions
    into several files per frame (tagged with proton gyroperiod time).
    Input:
    filename - the relative path to the file containing the result of hybrid
    kinetic simulations
    outfile_stamp - the stamp of the output file which precedes the gyroperiod time
    in the file header
    foldername - the path of the folder to save the files into
    framestart (optional) - the frame to start from; the default is zero
    frameend (optional) - the frame to end at; the default is -1. The frame equal to -1
    means that the file will be read until the end.
    kspi - number of species; default number is 4
    Output:
    No output. However, the Python pickle files are going to be saved.
    """
    # reading the file header
    with open(filename, 'rb') as file:
        print('Reading the file...')
        # setting up some parameters
        ihparm = 19     # number of parameters stored in the header
        ihdr = 19 + 8*kspi
        # reading the header information
        header = np.empty([ihdr], dtype=np.float32)
        header = np.fromfile(file, dtype=np.float32, count=len(header))
        # understanding how many frames are in the file
        nis = int(header[6])
        npts = np.empty([nis], dtype=np.int32)
        for _is in range (0, nis, 1):
            ihs = ihparm + 8*_is
            npts[_is] = np.int32(header[ihs + 7])
        nptsm=np.int32(np.amax(npts))
        xp = np.empty([nptsm,nis], dtype=np.float32)
        yp = np.empty([nptsm,nis], dtype=np.float32)
        vxp = np.empty([nptsm,nis], dtype=np.float32)
        vyp = np.empty([nptsm,nis], dtype=np.float32)
        vzp = np.empty([nptsm,nis], dtype=np.float32)
        framecur = 1
        while (True):
            timfrm = np.fromfile(file, dtype=np.float32, count=1)
            timep = np.fromfile(file, dtype=np.float32, count=1)
            for i in range (0, nis, 1):
                var = np.empty([5*npts[i]], dtype=np.float32)
                var = np.fromfile(file, dtype=np.float32, count=len(var))
                xp[0:npts[i],i] = np.copy(var[0::5])
                yp[0:npts[i],i] = np.copy(var[1::5])
                vxp[0:npts[i],i] = np.copy(var[2::5])
                vyp[0:npts[i],i] = np.copy(var[3::5])
                vzp[0:npts[i],i] = np.copy(var[4::5])
            if ( (framecur >= framestart) and ( (framecur <= frameend) or (frameend == -1) ) ):
                pklfile = foldername + outfile_stamp + str(timep[0]) + '.pkl'
                data = {'timeframe':timfrm, 'timep':timep, 'xp':xp, 'yp':yp, \
                        'vxp':vxp, 'vyp':vyp, 'vzp':vzp}
                with open(pklfile, "wb") as outfile:
                    pickle.dump(data, outfile)
                print("FILE WRITTEN: " + str(pklfile))
            framecur += 1
            if ( (timfrm[0] + header[11] >= header[3]) or ( (framecur > frameend) and (frameend != -1) ) ):
                break
    return


def retrieve_simulationframe(filename):
    """
    The routine retrieves the information in pickle file (see previous) back into arrays.
    Input:
    filename - the pickle file to retrieve the arrays from
    Output:
    timfrm - frame number
    timep - timing of the frame
    xp - X-positions of particles
    yp - Y-positions of particles
    vxp - X-velocities of particles
    vyp - Y-velocities of particles
    vzp - Z-velocities of particles
    """
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    timfrm = data['timeframe']
    timep = data['timep']
    xp = data['xp']
    yp = data['yp']
    vxp = data['vxp']
    vyp = data['vyp']
    vzp = data['vzp']
    return timfrm, timep, xp, yp, vxp, vyp, vzp


def generate_histogram(vxp, vyp, resl=0.05, xlim=[-3.5,3.5], ylim=[-3.5,3.5], kspi_indexes = [0]):
    """
    The routine constructs a histogram based on the velocity distributions and
    resolution provided.
    Input:
    vxp - X-velocities of particles (along the external magnetic field)
    vyp - Y-velocities of particles
    resl - resolution (in Alfven speed units) of the histogram
    kspi_indexes - indexes of species to generate the histogram for. Default is 0 only.
    Output:
    hist - two-dimensional histogram of particle distribution
    vx_edges, vy_edges - bins used for the histogram
    """
    vx_edges = np.arange(xlim[0], xlim[1] + resl, resl)
    vy_edges = np.arange(ylim[0], ylim[1] + resl, resl)
    hist, vx_edges, vy_edges = np.histogram2d(vxp[:,kspi_indexes].flatten(), vyp[:,kspi_indexes].flatten(), \
                                              bins=(vx_edges, vy_edges))
    return hist, vx_edges, vy_edges


def visualize_histogram(hist, vx_edges, vy_edges, title, to_image = False, imfile = 'test.png', logflag=False):
    """
    The routine visualizes the histogram.
    Input:
    hist - two-dimensional histogram of particle distribution
    vx_edges, vy_edges - bins used for the histogram
    title - the title to use for the histogram
    Output:
    None
    """
    matplotlib.rcParams.update({'font.size':10})
    im, ax = plt.subplots(1, 1, figsize = (7, 5), dpi=150)
    if (logflag == False):
        fig = ax.imshow(hist.T, interpolation='nearest', origin='lower', \
                  extent=[vx_edges[0], vx_edges[-1], vy_edges[0], vy_edges[-1]], cmap='rainbow')
        ax.set(xlabel = r'Vx / V$_{A}$', ylabel = r'Vy / V$_{A}$', title = title)
        cbar = im.colorbar(fig, ax=ax, label='# of particles')
    else:
        fig = ax.imshow(np.log10(hist.T+1), interpolation='nearest', origin='lower', \
                  extent=[vx_edges[0], vx_edges[-1], vy_edges[0], vy_edges[-1]], cmap='rainbow')
        ax.set(xlabel = r'Vx / V$_{A}$', ylabel = r'Vy / V$_{A}$', title = title)
        cbar = im.colorbar(fig, ax=ax, label='log10(# of particles + 1)')
    plt.tight_layout()
    if (to_image):
        plt.savefig(imfile)
        plt.close()
    else:
        plt.show()
        plt.close()
    return


def calculate_anisotropy_moments(vxp, vyp, vzp, kspi_indexes = [0]):
    """
    The routine calculates anisotropy of the VDF.
    Input:
    vxp - X-velocities of particles (along the external magnetic field)
    vyp - Y-velocities of particles
    vzp - Z-velocities of particles
    kspi_indexes - indexes of species to generate the histogram for. Default is 0 only.
    Output:
    anisotropy - anisotropy of the VDF
    moments[mn,i] - an array of moments where mn is a moment number and i is
    a dimension (x, y, z)
    """
    # computing moments
    moments = np.zeros([4,3], dtype=float)
    # 1st moment will be computed with respect to 0th point as reference instead of the mean
    mn = 0
    moments[mn,0] = np.mean(vxp[:,kspi_indexes].flatten())
    moments[mn,1] = np.mean(vyp[:,kspi_indexes].flatten())
    moments[mn,2] = np.mean(vzp[:,kspi_indexes].flatten())
    for mn in range (1, 4, 1):
        moments[mn,0] = moment(vxp[:,kspi_indexes].flatten(), moment=mn+1)
        moments[mn,1] = moment(vyp[:,kspi_indexes].flatten(), moment=mn+1)
        moments[mn,2] = moment(vzp[:,kspi_indexes].flatten(), moment=mn+1)
    anisotropy = (moments[1,1] + moments[1,2])/moments[1,0]/2.0
    return anisotropy, moments


def calculate_anisotropies_moments_selectedframes(filename, framestart = 0, frameend = -1, kspi=4, \
                                                  kspi_indexes_protons = [0,1], kspi_indexes_he = [2,3]):
    """
    The routine calculates anisotroies and moments for both proton and He populations
    for the selected range of frames
    Input:
    filename - the relative path to the file containing the result of hybrid
    framestart (optional) - the frame to start from; the default is zero
    frameend (optional) - the frame to end at; the default is -1. The frame equal to -1
    means that the file will be read until the end.
    kspi - number of species; default number is 4
    kspi_indexes_protons - indexes of proton species. Default is [0,1]
    kspi_indexes_he - indexes of He species. Default is [2,3]
    Output:
    anisotropies_p[t] - anisotropy of the proton VDFs. t is a frame number
    moments_p[t,mn,i] - an array of proton moments where mn is a moment number, i is
    a spatial dimension (x, y, or z), and t is a frame number
    anisotropies_he[t] - anisotropy of the proton VDFs. t is a frame number
    moments_he[t,mn,i] - an array of proton moments where mn is a moment number, i is
    a spatial dimension (x, y, or z), and t is a frame number
    timing[t] - times (in proton gyroperiods) of the saved frames
    """
    # reading the file header
    with open(filename, 'rb') as file:
        print('Reading the file...')
        # setting up some parameters
        ihparm = 19     # number of parameters stored in the header
        ihdr = 19 + 8*kspi
        # reading the header information
        header = np.empty([ihdr], dtype=np.float32)
        header = np.fromfile(file, dtype=np.float32, count=len(header))
        # understanding how many frames are in the file
        nis = int(header[6])
        npts = np.empty([nis], dtype=np.int32)
        for _is in range (0, nis, 1):
            ihs = ihparm + 8*_is
            npts[_is] = np.int32(header[ihs + 7])
        nptsm=np.int32(np.amax(npts))
        xp = np.empty([nptsm,nis], dtype=np.float32)
        yp = np.empty([nptsm,nis], dtype=np.float32)
        vxp = np.empty([nptsm,nis], dtype=np.float32)
        vyp = np.empty([nptsm,nis], dtype=np.float32)
        vzp = np.empty([nptsm,nis], dtype=np.float32)
        framecur = 1
        anisotropies_p = []
        moments_p = []
        anisotropies_he = []
        moments_he = []
        timing = []
        while (True):
            timfrm = np.fromfile(file, dtype=np.float32, count=1)
            timep = np.fromfile(file, dtype=np.float32, count=1)
            for i in range (0, nis, 1):
                var = np.empty([5*npts[i]], dtype=np.float32)
                var = np.fromfile(file, dtype=np.float32, count=len(var))
                xp[0:npts[i],i] = np.copy(var[0::5])
                yp[0:npts[i],i] = np.copy(var[1::5])
                vxp[0:npts[i],i] = np.copy(var[2::5])
                vyp[0:npts[i],i] = np.copy(var[3::5])
                vzp[0:npts[i],i] = np.copy(var[4::5])
            if ( (framecur >= framestart) and ( (framecur <= frameend) or (frameend == -1) ) ):
                _anisotropies_p, _moments_p = calculate_anisotropy_moments(vxp, vyp, vzp, kspi_indexes = kspi_indexes_protons)
                _anisotropies_he, _moments_he = calculate_anisotropy_moments(vxp, vyp, vzp, kspi_indexes = kspi_indexes_he)
                anisotropies_p.append(_anisotropies_p)
                moments_p.append(_moments_p)
                anisotropies_he.append(_anisotropies_he)
                moments_he.append(_moments_he)
                timing.append(timep)
                print("MOMENTS READ FOR FRAME: " + str(framecur))
            framecur += 1
            if ( (timfrm[0] + header[11] >= header[3]) or ( (framecur > frameend) and (frameend != -1) ) ):
                break
    anisotropies_p = np.array(anisotropies_p, dtype=float)
    moments_p = np.array(moments_p, dtype=float)
    anisotropies_he = np.array(anisotropies_he, dtype=float)
    moments_he = np.array(moments_he, dtype=float)
    timing = np.array(timing, dtype=float)
    return anisotropies_p, moments_p, anisotropies_he, moments_he, timing
    
    
