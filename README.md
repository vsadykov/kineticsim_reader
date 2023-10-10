# kineticsim_reader
A set of routines for hybrid kinetic simulations.

Notes for Jupyter Notebook 'demo_simulationreader.ipynb':
This notebook demonstrates some capabilities of the kinetic simulation reader developed. In order to run it for the particular file, please change the following variables in the second cell of the notebook:
- simulationfile: please indicate the path on your machine to the simulation file of interest
- kspi: number of species for the simulation of interest
- kspi_pr: list of indexes for the proton species. Typically, kspi_pr = [0,1] in case kspi = 4, and kspi_pr = [0] in case kspi = 2
- kspi_he: list of indexes for the He species. Typically, kspi_he = [2,3] in case kspi = 4, and kspi_he = [1] in case kspi = 2


Other notes:
The 'processing_pipeline' folder contains the rountines which are used to process the simulations and produce and visualize the VDFs and calculate the momenta on the automated regime.
