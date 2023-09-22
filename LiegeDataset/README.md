The ALPS software package[1] is a series of Python routines that can be used to study polaronic properties of materials. It does so by extracting already existing information from databases and by employing different Frohlich models to compute the Zero-point renormalisation (ZPR) and the electron-phonon coupling constant. In case any necessary parcel of data is missing (e.g. effective masses, dielectric constants, phonon displacements), the ALPS software package can interface with Abinit and its tools to compute and add them to its own local database. 

This package is distributed through the following paths:
src/                  - python scripts for the polaronic properties calculations using the Frohlich models
Repository/           - folders containing the files used as inputs for the calculations
abinit_effmasses/     - files containing ABINIT outputs used as the effective masses inputs for the calculations
eff_masses/           - backup files with effective masses using the Materials Project calculations
phonon/               - files with the data obtained by P. Guido et al. with the dielectric constant and phonon properties for the standard Frohlich model
Results/              - data of the calculations made for the paper (?)
q-random-1000/        - data for the average using the generalized Frohlich model 
atoms/                - data about atoms necessary for calculations


This package was developed by (Last name, First name, current affiliation, contact email):
Melo, Pedro (VASP), pedro.melo@vasp.at
Abreu, Joao (Université de Liège), jcabreu@uliege.be 

Future maintenance by (Last name, First name, current affiliation, contact email):
Verstraete, Matthieu (Université de Liège), matthieu.verstraete@uliege.be


NOTE

The effective masses were calculated using ABINIT[2] and the script in this software to prepare the calculations. As a backup the user can use the ones from Materials Project[3]
The electronic dielectric constant and all phonon properties used in the standard Frohlich model were obtained from database [4] by Guido Petretto et al.

The rights and permissions for the python scripts are noted in src/LICENSE and the rest of the data is attributed the license CC BY 4.0

Citations

[1] Melo, P. M. M. C. et al, High-throughput analysis of Fröhlich-type polaron models, arXiv:2207.00364 (2022)
[2] Gonze, X. et al., The Abinit project: Impact, environment and recent developments, Comput. Phys. Commun. 248, 107042 (2020)
[3] Ricci, F. et al., An ab initio electronic transport database for inorganic materials, doi: 10.1038/sdata.2017.85 (2017)
[4] Petretto, G. et al., High-throughput Density-Functional Perturbation Theory phonons for inorganic materials, doi: 10.6084/m9.figshare.c.3938023.v1 (2018)

