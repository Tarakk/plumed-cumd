# plumed-cumd (Spherical model)
Constant Chemical Potential Molecular Dynamics Simulations

This repository contains codes and input files required to run CmuMD nucleation simulations. 
The details of the method can be found in the following paper: 

Molecular Dynamics Simulations of Crystal Nucleation from Solution at Constant Chemical Potential
Tarak Karmakar, Pablo M. Piaggi, and Michele Parrinello*

Cite this: J. Chem. Theory Comput. 2019, 15, 12, 6923–6930
Publication Date:October 28, 2019
https://doi.org/10.1021/acs.jctc.9b00795

Route CmuMD method:
Perego, Salvalaglio, Parrinello J. Chem. Phys. 142 144113 (2015)
http://scitation.aip.org/content/aip/journal/jcp/142/14/10.1063/1.4917200


The codes can be found in the 'codes' folder, while the 'Examples' folder contains GROMACS-PLUMED input files.

To install GROMACS and PLUMED please visit the official websites,
GROMACS: http://www.gromacs.org/
PLUMED: https://github.com/plumed/plumed2

To install PLUMED with CmuMD, copy NshellSp.cpp to the plumed2/src/colvar/

To use the RefCV, copy RefCV.cpp to the plumed2/src/multicolvar/
Reference for RefCV: https://github.com/PabloPiaggi/JCP-2019.git


(It is to be noted that, the current version is a development version. In order to improve the method, changes will be made in the future release of the code.)
