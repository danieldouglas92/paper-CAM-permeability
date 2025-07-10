This repository is assosciated with the publication

*Constraining the Permeability and Outer-Rise Hydration at the Central America Margin*

by

Douglas, D.,
Aagaard, B.,
Naliboff, J.,
Naif, S.,

Our numerical simulations were run using the open source geodynamics software PyLith ([[https://aspect.geodynamics.org/](https://pylith.readthedocs.io/en/latest/)]([https://aspect.geodynamics.org/](https://pylith.readthedocs.io/en/latest/))) 
, as well as the initial conditions software the Geodynamic World Builder ([https://gwb.readthedocs.io/en/latest/](https://gwb.readthedocs.io/en/latest/)). 
Specifically, we utilize the PyLith version 4.2.0 and the Geodynamic World Builder version 1.0.0.

# Overview
All scripts required to reproduce the results and the figures in Douglas et al. (2025), "Constraining the Permeability and Outer-Rise Hydration at the Central America Margin", submitted to JGR: Solid Earth. 
The *.pvsm files can be run by specifying the path to the solution files on your local machine upon opening the state files. The summary of each directory is outlined below.

## data_files
Contains various data files which were used to construct the initial conditions of the PyLith model, as well as data files used to constrain the PyLith models from Naif et al., 2015. Additionally, within the directory Slab2_CAM contains
the data files from Slab2.0 (Hayes et al., 2018) used to construct the geometry in the Geodynamic World Builder.

## files_for_conversion
Contains various data files before they were converted into a format that were compatible with PyLith.

## model_data_files
Contains the PyLith compatible data files which specify the boundary conditions and the initial conditions within the PyLith models. 
### pressure
Specifies the pressure boundary condition applied to the top boundary of the PyLith models.
### velocity
Specifies the velocity boundary condition applied to the top, right, and bottom boundary of the PyLith models.
### spatialgrid
Specifies the initial conditions for the PyLith models for each of the four reference peremeability models.

## model_parameter_files
Contains the script used to generate the mesh for the PyLith models, and also the input files used to run the PyLith models. Each of the four sub-directories contain the input files used for each of the four reference permeability models.

## paraview_state_files
Contains paraview state files used to generate some of the figures within the manuscript.

## scripts
Contains various python scripts used to either create the figures within the manuscript (within `scripts_for_figures`), to generate the worldbuilder file (within `scripts_for_wb`), and to generate the data files for the PyLith model.

## worldbuilder_output
Contains the worldbuilder file used for the manuscript, and also the output from the worldbuilder used to determine the Darcy velocity for the fluid flux calculation.
