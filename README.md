This repository is assosciated with the publication:

*Constraining the Permeability and Outer-Rise Hydration at the Central America Margin*

by

Douglas, D.,
Aagaard, B.,
Naliboff, J.,
Naif, S.,

which is currently in review. Below is a broad summary of the contents of this repo.

# data_files
Contains various data files used for generating the model and constraining the models.

# files_for_conversion
Contains files which prescribe what the values of the velocity and the pressure should be on the model boundaries, but they are not in the correct format to be appled to the PyLith models.

# model_data_files
Contains the spatialdb and spatialgrid files which prescribe the pressure, the velocity, and the material parameters of the PyLith models. 
## pressure
Contains spatialdb files which prescribe the pressure on the PyLith model boundaries.
## velocity
Contains spatialdb files which prescribe the velocity on the PyLith model boundaries.
## spatialgrid
Contains the spatialgrid files which prescribe the material values at each point within the PyLith model domain. Each sub-directory within spatialgrid/ prescribe the material parameters for each of the four reference permeability models tested in the publication.

# model_parameter_files
Contains the .cfg files used to run the PyLith models, as well as the mesh file used to run the PyLith models. Each sub-directory corresponds to the four reference permeability models tested in this publication.

# paraview_state_files
Contains various state files used to generate the figures within the publication.

# scripts
Contains various python scripts used to generate the spatialdbfiles, spatialgridfiles, figures, and the world builder file used in this publication.
## scripts_for_figures
Various python notebooks which were used to create some of the figures in the publication.
## scripts_for_wb
Scripts for generating the worldbuilder file used in the publication.

# worldbuilder_output
Contains a .vtu file which contains the output used to determine the pore fluid flux within the publication.
