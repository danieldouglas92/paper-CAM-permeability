import numpy as np
import spatialdb_writer
import cfg_writer
import os
from os import system as sys

########################################### WRITE THE SPATIALDB FILES FOR EACH OF THE FAULTS ###########################################
number_of_faults = 16 # Specify the number of faults
file_directory = 'spatialdb_files/' # Where the store the spatialdb files
sys('rm -rf ' + file_directory) # remove and remake directory for clean file output
sys('mkdir ' + file_directory)

# setup the spatial arrays for defining a 1D profile of how slip
# should be generated along the fault
y_values = np.linspace(0, -100e3, 200)
y_values = np.insert(y_values, 0, 1e3)
x_values = np.zeros(len(y_values))
min_slip_depth = 0

# Create vertices array for spatialdb files, as well as rupture_times and opening_slips.
# Don't want to include time-dependence or opening_angles, so set to 0.
vertices = np.array([x_values, y_values]).T
rupture_times = np.zeros(len(x_values))
opening_slips = np.zeros(len(x_values))

# loop through the number of faults, creating a spatialdb file for each
# fault. Here, the specified_slip and max depths change based on the index
# of the fault, but this can be generalized by creating arrays of max_slip_depth
# and specified_slip_points and indexing these arrays in the loop. left_lateral_slips
# is reset for each fault as a quadratic fitting 3 points.
for i in range(number_of_faults):
    filename = 'outer_rise_fault_' + str(i) + '.spatialdb'
    max_slip_depth = -15e3 - 1e3*(2*i)
    specified_slip_points = np.array([[min_slip_depth, 1.0], \
                                      [-5e3, 2.0 + i/10], \
                                      [max_slip_depth, 0.0]])
    left_lateral_slips = spatialdb_writer.fault_slip(x_values, y_values, max_slip_depth, min_slip_depth, specified_slip_points, polynomial_degree=2)
    
    file_path = file_directory + filename
    spatialdb_writer.test_writer(vertices, rupture_times, left_lateral_slips, opening_slips, file_path)



########################################### NOW WE WRITE THE .CFG FILE ############################################
file_name = 'step01_outer_rise_faulting.cfg' # name of .cfg file
sys('rm ' + file_name) # remove file, if it exists, for clean file write
cfg_file = open(file_name, 'a+')

# Write the header part of the fule
cfg_writer.general_definitions(cfg_file, file_name)

# Specify features
features = ['Static simulation', \
           'pylith.faults.FaultCohesiveKin', \
           'pylith.bc.DirichletTimeDependent', \
           'spatialdata.spatialdb.UniformDB', \
           'pylith.faults.KinSrcStep', \
           'pylith.bc.ZeroDB']

# Join features into a single string and write to file
feature_line = ", \n            ".join(features)
cfg_file.write('features = [' + feature_line + ']\n \n')

# Write output
cfg_file.write('[pylithapp] \n')
cfg_file.write('dump_parameters.filename = output/step01_outerrise-parameters.json \n')
cfg_file.write('problem.progress_monitor.filename = output/step01_outerrise-progress.txt \n \n')

cfg_file.write('problem.defaults.name = step01_outerrise \n \n')

# Create an empty array for fault interfaces, then fill it
interfaces = np.empty(0)
for i in range(number_of_faults):
    interfaces = np.append(interfaces, 'outer_rise_fault_' + str(i))

# Create array of label_values and edge_values, based on the gmsh script.
label_values = np.arange(201, 201 + number_of_faults, 1)
edge_values = np.arange(301, 301 + number_of_faults, 1)

# Write the faults with spatialDB to the .cfg file
cfg_writer.fault_spatialDB(cfg_file, interfaces, label_values, edge_values, file_directory)

# Specify the boundary information.
boundaries = ['bc_east', 'bc_west', 'bc_bottom']
labels = ['bndry_east', 'bndry_west', 'bndry_bot']
label_values = [12, 11, 14]
boundary_conditions = ['ZeroDB', 'ZeroDB', 'ZeroDB']
constraints = [ [0], [0], [1] ]

# Write boudnary conditions to .cfg file
cfg_writer.boundary_condition(cfg_file, boundaries, labels, label_values, boundary_conditions, constraints)

# Close .cfg file
cfg_file.write('# End of file \n')
cfg_file.close()

# End of file
