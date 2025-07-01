import numpy as np
import sys
import os
os.chdir("/Users/danieldouglas/src/cig/Outer-Rise-Faulting-PyLith/CAM_outerrise/PUBLICATION_MODELS/scripts/")
import spatialdb_writer
os.chdir("/Users/danieldouglas/src/cig/Outer-Rise-Faulting-PyLith/CAM_outerrise/PUBLICATION_MODELS/")

time = 50

velocity_file = 'TEXT_FILES_FOR_DB_CONVERSION/CAM_velocity_' + str(time) + '.txt'
x_values = np.loadtxt(fname=velocity_file, usecols=0)
y_values = np.zeros(len(x_values))

########################################## TOP BOUNDARY VELOCITY ############################################
#############################################################################################################
velocity_y = np.loadtxt(fname=velocity_file, usecols=1)
velocity_x = np.zeros(len(velocity_y))
start_times = np.zeros(len(velocity_y))

vertices = np.array([x_values, y_values]).T

field_names = ["rate_amplitude_x", "rate_amplitude_y", "rate_start_time"]

field_units = ["m/year", "m/year", "year"]

field_data = np.array([[velocity_x], [velocity_y], [start_times]])

filename = "PYLITH_MODEL_FILES/velocity/top_boundary_" + str(time) + ".txt"

spatialdb_writer.write_spatialdb_file(vertices, field_names, field_units, field_data, filename, dimension=2)

######################################### BOTTOM BOUNDARY VELOCITY ##########################################
#############################################################################################################
bottom_boundary_y = 22.5e+3
vertices = np.array([x_values, y_values - bottom_boundary_y]).T

field_names = ["rate_amplitude_x", "rate_amplitude_y", "rate_start_time"]

field_units = ["m/year", "m/year", "year"]

field_data = np.array([[velocity_x], [velocity_y], [start_times]])

filename = "PYLITH_MODEL_FILES/velocity/bot_boundary_" + str(time) + ".txt"

spatialdb_writer.write_spatialdb_file(vertices, field_names, field_units, field_data, filename, dimension=2)


########################################## RIGHT BOUNDARY VELOCITY ##########################################
#############################################################################################################
x_values = np.loadtxt(fname=velocity_file, usecols=0)
x_boundary_location = 60.0e+3

x_index = abs(x_values - x_boundary_location).argmin()
x_values[x_values != x_values[x_index]] = x_values[x_index]
y_values = np.linspace(-22.6e+3, 0.1e+3, len(x_values))

velocity_y = np.loadtxt(fname=velocity_file, usecols=1)
velocity_y[velocity_y != velocity_y[x_index]] = velocity_y[x_index]
velocity_x = np.zeros(len(velocity_y))
start_times = np.zeros(len(velocity_y))

vertices = np.array([x_values, y_values]).T

field_names = ["rate_amplitude_x", "rate_amplitude_y", "rate_start_time"]

field_units = ["m/year", "m/year", "year"]

field_data = np.array([[velocity_x], [velocity_y], [start_times]])

filename = "PYLITH_MODEL_FILES/velocity/east_boundary_" + str(time) + ".txt"

spatialdb_writer.write_spatialdb_file(vertices, field_names, field_units, field_data, filename, dimension=2)

# End file
