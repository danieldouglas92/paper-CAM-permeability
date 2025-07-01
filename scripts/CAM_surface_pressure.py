import numpy as np
import sys
import os
import os
os.chdir("/Users/danieldouglas/src/cig/Outer-Rise-Faulting-PyLith/CAM_outerrise/PUBLICATION_MODELS/scripts/")
import spatialdb_writer
os.chdir("/Users/danieldouglas/src/cig/Outer-Rise-Faulting-PyLith/CAM_outerrise/PUBLICATION_MODELS/")

pressure_file = 'TEXT_FILES_FOR_DB_CONVERSION/CAM_fluid_pressure.txt'
x_values = np.loadtxt(fname=pressure_file, usecols=0)
y_values = np.zeros(len(x_values))

pressure = np.loadtxt(fname=pressure_file, usecols=1)

vertices = np.array([x_values, y_values]).T

field_names = ["initial_amplitude"]

field_units = ["Pa"]

field_data = np.array([[pressure]])

filename = "PYLITH_MODEL_FILES/pressure/surface_fluid_pressure.txt"

spatialdb_writer.write_spatialdb_file(vertices, field_names, field_units, field_data, filename, dimension=2)

# End file
