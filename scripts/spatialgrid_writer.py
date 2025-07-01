from spatialdata.spatialdb.SimpleGridAscii import createWriter
from spatialdata.geocoords.CSCart import CSCart
import numpy as np



def write_spatialgrid_file(vertices, x_points, y_points, z_points, field_names, field_units, field_data, filename, dimension=2):
    '''
    Writes the spatialDB file. vertices is an nxdimension array containing the spatial coordinates of the spatialDB file.
    rupture_times, left_lateral_slips, opening_angles correspond to when the fault will rupture, how much it will slip,
    and how much it will open, respectively. filename is the name of the spatialDB file.
    '''
    cs = CSCart()
    cs._configure()
    cs.setSpaceDim(dimension)

    spatialdb_string = []
    for i in range(len(field_names)):
        spatialdb_string.append({'name':field_names[i], 'units':field_units[i], 'data':field_data[i]})
   #  print(spatialdb_string)
    writer = createWriter(filename)
    writer.write({'points': vertices,
                  'x': x_points,
                  'y': y_points,
                  'z': z_points,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': spatialdb_string})
# End File