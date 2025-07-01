import numpy as np
import os

def general_definitions(cfg_file, file_name, description='', authors='', keywords=''):
    '''
    Writes the header (usually not changing) part of the .cfg file. Optional inputs are the 
    description of what the .cfg file does, keywords, and authors
    '''
    cfg_file.write('[pylithapp.metadata] \n \n')
    cfg_file.write('base = [pylithapp.cfg] \n')
    cfg_file.write('description = [' + description +  ']\n')
    cfg_file.write('authors = [' + authors + ']\n')
    cfg_file.write('keywords = [' + keywords + ']\n')
    cfg_file.write('arguments = [' + file_name + ']\n')
    cfg_file.write('version = 2.0.0 \n')
    cfg_file.write('pylith_version = [>=3.0, <4.0] \n \n')
    return



def fault_uniformDB(cfg_file, interfaces, label_values, edge_values, fault_parameters):
    '''
    Writes out the fault section of the .cfg file as a uniformDB. Takes the cfg_file path, an interfaces array
    containing the names of the faults, label_values and edge_values which contain the tags for the fault and the 
    buried edges, respectively, and fault_parameters, which is an nx3 dimensional array containing the slip parameters
    for each fault.
    '''
    cfg_file.write('# ---------------------------------------------------------------------- \n')
    cfg_file.write('# faults \n')
    cfg_file.write('# ---------------------------------------------------------------------- \n')
    cfg_file.write('[pylithapp.problem]\n \n')
    interface_line = ", ".join(interfaces)
    cfg_file.write('interfaces = [' + interface_line + ']\n \n')
    for i in range(len(interfaces)):
        cfg_file.write('[pylithapp.problem.interfaces.' + interfaces[i] + ']\n')
        cfg_file.write('label = ' + interfaces[i] + '\n')
        cfg_file.write('label_value = ' + str(label_values[i]) + '\n')
        cfg_file.write('edge = edge_' + interfaces[i] + '\n')
        cfg_file.write('edge_value = ' + str(edge_values[i]) + '\n')

        cfg_file.write('observers.observer.data_fields = [slip] \n \n')

        cfg_file.write('[pylithapp.problem.interfaces.' + interfaces[i] + '.eq_ruptures.rupture] \n')

        cfg_file.write('db_auxiliary_field = spatialdata.spatialdb.UniformDB\n')
        cfg_file.write('db_auxiliary_field.description = Test \n')
        cfg_file.write('db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening] \n')
        cfg_file.write('db_auxiliary_field.data = ' + str(fault_parameters[i]) + '\n \n \n')

    return



def fault_spatialDB(cfg_file, interfaces, label_values, edge_values, spatialdb_file_directory, single_file=False):
    '''
    Writes out the fault section of the .cfg file as a uniformDB. Takes the cfg_file path, an interfaces array
    containing the names of the faults, label_values and edge_values which contain the tags for the fault and the 
    buried edges, respectively, spatialdb_file_directory, which is the path to the directory where the spatial_db files
    are written to using spatialdb_writer.py. Lastly, if single_file=True, then a single spatialdb file will be applied to
    all faults, otherwise there must be n spatialdb_files for n faults.
    '''
    cfg_file.write('# ---------------------------------------------------------------------- \n')
    cfg_file.write('# faults \n')
    cfg_file.write('# ---------------------------------------------------------------------- \n')
    cfg_file.write('[pylithapp.problem]\n \n')
    interface_line = ", ".join(interfaces)
    cfg_file.write('interfaces = [' + interface_line + ']\n \n')

    filelist = os.listdir(spatialdb_file_directory)
    for i in range(len(interfaces)):
        cfg_file.write('[pylithapp.problem.interfaces.' + interfaces[i] + ']\n')
        cfg_file.write('label = ' + interfaces[i] + '\n')
        cfg_file.write('label_value = ' + str(label_values[i]) + '\n')
        cfg_file.write('edge = edge_' + interfaces[i] + '\n')
        cfg_file.write('edge_value = ' + str(edge_values[i]) + '\n')

        cfg_file.write('observers.observer.data_fields = [slip] \n \n')

        cfg_file.write('[pylithapp.problem.interfaces.' + interfaces[i] + '.eq_ruptures.rupture] \n')

        cfg_file.write('db_auxiliary_field = spatialdata.spatialdb.SimpleDB \n')
        cfg_file.write('db_auxiliary_field.description = Test \n')
        if single_file:
            spatial_db_file_path = spatialdb_file_directory + filelist
            cfg_file.write('db_auxiliary_field.iohandler.filename = ' + spatial_db_file_path + '\n')

        else:
            spatial_db_file_path = spatialdb_file_directory + filelist[i]
            cfg_file.write('db_auxiliary_field.iohandler.filename = ' + spatial_db_file_path + '\n')
        cfg_file.write('db_auxiliary_field.query_type = linear \n \n \n')

    return


def boundary_condition(cfg_file, boundaries, labels, label_values, boundary_conditions, constraints):
    '''
    Writes out the boundary conditions part of the .cfg file. cfg_file is the path to the .cfg file, boundaries
    is an array containing the names of the boundaries, labels is an array containing the names of the boundaries 
    defined in gmsh, label_values is an array containing the tags of the boundaries, boundary_conditions is an array
    of conditions for each boundary (Neumann, Dirichlet, etc.), and constraints is an array of which spatial variable
    the boundary condition applies to.
    '''
    cfg_file.write('# ---------------------------------------------------------------------- \n')
    cfg_file.write('# boundary conditions \n')
    cfg_file.write('# ---------------------------------------------------------------------- \n')
    cfg_file.write('[pylithapp.problem] \n \n')
    boundary_line = ", ".join(boundaries)

    cfg_file.write('bc = [' + str(boundary_line) + ']\n \n')
    for i in range(len(boundaries)):
        cfg_file.write('[pylithapp.problem.bc.' + boundaries[i] + ']\n')
        cfg_file.write('label = ' + str(labels[i]) + '\n')
        cfg_file.write('label_value = ' + str(label_values[i]) + '\n')
        cfg_file.write('constrained_dof = ' + str(constraints[i]) + '\n')
        cfg_file.write('db_auxiliary_field = pylith.bc.' + boundary_conditions[i] + '\n')
        cfg_file.write('db_auxiliary_field.description = Testing automated BC \n \n \n')
    return
