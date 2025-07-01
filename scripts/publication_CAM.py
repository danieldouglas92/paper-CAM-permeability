import numpy as np
import os
os.chdir("/Users/danieldouglas/src/cig/Outer-Rise-Faulting-PyLith/CAM_outerrise/PUBLICATION_MODELS/scripts/")
import initial_permeability
import spatialgrid_writer
import initial_porosity
os.chdir("/Users/danieldouglas/src/cig/Outer-Rise-Faulting-PyLith/CAM_outerrise/PUBLICATION_MODELS/")

x_max = 75.2e+3
x_min = -110.2e+3
dx = 0.1e+3

# x_max = 0.5e+3
# x_min = -200.5e+3
# dx = 0.25e+3

trench_x = 0.0e+3

y_max = 0.2e+3
y_min = -15.2e+3
dy = 0.1e+3

# y_max = 0.5e+3
# y_min = -200.5e+3
# dy = 0.25e+3

x_array = np.arange(x_min, x_max + dx, dx)
y_array = np.arange(y_min, y_max + dy, dy)
z_array = np.array([0.0])
X_matrix, Y_matrix = np.meshgrid(x_array, y_array)
X_flat, Y_flat = X_matrix.flatten(), Y_matrix.flatten()

vs = np.ones(len(X_flat)) * 3000 # 3000
vp = np.ones(len(X_flat)) * 5300 # 5300
solid_density = np.ones(len(X_flat)) * 2900 # 2950 # Turcotte & Schubert Appendix
fluid_density = np.ones(len(X_flat)) * 1030 # 1000
fluid_viscosity = np.ones(len(X_flat)) * 1e-3 # 1e-3
shear_modulus = np.ones(len(X_flat)) * 50e9 # 50e9 # Turcotte & Schubert Appendix
drained_bulk_modulus = np.ones(len(X_flat)) * 10.0e9 # * 10e9
fluid_bulk_modulus = np.ones(len(X_flat)) * 2.0e9 # * 2.0e9
solid_bulk_modulus = np.ones(len(X_flat)) * 80.0e9 # 80e9 # Turcotte & Schubert Appendix (Calculated)
# biot_coefficient = np.ones(len(X_flat)) * 1 - drained_bulk_modulus / solid_bulk_modulus

###################################### HERE WE SPECIFY THE LAYERED MATERIAL PARAMETERS ######################################
MORB_depth   = -1.6e3
gabbro_depth = -5.4e3
transition_thickness = -300

# Densities, shear moduli, and Young's moduli are taken from Table B5 from
# Turcotte and Schubert. Bulk Moduli are determined from the following wikilink
# https://en.wikipedia.org/wiki/Elastic_modulus

MORB_density   = 2900
MORB_shear_modulus = 30e9
MORB_young_modulus = 70e9
MORB_drain_modulus = 10e9
# MORB_drain_modulus = 75e9
MORB_bulk_modulus = MORB_shear_modulus * MORB_young_modulus / (3*(3*MORB_shear_modulus - MORB_young_modulus))
solid_density[Y_flat >= MORB_depth] = MORB_density
shear_modulus[Y_flat >= MORB_depth] = MORB_shear_modulus
solid_bulk_modulus[Y_flat >= MORB_depth] = MORB_bulk_modulus
drained_bulk_modulus[Y_flat >= MORB_depth] = MORB_drain_modulus

gabbro_density = 3100
gabbro_shear_modulus = 35e9
gabbro_young_modulus = 80e9
gabbro_drain_modulus = 10e9
# gabbro_drain_modulus = 75e9
gabbro_bulk_modulus = gabbro_shear_modulus * gabbro_young_modulus / (3*(3*gabbro_shear_modulus - gabbro_young_modulus))
solid_density[np.where( (Y_flat < MORB_depth) & (Y_flat >= gabbro_depth) )] = gabbro_density
shear_modulus[np.where( (Y_flat < MORB_depth) & (Y_flat >= gabbro_depth) )] = gabbro_shear_modulus
solid_bulk_modulus[np.where( (Y_flat < MORB_depth) & (Y_flat >= gabbro_depth) )] = gabbro_bulk_modulus
drained_bulk_modulus[np.where( (Y_flat < MORB_depth) & (Y_flat >= gabbro_depth) )] = gabbro_drain_modulus

dunite_density = 3300
dunite_shear_modulus = 65e9
dunite_young_modulus = 150e9
dunite_drain_modulus = 60e9
# dunite_drain_modulus = 75e9
dunite_bulk_modulus = dunite_shear_modulus * dunite_young_modulus / (3*(3*dunite_shear_modulus - dunite_young_modulus))
solid_density[Y_flat < gabbro_depth] = dunite_density
shear_modulus[Y_flat < gabbro_depth] = dunite_shear_modulus
solid_bulk_modulus[Y_flat < gabbro_depth] = dunite_bulk_modulus
drained_bulk_modulus[Y_flat < gabbro_depth] = dunite_drain_modulus

print(MORB_bulk_modulus/1e9, gabbro_bulk_modulus/1e9, dunite_bulk_modulus/1e9)

MORB_trans_inds  = np.where((Y_flat >= MORB_depth - transition_thickness) & (Y_flat < MORB_depth))
MORB_trans_dense = gabbro_density - MORB_density
MORB_trans_shear = gabbro_shear_modulus - MORB_shear_modulus
MORB_trans_bulk  = gabbro_bulk_modulus - MORB_bulk_modulus
MORB_trans_drain = gabbro_drain_modulus - MORB_drain_modulus

solid_density[MORB_trans_inds] = MORB_density + (np.abs(Y_flat[MORB_trans_inds] - MORB_depth)) / transition_thickness * MORB_trans_dense
shear_modulus[MORB_trans_inds] = MORB_shear_modulus + (np.abs(Y_flat[MORB_trans_inds] - MORB_depth)) / transition_thickness * MORB_trans_shear
solid_bulk_modulus[MORB_trans_inds] = MORB_bulk_modulus + (np.abs(Y_flat[MORB_trans_inds] - MORB_depth)) / transition_thickness * MORB_trans_bulk
drained_bulk_modulus[MORB_trans_inds] = MORB_drain_modulus + (np.abs(Y_flat[MORB_trans_inds] - MORB_depth)) / transition_thickness * MORB_trans_drain

gabbro_trans_inds  = np.where((Y_flat >= gabbro_depth - transition_thickness) & (Y_flat < gabbro_depth))
gabbro_trans_dense = dunite_density - gabbro_density
gabbro_trans_shear = dunite_shear_modulus - gabbro_shear_modulus
gabbro_trans_bulk  = dunite_bulk_modulus - gabbro_bulk_modulus
gabbro_trans_drain = dunite_drain_modulus - gabbro_drain_modulus

solid_density[gabbro_trans_inds] = gabbro_density + (np.abs(Y_flat[gabbro_trans_inds] - gabbro_depth)) / transition_thickness * gabbro_trans_dense
shear_modulus[gabbro_trans_inds] = gabbro_shear_modulus + (np.abs(Y_flat[gabbro_trans_inds] - gabbro_depth)) / transition_thickness * gabbro_trans_shear
solid_bulk_modulus[gabbro_trans_inds] = gabbro_bulk_modulus + (np.abs(Y_flat[gabbro_trans_inds] - gabbro_depth)) / transition_thickness * gabbro_trans_bulk
drained_bulk_modulus[gabbro_trans_inds] = gabbro_drain_modulus + (np.abs(Y_flat[gabbro_trans_inds] - gabbro_depth)) / transition_thickness * gabbro_trans_drain

biot_coefficient = 1 - drained_bulk_modulus / solid_bulk_modulus
biot_coefficient[np.where(biot_coefficient <= 0.101)] = 0.101
###################################### HERE WE SPECIFY THE BACKGROUND PERMEABILITY ######################################
perm_depths = np.array([0, -100, -200, -300, -400, -500, -600, -700, -800, -900, -1000, -1100, -1200])
perms = np.array([1e-11, 1e-13, 5e-14, 1e-15, 5e-16, 1e-16, 2e-17, 1e-17, 1e-17, 1e-17, 8e-18, 6e-18, 5e-18])

background_permeability = np.empty(0)

permeability_method = "hatakeyama_mid_perm"

if permeability_method == "power_law":
    a = -16.4
    b = 3
    background_permeability = initial_permeability.powerlaw_permeability(Y_flat, a, b)

elif permeability_method == "kuang_jiao":
    log_kr = -19
    log_ks = -12
    alpha = 1.8
    background_permeability = initial_permeability.kuang_jiao_permeability(Y_flat, log_kr, log_ks, alpha)

elif permeability_method == "hatakeyama":
    P_0 = 100
    gamma = 2.9e-2
    k_0 = 2.3e-21
    background_permeability = initial_permeability.my_permeability(Y_flat, perm_depths, perms, P_0, k_0, gamma, solid_density)

elif permeability_method == "hatakeyama_high_perm":
    P_0 = 100
    gamma = 3.3e-2
    k_0 = 3.5e-19
    background_permeability = initial_permeability.my_permeability(Y_flat, perm_depths, perms, P_0, k_0, gamma, solid_density)

elif permeability_method == "hatakeyama_mid_perm":
    P_0 = 100
    gamma = 2.7e-2
    k_0 = 6.1e-20
    background_permeability = initial_permeability.my_permeability(Y_flat, perm_depths, perms, P_0, k_0, gamma, solid_density)

elif permeability_method == "my_method":
    P_0 = 100
    gamma = 3.6e-2
    k_0 = 2e-23
    background_permeability = initial_permeability.my_permeability(Y_flat, perm_depths, perms, P_0, k_0, gamma, solid_density)

else:
    print("Permeability method not recognized. Exiting.")
    exit()
total_permeability = np.copy(background_permeability)

# Load in the depth-porosity distribution from Naif at the farthest point from the trench (x=100 km)
base_dir = "/Users/danieldouglas/src/cig/Outer-Rise-Faulting-PyLith/CAM_outerrise/PUBLICATION_MODELS/data_files/"
# porosity_file = np.loadtxt(fname=base_dir + "NAIF_porosity.txt")
porosity_file = np.loadtxt(fname=base_dir + "ave_NAIF_porosity.txt")
porosity = np.zeros(len(Y_flat))

for i in range(len(porosity)):
    depth_index = np.abs(Y_flat[i] + porosity_file[:, 0] * 1e3).argmin()
    porosity[i] = porosity_file[:, 1][depth_index]

# porosity[Y_flat < -5.4e3] = 0.0
# porosity *= 0
# porosity += 0.1
###################################### HERE WE SPECIFY THE FAULT ZONE PERMEABILITY ######################################
# Fault locations taken from Naif et al., 2015 
OUTER_RISE_FAULT_SURFACES_X = trench_x - np.array([3.5, 7, 9.5, 12.5, 15, 19.5, 22, 25.5, 27, 29, 32, 36, 39, 40.5, \
                                                   46, 49, 55.5, 58.5, 62, 67, 71, 78]) * 1e3
OUTER_RISE_FAULT_SURFACES_Y = np.zeros(len(OUTER_RISE_FAULT_SURFACES_X))
surface_points = np.array([OUTER_RISE_FAULT_SURFACES_X, OUTER_RISE_FAULT_SURFACES_Y]).T
fault_dips = np.ones(len(OUTER_RISE_FAULT_SURFACES_X)) * np.deg2rad(60)
buried_edge_Y = np.ones(len(OUTER_RISE_FAULT_SURFACES_X)) * y_min / 2
buried_edge_X = surface_points[:, 0] - buried_edge_Y / np.tan(fault_dips)
buried_points = np.array([buried_edge_X, buried_edge_Y]).T
fault_perm_factor = 100
fault_zone_width  = 1000

for i in range(len(surface_points)):
    distance_from_fault = np.abs( (surface_points[i][0] - buried_points[i][0]) * (buried_points[i][1] - Y_flat) - \
                                (buried_points[i][0] - X_flat) * (surface_points[i][1] - buried_points[i][1]) ) / \
                          np.sqrt( (surface_points[i][0] - buried_points[i][0])**2 + (surface_points[i][1] - buried_points[i][1])**2 )
    
    within_fault_zone_index = np.where( (distance_from_fault < (fault_zone_width/2)) & (Y_flat >= buried_points[i][1]) )
    outside_fault_zone_index = np.where( (distance_from_fault >= (fault_zone_width/2)) & (Y_flat < buried_points[i][1]) )

    distance_from_fault[outside_fault_zone_index] = fault_zone_width/2

    total_permeability[within_fault_zone_index] = total_permeability[within_fault_zone_index] * fault_perm_factor - \
                                                  (distance_from_fault[within_fault_zone_index] / (fault_zone_width/2)) * \
                                                   total_permeability[within_fault_zone_index] * fault_perm_factor

background_permeability[np.where(background_permeability > 1e-12)] = 1e-12
background_permeability[np.where(background_permeability < 1e-30)] = 1e-30

total_permeability[np.where(total_permeability > 1e-12)] = 1e-12
total_permeability[np.where(total_permeability < 1e-30)] = 1e-30

vertices = np.array([X_flat, Y_flat]).T

reference_stress_xx = solid_density * 9.81 * Y_flat
reference_stress_yy = np.copy(reference_stress_xx)
reference_stress_zz = np.copy(reference_stress_xx)
reference_stress_xy = np.zeros(len(reference_stress_xx))

reference_strain_xx = np.zeros(len(Y_flat))
reference_strain_yy = np.copy(reference_strain_xx)
reference_strain_zz = np.copy(reference_strain_xx)
reference_strain_xy = np.copy(reference_strain_xx)

body_force_x        = np.zeros(len(Y_flat))
body_force_y        = fluid_density * -9.81

field_names = ["solid_density", "fluid_density", "fluid_viscosity", \
               "porosity", "shear_modulus", "drained_bulk_modulus", "biot_coefficient", 
               "fluid_bulk_modulus", "solid_bulk_modulus", "isotropic_permeability", \
               "reference_stress_xx", "reference_stress_yy", "reference_stress_zz", \
               "reference_stress_xy", "reference_strain_xx", "reference_strain_yy", \
               "reference_strain_zz", "reference_strain_xy", "body_force_x", "body_force_y"]

field_units = ["kg/m**3", "kg/m**3", "Pa*s", \
               "none", "Pa", "Pa", "none", \
               "Pa", "Pa", "m**2", "Pa", "Pa", "Pa", \
               "Pa", "none", "none", "none", "none", "Pa", "Pa"]

field_data = np.array([[solid_density], [fluid_density], [fluid_viscosity],\
                       [porosity], [shear_modulus], [drained_bulk_modulus], [biot_coefficient], \
                       [fluid_bulk_modulus], [solid_bulk_modulus], [background_permeability],
                       [reference_stress_xx], [reference_stress_yy], [reference_stress_zz], \
                       [reference_stress_xy], [reference_strain_xx], [reference_strain_yy], \
                       [reference_strain_zz], [reference_strain_xy], [body_force_x], [body_force_y]])

field_data_faults= np.array([[solid_density], [fluid_density], [fluid_viscosity],\
                             [porosity], [shear_modulus], [drained_bulk_modulus], [biot_coefficient], \
                             [fluid_bulk_modulus], [solid_bulk_modulus], [total_permeability],
                             [reference_stress_xx], [reference_stress_yy], [reference_stress_zz], \
                             [reference_stress_xy], [reference_strain_xx], [reference_strain_yy], \
                             [reference_strain_zz], [reference_strain_xy], [body_force_x], [body_force_y]])

faults_string =  str(int(fault_zone_width)) + "m_width_" + str(fault_perm_factor) + "_reduction.txt"
no_fault_file = "PYLITH_MODEL_FILES/spatialgrid/" + permeability_method + "/no_faults.txt"
fault_file = "PYLITH_MODEL_FILES/spatialgrid/" + permeability_method + '/' + faults_string
# no_fault_file = "PYLITH_MODEL_FILES/spatialgrid/junk.txt"
# fault_file = "PYLITH_MODEL_FILES/spatialgrid/body_force.txt"
spatialgrid_writer.write_spatialgrid_file(vertices, x_array, y_array, z_array, field_names, field_units, field_data, no_fault_file, dimension=2)
spatialgrid_writer.write_spatialgrid_file(vertices, x_array, y_array, z_array, field_names, field_units, field_data_faults, fault_file, dimension=2)
