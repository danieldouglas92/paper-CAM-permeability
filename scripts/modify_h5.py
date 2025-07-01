import h5py
import numpy as np
import os
from os import system as sys

base_dir = '/Users/danieldouglas/src/cig/Outer-Rise-Faulting-PyLith/CAM_outerrise/CDM_presentation_models/'
h5_info_file_path = base_dir + 'my_method/1000m_fault_width/TEST/step03_faults_flexure-subducting_info.h5'
h5_full_file_path = base_dir + 'my_method/1000m_fault_width/TEST/step03_faults_flexure-subducting.h5'

h5_initial = h5py.File(h5_info_file_path)
initial_porosity = h5_initial["vertex_fields"]["porosity"][:]

h5_total = h5py.File(h5_full_file_path, 'r+')
total_porosity = h5_total["vertex_fields"]["porosity"][:]
porosity_change = total_porosity - initial_porosity

vertices = h5_total["geometry/vertices"]
x_vals = vertices[:, 0]
y_vals = vertices[:, 1]

########################### TEMPERATURE CALC FROM PS 77 FOR WATER RESISTIVITY CALCULATION ###########################
T0 = 273
T1 = 1623
litho_thickness = 125e3
kappa = 1e-6
s2yr = 60 * 60 * 24 * 365.25
plate_age = 24e6 * s2yr # 24 Myr plate age for CAM
plate_sum = np.zeros(len(y_vals))

n = 1
while n < 20:
    plate_sum += 1/n * np.exp(-kappa * n**2 * np.pi**2 * plate_age / litho_thickness**2) * \
                       np.sin(n * np.pi * np.abs(y_vals) / litho_thickness)
    n += 1

T = T0 + (T1 - T0) * (np.abs(y_vals) / litho_thickness + 2/np.pi * plate_sum) - 273

# Values from Naif et al., 2015
prefac = 2.903916
a = 2.97175e-2
b = 1.5551e-4
c = -6.7e-7

seawater_resistivity = 1 / (prefac * (1 + a * T + b * T**2 + c * T**3))
seawater_resistivity[np.where(seawater_resistivity <= 0)] = 0.4
seawater_resistivity[np.where(seawater_resistivity >= 0.4)] = 0.4

# seawater_resistivity = 0.2
bulk_resistivity = np.zeros(np.shape(total_porosity))

for t in range(len(bulk_resistivity)):
    for i in range(len(y_vals)):
        test = np.array([total_porosity[t][i][0], 1e-6])
        bulk_resistivity[t][i] = seawater_resistivity[i] / (np.max(test))**2


h5_total.create_dataset("vertex_fields/porosity_change", data=porosity_change)

# h5_total.create_dataset("vertex_fields/P_wave_velocity", data=P_wave_model)

h5_total.create_dataset("vertex_fields/bulk_resistivity", data=bulk_resistivity)

sys("pylith_genxdmf " + base_dir)

#################### CALCULATION FOR DETERMINING THE SEISMIC VELOCITY BASED ON THE POROSITY
# P_wave_crust_velocity = np.loadtxt(fname='/Users/danieldouglas/src/cig/Outer-Rise-Faulting-PyLith/scripts/data_files/P_wave_crust.txt', \
#                                    comments='#', delimiter=',')
# shallow_velocity = 5.5e3 # m/s
# P_wave_water_velocity = 1.5e3 # m/s
# P_wave_model = np.zeros(np.shape(total_porosity))
# # I need to modify this loop because the other variables here are for ALL 50 time steps, not just the last time step.

# for t in range(len(P_wave_model)):
#     for i in range(len(P_wave_model[t])):
#         if y_vals[i] >= P_wave_crust_velocity[:, 0][0]:
#             P_wave_model[t][i] = P_wave_water_velocity * shallow_velocity / \
#                             (total_porosity[t][i] * (shallow_velocity - P_wave_water_velocity) + P_wave_water_velocity)

#         elif P_wave_crust_velocity[:, 0][0] >= y_vals[i] >= P_wave_crust_velocity[:, 0][1]:
#             P_wave_model[t][i] = P_wave_water_velocity * P_wave_crust_velocity[:, 1][0] / \
#                             (total_porosity[t][i] * (P_wave_crust_velocity[:, 1][0] - P_wave_water_velocity) + P_wave_water_velocity)
            
#         elif P_wave_crust_velocity[:, 0][1] >= y_vals[i] >= P_wave_crust_velocity[:, 0][2]:
#             P_wave_model[t][i] = P_wave_water_velocity * P_wave_crust_velocity[:, 1][1] / \
#                             (total_porosity[t][i] * (P_wave_crust_velocity[:, 1][1] - P_wave_water_velocity) + P_wave_water_velocity)
            
#         elif P_wave_crust_velocity[:, 0][2] >= y_vals[i] >= P_wave_crust_velocity[:, 0][3]:
#             P_wave_model[t][i] = P_wave_water_velocity * P_wave_crust_velocity[:, 1][2] / \
#                             (total_porosity[t][i] * (P_wave_crust_velocity[:, 1][2] - P_wave_water_velocity) + P_wave_water_velocity)
            
#         elif P_wave_crust_velocity[:, 0][3] >= y_vals[i] >= P_wave_crust_velocity[:, 0][4]:
#             P_wave_model[t][i] = P_wave_water_velocity * P_wave_crust_velocity[:, 1][3] / \
#                             (total_porosity[t][i] * (P_wave_crust_velocity[:, 1][3] - P_wave_water_velocity) + P_wave_water_velocity)
            
#         elif y_vals[i] <= P_wave_crust_velocity[:, 0][4]:
#             P_wave_model[t][i] = P_wave_water_velocity * P_wave_crust_velocity[:, 1][4] / \
#                             (total_porosity[t][i] * (P_wave_crust_velocity[:, 1][4] - P_wave_water_velocity) + P_wave_water_velocity)
