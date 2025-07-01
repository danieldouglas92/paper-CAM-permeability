import numpy as np

poly_order = 4
if poly_order == 2:
    extrusive_porosity = np.log10(2.2) # 12.2 %
    dikes_porosity = np.log10(0.055) # 5 %
    gabbro_porosity = np.log10(0.0105) # 1.05 %
    moho_porosity = np.log10(0.0045) # 0.45 %

if poly_order == 3:
    extrusive_porosity = np.log10(0.65) # 12.2 %
    dikes_porosity = np.log10(0.05) # 5 %
    gabbro_porosity = np.log10(0.0105) # 1.05 %
    moho_porosity = np.log10(0.0045) # 0.45 %

if poly_order == 4:
    extrusive_porosity = np.log10(0.11) # 11 %
    dikes_porosity = np.log10(0.055) # 5 %
    gabbro_porosity = np.log10(0.0105) # 1.05 %
    moho_porosity = np.log10(0.0045) # 0.45 %

extrusive_thickness = 1 # km
dike_thickness = 1 # km
gabbro_thickness = 4 # km

porosity_points = np.array([extrusive_porosity, \
                            dikes_porosity, \
                            gabbro_porosity, \
                            moho_porosity, \
                            np.log10(0.001), \
                            np.log10(0.0001), \
                            np.log10(0.00001)])

# extrusive_porosity = np.log10(0.182) # 12.2 %
# dikes_porosity = np.log10(0.1) # 5 %
# gabbro_porosity = np.log10(0.0205) # 1.05 %
# moho_porosity = np.log10(0.009) # 0.45 %

# extrusive_thickness = 1 # km
# dike_thickness = 1 # km
# gabbro_thickness = 4 # km

# porosity_points = np.array([extrusive_porosity, \
#                             dikes_porosity, \
#                             gabbro_porosity, \
#                             moho_porosity, \
#                             np.log10(0.002), \
#                             np.log10(0.0002), \
#                             np.log10(0.00002)])

depth_points = np.array([0, \
                         extrusive_thickness, \
                         extrusive_thickness + dike_thickness, \
                         extrusive_thickness + dike_thickness + gabbro_thickness, \
                         10, \
                         20, \
                         30])

def initial_naif_porosity(y_array):
    poro_poly_coeff = np.polyfit(depth_points, porosity_points, deg=poly_order)
    poro_fit = np.zeros(len(y_array))
    for j in range(len(y_array)):
        for i in range(len(poro_poly_coeff)):
            poro_fit[j] += poro_poly_coeff[i] * np.abs(y_array[j]/1e3)**(len(poro_poly_coeff) - 1 - i)

    return 10**poro_fit
