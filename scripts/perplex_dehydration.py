import numpy as np
import matplotlib.pyplot as plt
import pygmt

def perplex_dehydration(profile_P, profile_T, \
                        line_spacing, profile_hydration, profile_dehydration, \
                        Perp_P_mat, Perp_T_mat, \
                        sediment_thickness, sediment_perp, sediment_density, \
                        MORB_thickness, MORB_perp, MORB_density, \
                        gabbro_thickness, gabbro_perp, gabbro_density, \
                        peridotite_thickness, peridotite_perp, peridotite_density):
                            
    Perp_T_vals = Perp_T_mat[0, :]
    Perp_P_vals = Perp_P_mat[:, 0]
    profile_depth = 0.0
    for i in range(len(profile_hydration)):
        profile_hydration[i] = np.zeros(len(profile_T[i]))
        profile_dehydration[i] = np.zeros(len(profile_T[i]))
        for j in range(len(profile_hydration[i])):
            point_P = profile_P[i][j]
            point_T = profile_T[i][j]
    
            T_index = abs(point_T - Perp_T_vals).argmin()
            P_index = abs(point_P - Perp_P_vals).argmin()

            if profile_depth < sediment_thickness:
                Perp_water = sediment_perp[P_index][T_index] / 100 # convert from wt% to mass fraction
                density = sediment_density

            elif profile_depth < sediment_thickness + MORB_thickness:
                Perp_water = MORB_perp[P_index][T_index] / 100 # convert from wt% to mass fraction
                density = MORB_density

            elif profile_depth < sediment_thickness + MORB_thickness + gabbro_thickness:
                Perp_water = gabbro_perp[P_index][T_index] / 100 # convert from wt% to mass fraction
                density = gabbro_density

            elif profile_depth < sediment_thickness + MORB_thickness + gabbro_thickness + peridotite_thickness:
                Perp_water = peridotite_perp[P_index][T_index] / 100 # convert from wt% to mass fraction
                density = peridotite_density
            
            if j == 0:
                profile_hydration[i][j] = Perp_water
            else:
                if Perp_water < profile_hydration[i][j - 1]:
                    profile_hydration[i][j] = Perp_water
                    profile_dehydration[i][j] = (profile_hydration[i][j - 1] - Perp_water) * density
                else:
                    profile_hydration[i][j] = profile_hydration[i][j - 1]
        profile_depth += line_spacing

    return profile_hydration, profile_dehydration


def vertical_distance(slab_surface_x, slab_surface_y, profile_x, profile_y):
    closest_surface_value = np.abs(profile_x - slab_surface_x).argmin()
    y_distance = slab_surface_y[closest_surface_value] - profile_y
    return y_distance

def total_arc_length(x_input, y_input):
    arclength = 0
    for i in range(len(x_input) - 1):
        arclength += np.sqrt((x_input[i + 1] - x_input[i])**2 + (y_input[i + 1] - y_input[i])**2)
    return arclength

def incremential_arc_length(x_input, y_input):
    arclength = 0
    inc_arclength = np.zeros(len(x_input))
    for i in range(len(x_input) - 1):
        arclength += np.sqrt((x_input[i + 1] - x_input[i])**2 + (y_input[i + 1] - y_input[i])**2)
        inc_arclength[i + 1] = arclength
    return inc_arclength

def surface_flux_calculation(parallel_slab_arclengths, fluid_advection_in_solid, dehydration_container):

    surface_flux = np.zeros(len(parallel_slab_arclengths[0]))
    for i in range(len(parallel_slab_arclengths)):
        for j in range(len(parallel_slab_arclengths[i])):
            arclength_at_slab_surface = parallel_slab_arclengths[i][j] + fluid_advection_in_solid[i][j]
            indices_of_slab_surface = np.argsort(np.abs(parallel_slab_arclengths[0] - arclength_at_slab_surface))
            # the advected fluid will land somewhere between two points on the slab surface, choose the larger one to bin the flux
            index_for_flux = np.max(indices_of_slab_surface[:2])
            if arclength_at_slab_surface < parallel_slab_arclengths[0][-1]:  
                surface_flux[index_for_flux] += dehydration_container[i][j]

    return surface_flux


def surface_flux_calculation_no_v(parallel_slab_x, dehydration_container):

    surface_flux_no_v = np.zeros(len(parallel_slab_x[0]))
    for i in range(len(parallel_slab_x)):
        for j in range(len(parallel_slab_x[i])):
            index_at_slab_surface = np.abs(parallel_slab_x[0] - parallel_slab_x[i][j]).argmin()
            surface_flux_no_v[index_at_slab_surface] += dehydration_container[i][j]

    return surface_flux_no_v


def hatakeyama_et_al_permeability(Y_flat, P_0, k_0, gamma, solid_density):
    P_vals = np.abs(9.81 * solid_density * Y_flat) / 1e6
    k_vals = k_0 * np.exp(-gamma * (P_vals - P_0))
    return k_vals

def powerlaw_permeability(Y_flat, a, b):
    log_permeability = a - b * np.log10(abs(Y_flat/1e3) + 1)
    return 10**log_permeability

def kuang_jiao_permeability(Y_flat, log_kr, log_ks, alpha):
    log_permeability = log_kr + (log_ks - log_kr) * np.power((1 + abs(Y_flat/1e3)), -alpha)
    return 10**log_permeability

def my_permeability(Y_flat, depth_data, perm_data, P_0, k_0, gamma, solid_density):

    p = np.polyfit(abs(depth_data), np.log10(perm_data), deg=1)
    linear_perm_fit = np.zeros(len(Y_flat))
    print(p)

    for j in range(len(Y_flat)):
        for i in range(len(p)):
            linear_perm_fit[j] += p[i] * np.abs(Y_flat[j])**(len(p) - 1 - i)
    
    hatakeyama_perm = hatakeyama_et_al_permeability(Y_flat, P_0, k_0, gamma, solid_density)
    for m in range(len(hatakeyama_perm)):
        if linear_perm_fit[m] <= np.log10(hatakeyama_perm[m]):
            linear_perm_fit[m] = np.log10(hatakeyama_perm[m])

    return 10**linear_perm_fit