import numpy as np

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
    
    for j in range(len(Y_flat)):
        for i in range(len(p)):
            linear_perm_fit[j] += p[i] * np.abs(Y_flat[j])**(len(p) - 1 - i)

    hatakeyama_perm = hatakeyama_et_al_permeability(Y_flat, P_0, k_0, gamma, solid_density)
    for m in range(len(hatakeyama_perm)):
        if linear_perm_fit[m] <= np.log10(hatakeyama_perm[m]):
            linear_perm_fit[m] = np.log10(hatakeyama_perm[m])

    return 10**linear_perm_fit

def my_old_permeability(Y_flat, depth_data, perm_data, additional_depths, additional_perms):
    
    depth_for_fit = Y_flat[np.where(Y_flat >= depth_data[-1])]
    
    p = np.polyfit(abs(depth_data), np.log10(perm_data), deg=7)
    data_depth_perm = np.zeros(len(depth_for_fit))
    
    for j in range(len(depth_for_fit)):
        for i in range(len(p)):
            data_depth_perm[j] += p[i] * np.abs(depth_for_fit[j])**(len(p) - 1 - i)
        
    background_depths = depth_for_fit
    background_permeability = data_depth_perm
    
    for k in range(len(additional_depths) - 1):
        
        p = np.polyfit(abs(np.array([additional_depths[k], additional_depths[k + 1]])), \
                       np.log10(np.array([additional_perms[k], additional_perms[k + 1]])), deg=1)

        depth_for_fit = Y_flat[np.where( (Y_flat < additional_depths[k]) & \
                                         (Y_flat >= additional_depths[k + 1]) )]

        perm_at_depths = np.zeros(len(depth_for_fit))
        
        for j in range(len(depth_for_fit)):
            for i in range(len(p)):
                perm_at_depths[j] += p[i] * np.abs(depth_for_fit[j])**(len(p) - 1 - i)
            
        background_depths = np.insert(background_depths, 0, depth_for_fit)
        background_permeability = np.insert(background_permeability, 0, perm_at_depths)
        
    return background_permeability

def my_old_permeability(Y_flat, depth_data, perm_data, additional_depths, additional_perms):
    
    depth_for_fit = Y_flat[np.where(Y_flat >= depth_data[-1])]
    
    p = np.polyfit(abs(depth_data), np.log10(perm_data), deg=7)
    data_depth_perm = np.zeros(len(depth_for_fit))
    
    for j in range(len(depth_for_fit)):
        for i in range(len(p)):
            data_depth_perm[j] += p[i] * np.abs(depth_for_fit[j])**(len(p) - 1 - i)
        
    background_depths = depth_for_fit
    background_permeability = data_depth_perm
    
    for k in range(len(additional_depths) - 1):
        
        p = np.polyfit(abs(np.array([additional_depths[k], additional_depths[k + 1]])), \
                       np.log10(np.array([additional_perms[k], additional_perms[k + 1]])), deg=1)

        depth_for_fit = Y_flat[np.where( (Y_flat < additional_depths[k]) & \
                                         (Y_flat >= additional_depths[k + 1]) )]

        perm_at_depths = np.zeros(len(depth_for_fit))
        
        for j in range(len(depth_for_fit)):
            for i in range(len(p)):
                perm_at_depths[j] += p[i] * np.abs(depth_for_fit[j])**(len(p) - 1 - i)
            
        background_depths = np.insert(background_depths, 0, depth_for_fit)
        background_permeability = np.insert(background_permeability, 0, perm_at_depths)
        
    return background_permeability