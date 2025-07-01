import numpy as np

T0 = 273
T1 = 1623
litho_thickness = 125e3
kappa = 1e-6 
s2yr = 60 * 60 * 24 * 365.25
plate_age = 24e6 * s2yr

y_vals = np.linspace(0, litho_thickness, 5000)
plate_sum = np.zeros(len(y_vals))

n = 1
while n < 20:
    plate_sum += 1/n * np.exp(-kappa * n**2 * np.pi**2 * plate_age / litho_thickness**2) * \
                       np.sin(n * np.pi * y_vals / litho_thickness)
    n += 1
T = T0 + (T1 - T0) * (y_vals / litho_thickness + 2/np.pi * plate_sum) - 273

prefac = 2.903916
a = 2.97175e-2
b = 1.5551e-4
c = -6.7e-7

seawater_resistivity = 1 / (prefac * (1 + a * T + b * T**2 + c * T**3))
seawater_resistivity[np.where(seawater_resistivity <= 0)] = 0
into_text_file = np.array([-y_vals, seawater_resistivity]).T

np.savetxt(fname='water_resistivity_temp_profile.txt', X=into_text_file)
