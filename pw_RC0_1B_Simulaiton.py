# -*- coding: utf-8 -*-
"""

This code creates a multilayer to look at
1) Ellipsometric variables
2) TRA spectra
3) Refractive index

The purpose of this code is to verify that that the model for the material as 
well as the TMM and effective index code can reproduce measurement results.

@author: parkerwray
"""



from wptherml.wpml import multilayer
from wptherml.datalib import datalib
from matplotlib import pyplot as plt
import numpy as np

d_np_film1 = 1.0e-6
d_np_film2 = 1.0e-6
d_film1 = 3e-6
d_film2 = 650e-9
d_film3 = 0 #200e-9
d_reflector = 200e-9
d_substrate = 200e-6 # Note, very large substrates cause an overflow effect.
 
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        ### temperature of the structure - relevant for all thermal applications
        ### value is stored in attribute self.T
        'Temperature': 300,
        ### actual materials the structure is made from
        ### values are stored in the attribute self.n
        'Material_List': ['Air', 'SiO2', 'Si', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 1e-6,1e-6,  0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [250e-9, 14000e-9, 1000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

# Initialize the layer slab
np_slab = multilayer(structure)

#np_slab.plot_index(2)




#%%
# Change one of the layers to an effective index
fill_fraction = 0.3
layer = 1
np_slab.layer_alloy(layer,fill_fraction,'Air','SiO2','Bruggeman', plot = True)
#np_slab.layer_alloy(layer,fill_fraction,'Air','SiO2','MG', plot = False)
np_slab.thermal_emission_ea()
np_slab.cooling_power()

#%%
# Calculate standard spectra related to radiative cooling
AM = datalib.AM(np_slab.lambda_array)
T_atm = datalib.ATData(np_slab.lambda_array)
BB = datalib.BB(np_slab.lambda_array, np_slab.T_ml)

### plot results!
plt.figure()
mask = (np_slab.lambda_array >= 3000e-9) & (np_slab.lambda_array <= 30000e-9)
plt.plot(np_slab.lambda_array[mask]*1e6, T_atm[mask], 'cyan')
plt.plot(np_slab.lambda_array[mask]*1e6, np_slab.emissivity_array[mask], 'red')
#plt.plot(np_slab.lambda_array[mask]*1e6, np_slab.thermal_emission_array[mask], 'red')

plt.figure()
mask = (np_slab.lambda_array >= 250e-9) & (np_slab.lambda_array <= 3000e-9)
plt.plot(np_slab.lambda_array[mask]*1e6, np_slab.emissivity_array[mask], 'blue')
plt.plot(np_slab.lambda_array[mask]*1e6, AM[mask]/(1.4*1e9), 'red')


print("Radiative Power (cooling) is ",np_slab.radiative_power_val, "W/m^2")
print("Absorbed Solar Power (warming) is ",np_slab.solar_power_val, "W/m^2")
print("Absorbed Atmospheric Radiation (warming) is ",np_slab.atmospheric_power_val, "W/m^2")
print("Net Power flux out of the structure is ",np_slab.cooling_power_val, "W/m^2")


#np_slab.layer_alloy(layer,fill_fraction,'Air','Au','Maxwell-Garnett', plot = True)
#np_slab.fresnel()
##plt.plot(np_slab.lambda_array*1e+9, np_slab.transmissivity_array )



# https://www2.pvlighthouse.com.au/resources/photovoltaic%20materials/refractive%20index/refractive%20index.aspx