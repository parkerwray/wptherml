# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 13:43:32 2019

This code creates a multilayer to look at
1) Ellipsometric variables
2) TRA spectra
3) Refractive index

The purpose of this code is to verify that that the model for the material as 
well as the TMM and effective index code can reproduce measurement results.

@author: parkerwray
"""

um = 1e-6
nm = 1e-9

from wptherml.wptherml.wpml_v3 import multilayer
from wptherml.wptherml.datalib import datalib
from matplotlib import pyplot as plt
import numpy as np


 
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        ### temperature of the structure - relevant for all thermal applications
        ### value is stored in attribute self.T
        'Temperature': 300,
        ### actual materials the structure is made from
        ### values are stored in the attribute self.n
        'Material_List': ['Air', 'RC0_1B_SiO2', 'Si3N4', 'Ag', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 1.7*um, 900*nm, 200*nm, 0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [200*nm, 40*um, 10000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

# Initialize the layer slab
np_slab = multilayer(structure)

#np_slab.plot_index(2)




#%
# Change one of the layers to an effective index
fill_fraction = 0.35
layer = 1
np_slab.layer_alloy(layer,fill_fraction,'Air','Si3N4','Bruggeman', plot = False)
#np_slab.layer_alloy(layer,fill_fraction,'Air','Si3N4','MG', plot = False)
#layer = 2
#np_slab.layer_alloy(layer,fill_fraction,'Air','SiO2','Bruggeman', plot = True)
#np_slab.layer_alloy(layer,fill_fraction,'Air','SiO2','MG', plot = False)
np_slab.tmm()
np_slab.update()
np_slab.thermal_emission_ea()
np_slab.cooling_power(radiative = True, atmospheric = True, solar = True, total = True)

print("Radiative Power (cooling) is ",np_slab.radiative_power_val, "W/m^2")
print("Absorbed Solar Power (warming) is ",np_slab.solar_power_val, "W/m^2")
print("Absorbed Atmospheric Radiation (warming) is ",np_slab.atmospheric_power_val, "W/m^2")
print("Net Power flux out of the structure is ",np_slab.cooling_power_val, "W/m^2")


#%%

CP=[]
thickness = np.arange(1.7-1,1.8+2, 0.05)*um
for thic in thickness:
    np_slab.d[1] = thic
    print("NP Thickness: ", np_slab.d[1])
    np_slab.tmm()
    np_slab.update()
    np_slab.cooling_power(radiative = True, atmospheric = True, solar = True, total = True)
    print("Net Power Flux: ", np_slab.cooling_power_val)
    CP.append(np_slab.cooling_power_val)
    
plt.figure()
plt.plot(thickness/um, CP)  
plt.ylabel('Net power flux (W/m^2)')
plt.xlabel('SiO2 nanoparticle film thickness (um)')
plt.title('30% F.F. NP SiO2 Film on 900nm Si3N4 on 200nm Ag')    

#%%
    
plt.figure()
plt.plot(thickness/um, CP)  
plt.ylabel('Net power flux (W/m^2)')
plt.xlabel('SiO2 nanoparticle film thickness (um)')
plt.title('35% F.F. NP SiO2 Film on 900nm Si3N4 on 200nm Ag')    
   


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