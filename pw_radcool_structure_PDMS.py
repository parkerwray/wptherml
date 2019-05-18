# -*- coding: utf-8 -*-
"""
Sweep thickness parameters to determine the most efficient radiative cooler
"""



from wptherml.wptherml.wpml import multilayer
from wptherml.wptherml.datalib import datalib
from matplotlib import pyplot as plt
import numpy as np

nm = 1e-9
um = 1e-6 
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        ### temperature of the structure - relevant for all thermal applications
        ### value is stored in attribute self.T
        'Temperature': 275,
        ### actual materials the structure is made from
        ### values are stored in the attribute self.n
        #'Material_List': ['Air','SiO2', 'SiO2','Si3N4','Ag', 'Air'],
        'Material_List': ['Air','PDMS','PDMS','Air', 'Ag', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 500*nm, 600*nm, 900*nm, 200*nm, 0], # You can not have the back reflector as the last layer!!!
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [350*nm, 30*um, 1000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

# Initialize the layer slab
np_slab = multilayer(structure)

# Change one of the layers to an effective index
fill_fraction = 0.3
layer = 1
np_slab.layer_alloy(layer,fill_fraction,'PDMS','RC0_1B_SiO2','Bruggeman', plot = False)
layer = 2
np_slab.layer_alloy(layer,fill_fraction,'PDMS','Si3N4','Bruggeman', plot = False)
layer = 3
np_slab.layer_static_ri(layer, 2+1j*0)


#%%
# Run simulation 
np_slab.tmm()
Tss = np_slab.steady_state_temperature()

#%%
#### Plot the structures spectra
#AM = datalib.AM(np_slab.lambda_array)
#T_atm = datalib.ATData(np_slab.lambda_array)
#
#plt.figure()
#mask = (np_slab.lambda_array >= 3000e-9) & (np_slab.lambda_array <= 30000e-9)
#plt.plot(np_slab.lambda_array[mask]*1e6, T_atm[mask]*100, 'k', alpha = 0.1, label = 'AM1.5 or \n Atmospheric \n transmittance')
#plt.plot(np_slab.lambda_array[mask]*1e6, np_slab.emissivity_array[mask]*100, 'r', label = 'Structure \n absorption')
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('Absorption (%)')
#plt.tight_layout(rect=[-0.10,0,0.75,1])
#plt.legend(bbox_to_anchor=(1.04, 1))
#plt.show() 
#
#plt.figure()
#mask = (np_slab.lambda_array >= 250e-9) & (np_slab.lambda_array <= 3000e-9)
#plt.plot(np_slab.lambda_array[mask]*1e6, np_slab.emissivity_array[mask]*100, 'r', label = 'Structure \n absorption')
#plt.plot(np_slab.lambda_array[mask]*1e6, 100*AM[mask]/(1.4*1e9), 'k', alpha = 0.1, label = 'AM1.5 or \n Atmospheric \n transmittance')
#plt.xlabel('Wavelength (um)')
#plt.ylabel('Absorption (%)')
#plt.tight_layout(rect=[-0.10,0,0.75,1])
#plt.legend(bbox_to_anchor=(1.04, 1))
#plt.show() 

#%%
### Plot the emission properties
T_atm = datalib.ATData(np_slab.lambda_array)
BBamb = datalib.BB(np_slab.lambda_array, np_slab.T_amb)
BBml = datalib.BB(np_slab.lambda_array, np_slab.T_ml)

plt.figure()
mask = (np_slab.lambda_array >= 3000e-9) & (np_slab.lambda_array <= 30000e-9)
plt.plot(np_slab.lambda_array[mask]*1e6, BBamb[mask]*(1-T_atm[mask]))
plt.plot(np_slab.lambda_array[mask]*1e6, BBml[mask]*np_slab.emissivity_array[mask])
plt.show()

print("Radiative Power (cooling) is ",np_slab.radiative_power_val, "W/m^2")
print("Absorbed Solar Power (warming) is ",np_slab.solar_power_val, "W/m^2")
print("Absorbed Atmospheric Radiation (warming) is ",np_slab.atmospheric_power_val, "W/m^2")
print("Net Power flux out of the structure is ",np_slab.cooling_power_val, "W/m^2")




