# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 17:55:30 2019

@author: parkerwray

This code is for comparing the spectra from uv-vis and ftir to the simulated 
spectra using material properties (oscillators) derived from ellipsometery.


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
d_reflector = 2e-6
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
        #'Material_List': ['Air','SiO2', 'SiO2','Si3N4','Ag', 'Air'],
        'Material_List': ['Air','Si3N4','SiO2','SiO2','Si3N4', 'Ag', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 1.0e-6, 1.0e-6, 3.0e-6, 650e-9, 2.0e-6, 0], # You can not have the back reflector as the last layer!!!
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [250e-9, 15000e-9, 5000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

# Initialize the layer slab
np_slab = multilayer(structure)

#np_slab.plot_index(1)




#%%
## Change one of the layers to an effective index

fill_fraction = 0.3
layer = 1
np_slab.layer_alloy(layer,fill_fraction,'Air','Si3N4','Bruggeman', plot = False)
#np_slab.layer_alloy(layer,fill_fraction,'Air','Si3N4','MG', plot = False)
layer = 2
np_slab.layer_alloy(layer,fill_fraction,'Air','SiO2','Bruggeman', plot = False)
#np_slab.layer_alloy(layer,fill_fraction,'Air','SiO2', 'MG', plot = False)
#np_slab.plot_index(1)
np_slab.fresnel() # You need to update the fresnel Quantities to reflect the effective index change. 
np_slab.fresnel_ea()
#np_slab.thermal_emission()
np_slab.thermal_emission_ea()
np_slab.cooling_power()
