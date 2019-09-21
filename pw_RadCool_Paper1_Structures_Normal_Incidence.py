# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 14:21:30 2019

@author: parkerwray
"""

from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
from wptherml.wptherml.wpml import multilayer
import numpy as np
from numpy import linspace, inf, pi, stack, array, real, imag
import matplotlib.pyplot as plt
import matplotlib as mplib
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
import scipy.io as sio
import multiprocessing as mp
print("Number of processors: ", mp.cpu_count())


import time


#%%

def plot_cooling(slab):
    # Look at cooling properties
    T_atm = datalib.ATData(slab.lambda_array)
    AM = datalib.AM(slab.lambda_array)
#    BBamb = datalib.BB(slab.lambda_array, slab.T_amb)
#    BBml = datalib.BB(slab.lambda_array, slab.T_ml)
    
#    plt.figure()
#    mask = (slab.lambda_array >= 3000e-9) & (slab.lambda_array <= 30000e-9)
#    plt.plot(slab.lambda_array[mask]*1e6, T_atm[mask]*100, 'k', alpha = 0.1, label = 'AM1.5 or \n Atmospheric \n transmittance')
#    plt.plot(slab.lambda_array[mask]*1e6, slab.emissivity_array[mask]*100, 'r', label = 'Structure \n absorption')
#    #plt.plot(np_slab.lambda_array[mask]*1e6, np_slab.thermal_emission_array[mask], 'red')
#    plt.xlabel('Wavelength (nm)')
#    plt.ylabel('Absorption (%)')
#    plt.tight_layout(rect=[-0.10,0,0.75,1])
#    plt.legend(bbox_to_anchor=(1.04, 1))
#    plt.show() 
#    
#    plt.figure()
#    mask = (slab.lambda_array >= 250e-9) & (slab.lambda_array <= 3000e-9)
#    plt.plot(slab.lambda_array[mask]*1e6, slab.emissivity_array[mask]*100, 'r', label = 'Structure \n absorption')
#    plt.plot(slab.lambda_array[mask]*1e6, 100*AM[mask]/(1.4*1e9), 'k', alpha = 0.1, label = 'AM1.5 or \n Atmospheric \n transmittance')
#    plt.xlabel('Wavelength (um)')
#    plt.ylabel('Absorption (%)')
#    plt.tight_layout(rect=[-0.10,0,0.75,1])
#    plt.legend(bbox_to_anchor=(1.04, 1))
#    plt.show() 

    print("Radiative Power (cooling) is ",slab.radiative_power_val, "W/m^2")
    print("Absorbed Solar Power (warming) is ",slab.solar_power_val, "W/m^2")
    print("Absorbed Atmospheric Radiation (warming) is ",slab.atmospheric_power_val, "W/m^2")
    print("Net Power flux out of the structure is ",slab.cooling_power_val, "W/m^2")

#%%

#def main():
start_time = time.time() 
nm = 1e-9
um = 1e-6
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        ### temperature of the structure - relevant for all thermal applications
        ### value is stored in attribute self.T
        'Temperature':280,
        ### actual materials the structure is made from
        ### values are stored in the attribute self.n
        #'Material_List': ['Air','SiO2', 'SiO2','Si3N4','Ag', 'Air'],
        'Material_List': ['Air', 'SiO2', 'Si3N4', 'Ag', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 200*nm, 2500*nm, 200*nm, 0], # You can not have the back reflector as the last layer!!!
        ### range of wavelengths optical properties will be calculated for
        ### values are stored in the array self.lam
        'Lambda_List': [250*nm, 30*um, 2000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Number of degrees you calc in cooling power
        'DEG':7,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }


slab = multilayer(structure)
mat_list = structure['Material_List']
ff_top =0.2
ff_bottom = 0.25
slab.layer_alloy(1,ff_top,'Air',mat_list[1],'Bruggeman', plot = False)
slab.layer_alloy(2,ff_bottom,'Air',mat_list[2],'Bruggeman', plot = False)
slab.tmm()
plot_cooling(slab)

print("--- %s seconds ---" % (time.time() - start_time))
#print("%d" % slab.emissivity_array_p-slab.emissivity_array_s)

sio.savemat(
  'NPNP_SiO2_200nm_FF20_SiN_2500nm_FF25_280K.mat', {
                             'T_s': slab.transmissivity_array_s, 
                             'T_p': slab.transmissivity_array_p,
                             'R_s': slab.reflectivity_array_s,
                             'R_p': slab.reflectivity_array_p,    
                             'e_s': slab.emissivity_array_s, 
                             'e_p': slab.emissivity_array_p, 
                             'e_normal': slab.emissivity_array,                                 
                             'Te_s': slab.thermal_emission_array_s,
                             'Te_p': slab.thermal_emission_array_p,
                             'Te_normal': slab.thermal_emission_array,                            
                             'lda': slab.lambda_array})
    
    
    

#%%
#if __name__ == "__main__":
#    main()



