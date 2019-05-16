


from wptherml.wptherml.wpml import multilayer
from wptherml.wptherml.datalib import datalib
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
        'Material_List': ['Air','a-Al2O3','SiO2','SiO2', 'Ag', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 100.0e-9, 1.5e-6, 3.0e-6, 200.0e-9, 0], # You can not have the back reflector as the last layer!!!
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [250e-9, 25000e-9, 2000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

# Initialize the layer slab
np_slab = multilayer(structure)



## Change one of the layers to an effective index
layer = 1
np_slab.layer_alloy(layer,0.15,'Air','a-Al2O3','Bruggeman', plot = False)
#np_slab.layer_alloy(layer,fill_fraction,'Air','Si3N4','MG', plot = False)
layer = 2
np_slab.layer_alloy(layer,0.3,'Air','SiO2','Bruggeman', plot = False)
np_slab.fresnel() # You need to update the fresnel Quantities to reflect the effective index change. 
np_slab.fresnel_ea()
#np_slab.thermal_emission_ea()
#np_slab.cooling_power()

#print("Radiative Power (cooling) is ",np_slab.radiative_power_val, "W/m^2")
#print("Absorbed Solar Power (warming) is ",np_slab.solar_power_val, "W/m^2")
#print("Absorbed Atmospheric Radiation (warming) is ",np_slab.atmospheric_power_val, "W/m^2")
#print("Net Power flux out of the structure is ",np_slab.cooling_power_val, "W/m^2")
#


#%%
#elements = 80
#temp =  np.linspace(219,450,elements)
#rad_pow = np.zeros([elements,elements])
#sol_pow = np.zeros([elements,elements])
#at_pow = np.zeros([elements,elements])
#cool_pow = np.zeros([elements,elements])
#
#for idx0 in range(0,elements):
#    for idx1 in range(0,elements):
#        np_slab.T_ml = temp[idx0]
#        np_slab.T_amb = temp[idx1]
#        
#        #np_slab.thermal_emission()
#        np_slab.thermal_emission_ea()
#        np_slab.cooling_power()
#        BB = datalib.BB(np_slab.lambda_array, np_slab.T_ml)
#        
#        rad_pow[idx0][idx1] = np_slab.radiative_power_val
#        sol_pow[idx0][idx1] = np_slab.solar_power_val
#        at_pow[idx0][idx1] = np_slab.atmospheric_power_val
#        cool_pow[idx0][idx1] = np_slab.cooling_power_val
#        
#
#    
#    
## Calculate standard spectra related to radiative cooling
#AM = datalib.AM(np_slab.lambda_array)
#T_atm = datalib.ATData(np_slab.lambda_array)

#%%

plt.figure()
tempf = (temp-273.15)*(9/5)+32
X,Y = np.meshgrid(tempf,tempf)
#plt.contourf(X, Y, cool_pow, levels = 100, cmap = 'jet')
plt.contourf(X, Y, cool_pow)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Net Outward power flux ($W/cm^2$)')
plt.xlabel('Ambient temperature ($^{\circ}F$)')
plt.ylabel('Structure temperature ($^{\circ}F$)')



#%%

#plt.figure()
#plt.plot(temp,cool_pow, label ='Net outward power flux')
#plt.plot(temp,rad_pow, label = 'Radiative power flux \n (cooling)')
#plt.plot(temp,sol_pow, label = 'Solar power flux \n (warming)')
#plt.plot(temp,sol_pow, label = 'Atmospheric power flux \n (warming)')
#plt.xlabel('Temperature (K)')
#plt.ylabel('Power flux ( $ W/cm^2 $ )')
#plt.tight_layout(rect=[-0.10,0,0.75,1])
#plt.legend(bbox_to_anchor=(1.04, 1))
#plt.show() 
#
#





#%% plot results!
#plt.figure()
#mask = (np_slab.lambda_array >= 3000e-9) & (np_slab.lambda_array <= 30000e-9)
##plt.plot(np_slab.lambda_array[mask]*1e6, T_atm[mask]*100, 'k', alpha = 0.1, label = 'AM1.5 or \n Atmospheric \n transmittance')
#plt.plot(np_slab.lambda_array[mask]*1e6, np_slab.emissivity_array[mask]*100, 'r', label = 'Structure \n absorption')
##plt.plot(np_slab.lambda_array[mask]*1e6, np_slab.thermal_emission_array[mask], 'red')
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('Absorption (%)')
#plt.tight_layout(rect=[-0.10,0,0.75,1])
#plt.legend(bbox_to_anchor=(1.04, 1))
#plt.show() 
#
#plt.figure()
#mask = (np_slab.lambda_array >= 250e-9) & (np_slab.lambda_array <= 3000e-9)
#plt.plot(np_slab.lambda_array[mask]*1e6, np_slab.emissivity_array[mask]*100, 'r', label = 'Structure \n absorption')
##plt.plot(np_slab.lambda_array[mask]*1e6, 100*AM[mask]/(1.4*1e9), 'k', alpha = 0.1, label = 'AM1.5 or \n Atmospheric \n transmittance')
#plt.xlabel('Wavelength (um)')
#plt.ylabel('Absorption (%)')
#plt.tight_layout(rect=[-0.10,0,0.75,1])
#plt.legend(bbox_to_anchor=(1.04, 1))
#plt.show() 
#



#
#print("Radiative Power (cooling) is ",np_slab.radiative_power_val, "W/m^2")
#print("Absorbed Solar Power (warming) is ",np_slab.solar_power_val, "W/m^2")
#print("Absorbed Atmospheric Radiation (warming) is ",np_slab.atmospheric_power_val, "W/m^2")
#print("Net Power flux out of the structure is ",np_slab.cooling_power_val, "W/m^2")
#
































