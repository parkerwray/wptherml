


from __future__ import division, print_function, absolute_import

import os
os.chdir('C:/Users/parkerwray/Documents/GitHub')

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
#print("Number of processors: ", mp.cpu_count())


import time


# This code makes a two large arrays (1) sio2 on sin and (2) sin on sio2
# The code sweeps temperature
# The code sweeps thickness of top and bottom layer
# The code sweeps fill fraction (including thin film, 100% !)
# In doing this we get the cooling data for all nanoparticle films, nanoparticle
# films + thin film composites, and all thin film films. 
# The maximum achieved cooling power as well as the operating paramaters to 
# achieve this is saved as a function of operating temperature. 
# There are 8 "Opt_P_cool" variables for this. 




#%%
###############################################################################

def plot_cooling(slab):
    # Look at cooling properties
    T_atm = datalib.ATData(slab.lambda_array)
    AM = datalib.AM(slab.lambda_array)
#    BBamb = datalib.BB(slab.lambda_array, slab.T_amb)
#    BBml = datalib.BB(slab.lambda_array, slab.T_ml)
    
    plt.figure()
    mask = (slab.lambda_array >= 3000e-9) & (slab.lambda_array <= 30000e-9)
    plt.plot(slab.lambda_array[mask]*1e6, T_atm[mask]*100, 'k', alpha = 0.1, label = 'AM1.5 or \n Atmospheric \n transmittance')
    plt.plot(slab.lambda_array[mask]*1e6, slab.emissivity_array[mask]*100, 'r', label = 'Structure \n absorption')
    #plt.plot(np_slab.lambda_array[mask]*1e6, np_slab.thermal_emission_array[mask], 'red')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Absorption (%)')
    plt.tight_layout(rect=[-0.10,0,0.75,1])
    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.show() 
    
    plt.figure()
    mask = (slab.lambda_array >= 250e-9) & (slab.lambda_array <= 3000e-9)
    plt.plot(slab.lambda_array[mask]*1e6, slab.emissivity_array[mask]*100, 'r', label = 'Structure \n absorption')
    plt.plot(slab.lambda_array[mask]*1e6, 100*AM[mask]/(1.4*1e9), 'k', alpha = 0.1, label = 'AM1.5 or \n Atmospheric \n transmittance')
    plt.xlabel('Wavelength (um)')
    plt.ylabel('Absorption (%)')
    plt.tight_layout(rect=[-0.10,0,0.75,1])
    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.show() 

    print("Radiative Power (cooling) is ",slab.radiative_power_val, "W/m^2")
    print("Absorbed Solar Power (warming) is ",slab.solar_power_val, "W/m^2")
    print("Absorbed Atmospheric Radiation (warming) is ",slab.atmospheric_power_val, "W/m^2")
    print("Net Power flux out of the structure is ",slab.cooling_power_val, "W/m^2")

def get_cp(slab, mat_list, T, HT, HB, FFT, FFB):   
    slab.T_ml = T
    slab.d[1] = HT
    slab.d[2] = HB
    slab.layer_alloy(1,FFT,'Air',mat_list[1],'MG', plot = False)
    slab.layer_alloy(2,FFB,'Air',mat_list[2],'MG', plot = False)
    #slab.fresnel()
    slab.tmm()

   # plot_cooling(slab)
    return slab.cooling_power_val


def do_loop(HT):
    start_time = time.time() 
    nm = 1e-9
    um = 1e-6
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
            'Material_List': ['Air', 'SiO2', 'Si3N4', 'Ag', 'Air'],
            ### thickness of each layer... terminal layers must be set to zero
            ### values are stored in attribute self.d
            'Thickness_List': [0, 100*nm, 100*nm, 200*nm, 0], # You can not have the back reflector as the last layer!!!
            ### range of wavelengths optical properties will be calculated for
            ### values are stored in the array self.lam
            'Lambda_List': [250*nm, 30*um, 2000],
            ## Calculate for explicit angular dependence
            'EXPLICIT_ANGLE': 1,
            ## Calculate quantities related to radiative cooling
            'COOLING': 1
            }
    
    
    slab = multilayer(structure)
    mat_list = structure['Material_List']
    
#    slab.lambda_array = lda*nm
    
    FF = np.array([0,20,25,30,35,40])/100;    
    T = np.array([300, 290, 280, 270])  
    H_fine = np.linspace(0,75,num = 4) #25nm resolution
    H_med = np.linspace(100,1900,num = 19) #100nm resolution
    H_coarse = np.linspace(2000,3000, num = 5) #250nm resolution
    H = np.concatenate((H_fine,H_med, H_coarse), axis=0)*nm
    
    cp = np.zeros((len(T),len(H),len(FF),len(FF)))
    for idx_T in range(0,len(T)):  # Change temperature
        print("--- %s seconds ---" % (time.time() - start_time))
        print("Temp: %s" % (T[idx_T]))
        for idx_HB in range(0,len(H)):  # Change the top layer thickness
            for idx_FFT in range(0,len(FF)):   # Change the bottom layer thickness
                for idx_FFB in range(0,len(FF)):  # Change the top layer fill fraction   
                    cp[idx_T,idx_HB,idx_FFT, idx_FFB] = get_cp(
                            slab, mat_list, T[idx_T], HT, H[idx_HB], FF[idx_FFT], FF[idx_FFB]) 

    return cp


##############################################################################

# Serial Version
#start_time = time.time()     
#slab = slab_sio2_sin
#mat_list = structure_sio2_sin['Material_List']
#cp = np.zeros((len(T),len(H),len(H),len(FF),len(FF)))
#for idx_T in range(0,len(T)):  # Change temperature
#    print("--- %s seconds ---" % (time.time() - start_time))
#    print("Temp: %s" % (T[idx_T]))
#    for idx_HT in range(0,len(H)):  # Change the top layer thickness
#        for idx_HB in range(0,len(H)):   # Change the bottom layer thickness
#            for idx_FFT in range(0,len(FF)):  # Change the top layer fill fraction
#                for idx_FFB in range(0,len(FF)):  # Change the bottom layer fill fraction
#                    cp[idx_T,idx_HT,idx_HB,idx_FFT,idx_FFB] = get_cp(
#                            slab, mat_list, T[idx_T], H[idx_HT], H[idx_HB], FF[idx_FFT], FF[idx_FFB]) 
#print("--- %s seconds ---" % (time.time() - start_time))
# 
#
#
#sio.savemat('P_SiO2_SiN_Ag.mat', {'sio2_sin_cp_array': cp, 'H':H, 'T':T, 'FF':FF})  


#%%
    
def main(): 
    print("Number of processors: ", mp.cpu_count())
    nm = 1e-9
#    H_fine = np.linspace(120,1000,num = 18) #50nm resolution
#    H_coarse = np.linspace(1000,5000, num = 40) #100nm resolution
#    H = np.concatenate((H_fine, H_coarse), axis=0)*nm
#    test = do_loop(1)

# 130 opps/s
    #FF = np.array([20,25,30,35,40,45,50,55,60,65,100])/100;  
    H_fine = np.linspace(0,75,num = 4) #25nm resolution
    H_med = np.linspace(100,1900,num = 19) #100nm resolution
    H_coarse = np.linspace(2000,3000, num = 5) #250nm resolution
    H = np.concatenate((H_fine,H_med, H_coarse), axis=0)*nm
    #cp2 = np.zeros((len(T),len(H),len(H),len(FF),len(FF)))
    pool = mp.Pool(10)
    start_time = time.time() 
#    cp2 = [pool.map(do_loop,ffo)for ffo in FF] 
    cp2 = pool.map(do_loop,H) 
    print("--- %s seconds ---" % (time.time() - start_time))
    pool.close()  
#
    
    
    FF = np.array([20,25,30,35,100])/100;    
    T = np.array([300, 290, 280, 270, 260])  
    
    sio.savemat('Parallel_Sweep_MG_SiO2_On_Si3N4.mat', {'cp_parallel_array': cp2, 'T':T, 'H':H, 'FF':FF})                  
             
if __name__ == "__main__":
    main()





#operations/sec = 4

## Parallel Version 
#pool = mp.Pool(len(FF))
#start_time = time.time()                     
#slab = slab_sio2_sin
#mat_list = structure_sio2_sin['Material_List']
#for idx_T in range(0,len(T)):  # Change temperature
#    for idx_HT in range(0,len(H)):  # Change the top layer thickness
#        for idx_HB in range(0,len(H)):   # Change the bottom layer thickness
#            for idx_FFT in range(0,len(FF)):  # Change the top layer fill fraction
#                cp2 = [pool.apply(get_cp,
#                            args=(slab, 
#                                  mat_list,
#                                  T[idx_T],
#                                  H[idx_HT],
#                                  H[idx_HB],
#                                  FF[idx_FFT],
#                                  FF[idx_FFB]) )
#                    for idx_FFB in range(0,len(FF))] 
#                    
#print("--- %s seconds ---" % (time.time() - start_time))
#pool.close()                      
#                    
                










