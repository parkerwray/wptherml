
from __future__ import division, print_function, absolute_import

import os
os.chdir('C:/Users/parkerwray/Documents/GitHub')

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.datalib import datalib
from wptherml.wpml import multilayer
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


nm = 1e-9
um = 1e-6
structure_sio2_sin = {
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
        'Lambda_List': [250*nm, 30*um, 2],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

structure_sin_sio2 = structure_sio2_sin
structure_sin_sio2['Material_List']=['Air', 'Si3N4', 'SiO2', 'Ag', 'Air']

slab_sio2_sin = multilayer(structure_sio2_sin)
slab_sin_sio2 = multilayer(structure_sin_sio2)

lda_uv = np.linspace(251,370, num = 120) # 1nm resolution
lda_vis = np.linspace(373,1000,num = 210) # 3nm resolution
lda_ir = np.linspace(1100,30000,num = 290) # 100nm resolution
lda = np.concatenate((lda_uv,lda_vis,lda_ir), axis=0)
slab_sio2_sin.lam = lda*nm
slab_sin_sio2.lam = lda*nm


H_fine = np.linspace(120,1000,num = 18) #50nm resolution
H_coarse = np.linspace(1000,5000, num = 40) #100nm resolution
H = np.concatenate((H_fine, H_coarse), axis=0)*nm

T = np.array([300, 290]) #, 280, 270, 260, 250])  
H = np.linspace(100, 5000, num = 2)*nm
FF = np.array([5,10])/100 # np.array([5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])/100;

#%%
###############################################################################
def get_cp(slab, mat_list, T, HT, HB, FFT, FFB):   
    slab.T_ml = T
    slab.d[1] = HT
    slab.d[2] = HB
    slab.layer_alloy(1,FFT,'Air',mat_list[1],'Bruggeman', plot = False)
    slab.layer_alloy(2,FFB,'Air',mat_list[2],'Bruggeman', plot = False)
    slab.fresnel()
    slab.tmm()
    return slab.cooling_power_val

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

pool = mp.Pool(10)
start_time = time.time()  
cp2 = [pool.apply(get_cp,
            args=(slab, 
                  mat_list,
                  T[idx_T],
                  H[idx_HT],
                  H[idx_HB],
                  FF[idx_FFT],
                  FF[idx_FFB]))
    for idx_FFB in range(0,len(FF))] 
print("--- %s seconds ---" % (time.time() - start_time))
pool.close()                      
             


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
                    
                    































