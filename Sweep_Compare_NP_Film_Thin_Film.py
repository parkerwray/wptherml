
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




#%% 
# GET EFFECTIVE INDEX DATA FROM BRUGGEMAN APPROXIMATION 
# GET DATA FOR SIO2 AND SIN OVER DENSE WAVELENGTH RANGE

nm = 1e-9
lda = linspace(200,30000,10000) # list of wavelengths in nm

m = datalib.Material_RI(lda*nm, 'Si3N4') #convert lda to SI unit
fig1 = plt.figure()
plt.plot(lda, real(m),'k')
plt.plot(lda, imag(m),'r')

#%%
ff = np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])/100;

m_sio2 = np.zeros((len(lda),len(ff)+1), dtype = np.complex64);
for idx in range(0,len(ff)):
    m_sio2[:,idx] = datalib.alloy(lda*nm, ff[idx], 'Air','SiO2','Bruggeman')
m_sio2[:,-1] = lda

m_sin = np.zeros((len(lda),len(ff)+1), dtype = np.complex64);
for idx in range(0,len(ff)):
    m_sin[:,idx] = datalib.alloy(lda*nm, ff[idx], 'Air','SiN','Bruggeman')
m_sin[:,-1] = lda    
    
sio.savemat('SiO2_Brugg_FF_0_5_100_lda.mat', {'m_sio2': m_sio2})
sio.savemat('SiN_Brugg_FF_0_5_100_lda.mat', {'m_sin': m_sin})
   
    
#%%

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
        'Thickness_List': [0, 100*nm, 100*nm, 1*um, 0], # You can not have the back reflector as the last layer!!!
        ### range of wavelengths optical properties will be calculated for
        ### values are stored in the array self.lam
        'Lambda_List': [250*nm, 30*um, 5000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

structure_sin_sio2 = structure_sio2_sin
structure_sin_sio2['Material_List']=['Air', 'Si3N4', 'SiO2', 'Ag', 'Air']

slab_sio2_sin_np = multilayer(structure_sio2_sin)
slab_sin_sio2_np = multilayer(structure_sin_sio2)

H = np.linspace(100, 5000, num = 50)*nm
T = np.array([300, 290, 280, 270, 260, 250])
FF = np.array([5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])/100;


#H = np.linspace(500, 2000, num = 2)*nm
#T = np.array([300, 280])
#FF = np.array([30,100])/100;

P_cool_sio2_sin_np = np.zeros((len(T),len(H),len(H),len(FF),len(FF)))
P_cool_sin_sio2_np = np.zeros((len(T),len(H),len(H),len(FF),len(FF)))

Opt_P_cool_sio2_sin = np.zeros((len(T),5))+(-1000)
Opt_P_cool_sio2_sin_np = np.zeros((len(T),5))+(-1000)
Opt_P_cool_sio2_np_sin = np.zeros((len(T),5))+(-1000)
Opt_P_cool_sio2_np_sin_np = np.zeros((len(T),5))+(-1000)

Opt_P_cool_sin_sio2 = np.zeros((len(T),5))+(-1000)
Opt_P_cool_sin_sio2_np = np.zeros((len(T),5))+(-1000)
Opt_P_cool_sin_np_sio2 = np.zeros((len(T),5))+(-1000)
Opt_P_cool_sin_np_sio2_np = np.zeros((len(T),5))+(-1000)

#%%

for idx_T in range(0,len(T)):  # Change the operating temp of structure
    slab_sin_sio2_np.T_ml = T[idx_T]
    slab_sio2_sin_np.T_ml = T[idx_T]
    
    for idx_HT in range(0,len(H)):  # Change the top layer thickness
        slab_sin_sio2_np.d[1] = H[idx_HT]
        slab_sio2_sin_np.d[1] = H[idx_HT]
        
        for idx_HB in range(0,len(H)):   # Change the bottom layer thickness
            slab_sin_sio2_np.d[2] = H[idx_HT]
            slab_sio2_sin_np.d[2] = H[idx_HT]
                    
            for idx_ff_T in range(0,len(FF)):  # Change the top layer fill fraction
                slab_sin_sio2_np.layer_alloy(1,FF[idx_ff_T],'Air','Si3N4','Bruggeman', plot = False)
                slab_sin_sio2_np.fresnel()
                slab_sio2_sin_np.layer_alloy(1,FF[idx_ff_T],'Air','SiO2','Bruggeman', plot = False)          
                slab_sin_sio2_np.fresnel()

                for idx_ff_B in range(0,len(FF)):  # Change the bottom layer fill fraction
                    slab_sin_sio2_np.layer_alloy(2,FF[idx_ff_B],'Air','SiO2','Bruggeman', plot = False)
                    slab_sin_sio2_np.fresnel()
                    slab_sio2_sin_np.layer_alloy(2,FF[idx_ff_B],'Air','Si3N4','Bruggeman', plot = False)    
                    slab_sin_sio2_np.fresnel()
                    
                    slab_sin_sio2_np.tmm()
                    slab_sio2_sin_np.tmm()
                    
                    P_cool_sin_sio2_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B] = slab_sin_sio2_np.cooling_power_val                    
                    P_cool_sio2_sin_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B] = slab_sio2_sin_np.cooling_power_val
                    
                    # Update data for top and bottom both thin film cases
                    if (idx_ff_B == len(FF)-1 and idx_ff_T == len(FF)-1):
                            if P_cool_sin_sio2_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B] > Opt_P_cool_sin_sio2[idx_T][0]:
                                Opt_P_cool_sin_sio2[idx_T,0] = P_cool_sin_sio2_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B]
                                Opt_P_cool_sin_sio2[idx_T,1:] = [H[idx_HT],H[idx_HB],FF[idx_ff_T],FF[idx_ff_B]]
                                
                            if P_cool_sio2_sin_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B] > Opt_P_cool_sio2_sin[idx_T][0]:
                                Opt_P_cool_sio2_sin[idx_T,0] = P_cool_sio2_sin_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B]
                                Opt_P_cool_sio2_sin[idx_T,1:] = [H[idx_HT],H[idx_HB],FF[idx_ff_T],FF[idx_ff_B]]  

                    # Update data for top is NP and bottom is thin film cases
                    if (idx_ff_B == len(FF)-1 and idx_ff_T != len(FF)-1):
                            if P_cool_sin_sio2_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B] > Opt_P_cool_sin_np_sio2[idx_T][0]:
                                Opt_P_cool_sin_np_sio2[idx_T,0] = P_cool_sin_sio2_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B]
                                Opt_P_cool_sin_np_sio2[idx_T,1:] = [H[idx_HT],H[idx_HB],FF[idx_ff_T],FF[idx_ff_B]]
                                
                            if P_cool_sio2_sin_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B] > Opt_P_cool_sio2_np_sin[idx_T][0]:
                                Opt_P_cool_sio2_np_sin[idx_T,0] = P_cool_sio2_sin_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B]
                                Opt_P_cool_sio2_np_sin[idx_T,1:] = [H[idx_HT],H[idx_HB],FF[idx_ff_T],FF[idx_ff_B]] 
                
                    # Update data for bottom is NP and top is thin film cases
                    if (idx_ff_B != len(FF)-1 and idx_ff_T == len(FF)-1):
                            if P_cool_sin_sio2_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B] > Opt_P_cool_sin_sio2_np[idx_T][0]:
                                Opt_P_cool_sin_sio2_np[idx_T,0] = P_cool_sin_sio2_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B]
                                Opt_P_cool_sin_sio2_np[idx_T,1:] = [H[idx_HT],H[idx_HB],FF[idx_ff_T],FF[idx_ff_B]]
                                
                            if P_cool_sio2_sin_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B] > Opt_P_cool_sio2_sin_np[idx_T][0]:
                                Opt_P_cool_sio2_sin_np[idx_T,0] = P_cool_sio2_sin_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B]
                                Opt_P_cool_sio2_sin_np[idx_T,1:] = [H[idx_HT],H[idx_HB],FF[idx_ff_T],FF[idx_ff_B]] 

                    # Update data for top and bottom both being NP films case
                    if (idx_ff_B != len(FF)-1 and idx_ff_T != len(FF)-1):
                            if P_cool_sin_sio2_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B] > Opt_P_cool_sin_np_sio2_np[idx_T][0]:
                                Opt_P_cool_sin_np_sio2_np[idx_T,0] = P_cool_sin_sio2_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B]
                                Opt_P_cool_sin_np_sio2_np[idx_T,1:] = [H[idx_HT],H[idx_HB],FF[idx_ff_T],FF[idx_ff_B]]
                                
                            if P_cool_sio2_sin_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B] > Opt_P_cool_sio2_np_sin_np[idx_T][0]:
                                Opt_P_cool_sio2_np_sin_np[idx_T,0] = P_cool_sio2_sin_np[idx_T][idx_HT][idx_HB][idx_ff_T][idx_ff_B]
                                Opt_P_cool_sio2_np_sin_np[idx_T,1:] = [H[idx_HT],H[idx_HB],FF[idx_ff_T],FF[idx_ff_B]]


#%%
                                
sio.savemat('P_SiO2_SiN_Ag.mat', {'P_cool_sio2_sin_np': P_cool_sio2_sin_np})                                
sio.savemat('Optimal_SiO2_SiN_Ag.mat', {'Opt_P_cool_sio2_sin': Opt_P_cool_sio2_sin})
sio.savemat('Optimal_SiO2_NP_SiN_Ag.mat', {'Opt_P_cool_sio2_np_sin': Opt_P_cool_sio2_np_sin})
sio.savemat('Optimal_SiO2_SiN_NP_Ag.mat', {'Opt_P_cool_sio2_sin_np': Opt_P_cool_sio2_sin_np})
sio.savemat('Optimal_SiO2_NP_SiN_NP_Ag.mat', {'Opt_P_cool_sio2_np_sin_np': Opt_P_cool_sio2_np_sin_np})

sio.savemat('P_SiN_SiO2_Ag.mat', {'P_cool_sin_sio2_np': P_cool_sin_sio2_np})                                
sio.savemat('Optimal_SiN_SiO2_Ag.mat', {'Opt_P_cool_sin_sio2': Opt_P_cool_sin_sio2})
sio.savemat('Optimal_SiN_NP_SiO2_Ag.mat', {'Opt_P_cool_sin_np_sio2': Opt_P_cool_sin_np_sio2})
sio.savemat('Optimal_SiN_SiO2_NP_Ag.mat', {'Opt_P_cool_sin_sio2_np': Opt_P_cool_sin_sio2_np})
sio.savemat('Optimal_SiN_NP_SiO2_NP_Ag.mat', {'Opt_P_cool_sin_np_sio2_np': Opt_P_cool_sin_np_sio2_np})













