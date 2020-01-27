
from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
import numpy as np
import matplotlib.pyplot as plt
import wptherml.Feedback_Thermal_Emitter.pw_SpaceRadCool_v1 as rcs

import itertools as it
import matplotlib as mplib
#import random 
#import matplotlib.animation as animation
#from wptherml.wptherml import coolinglib as cool
import time

# DEFINE PLOT PARAMETERS
mplib.rcParams['lines.linewidth'] = 8
mplib.rcParams['lines.markersize'] = 6
mplib.rcParams['axes.titlesize'] = 30
mplib.rcParams['axes.labelsize'] = 26
mplib.rcParams['xtick.labelsize'] = 24
mplib.rcParams['ytick.labelsize'] = 24
mplib.rcParams['font.size'] = 26

# DEFINE SOME UNITS.        
nm = 1e-9
um = 1e-6 
minute = 60
hour = 60*minute

def ktoc(T):
    #CONVERT K TO C TEMP.
    C = T-273.15
    return C

def ctok(T):
    #CONVERT C TO K TEMP.
    K = T+273.15
    return K


# DEFINE SWEEP PARAMETERS
e_values = np.linspace(0, 1, 2)
lda = np.linspace(1, 40, 10)


# DEFINE STRUCTURE.
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        
        'Structure_Temperature': rcs.ctok(40),
        'Ambient_Temperature': 2,
        'Lambda_List': [np.amin(lda)*um, np.amax(lda)*um, len(lda)], #SI unit of m
        
        # Structure depenednt properties
        'Internal_Load_Flux': 40,
        
        'Material_List': ['Air', 'Air', 'Air'],
        'Thickness_List': [0, 0, 0], # You can not have the back reflector as the last layer!!!
        
        ## Flags
        'EXPLICIT_ANGLE': 1,
        'COOLING': 1
        }

e_list = list(it.product(e_values, repeat = len(lda))) # size is e_values ^ lda

#%%
slab = rcs.radiator(structure)

CPS = []
Es = []
start = time.time()
elem = range(len(e_list))
for i in elem:
    slab.update_emissivity(e_list[i])
    slab.CPvT()
#    CPS.append(slab.Q_rad_array)
#    Es.append(slab.total_emissivity_array)
    #print('loop time: '+str(time.time()-start)+'sec')
CPS = np.reshape(np.array(slab.Q_rad_array),(len(elem),len(slab.temp_array)))
Es = np.reshape(np.array(slab.total_emissivity_array),(len(elem),len(slab.temp_array)))

#%%
CPSmin = np.amin(CPS, axis = 1)
CPSmax = np.amax(CPS, axis = 1)
CPSdelta = CPSmax-CPSmin;
best_spectra = np.argmax(CPSdelta)



slab2 = rcs.radiator(structure)
slab2.step_emissivity(3*um, 18*um, plot = False)
slab2.CPvT()
CPSdelta2 = (np.max(slab2.Q_rad_array)-np.min(slab2.Q_rad_array))
#%% Get Space
t = slab.temp_array
e = np.arange(0.1,1.001,0.01)
E, T = np.meshgrid(e, t)
sigma =  5.670373e-8
Q = sigma*E*(T**4)


#%%

plt.figure()
plt.contourf(ktoc(T),
             Q,
             E, 
             100) 
plt.colorbar(ticks = (np.linspace(np.amin(ktoc(T)),np.amax(ktoc(T)),10)))   
plt.plot(ktoc(t), CPS[best_spectra])
plt.plot(ktoc(t), slab2.Q_rad_array)

#for i in range(len(elem)):
#    plt.plot(ktoc(t),CPS[i])    
   
#%%

#plt.figure()
#plt.plot(ktoc(slab.temp_array), CPS[0])
#plt.plot(ktoc(slab.temp_array), CPS[1])   









    