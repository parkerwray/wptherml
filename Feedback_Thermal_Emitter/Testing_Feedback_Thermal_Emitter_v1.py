

from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
from wptherml.wptherml.wpml import multilayer
import tmm.tmm_core as tmm
from numpy import linspace, inf, pi, stack, array
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mplib
#import random 
#import matplotlib.animation as animation
#from wptherml.wptherml import coolinglib as cool
import time

def ktoc(T):
    #CONVERT K TO C TEMP.
    C = T-273.15
    return C

def ctok(T):
    #CONVERT C TO K TEMP.
    K = T+273.15
    return K

class radiator(multilayer):

    def __init__(self, args):
       
        if 'Thermal_Capacitance' in args:
            self.C = args['Thermal_Capacitance'] # J/(kg*K)
        else:
            self.C = 1000
            

        self.C = 1500  
        self.dt = 1*60 #1 minute
        self.Q = 0
        self.tmm()
    
    
    def step_emissivity(self, lda0 = 0, ldaf = 0, plots = False):
        # MAKE BOXCAR FUNCTION EMISSIVITY
        
        for i in range(0,len(self.t)):
            for j in range(0,len(self.lambda_array)):
                if (self.lambda_array[j] > lda0 and self.lambda_array[j] < ldaf):
                    self.emissivity_array_p[i,j] = self.emissivity_array_s[i,j] =  1
                    self.emissivity_array[j] = 1            
                else:
                    self.emissivity_array_p[i,j] = self.emissivity_array_s[i,j] =  0
                    self.emissivity_array[j] = 0     
           
            self.update()
            
        if (plots):
            self.plot_emissivity
            self.plot_te
            
            
    def get_R(self): 
        # CALCULATE EFFECTIVE RESISTANCE OF RADCOOLER **CURRENTLY INCORRECT**
         self.cooling_power()
         self.R = self.T_ml/self.radiative_power_val
         
         
    def inc_T(self, sun = False):   
        # EULER METHOD TIME DEPENDENT INCREMENT OF THE STRUCTURES TEMPERATURE
        # Tn = Tn-1+(Q-Qrad)dt/C (FROM CIRCUIT DIAGRAM)
        # USER DEFINED INPUTS:
        # 1) IN OR OUR OF THE SUN
        # 2) TEMPERATURE OF THE STRUCTURE AT T-1
        # 3) SPECTRAL EMISSIVITY PROFILE
        # CURRENT HARD WIRED INPUTS
        # 1) INTERNAL HEATING
        # 2) % SOLAR ABSORPTION
        

         Ti = self.T_ml
         Qrad = self.radiative_power_val
         C = self.C
         dt = self.dt
         Q = 50 # Some number to represent internal heating from ship
         
         if (sun == True):
             dl = self.lambda_array[1]-self.lambda_array[0]
             am15 = datalib.AM(self.lambda_array)  #*self.emissivity_array
             am15_integral = am15.sum()*dl
             Q = Q+am15_integral*0.2
             #Q = Q+self.solar_power_val
         
         
         Tf = Ti +(Q-Qrad)*dt/C
         self.Q = Q;
         self.T_ml = Tf
         self.update()
         self.cooling_power()
         
    def CPvT(self, To, Tf):
        # CALCULATE THE CHANGE IN COOLING POWER AS A FUNCTION OF TEMPERATURE 
        # THIS IS BASED ON THE EMISSIVITY PROFILE OF THE STRUCTURE
        
        dT = (Tf-To)/50
        self.TCP = []
        self.T_array = np.arange(To,Tf,dT)
        for T in self.T_array:
            self.T_ml = T
            self.update()
            self.TCP.append(self.radiative_power_val) 
             
    def plot_emissivity(self):
        # PLOT THE SPECTRAL EMISSIVITY AND THE THERMAL EMISSION (BB*EMISSIVITY)
        
        um = 1e-6
        fig, ax1 = plt.subplots()
        ax1.plot(self.lambda_array/um, self.emissivity_array, 'b.-')
        ax1.set_xlabel('Wavelength (um)')
        # Make the y-axis label, ticks and tick labels match the line color.
        ax1.set_ylabel('Emissivity', color='b')
        ax1.tick_params('y', colors='b')
        
        ax2 = ax1.twinx()
        ax2.plot(self.lambda_array/um, self.thermal_emission_array*um)
        #ax2.plot(slab.sim_time, Q, 'k' )
        ax2.set_ylabel('Thermal Emission ($W/m^2/ \mu m$)', color='r')
        ax2.tick_params('y', colors='r')
        
        fig.tight_layout()
        plt.show()   
 
     
         


                    
#%%                    
######################################################################################################                    
                    
# DEFINE SOME UNITS.        
nm = 1e-9
um = 1e-6 
minute = 60
hour = 60*minute

# DEFINE STRUCTURE.
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        
        'Structure_Temperature': 300,
        'Ambient_Temperature': 2,
        'Lambda_List': [2*um, 10*um, 500], #SI unit of m
        'Time_List':[0,2*hour, 20], #SI unit of s   
        
        # Structure depenednt properties
        'Thermal_Capacitance': 1500,  #SI unit of J/(kg*K)
        'Material_List': ['Air', 'Air', 'Air'],
        'Thickness_List': [0, 0, 0], # You can not have the back reflector as the last layer!!!
        
        ## Flags
        'EXPLICIT_ANGLE': 0,
        'COOLING': 1
        }

slab = radiator(structure)
slab.step_emissivity(5*um, 10*um, True)


#%%
T = []
CP = []
Q = []
emission_spectras = []
print(slab.T_ml)
slab.sim_time = np.array(range(0,2*hour, slab.dt))
start_time = time.time()
for i in slab.sim_time:
    CP.append(slab.radiative_power_val)
    T.append(ktoc(slab.T_ml))
    Q.append(slab.Q)
    emission_spectras.append(slab.thermal_emission_array*um)
    
    if (np.floor(i/hour) % 2) == 0:
        slab.inc_T(sun=False)
    else:
        slab.inc_T(sun=True)
        
    #print('Time since loop start: ',time.time()-start_time, ' seconds')
    print('Simulation time: ', i/hour, ' hours')

#%%
    
fig, ax1 = plt.subplots()
ax1.plot(slab.sim_time/hour, T, 'b-')
ax1.set_xlabel('Time (hour)')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('Temperature (C)', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(slab.sim_time/hour, CP, 'r-.')
#ax2.plot(slab.sim_time, Q, 'k' )
ax2.set_ylabel('Radiative Power ($W/m^2$)', color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()
plt.show()   
    
    
#%%    
fig, ax1 = plt.subplots()
ax1.plot(slab.lambda_array/um, slab.emissivity_array, 'b-')
ax1.set_xlabel('Wavelength (um)')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('Emissivity', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(slab.lambda_array/um, np.transpose(emission_spectras))
#ax2.plot(slab.sim_time, Q, 'k' )
ax2.set_ylabel('Thermal Emission ($W/m^2/ \mu m$)', color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()
plt.show()   


#%%
sorted_idx = np.argsort(T)



    
plt.figure()
plt.plot(np.array(T)[sorted_idx], np.array(CP)[sorted_idx])
#
#plt.figure()
#plt.plot(slab.sim_time, CP)

#plt.figure()
#plt.plot(Q)


#%% Get Cooling Power as a Function of Temperature for different Emission Spectras

slab = []
slab.append(radiator(structure))
slab[0].T_ml = ctok(0)
slab[0].T_amb = 2
slab.append(radiator(structure))
slab[1].T_ml = ctok(0)
slab[1].T_amb = 2
#slab[0].tmm()
slab[0].step_emissivity(0.2*um, 20*um, False)
slab[1].step_emissivity(5*um, 6*um, False)

#%%
#plt.figure()
slab[1].plot_emissivity()
#plt.show()



#%%
slab[0].CPvT(ctok(-200), ctok(200))
slab[1].CPvT(ctok(-200), ctok(200))

#%%
plt.figure()
plt.plot(ktoc(slab[0].T_array), slab[0].TCP)
plt.plot(ktoc(slab[1].T_array), slab[1].TCP)






#%%




















#%% 
#plt.figure()
#plt.plot(slab.BBs*um)
#plt.ylabel('$W / m^2 \cdot um$')




#slab.update()
#slab.T_ml_0 = slab.T_ml
#slab.cooling_power

#%%
#plt.figure() 
#plt.plot(slab.lda*1e9, slab.emissivity_array, 'red')
#string = "Thermal Emission at " + str(slab.T_ml) + " K"
#plt.legend(string)
#plt.show()
#
#plt.figure(), 
#plt.plot(slab.lda*1e9, slab.thermal_emission_array, 'red')
#string = "Thermal Emission at " + str(slab.T_ml) + " K"
#plt.legend(string)
#plt.show()
#
#
#
#slab.T_ml = 900
#slab.update()
#
#
#plt.figure(), 
#plt.plot(slab.lda*1e9, slab.thermal_emission_array, 'red')
#string = "Thermal Emission at " + str(slab.T_ml) + " K"
#plt.legend(string)
#plt.show()


#slab.T_ml = ctok(40)
#slab.T_amb = ctok(-100)
#
#
#
#
#
#    slab.Tml = T
#    T_atm = datalib.ATData(slab.lambda_array)
#    BBamb = datalib.BB(slab.lambda_array, slab.T_amb)
#    BBml = datalib.BB(slab.lambda_array, slab.T_ml)
#    e_amb = BBamb*(1-T_atm)
#    
#      
#    
#    slab.thermal_emission_ea()    
#    slab.thermal_emission()        
#    slab.cooling_power()























































