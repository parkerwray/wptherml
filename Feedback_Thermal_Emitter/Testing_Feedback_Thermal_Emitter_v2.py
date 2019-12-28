

from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
from wptherml.wptherml.wpml_v2 import multilayer
import tmm.tmm_core as tmm
from numpy import linspace, inf, pi, stack, array
import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
import copy
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
            
        self.solar_abs = []
        self.T_ideal = ctok(40)
        self.BB_array = []
        self.BB_total_array = []
        self.temp_array = []
        self.Q_rad_array = []
        self.total_emissivity_array = []
        
        multilayer.__init__(self,args)
        self.tmm()
        self.update()
    
    def step_emissivity(self, lda0 = 0, ldaf = 0, plot = False):
        # MAKE BOXCAR FUNCTION EMISSIVITY
        
        for i in range(0,len(self.t)):
            for j in range(0,len(self.lambda_array)):
                if (self.lambda_array[j] > lda0 and self.lambda_array[j] < ldaf):
                    self.emissivity_array_p[i,j] = self.emissivity_array_s[i,j] =  1
                    self.emissivity_array[j] = 1            
                else:
                    self.emissivity_array_p[i,j] = self.emissivity_array_s[i,j] =  0
                    self.emissivity_array[j] = 0     
            
            self.cooling_calc = False
            self.update()
            self.cooling_calc = True
            
        if (plot):
            self.plot_emissivity()

    def get_e_total(self):
        # GET THE TOTAL EMISSIVITY
        self.integrate_BB()
        self.e_total = self.radiative_power_val/self.BB_total
        
    def integrate_BB(self):
        # INTEGRATE THE BB TO GET W/M^2
        dlda=(self.lambda_array[1]-self.lambda_array[0])
        self.BB_total = np.sum(np.array(self.BBs)*dlda)*np.pi
        
    def get_R(self): 
        # CALCULATE EFFECTIVE RESISTANCE OF RADCOOLER **CURRENTLY INCORRECT**
         self.cooling_power()
         self.R = self.T_ml/self.radiative_power_val
              
    def CPvT(self, plot = False, temp_difference = True):
        # CALCULATE THE CHANGE IN COOLING POWER AS A FUNCTION OF TEMPERATURE 
        # THIS IS BASED ON THE EMISSIVITY PROFILE OF THE STRUCTURE AND IS NOT 
        # TIME DEPENDENT
        Tf = ctok(200)
        To = ctok(-200)
        dT = 10
        self.temp_array = np.arange(To,Tf,dT)
        for T in self.temp_array:
            self.T_ml = T
            self.update()
            self.Q_rad_array.append(self.radiative_power_val) 
            self.BB_array.append(self.BBs)
            
            self.get_e_total()
            self.BB_total_array.append(self.BB_total)   
            self.total_emissivity_array.append(self.e_total)            

        
        self.Boltzmann_law()
        
        if (plot):
            self.plot_CPvT(temp_difference)

            
    def Boltzmann_law(self):
        sigma = 5.670373e-8 #W/m^2K^4
        self.Q_Boltzmann_array = sigma*np.array(self.temp_array)**4   #W/m^2
    #np.array(self.total_emissivity_array)*
    
    
    def plot_CPvT(self, temp_difference = True):
        
        if temp_difference == True:
            x = (self.temp_array-self.T_ideal)
        else:
            x = ktoc(self.temp_array)
        
        fig, ax1 = plt.subplots()
        ax1.plot(x, self.total_emissivity_array, 'b')
        ax1.set_xlabel('$\Delta$ Temperature (C)')
        # Make the y-axis label, ticks and tick labels match the line color.
        ax1.set_ylabel('Total emissivity', color='b')
        ax1.tick_params('y', colors='b')
        
        ax2 = ax1.twinx()
        ax2.plot(x, np.array(self.Q_rad_array),'r', label = 'Strucutre Emission')
        ax2.plot(x, np.array(self.BB_total_array),'g-.', label = 'Blackbody Emission')
        ax2.plot(x, np.array(self.Q_Boltzmann_array),'k.', label = 'Boltzmann Emission')
        #ax2.plot(slab.sim_time, Q, 'k' )
        ax2.set_ylabel('Radiative emission ($W/m^2$)', color='r')
        ax2.tick_params('y', colors='r')
        ax2.legend()
        fig.tight_layout()
        plt.show()          
        
             
    def plot_emissivity(self):
        # PLOT THE SPECTRAL EMISSIVITY AND THE THERMAL EMISSION (BB*EMISSIVITY)
        
        um = 1e-6
        fig, ax1 = plt.subplots()
        ax1.plot(self.lambda_array/um, self.emissivity_array, 'b:')
        ax1.set_xlabel('Wavelength (um)')
        # Make the y-axis label, ticks and tick labels match the line color.
        ax1.set_ylabel('Emissivity', color='b')
        ax1.tick_params('y', colors='b')
        
        ax2 = ax1.twinx()
        ax2.plot(self.lambda_array/um, self.thermal_emission_array*um, 'r')
        #ax2.plot(slab.sim_time, Q, 'k' )
        ax2.set_ylabel('Thermal Emission ($W/m^2/ \mu m$)', color='r')
        ax2.tick_params('y', colors='r')
        
        fig.tight_layout()
        plt.show()   
 
     
class time_dependent_radiator(radiator):
    
    def __init__(self, args):
        
        if 'Time_List' in args:
            timelist = args['Time_List']
            self.time_array = np.linspace(timelist[0],timelist[1],int(timelist[2]), endpoint = False)
        else:
            self.time_array = np.linspace(0,10,int(10), endpoint = False)
            
        if 'Sun_Period' in args:
            period = 2*args['Sun_Period']
            #time_range = (timelist[1]-timelist[0])
            t_norm = 2*np.pi*self.time_array/period
            self.sun_mask = 0.5+0.5*sci.signal.square(t_norm, duty = 0.5)
        else:
            self.sun_mask = 0*self.time_array
            
        if 'Internal_Load_Flux' in args:
            self.Qload = args['Internal_Load_Flux'] # W/m^2
        else:
            self.Qload = 10   

        radiator.__init__(self,args)
        self.define_Qinput()    
                
        
    def define_Qinput(self):
        # DEFINE THE PULSE SHAPE FOR THE INPUT EXCITATION. THIS IS ASSUMED TO BE 
        # A CONSTANT "CURRENT = POWER" SOURCE. I.E., dT/dt. 
        # THE PURPOSE OF THIS FUNCTION IS TO EXPAND IT TO ENABLE ARBITRATY INPUT
        # WAVEFORMS THAT CAN REPRESENT ORBITAL DYNAMICS.
        
        dl = self.lambda_array[1]-self.lambda_array[0]
        am15 = datalib.AM(self.lambda_array)  #*self.emissivity_array
        am15_integral = am15.sum()*dl
        self.Qsolar = am15_integral*0.2
        
        self.Qinput_array = self.Qload + self.Qsolar*self.sun_mask
            
    
    
    def inc_T(self, Qinput):   
        # EULER METHOD TIME DEPENDENT INCREMENT OF THE STRUCTURES TEMPERATURE
        # Tn = Tn-1+(Q-Qrad)dt/C (FROM CIRCUIT DIAGRAM)
        # USER DEFINED INPUTS:
        # 1) IN OR OUR OF THE SUN
        # 2) TEMPERATURE OF THE STRUCTURE AT T-1
        # 3) SPECTRAL EMISSIVITY PROFILE
        # CURRENT HARD WIRED INPUTS
        # 1) INTERNAL HEATING
        # 2) % SOLAR ABSORPTION
        
         dt = self.time_array[1]-self.time_array[0]                  
         Tf = self.T_ml +(Qinput-self.radiative_power_val)*dt/self.C
         
         self.T_ml = Tf
         self.update()
         self.cooling_power()

         
    def time_dependent_temp(self):  
        # RUN TIME DEPENDENT HEATING EQUATION FOR A GIVEN 
        # 1) TIME STEP AND TIME WINDOW (GIVEN BY TIME_ARRAY)
        # 2) SUN/NO SUN INTERVALS (GIVEN BY SUN MASK)
        
        #emission_spectras = []
        for i in range(len(self.time_array)):
            self.Q_rad_array.append(self.radiative_power_val) # Record Qrad condition before time/Temp increment (array is always n-1 behind)
            self.temp_array.append(ktoc(slab.T_ml)) # Record Structure Temp condition before time/Temp increment (array is always n-1 behind)
            self.inc_T(self.Qinput_array[i]) # Increment structure temp based on Qinput and time step
            print('Simulation time: ', (self.time_array[i]/(60*60)), ' hours')                 
            #emission_spectras.append(slab.thermal_emission_array*um)


    #def analytical_time_dependent_temp(self):
        

    
    def plot_time_dependent_temp(self):
        fig, ax1 = plt.subplots()
        ax1.plot(self.time_array/(60*60), self.temp_array, 'b-')
        ax1.set_xlabel('Time (hour)')
        # Make the y-axis label, ticks and tick labels match the line color.
        ax1.set_ylabel('Temperature (C)', color='b')
        ax1.tick_params('y', colors='b')
        
        ax2 = ax1.twinx()
        ax2.plot(self.time_array/(60*60), self.Q_rad_array, 'r-.')
        #ax2.plot(slab.sim_time, Q, 'k' )
        ax2.set_ylabel('Radiative Power ($W/m^2$)', color='r')
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
        'Lambda_List': [0.2*um, 40*um, 2000], #SI unit of m
        'Time_List':[0,2*hour, 20], #SI unit of s  
        'Sun_Period': 1*hour,
        
        # Structure depenednt properties
        'Internal_Load_Flux': 10,
        'Thermal_Capacitance': 1500,  #SI unit of J/(kg*K)
        'Material_List': ['Air', 'Air', 'Air'],
        'Thickness_List': [0, 0, 0], # You can not have the back reflector as the last layer!!!
        
        ## Flags
        'EXPLICIT_ANGLE': 1,
        'COOLING': 1
        }


#slab = radiator(structure)
#slab.step_emissivity(2*um, 5*um, plot = False)
#slab.CPvT(plot=True, temp_difference = False)
#slab.plot_CPvT
#%%

slab = []
for i in range(30):
    start_time = time.time()
    
    slab.append(radiator(structure))
    make_structure_time = time.time()
    
    slab[i].step_emissivity((1+i)*um, (6+i)*um, plot = False)
    make_emissivity_time = time.time()
    
    slab[i].CPvT(plot=False, temp_difference = False)
    CPvT_time = time.time()
    
    print('Make Structure time: ', make_structure_time-start_time, ' s')
    print('Make Emissivity time: ', make_emissivity_time-make_structure_time, ' s')
    print('Make CPV time: ', CPvT_time-make_emissivity_time, ' s')




#%%

plt.figure()
for i in range(len(slab)):
    plt.plot(ktoc(slab[i].temp_array), slab[i].total_emissivity_array,
             label = 'Spectral Window: ', (1+i), '-', (6+i),'um')
plt.set_xlabel('Temperature (C)')
plt.set_ylabel('Total Emissivity')
plt.legend()
fig.tight_layout()
plt.show()

plt.figure()
for i in range(len(slab)):
    plt.plot(ktoc(slab[i].temp_array), np.array(slab[i].Q_rad_array),
             label = 'Spectral Window: ', (1+i), '-', (6+i),'um')
plt.plot(ktoc(slab[i].temp_array),np.array(self.BB_total_array)
    ,':', label = 'Blackbody Emission')
plt.set_xlabel('Temperature (C)')
plt.set_ylabel('Radiative Emission ($W/m^2$)')
plt.legend()
fig.tight_layout()
plt.show()










#%%
#plt.figure()
#for i in range(len(slab.BB_array)):
#    plt.plot(slab.lambda_array/um, slab.BB_array[i]*um)

#slab = time_dependent_radiator(structure)
#slab.step_emissivity(5*um, 10*um, False)
#slab.time_dependent_temp()
#slab.plot_time_dependent_temp()
#plt.figure()
#plt.plot(slab.time_array)
#plt.plot(slab.time_array/hour, slab.sun_mask)

#%%
#


#plt.figure()
#plt.plot(ktoc(slab.T_array), slab.Q_rad_array)
#plt.plot(ktoc(slab.T_array), slab.e_total_array)


#%%
    

    
#%%    
#fig, ax1 = plt.subplots()
#ax1.plot(slab.lambda_array/um, slab.emissivity_array, 'b-')
#ax1.set_xlabel('Wavelength (um)')
## Make the y-axis label, ticks and tick labels match the line color.
#ax1.set_ylabel('Emissivity', color='b')
#ax1.tick_params('y', colors='b')
#
#ax2 = ax1.twinx()
#ax2.plot(slab.lambda_array/um, np.transpose(emission_spectras))
##ax2.plot(slab.sim_time, Q, 'k' )
#ax2.set_ylabel('Thermal Emission ($W/m^2/ \mu m$)', color='r')
#ax2.tick_params('y', colors='r')
#
#fig.tight_layout()
#plt.show()   


#%%
#sorted_idx = np.argsort(T)
#plt.figure()
#plt.plot(np.array(T)[sorted_idx], np.array(CP)[sorted_idx])
#
#plt.figure()
#plt.plot(slab.sim_time, CP)

#plt.figure()
#plt.plot(Q)


#%% Get Cooling Power as a Function of Temperature for different Emission Spectras

#slab = []
#slab.append(radiator(structure))
#slab[0].T_ml = ctok(0)
#slab[0].T_amb = 2
#slab.append(radiator(structure))
#slab[1].T_ml = ctok(0)
#slab[1].T_amb = 2
##slab[0].tmm()
#slab[0].step_emissivity(0.2*um, 20*um, False)
#slab[1].step_emissivity(5*um, 6*um, False)
#
##%%
##plt.figure()
#slab[1].plot_emissivity()
##plt.show()
#
#
#
##%%
#slab[0].CPvT(ctok(-200), ctok(200))
#slab[1].CPvT(ctok(-200), ctok(200))
#
##%%
#plt.figure()
#plt.plot(ktoc(slab[0].T_array), slab[0].TCP)
#plt.plot(ktoc(slab[1].T_array), slab[1].TCP)






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























































