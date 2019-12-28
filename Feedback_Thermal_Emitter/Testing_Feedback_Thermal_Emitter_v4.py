

from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
from wptherml.wptherml.wpml_v3 import multilayer
import tmm.tmm_core as tmm
from numpy import linspace, inf, pi, stack, array
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy as sci


import copy
import matplotlib as mplib
#import random 
#import matplotlib.animation as animation
#from wptherml.wptherml import coolinglib as cool
import time


mplib.rcParams['lines.linewidth'] = 8
mplib.rcParams['lines.markersize'] = 6
mplib.rcParams['axes.titlesize'] = 30
mplib.rcParams['axes.labelsize'] = 26
mplib.rcParams['xtick.labelsize'] = 24
mplib.rcParams['ytick.labelsize'] = 24
mplib.rcParams['font.size'] = 26


#%%



def ktoc(T):
    #CONVERT K TO C TEMP.
    C = T-273.15
    return C

def ctok(T):
    #CONVERT C TO K TEMP.
    K = T+273.15
    return K


# The radiator class is used to simulate non-time evolving relationships 
# between object temperature and spectral emissivity window.
# The following can be calculated:
# 1) an arbitrary wavelength square function emissivity profile using the 
#   "step_emissivity" function.
# 2) The wavelength integrated emissivity at a given temperature, using the 
#   "get_e_total" function.     
# 3) The wavelength integrated black body emission, using the 
#   "integrate_BB" function.
# 4) The relationship between a structures BB and total emissivity as a 
#   function of temperature, using the "CPvT" function.
# Note: The CPvT function will populate the following arrays
# BB_total_array = integrated BB as a function of temp.
# Q_rad_array = structures radiative emission as a function of temp.
# total_emissivity_array = total emissivity as a function of temp.





class radiator(multilayer):

    def __init__(self, args):
       
        if 'Thermal_Capacitance' in args:
            self.C = args['Thermal_Capacitance'] # J/(kg*K)
        else:
            self.C = 1000
            
        if 'Internal_Load_Flux' in args:
            self.Qload = args['Internal_Load_Flux'] # W/m^2
        else:
            self.Qload = 1000       
            
        self.Qinput_array = []    
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
            
            self.thermal_emission_ea()
            self.thermal_emission()
#            self.cooling_calc = False
#            self.update()
#            self.cooling_calc = True
            
        if (plot):
            self.plot_emissivity()

    def sawtooth_emissivity(self, lda0 = 0, ldaf = 0, plot = False):
        # MAKE SAWTOOTH FUNCTION EMISSIVITY
        LY = 1
        LX = ldaf-lda0

        theta = np.arctan(LY/LX)
        
        for i in range(0,len(self.t)):
            for j in range(0,len(self.lambda_array)):
                if (self.lambda_array[j] > lda0 and self.lambda_array[j] < ldaf):
                    LX = ldaf-self.lambda_array[j]
                    #H = np.sqrt(LX**2+LY**2)
                    value = LX*np.tan(theta)
                    
                    self.emissivity_array_p[i,j] = self.emissivity_array_s[i,j] =  value
                    self.emissivity_array[j] = value            
                else:
                    self.emissivity_array_p[i,j] = self.emissivity_array_s[i,j] =  0
                    self.emissivity_array[j] = 0     
            
            self.thermal_emission_ea()
            self.thermal_emission()
        
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
        Tf = 500 # ctok(100.1)
        To = 100 #ctok(-100)
        dT = 10
        self.temp_array = np.arange(To,Tf,dT)
        for T in self.temp_array:
            self.T_ml = T
            self.thermal_emission_ea()
            self.cooling_power(radiative = True, atmospheric=False, solar=False, total = False )
            #self.update()
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
        self.Q_Boltzmann = sigma*np.array(self.T_ml)**4 
    #np.array(self.total_emissivity_array)*
    
    
    def optimal_total_emissivity(self):
        
        self.temp_array = []
        #self.temp_array.append(self.T_ml)
        self.total_emissivity_array = np.arange(0.01,1.01,0.01)
        dT = 0.5
    
        for e in self.total_emissivity_array:
           
            # Increment temperature until power flow is below threshold.
            # Since temperature starts at the ideal temp. The first time 
            # temp is below threshold should be the smallest dT at which 
            # the condition is satisfied, for that given emissivity.
            pflow = 100
            self.T_ml = self.T_ideal
            
            while (abs(pflow)>10):
                self.Boltzmann_law()
                pstructure = e*self.Q_Boltzmann
                pflow = self.Qload-pstructure
                if (pflow > 0):
                    self.T_ml = self.T_ml+dT
                    #self.temp_array.append(self.T_ml)
                if (pflow < 0):
                    self.T_ml = self.T_ml-dT
                    #self.temp_array.append(self.T_ml)
               # print(str(self.T_ml))  
               # print(str(e))
#                print(str(pflow))
#                print(str(pstructure))

            
            self.temp_array.append(self.T_ml)
    
        # Result are the closest temperatures to the operating temperature
        # that satisies the steady state power balance, for a given Qi,
        # as a function of total emissivity.

    
    def QvT(self):
        
        optimal_e = []
        min_temp = []
        dT_array = []
        self.Qinput_array = np.arange(10,1300.1,10)
        for q in self.Qinput_array:
            self.Qload = q
            self.optimal_total_emissivity()
            dT =(np.array(self.temp_array)-np.array(self.T_ideal))
            dT_array.append(dT)
            min_temp.append(min(abs(dT)))
            min_temp_idx = np.argmin(abs(dT))
            optimal_e.append(self.total_emissivity_array[min_temp_idx])
           # print('For an input power = ' + str(q) + 'W/m^2')
           # print('Minimum change in temp is ' + str(min(dT)))
           # print('From a total emissivity of ' + str(self.total_emissivity_array[min_temp_idx]))
        return min_temp, optimal_e, dT_array
    
    
    
    
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
 
    def plot_QvT(self, dT_array, X = 0, Y = 0):

        levels = np.arange(-100,100.5,1) 
        plt.figure()
        plt.contourf(slab.total_emissivity_array, 
                     slab.Qinput_array, 
                     dT_array,levels) 

        plt.xlabel('Total emissivity')
        plt.ylabel('Input power ($W/m^2$)')
        plt.ylim(0,1300)
        plt.xlim(0,1)    
        cbar = plt.colorbar(ticks = np.arange(-100,100.1,20))
        cbar.ax.set_ylabel('Temperature change \n for equilibrium ($\degree$K or $\degree$C)')
        plt.plot(optimal_e, slab.Qinput_array,'k-', linewidth = 6 , label = 'Ideal')
        
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        colors = [
                  'darkorange',
                  'purple',
                  'hotpink',
                  'cyan',
                  ]
        if (len(X)>1):
            for i in range(1,len(X)):
                plt.plot(X[i], Y[i], 
                         '-', 
                         linewidth = 6, 
                         label = 'Triangle (3-'+str(5*i+3)+'um)',
                         color = colors[i])
                
            for i in range(1,len(X2)):
                plt.plot(X2[i], Y2[i],
                         ':',
                         linewidth = 6,
                         label = 'Rectangle (3-'+str(5*i+3)+'um)',
                         color = colors[i])        
                
                
        plt.legend()
        
        
        
        
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
            


        radiator.__init__(self,args)
        self.define_Qinput()    
                
        
    def define_Qinput(self):
        # DEFINE THE PULSE SHAPE FOR THE INPUT EXCITATION. THIS IS ASSUMED TO BE 
        # A CONSTANT "CURRENT = POWER" SOURCE. I.E., dT/dt. 
        # THE PURPOSE OF THIS FUNCTION IS TO EXPAND IT TO ENABLE ARBITRATY INPUT
        # WAVEFORMS THAT CAN REPRESENT ORBITAL DYNAMICS.
        
        #dl = self.lambda_array[1]-self.lambda_array[0]
        #am15 = datalib.AM(self.lambda_array)  #*self.emissivity_array
        #am15_integral = am15.sum()*dl
        self.Qsolar = 1150 #am15_integral*0.6  #0.1 For 3-11 window
        
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
         self.thermal_emission_ea()
         self.thermal_emission()
         self.cooling_power(radiative = True, atmospheric=False, solar=False, total = False )         
         
         #self.update()
         #self.cooling_power()

         
    def time_dependent_temp(self, plot = False):  
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

        if plot:
            self.plot_time_dependent_temp()
    #def analytical_time_dependent_temp(self):
        

    
    def plot_time_dependent_temp(self):
        fig, ax1 = plt.subplots()
        ax1.plot(self.time_array/(60), self.temp_array, 'b-')
        ax1.set_xlabel('Hubble orbit time (min)')
        # Make the y-axis label, ticks and tick labels match the line color.
        ax1.set_ylabel('Structure temperature (C)', color='b')
        ax1.tick_params('y', colors='b')
        
        ax2 = ax1.twinx()
        ax2.plot(self.time_array/(60), self.Q_rad_array, 'r-.')
        #ax2.plot(slab.sim_time, Q, 'k' )
        ax2.set_ylabel('Structure radiative power ($W/m^2$)', color='r')
        ax2.tick_params('y', colors='r')
        
        fig.tight_layout()
        plt.show()   
    
        
def pw_plot(x):
    plt.figure()
    plt.plot(x)


       
#%            
######################################################################################################                    
                    
# DEFINE SOME UNITS.        
nm = 1e-9
um = 1e-6 
minute = 60
hour = 60*minute

sio2_specific_heat = (730 + 680)/2
sio2_density = 1000*(2.65 + 2.17)/2  #kg/m^3
ocradecane_specific_heat = 244*1e+3
#cooling_film_area = 0.06 # m^2
cooling_film_thickness = 1000*um # m
#cooling_film_volume = cooling_film_area*cooling_film_thickness
#cooling_film_mass = sio2_density*cooling_film_volume
thermal_capacitance = 3000 #cooling_film_thickness*sio2_density*sio2_specific_heat

# DEFINE STRUCTURE.
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        
        'Structure_Temperature': ctok(40),
        'Ambient_Temperature': 2,
        'Lambda_List': [0.2*um, 40*um, 2000], #SI unit of m
        'Time_List':[0,1.2*94*minute, 1.2*4*100], #SI unit of s  
        'Sun_Period': 47*minute,
        
        # Structure depenednt properties
        'Internal_Load_Flux': 40,
        
        # Thermal capacitance = specific heat * mass. 
        # Specific heat of SiO2 ~680-730 J/(kg*K)
        'Thermal_Capacitance': thermal_capacitance,  #SI unit of J/(K)
        'Material_List': ['Air', 'Air', 'Air'],
        'Thickness_List': [0, 0, 0], # You can not have the back reflector as the last layer!!!
        
        ## Flags
        'EXPLICIT_ANGLE': 1,
        'COOLING': 1
        }
QvT_flag = False


slab = time_dependent_radiator(structure)
slab.step_emissivity(3*um, 15*um, plot = False)
slab.time_dependent_temp(plot = True)


#%%
#slab.sawtooth_emissivity(3*um, 4*um, True)
slab.plot_time_dependent_temp()


if QvT_flag: 
    
    start = time.time()
    [min_temp, optimal_e, dT_array] = slab.QvT()
    print('Elappsed time: ' + str(time.time()-start))

    slab2 = radiator(structure)
    
    X =[]
    Y = []
    start = time.time()
    for ldaf in np.arange(5,20.1,5):
        slab2 = radiator(structure)
        slab2.sawtooth_emissivity(3*um, ldaf*um, plot = False)
        slab2.CPvT()
        X.append(slab2.total_emissivity_array)
        Y.append(slab2.Q_rad_array)
        del slab2
    print('Elappsed time: ' + str(time.time()-start))       
    #%
    
    X2 =[]
    Y2 = []
    for ldaf in np.arange(5,20.1,5):
        slab2 = radiator(structure)
        slab2.step_emissivity(3*um, ldaf*um, plot = False)
        slab2.CPvT()
        X2.append(slab2.total_emissivity_array)
        Y2.append(slab2.Q_rad_array)
        del slab2        
#%%
 


       
slab = radiator(structure)
slab.step_emissivity(3*um, 13*um, plot = False)
slab2 = radiator(structure)
slab2.sawtooth_emissivity(3*um, 13*um, plot = False)

#%%

mplib.rcParams['lines.linewidth'] = 12
mplib.rcParams['lines.markersize'] = 12
mplib.rcParams['axes.titlesize'] = 30
mplib.rcParams['axes.labelsize'] = 36
mplib.rcParams['xtick.labelsize'] = 26
mplib.rcParams['ytick.labelsize'] = 26
mplib.rcParams['font.size'] = 36

fig, ax = plt.subplots()
ax.plot(slab.lambda_array/um, slab.emissivity_array, 'k', label = 'Rectangle \n (3-13um)', linewidth = 12)
ax.plot(slab2.lambda_array/um, slab2.emissivity_array, 'c:', label = 'Triangle \n (3-13um)', linewidth = 12)
ax.set_xlim(1, 20)
ax.set_xlabel('Wavelength (um)')
ax.set_ylabel('Emissivity')
leg = plt.legend()
leg.draggable()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#%%
## Set up figure settings
fig, ax1 = plt.subplots(figsize=(20,10))
ax2 = ax1.twinx()
ax2.set_ylim(0, 1.1e7*1e-6)

#plt.xlim(min(slab.lambda_array)*1e6, max(slab.lambda_array)*1e6)
#plt.ylim(0, 1.1e7*1e-6)
#plt.xlabel('Wavelength ($um$)',fontsize=20)
#plt.ylabel('Power density ($W \cdot m^{-2} \cdot um^{-1}$)',fontsize=20)
#mng = plt.get_current_fig_manager()
#mng.full_screen_toggle()
#plt.title('Ideal Spectrum',fontsize=20)
ln1, = plt.plot([],[],'k', label = 'Structure emmisivity')
#ln2, = plt.plot([],[],'b', label = 'Atmospheric emmisivity')
#ln3, = plt.plot([],[], label = 'Structure BB')
#temp_text = ax1.text(0.5, 0.8, '', transform=ax1.transAxes)

def init():

    
    
    # temp_text.set_text('')
    return ln1
#, fill1, fill2
        
def animate(i):
    global slab, T, lda0, ldaf, ax1, ax2
    slab.T_ml = T[i]
    slab.update()
    slab.step_emissivity(lda0, ldaf, False)
    slab.get_e_total()
    slab.Boltzmann_law()
   
    ax1.clear()
    ax2.clear()

    ax1.plot(slab.lambda_array*1e6, slab.emissivity_array, 'k-')
    ax1.set_ylabel('Emissivity')
    ax1.set_xlabel('Wavelength (um)')
    
    ax2.plot(slab.lambda_array*1e6, slab.BBs*1e-6, 
             label = 'Black body at ' + str(T[i]) +'K')   
    ax2.legend(loc = 'center right')

    #ax1.title('Emissivity (' + str(round(lda0*1e6))+'-'+str(round(ldaf*1e6))+'um)')
    #ln1.color(colors[i+1])
    #ax1.collections.clear()
    ax2.fill_between(slab.lambda_array*1e6,0,slab.thermal_emission_array*1e-6,
                             color = 'r',
                             alpha=0.5)
    ax2.set_ylabel('Thermal emission ($W/m^2um$)')
    #ax2.set_ylim(0, 1.1e7)
    return 
    
lda0 = 5*um
ldaf = 15*um    
 
T = [500,400,300,250,200,150,100]
ani = animation.FuncAnimation(fig, animate, frames=len(T), init_func = init,  repeat=False)

plt.show()
ani.save('ATEST.mp4', fps = 1)  
    

    
    
    
#%%
    
 ## Set up figure settings
fig, ax1 = plt.subplots(figsize=(20,10))
ax2 = ax1.twinx()
ax2.set_ylim(0, 1.1e7*1e-6)

#plt.xlim(min(slab.lambda_array)*1e6, max(slab.lambda_array)*1e6)
#plt.ylim(0, 1.1e7*1e-6)
#plt.xlabel('Wavelength ($um$)',fontsize=20)
#plt.ylabel('Power density ($W \cdot m^{-2} \cdot um^{-1}$)',fontsize=20)
#mng = plt.get_current_fig_manager()
#mng.full_screen_toggle()
#plt.title('Ideal Spectrum',fontsize=20)
ln1, = plt.plot([],[],'k', label = 'Structure emmisivity')
#ln2, = plt.plot([],[],'b', label = 'Atmospheric emmisivity')
#ln3, = plt.plot([],[], label = 'Structure BB')
#temp_text = ax1.text(0.5, 0.8, '', transform=ax1.transAxes)


def init():

    
    
    # temp_text.set_text('')
    return ln1
#, fill1, fill2
        
def animate(i):
    global slab, T, lda0, ldaf, ax1, ax2
    slab.T_ml = T[i]
    slab.update()
    slab.step_emissivity(lda0, ldaf, False)
    slab.get_e_total()
    slab.Boltzmann_law()
   
    #ax1.clear()
    #ax2.clear()

    ax1.plot(slab.Tml, slab.e_total, 'k-')
    ax1.set_ylabel('Emissivity')
    ax1.set_xlabel('Wavelength (um)')
    
#    ax2.plot(slab.lambda_array*1e6, slab.BBs*1e-6, 
#             label = 'Black body at ' + str(T[i]) +'K')   
#    ax2.legend(loc = 'center right')
#
#    #ax1.title('Emissivity (' + str(round(lda0*1e6))+'-'+str(round(ldaf*1e6))+'um)')
#    #ln1.color(colors[i+1])
#    #ax1.collections.clear()
#    ax2.fill_between(slab.lambda_array*1e6,0,slab.thermal_emission_array*1e-6,
#                             color = 'r',
#                             alpha=0.5)
#    ax2.set_ylabel('Thermal emission ($W/m^2um$)')
#    #ax2.set_ylim(0, 1.1e7)
    return 
    
lda0 = 5*um
ldaf = 15*um    
 
T = [500,400,300,250,200,150,100]
ani = animation.FuncAnimation(fig, animate, frames=len(T), init_func = init,  repeat=False)

plt.show()
ani.save('ATEST.mp4', fps = 1)     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#plt.plot(slab3.total_emissivity_array, slab3.Q_rad_array)
#%%
#plt.figure()
#plt.plot(slab.Qinput_array,optimal_e)
#
#plt.figure()
#plt.plot(slab.Qinput_array,min_temp)
##plt.figure, 
#plt.plot(slab.total_emissivity_array, slab.temp_array-np.array(slab.T_ideal))




#slab = time_dependent_radiator(structure)
##slab.step_emissivity(3*um, 11*um, plot = False)  #3-11 is best
#slab.step_emissivity(5*um, 15*um, plot = False)  #3-11 is best
#slab.time_dependent_temp(plot=True)

#slab.CPvT(plot=True, temp_difference = False)
#slab.plot_CPvT



#
#
#
##%%

slab = radiator(structure)
lda0 = 5
dlda = 10
slab = []
for i in range(3):
    start_time = time.time()
    
    slab.append(radiator(structure))
    make_structure_time = time.time()
    
    slab[i].step_emissivity((lda0+i*5)*um, (lda0+dlda+i*5)*um, plot = False)
    make_emissivity_time = time.time()
    
    slab[i].CPvT(plot=False, temp_difference = False)
    CPvT_time = time.time()
    
    print('Make Structure time: ', make_structure_time-start_time, ' s')
    print('Make Emissivity time: ', make_emissivity_time-make_structure_time, ' s')
    print('Make CPV time: ', CPvT_time-make_emissivity_time, ' s')




#%%

plt.figure()
for i in range(len(slab)):
    plt.plot((slab[i].temp_array), slab[i].total_emissivity_array,
             label = 'Spectral Window: ' + str(lda0+i*5) + '-' + str(lda0+dlda+i*5) +'um')
plt.xlabel('Temperature (K)')
plt.ylabel('Total Emissivity')
leg = plt.legend()
leg.draggable()
#plt.tight_layout()
plt.show()

##%%
#plt.figure()
#for i in range(len(slab)):
#    plt.plot(ktoc(slab[i].temp_array), np.array(slab[i].Q_rad_array),
#             label = 'Spectral Window: ' + str(lda0+i) + '-' + str(lda0+dlda+i)+'um')
##plt.plot(ktoc(slab[i].temp_array),np.array(slab[i].BB_total_array)
##    ,':', label = 'Blackbody Emission')
#plt.xlabel('Temperature (C)')
#plt.ylabel('Radiative Emission ($W/m^2$)')
#leg = plt.legend()
#leg.draggable()
##plt.tight_layout()
#plt.show()

############################################################################################3
#%%
#np.save('Slabs_5umBoxcar',slab)
#io.savemat('Slabs_5umBoxcar.mat',{'slabs': slab}, long_field_names = True)





#%%
plt.figure()
for i in range(len(slab.BB_array)):
    plt.plot(slab.lambda_array/um, slab.BB_array[i]*um)

slab = time_dependent_radiator(structure)
slab.step_emissivity(5*um, 10*um, False)
slab.time_dependent_temp()
slab.plot_time_dependent_temp()
plt.figure()
plt.plot(slab.time_array)
plt.plot(slab.time_array/hour, slab.sun_mask)

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























































