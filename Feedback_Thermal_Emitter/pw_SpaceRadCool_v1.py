from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
from wptherml.wptherml.wpml_v3 import multilayer
from numpy import linspace, inf, pi, stack, array
import numpy as np
import matplotlib.pyplot as plt
import scipy as sci



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
    
    
    def update_emissivity(self, e_array, plot = False):
        
        for i in range(0,len(self.t)):
            self.emissivity_array_p[i,:] = self.emissivity_array_s[i,:] =  e_array
           
        self.emissivity_array = e_array     
        self.thermal_emission_ea()
        self.thermal_emission()
    
        if (plot):
            self.plot_emissivity()
    
    
    
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
            self.update()
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
              
    def CPvT(self, plot = False, temp_difference = False):
        # CALCULATE THE CHANGE IN COOLING POWER AS A FUNCTION OF TEMPERATURE 
        # THIS IS BASED ON THE EMISSIVITY PROFILE OF THE STRUCTURE AND IS NOT 
        # TIME DEPENDENT
        Tf = ctok(150.1) # ctok(100.1)
        To = ctok(-100) #ctok(-100)
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
 
    def plot_QET(self, dT_array, X = 0, Y = 0):

        levels = np.arange(-100,100.5,1) 
        plt.figure()
        plt.contourf(self.total_emissivity_array, 
                     self.Qinput_array, 
                     dT_array,levels) 

        plt.xlabel('Total emissivity')
        plt.ylabel('Input power ($W/m^2$)')
        plt.ylim(0,1300)
        plt.xlim(0,1)    
        cbar = plt.colorbar(ticks = np.arange(-100,100.1,20))
        cbar.ax.set_ylabel('Temperature change \n for equilibrium ($\degree$K or $\degree$C)')
        plt.plot(self.optimal_e, self.Qinput_array,'k-', linewidth = 6 , label = 'Ideal')
        
#        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#        colors = [
#                  'darkorange',
#                  'purple',
#                  'hotpink',
#                  'cyan',
#                  ]
#        if (len(X)>1):
#            for i in range(1,len(X)):
#                plt.plot(X[i], Y[i], 
#                         '-', 
#                         linewidth = 6, 
#                         label = 'Triangle (3-'+str(5*i+3)+'um)',
#                         color = colors[i])
#                
#            for i in range(1,len(X2)):
#                plt.plot(X2[i], Y2[i],
#                         ':',
#                         linewidth = 6,
#                         label = 'Rectangle (3-'+str(5*i+3)+'um)',
#                         color = colors[i])        
#                
                
        plt.legend()
        
    def plot_TEQ(self, dT_array):

        #levels = np.arange(-100,100.5,1) 
        plt.figure()
        plt.contourf(self.total_emissivity_array, 
                     dT_array,
                     self.Qinput_array) 

        plt.xlabel('Total emissivity')
        plt.ylabel('Temperature change \n for equilibrium ($\degree$K or $\degree$C)')
        #plt.ylim(0,1300)
        plt.xlim(0,1)    
        #cbar = plt.colorbar(ticks = np.arange(-100,100.1,20))
        #cbar.ax.set_ylabel('Input power ($W/m^2$)')
        #plt.plot(self.optimal_e, self.Qinput_array,'k-', linewidth = 6 , label = 'Ideal')
        
        
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
            self.temp_array.append(ktoc(self.T_ml)) # Record Structure Temp condition before time/Temp increment (array is always n-1 behind)
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


       







