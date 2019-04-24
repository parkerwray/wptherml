# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:12:53 2019
This code plots the index of material data stored in wptherml
@author: parkerwray
"""
# from [file location] import [.py file]
from wptherml.datalib import datalib
import matplotlib.pyplot as plt
import numpy as np


lam = np.linspace(150e-9,1000e-9,1000)
m = datalib.Material_RI(lam, 'SiO2NP_model')

c = 299792458
h = 4.135667516*10**(-15)
e = m**2
ev = h*c/lam

lam_min = np.min(lam)
lam_max = np.max(lam)

if (lam_min <= 1000e-9):
    nm = 1e+9  
    #  Plot the real part of the refractive index in the visible regime  
    mask = (lam>=lam_min) & (lam <= np.min([1000e-9, lam_max]))
    
    plt.fig, ax1 = plt.subplots()
    ax1.plot(lam[mask]*nm, m[mask].real, 'k-')
    ax1.autoscale()
    ax1.set_ylabel('n', color = 'k')  
    ax1.set_xlabel('Wavelength (nm)')    
    ax2 = ax1.twinx()
    ax2.plot(lam[mask]*nm, m[mask].imag, 'r:')
    ax2.autoscale()
    ax2.set_ylabel('k', color = 'r')        
    plt.title('Refractive index in the visible')
    plt.show()
    
    plt.fig, ax1 = plt.subplots()
    ax1.plot(ev[mask], e[mask].real, 'k-')
    ax1.autoscale()
    ax1.set_ylabel('er', color = 'k')  
    ax1.set_xlabel('Energy (ev)')    
    ax2 = ax1.twinx()
    ax2.plot(ev[mask], e[mask].imag, 'r:')
    ax2.autoscale()
    ax2.set_ylabel('ei', color = 'r')        
    plt.title('Epsilon in the visible')
    plt.show()
    
    
    
    
    
    

if (lam_max > 1000e-9):
    um = 1e+6  
    #  Plot the real part of the refractive index in the Infrared regime   
    mask = (lam<=lam_max) & (lam >= np.max([1000e-9, lam_min]))
    plt.fig, ax1 = plt.subplots()
    ax1.plot(lam[mask]*um, m[mask].real, 'k-')
    ax1.autoscale()
    ax1.set_ylabel('n', color = 'k')  
    ax1.set_xlabel('Wavelength (um)')    
    ax2 = ax1.twinx()
    ax2.plot(lam[mask]*um, m[mask].imag, 'r:')
    ax2.autoscale()
    ax2.set_ylabel('k', color = 'r')        
    plt.title('Refractive index in the IR')
    plt.show()


