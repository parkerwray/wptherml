
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert
from scipy.interpolate import InterpolatedUnivariateSpline

q = 1.60217662e-19
c=299792458
h=6.626e-34
k=1.38064852e-23


### functions for refractive index of different common dielectrics
def Material_RI(lam, arg):
    ## array of wavelengths should be in meters
    ## but we might need it in nm or microns
    l_nm = lam*1e9
    lmic = lam*1e6
    if (arg=='HfO2_model'):
        A = 187178
        B = 9993.46
        C = 0.0173801
        D = 1.939999
        n = A/(l_nm**4) + B/(l_nm**2) + C/l_nm + D + 0j/l_nm
        
    elif (arg=='Al2O3_model'):
        A = 187178
        B = 9993.46
        C = 0.0173801
        D = 1.69999
        n = A/(l_nm**4) + B/(l_nm**2) + C/l_nm + D + 0j/l_nm
    ### This model works well for glass into near IR
    
    elif (arg=='SiO2_model' ):
        print('This model for SiO2 works well for lambda < 5um (NIR)' )
        A = 187178
        B = 9993.46
        C = 0.0173801
        D = 1.45
        n = A/(l_nm**4) + B/(l_nm**2) + C/l_nm + D + 0j/l_nm
        
    elif (arg=='TiO2_model'):
        print('This model for TiO2 works well for lambda < 5um (NIR)' )
        A = 187178
        B = 9993.46
        C = 0.0173801
        D = 2.4
        n = A/(l_nm**4) + B/(l_nm**2) + C/l_nm + D + 0j/l_nm
        
    elif (arg=='AlN_model'):
        print('This model for AlN works well for lambda < 10um (NIR)' )
        A = 1.859
        B = 0.3401
        n = A + B/(lmic*lmic) + 0j/lmic
        
    elif (arg=='Air'):
        A = 0.
        n = A/lam + 0j/lam + 1
        
    elif (arg=='TiN_model'):
        n = TiN_Drude_Lorentz(lam)

    elif (arg=='a-Al2O3'
          or arg=='Ag'
          or arg=='AlN'
          or arg=='Au'
          or arg=='HfN' 
          or arg=='Pd'
          or arg=='Pt'
          or arg=='Re' 
          or arg=='Rh' 
          or arg=='Ru'
          or arg=='Si'
          or arg =='Si3N4'
          or arg=='SiC'
          or arg=='SiO2'
          or arg=='TiO2' 
          or arg=='W' 
          ):
        n = Read_RI_from_File(lam, arg)
        
    elif (arg=='RC0_1A_Si'
          or arg=='RC0_1A_SiO2nox'
          or arg=='RC0_1A_SiO2brugg'
          or arg=='RC0_1A_SiO2rough'):
        n = Read_RI_from_File(lam, arg)      
        
    elif (arg=='RC0_1B_Si'
          or arg=='RC0_1B_SiO2'
          or arg=='RC0_1B_SiO2nox'
          or arg=='RC0_1B_SiO2brugg'
          or arg=='RC0_1A_SiO2rough'):
        n = Read_RI_from_File(lam, arg)            

    ### default is air    
    else:
        A = 0.
        n = A/lam + 0j/lam + 1
    return n


def TiN_Drude_Lorentz(lam):
    ci = 0+1j
    epsinf = 3.59
    rho = 1.09e-6
    eps0 = 8.854e-12
    hbar = 6.582e-16
    tau = 7.82e-16
    amp = 4.658129
    br = 4.3827
    en = 5.778

    l_nm = lam*1e9
    E = 1240./l_nm
    eps = epsinf - hbar*hbar/(eps0*rho*(tau*E*E + ci*hbar*E))
    eps = eps + amp*br*en/(en*en - E*E - ci*E*br)
    return np.sqrt(eps)



def Read_RI_from_File(lam, matname):

    if (matname=='W'):
        a = np.loadtxt('wptherml/datalib/W_Palik_RI_f.txt')
        
    elif (matname=='a-Al2O3'): 
        print('a-Al2O3 is valid for lda = 210nm - 55um' )
        a = np.loadtxt('wptherml/wptherml/datalib/a-Al2O3_Querry.txt')
        for i in range(0,len(a)):
            a[i][0] = a[i][0]*1e-9 # lda data in nm, convert to m
            
    elif (matname=='TiO2'): 
        a = np.loadtxt('wptherml/datalib/TiO2_Siefke.txt')
        
    elif (matname=='Re'):
        a = np.loadtxt('wptherml/datalib/Re_Palik_RI_f.txt')
        
    elif (matname=='Ru'):
        a = np.loadtxt('wptherml/datalib/Ru_Palik_RI_f.txt')
        
    elif (matname=='Rh'):
        a = np.loadtxt('wptherml/datalib/Rh_Palik_RI_f.txt')
        
    elif (matname=='Ag' and lam[len(lam)-1]<=1000e-9):
        a = np.loadtxt('wptherml/datalib/Ag_JC_RI_f.txt')
        
    elif (matname=='Ag' and lam[len(lam)-1]>1000e-9):
        a = np.loadtxt('wptherml/datalib/Ag_Yang.txt')
        
    elif (matname=='Au' and lam[len(lam)-1]<=1000e-9):
        a = np.loadtxt('wptherml/datalib/Au_JC_RI_f.txt')
        
    elif (matname=='Au' and lam[len(lam)-1]>1000e-9):
        a = np.loadtxt('wptherml/wptherml/datalib/Au_IR.txt')
        
    elif (matname=='Pd'):
        a = np.loadtxt('wptherml/datalib/Pd_Palik_RI_f.txt')
        
    elif (matname=='Pt'):
        a = np.loadtxt('wptherml/datalib/Pt_Palik_RI_f.txt')
        
    elif (matname=='SiO2'):
        a = np.loadtxt('wptherml/wptherml/datalib/SiO2_IR.txt')
        
    elif (matname=='AlN'):
        a = np.loadtxt('wptherml/datalib/AlN_IR.txt')
        
    elif (matname=='Si'):
        a = np.loadtxt('wptherml/wptherml/datalib/Si_Schinke.txt')
        
    elif (matname=='SiC'):
        a = np.loadtxt('wptherml/datalib/SiC_Larruquert.txt')  
        
    elif (matname=='Si3N4'):
        a = np.loadtxt('wptherml/datalib/Si3N4_Philipp.txt')  
         
    elif (matname=='W_Al2O3_Alloy'):
        a = np.loadtxt('wptherml/datalib/W_Al2O3_Alloy.txt')
        
    elif (matname=='RC0_1B_SiO2'):
        # This is a custom bulk visible SiO2 properties with broadband absorption. 
        # The goal of the broadband abs. is to match spectroscopic data.
        # Still need bruggeman model 
        # wavelength is in nm, so we convert to SI
        #a = np.loadtxt('wptherml/wptherml/datalib/rc0_1b_sio2_1_nk.txt')
        
        # This is a custom bulk visible SiO2 properties with no broadband absorption. 
        # There is some UV absorption to match spectroscopic data.
        # Still need bruggeman model 
        # wavelength is in nm, so we convert to SI
        a = np.loadtxt('wptherml/wptherml/datalib/rc0_1b_sio2_3_nk.txt')        

        for i in range(0,len(a)):
            a[i][0] = a[i][0]*1e-9 # lda data in nm, convert to m
            
    elif (matname=='RC0_1B_Si'):
        # This is a custom bulk visible Si properties 
        # The purpose of this is because the measured Si substrate is more 
        # reflective in the visible compared to the bulk material model that 
        # was previously used. ALSO, the bulk SI data doesn't go far into the 
        # IR! The loaded n,k data is from ellipsometery, which is saved in nm. 
        # so we convert to SI units!
        
        a_vis = np.loadtxt('wptherml/wptherml/datalib/rc0_1b_si_1_nk.txt')        
        a_ir = np.loadtxt('wptherml/wptherml/datalib/rc0_1a_si_p-type_ir_substrate_nk.txt')    
        
        a = np.concatenate((a_vis,a_ir), axis=0) #Mat applicable for Vis and IR
        for i in range(0,len(a)):
            a[i][0] = a[i][0]*1e-9 # lda data in nm, convert to m


    else:
        a = np.loadtxt('wptherml/datalib/W_Palik_RI_f.txt')

        
    ### now that we have read in the text, interpolate/extrapolate RI
    datlam = np.zeros(len(a))
    datn   = np.zeros(len(a))
    datk   = np.zeros(len(a))
    for i in range(0,len(a)):
        datlam[i]  = a[i][0]
        datn[i] = a[i][1]
        datk[i] = a[i][2]
    ### in case of duplicate wavelength values...     
    datlam, unique_idx = np.unique(datlam, return_index=True)
    datn = datn[unique_idx]
    datk = datk[unique_idx]
        
    ### use linear interpolation/extrapolation
    order = 1
    ### form the interpolator/extrapolator object for datn
    sn = InterpolatedUnivariateSpline(datlam, datn, k=order)
    ### form the interpolator/extrapolator object for datk
    sk = InterpolatedUnivariateSpline(datlam, datk, k=order)
    ### compute the interpolated/extrapolated values for real part of RI
    yn = sn(lam)
    ### compute the interpolated/extrapolated values for imaginary part of R
    yk = sk(lam)
    ### for complex RI array for each value of lambda
    n = yn + 1j*yk
    return n


### This function reads the data files for reflection data collected by UV-Vis
### and FTIR, respectively. The UV-Vis data is already normalized to absolute 
### reflection (this is done in a matlab program that makes the files loaded
### this python program). Both the UV-Vis and FTIR data have been conditioned
### and filtered such that they are completely usable, without need for further
### processing. NOTE: The visible data is normalized to Teflon (and void). The 
### Teflon has absorption for lda > 300nm. This will cause spectra to have
### artificially higher refleciton in these spectral regions. It is encouraged 
### to normalize a reference sample to known data online (or from simulation) 
### and apply the necessary correction factors from there. 
### This correction factor was not performed     
def Read_spectra_from_File(matname):
    
    if (matname=='RC0_1B_SiO2'):
        vis = np.loadtxt('wptherml/wptherml/datalib/measured_spectras/RC0_1B_sample_vis_TR_FR_DR.txt')
        ir = np.loadtxt('wptherml/wptherml/datalib/measured_spectras/RC0_1B_sample_ir_R_and_T.txt')      
    elif (matname=='RC0_1B_Si'):  #Need to re-order HfN data
        vis = np.loadtxt('wptherml/wptherml/datalib/measured_spectras/RC0_1B_ref_vis_TR_FR_DR.txt')
        ir = np.loadtxt('wptherml/wptherml/datalib/measured_spectras/RC0_1B_ref_ir_R_and_T.txt')
    elif (matname=='RC0_1D_Al2O3'):
        vis = np.loadtxt('wptherml/wptherml/datalib/measured_spectras/RC0_1D_sample_vis_TR_FR_DR.txt')
        ir = np.loadtxt('wptherml/wptherml/datalib/measured_spectras/RC0_1D_sample_ir_R_and_T.txt')      
    elif (matname=='RC0_1D_Si'):  #Need to re-order HfN data
        vis = np.loadtxt('wptherml/wptherml/datalib/measured_spectras/RC0_1D_ref_vis_TR_FR_DR.txt')
        ir = np.loadtxt('wptherml/wptherml/datalib/measured_spectras/RC0_1D_ref_ir_R_and_T.txt')        
    else:
        print('No known spectra selected')

    for i in range(0,len(vis)):
        vis[i][0] = vis[i][0]*1e-9 # lda data in nm, convert to m
    for i in range(0,len(ir)):
        ir[i][0] = ir[i][0]*1e-9 # lda data in nm, convert to m

    return vis, ir 


### Define the RI of a specified layer to be an alloy
### between two specified materials, mat1 and mat2,
### using Bruggenmans approximation
def alloy(lambda_array, fraction, mat1, mat2, model):
    ### Bruggeman model must be specified
    n = np.zeros(len(lambda_array),dtype=complex)
    if (model=='Bruggeman'):
        ### Get RIs of two materials... mat1 can be 
        ### a string that codes a material name or 
        ### it can be a single number
        if(isinstance(mat1, str)):
            n_1 = Material_RI(lambda_array, mat1)
        else:
            n_1 = mat1
            
        n_2 = Material_RI(lambda_array, mat2)
        
        for i in range(0,len(lambda_array)):
            if(isinstance(mat1, str)):
                eps1 = n_1[i]*n_1[i]
            else:
                eps1 = n_1*n_1
                
            eps2 = n_2[i]*n_2[i]
            flag = 1
            f1 = (1-fraction)
            f2 = fraction
            b = (2*f1-f2)*eps1 + (2*f2 - f1)*eps2
            arg = 8*eps1*eps2 + b*b
            srarg = np.sqrt(arg)
            
            if (np.imag(arg)<0):
                flag = -1
            else:
                flag = 1
                
            epsBG = (b+flag*srarg)/4.
            n[i] = np.sqrt(epsBG)
    #### Default is Maxwell-Garnett        
    else:
        if(isinstance(mat1, str)):
            n_1 = Material_RI(lambda_array, mat1)
        else:
            n_1 = mat1
            
        n_2 = Material_RI(lambda_array, mat2)
        f = fraction
        

        for i in range(0,len(lambda_array)):
            ### eps1 == epsD and eps2 == epsM in MG notation
            if(isinstance(mat1, str)):
                epsD = n_1[i]*n_1[i]
            else:
                epsD = n_1*n_1

            epsM = n_2[i]*n_2[i]
            num = epsD*(2*f*(epsM-epsD) + epsM + 2*epsD)
            denom = 2*epsD + epsM + f*(epsD-epsM)
            n[i] = np.sqrt((num/denom))
            
    return n


### returns interpolated/extrapolated EQE of monocrystaline silicon PV cells
### given an input array of wavelengths
def SR_Si(lam):
    ### values of lambda along which EQE is experimeintally known
    datlam = np.linspace(260e-9, 1310e-9, 22)
    ### experimental values of EQE for monocrystaline Si... 
    dateqe = 0.01*np.array([0., 0., 11.5, 23., 33., 37., 41., 45., 49., 52., 56., 60., 64., 62.5, 51., 35., 27.5, 20., 12.5, 7.5, 0., 0.])
    ### order of the spline
    order = 1
    ### form the interpolator/extrapolator object
    s = InterpolatedUnivariateSpline(datlam, dateqe, k=order)
    ### compute the interpolated/extrapolated values
    y = s(lam)
    
    return y
    


### returns interpolated/extrapolated EQE of InGaAsSb PV cells given an input
### array of wavelengths
def EQE_InGaAsSb(lam):
    ### values of lambda along which EQE is experimentally known
    datlam = np.linspace(1000.0e-9,3000.0e-9,21)
    ### values of EQE that were experimentally measured
    dateqe = 0.01*np.array([48., 50., 52.5, 54., 58., 59., 60., 62., 61., 62., 61.5, 59., 54., 22., 2., 0., 0., 0., 0., 0., 0.])
    
    ### use linear interpolation/extrapolation
    order = 1
    ### form the interpolator/extrapolator object
    s = InterpolatedUnivariateSpline(datlam, dateqe, k=order)
    ### compute the interpolated/extrapolated values
    y = s(lam)

    ### uncomment to plot experimental values
    '''
    plt.figure()
    plt.plot(datlam, dateqe, 'o')
    ### plot interpolated/extrapolated values
    plt.plot(lam, y, 'blue')
    '''
    ### return extrapolated EQE
    return y

def SR_InGaAsSb(lam):

    ### values of lambda along which EQE is experimentally known
    #datlam = np.linspace(1000.0e-9,3000.0e-9,21)
    datlam = 1e-9*np.array([200,320,1000,1100,1200,1300,1400,1500,
                            1600,1700,1800,1900,2000,2100,2200,2300,2400,
                            2500,2600,2700,2800,2900,3000])
    ### values of EQE that were experimentally measured
    dateqe = 0.01*np.array([0.,0.,48., 50., 52.5, 54., 58., 59., 60., 62., 61., 62., 61.5, 59., 54., 22., 2., 0., 0., 0., 0., 0., 0.])
    datsr = q*dateqe*datlam/(h*c)
    ### use linear interpolation/extrapolation
    order = 1
    ### form the interpolator/extrapolator object
    s = InterpolatedUnivariateSpline(datlam, datsr, k=order)
    ### compute the interpolated/extrapolated values
    y = s(lam)

    ### uncomment to plot experimental values
    '''
    plt.figure()
    plt.plot(datlam, datsr, 'o')
    ### plot interpolated/extrapolated values
    plt.plot(lam, y, 'blue')
    '''
    ### return extrapolated EQE
    return y

### spectral response for Silicon 
    

### returns interpolated/extrapolated spectral response (A/W) of GaSb PV cells
### given an input array of wavelenghts
def SR_GaSb(lam):
    ### values of lambda along which SR is experimentally known
    datlam = np.linspace(400.0e-9,2000.0e-9,17)
    ### values of SR that were experimentally measured
    datsr = 0.01*np.array([0., 0., 17.,34.,49.,58.,61. ,68., 73., 80., 85., 88., 87., 70.,  2., 0.,  0.])
    ### use linear interpolation/extrapolation
    order = 1
    ### form the interpolator/extrapolator object
    s = InterpolatedUnivariateSpline(datlam, datsr, k=order)
    ### compute the interpolated/extrapolated values
    y = s(lam)
    
    ### uncomment to plot experimental values
    '''
    plt.figure()
    plt.plot(datlam, datsr, 'o')
    plt.plot(lam, y, '-')
    '''
    return y

def BB(lam, T):
    ### speed of light in SI
    c = 299792458
    ### plancks constant in SI
    h = 6.62607004e-34
    ### boltzmanns constant in SI
    kb = 1.38064852e-23
    
    rho = np.zeros_like(lam)
    
    for i in range(0,len(rho)):
        rho[i] = 2*h*c*c/lam[i]**5 * 1./(np.exp(h*c/(lam[i]*kb*T))-1)
     
    #plt.figure()
    #plt.plot(lam, rho, '-')
    return rho



### Example for how to call the EQE_InGaAsSb function with a custom range of wavelengths!
#lam  = np.linspace(10e-9,3000e-9,1000)
#qe = SR_InGaAsSb(lam)
#SR = SR_GaSb(lam)

#B = BB(lam, 3000)



###AM 1.5
  
def AM(lam):  ###lam is x SI is y
    a = np.loadtxt('wptherml/wptherml/datalib/scaled_AM_1_5.txt')
    x = np.zeros(len(a))
    y = np.zeros(len(a))
    ###  issue was just that the AM1.5 data had wavelength
    ###  in nanometers and we were assuming it was in meters!
    ###  now it is converted to meters at the time that it
    ###  is stored to the array called x
    for i in range(0,len(a)):
        x[i] = a[i][0]
        y[i] = a[i][1]
    datlam = x
    dateqe = y
    order = 1
    s = InterpolatedUnivariateSpline(datlam, dateqe, k=order)
    z = s(lam)
    #plt.plot(x,y,'blue')
    #plt.plot(lam, z, 'r--')
    #plt.show()
    return z

def ATData(lam):
    a = np.loadtxt('wptherml/wptherml/datalib/ATrans.txt')
    x = np.zeros(len(a))
    y = np.zeros(len(a))

    for i in range(0,len(a)):
        x[i] = a[i][0]*1e-6
        y[i] = a[i][1]
    #plt.plot(x,y)
    datlam = x
    dateqe = y
    order = 1
    s = InterpolatedUnivariateSpline(datlam, dateqe, k=order)
    z = s(lam)
    #plt.plot(lam,z)
    #plt.show()
    return z
    

''' something wrong with reading this pl.txt text file!
def PhLum(lam):
    a = np.loadtxt('wptherml/datalib/pl.txt')
    x = np.zeros(len(a))
    y = np.zeros(len(a))
    ###  issue was just that the AM1.5 data had wavelength
    ###  in nanometers and we were assuming it was in meters!
    ###  now it is converted to meters at the time that it
    ###  is stored to the array called x
    for i in range(0,len(a)):
        x[i] = a[i][0]*1e-9
        y[i] = a[i][1]
    datlam = x
    datph = y
    order = 1
    s = InterpolatedUnivariateSpline(datlam, datph, k=order)
    z = s(lam)
    plt.plot(x,y,'blue')
    plt.plot(lam, z, 'r--')
    plt.show()
    return z
    
'''
### Given an array of lambda values, this function 
### will evaluate a Gaussian fit to photopic luminosity function
### and return an array of values from that fit.
def PhLum(lam):
    ## if changed change light lib functions

    ### gaussian parameters determined from fit to Photopic luminosity 
    ### function data, where raw data was accessed from here: 
    ### http://www.cvrl.org/database/data/lum/linCIE2008v2e_5.htm
    ### and fit was performed with gnuplot with the function f(x) = a*exp(-b*(x-c)**2)
    a               = 1.02433
    b               = 2.59462e+14
    c               = 5.60186e-07
    ### should be able to evaluate function at all values of array in one line
    z = a*np.exp(-b*(lam-c)**2)
    return z

### Read in CIE color matching functions from data file and interpolate/extrapolate
### on lam
def CIE(lam):
    a = np.loadtxt('wptherml/datalib/cie-cmf.txt')
    l = np.zeros(len(a))
    x = np.zeros(len(a))
    y = np.zeros(len(a))
    z = np.zeros(len(a))

    for i in range(0,len(a)):
        l[i] = a[i][0]*1e-9
        x[i] = a[i][1] 
        y[i] = a[i][2]
        z[i] = a[i][3]
        
    ## now interpolate over lam range using linear interpolation
    order = 1
    ## red response
    xint = InterpolatedUnivariateSpline(l, x, k=order)
    xbar = xint(lam)
    ## green response
    yint = InterpolatedUnivariateSpline(l, y, k=order)
    ybar = yint(lam)
    ## blue response
    zint = InterpolatedUnivariateSpline(l, z, k=order)
    zbar = zint(lam)

    ### store all the interpolated color matching functions 
    ### in a dictionary and return
    cie = {"xbar": xbar, 
         "ybar": ybar, 
         "zbar": zbar } 

    return cie

### Methods for plotting index data
    
#def plot_index(self, lam, material):
#    lam = self.lambda_array
#    lam_min = np.min(lam)
#    lam_max = np.max(lam)
#    
#    if (lam_min < 1000e-9):
#        self.plot_index_vis(layer)
#    
#    if (lam_max > 1000e-9):
#        self.plot_index_ir(layer)
#    return 1
#
#            
#def plot_index_vis(self, layer):
#    lam = self.lambda_array
#    RI = self.layer_ri(layer)
#    lam_min = np.min(lam)
#    lam_max = np.max(lam)
#    nm = 1e+9  
#    #  Plot the real part of the refractive index in the visible regime  
#    mask = (lam>=lam_min) & (lam <= np.min([1000e-9, lam_max]))
#    plt.fig, ax1 = plt.subplots()
#    ax1.plot(lam[mask]*nm, RI[mask].real, 'k-')
#    ax1.autoscale()
#    ax1.set_ylabel('n', color = 'k')  
#    ax1.set_xlabel('Wavelength (nm)')    
#    ax2 = ax1.twinx()
#    ax2.plot(lam[mask]*nm, RI[mask].imag, 'r:')
#    ax2.autoscale()
#    ax2.set_ylabel('k', color = 'r')        
#    plt.title('Refractive index in the visible')
#    plt.show()
#    return 1
#
#def plot_index_ir(self,layer):
#    lam = self.lambda_array
#    RI = self.layer_ri(layer)
#    lam_min = np.min(lam)
#    lam_max = np.max(lam)
#    um = 1e+6  
#    #  Plot the real part of the refractive index in the Infrared regime   
#    mask = (lam<=lam_max) & (lam >= np.max([1000e-9, lam_min]))
#    plt.fig, ax1 = plt.subplots()
#    ax1.plot(lam[mask]*um, RI[mask].real, 'k-')
#    ax1.autoscale()
#    ax1.set_ylabel('n', color = 'k')  
#    ax1.set_xlabel('Wavelength (nm)')    
#    ax2 = ax1.twinx()
#    ax2.plot(lam[mask]*um, RI[mask].imag, 'r:')
#    ax2.autoscale()
#    ax2.set_ylabel('k', color = 'r')        
#    plt.title('Refractive index in the IR')
#    plt.show()
#    return 1
#    
#    
#    def plot_index_alloy(self, mat1, mat2, RI_alloy):
#        lam = self.lambda_array
#       # RI_alloy = self.layer_alloy(layer, fraction, mat1, mat2, model)
#        lam_min = np.min(lam)
#        lam_max = np.max(lam)
#        
#        if (lam_min < 1000e-9):
#            nm = 1e+9  
#            #  Plot the real part of the refractive index in the visible regime  
#            mask = (lam>=lam_min) & (lam <= np.min([1000e-9, lam_max]))
#            plt.fig, ax1 = plt.subplots()
#            ax1.plot(lam[mask]*nm, mat1[mask].real, 'k--')
#            ax1.plot(lam[mask]*nm, mat2[mask].real, 'k-.')
#            ax1.plot(lam[mask]*nm, RI_alloy[mask].real, 'k-')
#            ax1.autoscale()
#            ax1.set_ylabel('n', color = 'k')  
#            ax1.set_xlabel('Wavelength (nm)')    
#            ax2 = ax1.twinx()
#            ax2.plot(lam[mask]*nm, mat1[mask].imag, 'r--')
#            ax2.plot(lam[mask]*nm, mat2[mask].imag, 'r-.')            
#            ax2.plot(lam[mask]*nm, RI_alloy[mask].imag, 'r-')            
#            ax2.autoscale()
#            ax2.set_ylabel('k', color = 'r')        
#            ax1.legend(['Hoast', 'Inclusion', 'Effective'],bbox_to_anchor=(1.04,1))
#            plt.tight_layout(rect=[0,0,0.99,1])
#            plt.title('Refractive index in the visible')
#            plt.show()
#
#        if (lam_max > 1000e-9):
#            um = 1e+6  
#            #  Plot the real part of the refractive index in the Infrared regime   
#            mask = (lam<=lam_max) & (lam >= np.max([1000e-9, lam_min]))
#            plt.fig, ax1 = plt.subplots()
#            ax1.plot(lam[mask]*um, mat1[mask].real, 'k--')
#            ax1.plot(lam[mask]*um, mat2[mask].real, 'k-.')
#            ax1.plot(lam[mask]*um, RI_alloy[mask].real, 'k-')
#            ax1.autoscale()
#            ax1.set_ylabel('n', color = 'k')  
#            ax1.set_xlabel('Wavelength (nm)')    
#            ax2 = ax1.twinx()
#            ax2.plot(lam[mask]*um, mat1[mask].imag, 'r--')
#            ax2.plot(lam[mask]*um, mat2[mask].imag, 'r-.')            
#            ax2.plot(lam[mask]*um, RI_alloy[mask].imag, 'r-')            
#            ax2.autoscale()
#            ax2.set_ylabel('k', color = 'r')        
#            ax1.legend(['Hoast', 'Inclusion', 'Effective'],bbox_to_anchor=(1.04,1))
#            plt.tight_layout(rect=[0,0,0.99,1])            
#            plt.title('Refractive index in the IR')
#            plt.show()
#        return 1
           
