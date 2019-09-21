# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 11:33:57 2019

@author: parkerwray
"""

best_cp_sio2_sin = np.zeros((1,5))+(-1000)
best_cp_sio2_sin_np = np.zeros((1,5))+(-1000)
best_cp_sio2_np_sin = np.zeros((1,5))+(-1000)
best_cp_sio2_np_sin_np = np.zeros((1,5))+(-1000)

for idx_T in range(0,len(T)):  # Change structure temperature
    for idx_HT in range(0,len(H)):  # Change the top layer thickness
        for idx_HB in range(0,len(H)):   # Change the bottom layer thickness
            for idx_FFT in range(0,len(FF)):  # Change the top layer fill fraction
                for idx_FFB in range(0,len(FF)):  # Change the bottom layer fill fraction   
                    
                    if (idx_ff_B == len(FF)-1 and idx_ff_T == len(FF)-1):
                        best_cp = best_cp_sio2_sin
                        
                        
                        compile_results(
                                P_cool[idx_T][idx_HT][idx_HB][idx_FFT][idx_FFB],
                                )
                        
                        
                        
                    if (idx_ff_B == len(FF)-1 and idx_ff_T != len(FF)-1):
                    # Update data for top is NP and bottom is thin film cases

                    
                    if (idx_ff_B != len(FF)-1 and idx_ff_T == len(FF)-1):
                    # Update data for bottom is NP and top is thin film cases
                    
                    
                    
                    if (idx_ff_B != len(FF)-1 and idx_ff_T != len(FF)-1):
                    # Update data for top and bottom both being NP films case
                    