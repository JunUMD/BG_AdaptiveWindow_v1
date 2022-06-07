# -*- coding: utf-8 -*-
"""                 
Created on 10/29/2020
                    
@author: Jun@UMD             
"""
import matplotlib as mpl
#mpl.use('Agg')      

import calculate_BGcoef.calculateUVS_fortranlib.calculate_uvs

def calculateUVS_single(input1):
    [save_file,tgt_grd_dir,halfwin,atm_grd_dir,window_fname,ifr,
     numxgrd,numygrd] = input1

    calculate_BGcoef.calculateUVS_fortranlib. \
    calculate_uvs.calculateuvs_single_mod. \
    calculateuvs_single(save_file,tgt_grd_dir,halfwin,atm_grd_dir,window_fname,ifr,
                        numxgrd,numygrd)
    
