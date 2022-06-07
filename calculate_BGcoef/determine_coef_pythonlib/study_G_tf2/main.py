# -*- coding: utf-8 -*-
"""
Created on 9/16/2020

@author: Jun@UMD
"""
import matplotlib as mpl
mpl.use('Agg')
from pathlib import Path
import os,pickle
import numpy as np
from display_G import display_G_class

#output_path='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_Ch1_5p2_Ch1_3p3/BGcoef/'
output_path='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_ARTSprj_GeoCrr_Ch1_5p2_Ch1_3p3/BGcoef/'
window_style='window_out0p01'
dBs=[-9] #np.arange(-11,0,1)

#window_style='window_square3'
#dBs=[-999]

#window_style='window_in0'
#dBs=[-999]

search_style='noiselv5_cycleSearch'
#search_style='noiselv999_singleSearch'

G_dir=os.path.join(output_path,
                   window_style,
                   'WindowAdjust_'+search_style,
                   'atmch0')+'/'


nfov=96
delx,dely=3,3

for ifr  in [0,1,2,3,4,5,6,24,47,71,89,90,91,92,93,94,95]: #[47]: 
     
    print('fov: ',ifr)
    G_fv_dir=os.path.join(G_dir,'fov'+str(ifr))+'/'
    G_fv_figdir=os.path.join(G_fv_dir,'fig')+'/'
    Path(G_fv_figdir).mkdir(parents=True,exist_ok=True)

    display_G_obj=display_G_class(ifr,G_fv_dir,G_fv_figdir,delx,dely)
    display_G_obj.display_G_0('Gsrc_wgtave')
    display_G_obj.display_G_0('Gtgt_wgtave')
    
    for dB in dBs:
        print('dB: ',dB)
        dB_str=str(dB).replace('-','neg')+'dB'
        display_G_obj.display_G_0('Gsim_wgtave_'+dB_str)




    

