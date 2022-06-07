# -*- coding: utf-8 -*-
"""
BG remapping with adaptive window

"Improving ATMS Remapping Accuracy Using Adaptive Window and Noise-tuning Method in Backus-Gilbert Inversion"
Published on TGRS 2022

Created on 9/16/2020

@author: Jun Zhou @ UMD
"""
import matplotlib as mpl
mpl.use('Agg')
import os
import sys
import pickle
import numpy as np
from pathlib import Path

from match_data.match_ATMS import match_ATMS_class
from calculate_BGcoef.calculate_BGcoef_ATMS import calculate_BGcoef_ATMS_class
from apply_BG.apply_BG import apply_BG_class


main_path='/mars/jzhou/BG_AdaptiveWindow_v1/'
data_path=os.path.join(main_path,'data')+'/'

# NOAA-20 channel 1 remapping from 5.2 deg to 3.3 deg 
sensor='N20_ATMS'
src_ch=1
tgt_ch=1
src_nedt=0.222    
src_bm=5.2       
tgt_bm=3.3       
window_dB=-5
noise_lev=2.5

output_path=data_path+'output/Remap_'+sensor+\
                '_Ch'+str(src_ch)+'_'+str(src_bm).replace('.','p')+\
                '_Ch'+str(tgt_ch)+'_'+str(tgt_bm).replace('.','p')
Path(output_path).mkdir(parents=True,exist_ok=True)

input_path=os.path.join(data_path,'input',sensor)+'/'
Path(input_path).mkdir(parents=True,exist_ok=True)

control={'match_data':0,      
         'calculate_BGcoef':0,
         'apply_BG':1}


if control['match_data']>=0:

    winsize=41 
    isc_atm=26
    output_dir=os.path.join(output_path,'match_data')+'/'
    Path(output_dir).mkdir(parents=True,exist_ok=True)

    match_data_obj=match_ATMS_class(control['match_data'],input_path,output_dir,sensor,
                                    [src_ch,tgt_ch],[src_bm,tgt_bm],
                                    isc_atm,winsize)
    match_data_obj.match_data()


if control['calculate_BGcoef']>0:

    delx=3 
    dely=3
    output_dir=os.path.join(output_path,'BGcoef')+'/'
    Path(output_dir).mkdir(parents=True,exist_ok=True)

    calculate_BGcoef_obj=calculate_BGcoef_ATMS_class(control['calculate_BGcoef'],
                                                     output_dir,
                                                     delx,dely,
                                                     match_data_obj,
                                                     src_nedt,
                                                     window_dB,
                                                     noise_lev)
    calculate_BGcoef_obj.calculate_BGcoef()
   

if control['apply_BG']>0:

    datasource='simu' 
    ymd='20190831'
    case='Dorian'

    #datasource='sdr' 
    #ymd='20190831'
    #case='Dorian'

    #datasource='sdr'
    #ymd='20190831'
    #case='AmazonRiverDelta'

    #datasource='sdr'  
    #ymd='20210419'
    #case='Surigae'


    input_dir=os.path.join(input_path,datasource,ymd,case)+'/'

    if datasource=='sdr':

        if 'Dorian' in case:
            local={'latmin':20,
                   'latmax':31,
                   'lonmin':-95,
                   'lonmax':-50}
            orb='asc'
            orbit_num='non'
            ta_vmin,ta_vmax=210,260
            dif_vmin,dif_vmax=-4,4
        
        elif 'AmazonRiverDelta' in case:
            local={'latmin':-8,
                   'latmax':12,
                   'lonmin':-90,
                   'lonmax':-30}
            orb='des'
            orbit_num='non'
            ta_vmin,ta_vmax=210,300
            dif_vmin,dif_vmax=-6,6
            
        elif 'Surigae' in case:
            local={'latmin':4,
                   'latmax':21,
                   'lonmin':90,
                   'lonmax':150}
            orb='des'
            orbit_num='non'
            ta_vmin,ta_vmax=210,275
            dif_vmin,dif_vmax=-5,5

    elif datasource=='simu':

        if 'Dorian' in case:
            local={'latmin':20,
                   'latmax':31,
                   'lonmin':-95,
                   'lonmax':-50}
            orb='non'
            orbit_num='non'
            ta_vmin,ta_vmax=210,260
            dif_vmin,dif_vmax=-4,4

 
    coef_dir=os.path.join(output_path,'BGcoef','window_'+str(window_dB)+'dB','Coeff','atmch'+str(src_ch-1)+'/')

    output_case_dir=os.path.join(output_path,'apply_BG',case+'_'+ymd+'_'+datasource,'window_'+str(window_dB)+'dB/')
    Path(output_case_dir).mkdir(parents=True,exist_ok=True)
       
    apply_BG_obj=apply_BG_class(case,
                                datasource,
                                input_dir, 
                                local,orb,orbit_num,
                                src_ch,
                                coef_dir,
                                output_case_dir,
                                match_data_obj.nfov,
                                ta_vmin,ta_vmax,
                                dif_vmin,dif_vmax)
    apply_BG_obj.ingest()
    apply_BG_obj.prepcoef()
    apply_BG_obj.apply()
    apply_BG_obj.display()
 
