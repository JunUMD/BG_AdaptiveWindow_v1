# -*- coding: utf-8 -*-
"""                 
Created on 10/21/2020
                    
@author: Jun@UMD             
"""                          
import matplotlib as mpl     
#mpl.use('Agg')              
import os,pickle
import numpy

def determine_window_single(input1):

    [save_file,tgt_grd_dir,halfwin,atm_grd_dir,ifr,xfovsiz,yfovsiz,
     atm_nscan,atm_nfov,window_dB]=input1

    fname_grid=os.path.join(tgt_grd_dir,'isc'+str(halfwin)+'_ifv'+str(ifr)+'.bin')
    f=open(fname_grid,'rb')
    [min_nx_target,max_nx_target,min_ny_target,max_ny_target,
     nx_foot_target,ny_foot_target]=numpy.fromfile(f,dtype='>i4',count=6)[:]
    antwgt_target=numpy.fromfile(f,dtype='>f8',count=xfovsiz*yfovsiz).reshape(yfovsiz,xfovsiz).T
    antwgt_ave_target=numpy.fromfile(f,dtype='>f8',count=xfovsiz*yfovsiz).reshape(yfovsiz,xfovsiz).T
    antwgt_ave_target=antwgt_ave_target/numpy.max(antwgt_ave_target)
    f.close()
    antpos_target_nx=numpy.zeros((xfovsiz,yfovsiz))
    antpos_target_ny=numpy.zeros((xfovsiz,yfovsiz))
    for idim in range(xfovsiz):
        for jdim in range(yfovsiz):
            antpos_target_nx[idim,jdim] = idim + min_nx_target 
            antpos_target_ny[idim,jdim] = jdim + min_ny_target

    atm_windowindx=[]
    atm_windowsize=0

    WG=10**(window_dB/10)

    max_src=-999.0
    idx_src=-999   
         
    for isc in range(atm_nscan):
        #print(isc)
        for ifv in range(atm_nfov):
            fname_grid=os.path.join(atm_grd_dir,'isc'+str(isc)+'_ifv'+str(ifv)+'.bin')

            f=open(fname_grid,'rb')
            [min_nx_source,max_nx_source,min_ny_source,max_ny_source,
             nx_foot_source,ny_foot_source]=numpy.fromfile(f,dtype='>i4',count=6)[:]
            antwgt_source=numpy.fromfile(f,dtype='>f8',count=xfovsiz*yfovsiz).reshape(yfovsiz,xfovsiz).T
            antwgt_ave_source=numpy.fromfile(f,dtype='>f8',count=xfovsiz*yfovsiz).reshape(yfovsiz,xfovsiz).T
            antwgt_ave_source=antwgt_ave_source/numpy.max(antwgt_ave_source)
            f.close()

            Uo_src=0
            Uo_tgt=0 
            antpos_source_nx=numpy.zeros((xfovsiz,yfovsiz))
            antpos_source_ny=numpy.zeros((xfovsiz,yfovsiz))
            for idim in range(xfovsiz):
                for jdim in range(yfovsiz):
                    if (antwgt_ave_source[idim,jdim] <= 0): continue
                    antpos_source_nx[idim,jdim] = idim + min_nx_source 
                    antpos_source_ny[idim,jdim] = jdim + min_ny_source
                    if (antpos_source_nx[idim,jdim]<=max_nx_target and
                        antpos_source_nx[idim,jdim]>=min_nx_target and
                        antpos_source_ny[idim,jdim]<=max_ny_target and
                        antpos_source_ny[idim,jdim]>=min_ny_target):
                        idim_target=int(antpos_source_nx[idim,jdim]-min_nx_target)
                        jdim_target=int(antpos_source_ny[idim,jdim]-min_ny_target)
                        if(idim_target >= 0 and idim_target < xfovsiz and
                           jdim_target >= 0 and jdim_target < yfovsiz ):

                            if (antwgt_ave_target[idim_target,jdim_target]>0 and antwgt_ave_source[idim,jdim]>0):
                                if antwgt_ave_source[idim,jdim]>Uo_src: Uo_src=antwgt_ave_source[idim,jdim]
                                if antwgt_ave_target[idim_target,jdim_target]>Uo_tgt: Uo_tgt=antwgt_ave_target[idim_target,jdim_target]

            if Uo_src>WG and Uo_tgt>WG:
                atm_windowsize += 1
                atm_windowindx.append([isc,ifv])

                idim_target=int(nx_foot_source-min_nx_target) 
                jdim_target=int(ny_foot_source-min_ny_target)
                if(idim_target >= 0 and idim_target < xfovsiz and
                   jdim_target >= 0 and jdim_target < yfovsiz):  
                    if antwgt_ave_target[idim_target,jdim_target]>max_src:
                        max_src = antwgt_ave_target[idim_target,jdim_target]
                        idx_src = atm_windowsize 

      
                    
    f=open(save_file,'wb')
    numpy.array([atm_windowsize]).astype('>i4').tofile(f)
    numpy.array(atm_windowindx).T.astype('>i4').tofile(f)
    numpy.array([idx_src]).astype('>i2').tofile(f)
    f.close()
    
                
