# -*- coding: utf-8 -*-
"""                 
Created on 10/21/2020
                    
@author: Jun@UMD             
"""
import matplotlib as mpl
#mpl.use('Agg')              


def project2grid_single(input):
    import math,numpy,geopy,os,pickle
    import calculate_BGcoef.project2grid_fortranlib.ant_prj

    [fname_latlon,fname_grid,bw,lat,lon,satzen,satazi,satran,
     delx,dely,grid_origin_lat,grid_origin_lon,numx,numy,xfovsiz,yfovsiz,
     select_enlarge,print_latlon]=input
     
    radiom_half_beamwidth_rad = bw/2.0
    theta_prime_max           = select_enlarge*radiom_half_beamwidth_rad
    delta_theta_prime_max     = theta_prime_max/200.0
    delta_phi_prime_max       = delta_theta_prime_max/math.sin(theta_prime_max)

    delta_theta_prime_equiv=(delx/2.0)*1000*math.cos(satzen)/satran
    delta_phi_prime_equiv=delta_theta_prime_equiv/math.sin(theta_prime_max)

    dtheta=min(delta_theta_prime_max,delta_theta_prime_equiv)
    dphi=min(delta_phi_prime_max,delta_phi_prime_equiv)
    ntheta=math.floor(theta_prime_max/dtheta)+2
    nphi=math.floor(2*numpy.pi/dphi)+2
    dtheta=theta_prime_max/float(ntheta-1)
    dphi=2*numpy.pi/float(nphi-1)

    theta=numpy.zeros(ntheta)
    phi=numpy.zeros(nphi)
    ant_response=numpy.zeros(ntheta)
    wgt=numpy.zeros((ntheta,nphi))
    sigma = radiom_half_beamwidth_rad/numpy.sqrt(2.0*math.log(2.0))
    sum_wgt=0.0
    for itheta in range(ntheta):
        theta[itheta] = float(itheta)*dtheta
        ant_response[itheta] = math.exp(-theta[itheta]*theta[itheta]/ \
                               (2.0*sigma*sigma))
        for iphi in range(nphi-1): 
            phi[iphi] = float(iphi)*dphi
            if phi[iphi]>numpy.pi: phi[iphi]=phi[iphi]-2.0*numpy.pi
            wgt[itheta, iphi]   = ant_response[itheta]
            sum_wgt             = sum_wgt+ant_response[itheta]*math.sin(theta[itheta])* \
                                  dtheta*dphi
    calculate_BGcoef.project2grid_fortranlib.ant_prj.ant_prj_mod.ant_prj(
                            ntheta,nphi,dtheta,dphi, 
                            lat,                   
                            lon,                  
                            satzen,
                            satazi,                
                            satran)         

    lat_sub=calculate_BGcoef.project2grid_fortranlib.ant_prj.ant_prj_mod.lat 
    lon_sub=calculate_BGcoef.project2grid_fortranlib.ant_prj.ant_prj_mod.lon 

    dx=geopy.distance.geodesic((grid_origin_lat,numpy.degrees(lon)),
                               (grid_origin_lat,grid_origin_lon)).km
    dy=geopy.distance.geodesic((numpy.degrees(lat),grid_origin_lon),
                               (grid_origin_lat,grid_origin_lon)).km
    nx_foot=round(dx/delx) 
    ny_foot=round(dy/dely) 

    temp_ant_wgt=numpy.zeros((numx,numy))
    temp_ant_wgt_ave=numpy.zeros((numx,numy))
    temp_ant_wgt_num=numpy.zeros((numx,numy))
    ant_pos=numpy.zeros((ntheta,nphi,2),dtype=int)
    for itheta in range(ntheta):
        for iphi in range(nphi-1):
            dx=geopy.distance.geodesic((grid_origin_lat,numpy.degrees(lon_sub[itheta,iphi])),
                                       (grid_origin_lat,grid_origin_lon)).km
            dy=geopy.distance.geodesic((numpy.degrees(lat_sub[itheta,iphi]),grid_origin_lon),
                                       (grid_origin_lat,grid_origin_lon)).km
            idx=round(dx/delx) 
            idy=round(dy/dely) 

            ant_pos[itheta,iphi,0]=idx
            ant_pos[itheta,iphi,1]=idy
            temp_ant_wgt[idx,idy]+=wgt[itheta,iphi]*math.sin(theta[itheta])*dtheta*dphi
            temp_ant_wgt_ave[idx,idy]+=wgt[itheta,iphi]
            temp_ant_wgt_num[idx,idy]+=1
    for i in range(numx):
        for j in range(numy):
            if temp_ant_wgt_num[i,j]>0: temp_ant_wgt_ave[i,j]=temp_ant_wgt_ave[i,j]/float(temp_ant_wgt_num[i,j])


    max_nx = ant_pos[:,:,0][ant_pos[:,:,0]>0].max()
    min_nx = ant_pos[:,:,0][ant_pos[:,:,0]>0].min()
    max_ny = ant_pos[:,:,1][ant_pos[:,:,1]>0].max()
    min_ny = ant_pos[:,:,1][ant_pos[:,:,1]>0].min()

    ant_wgt=numpy.zeros((xfovsiz,yfovsiz))
    ant_wgt_ave=numpy.zeros((xfovsiz,yfovsiz))
    for itheta in range(ntheta):
        for iphi in range(nphi-1):
            nx   = ant_pos[itheta,iphi,0]
            ny   = ant_pos[itheta,iphi,1]
            idim = nx - min_nx   
            jdim = ny - min_ny

            ant_wgt[idim,jdim] = temp_ant_wgt[nx,ny]
            ant_wgt_ave[idim,jdim] = temp_ant_wgt_ave[nx,ny]

    ant_wgt = ant_wgt/sum_wgt


    f=open(fname_grid,'wb')
    ant_pos=numpy.array([min_nx,max_nx,min_ny,max_ny,nx_foot,ny_foot])
    ant_pos.astype('>i4').tofile(f)
    ant_wgt.T.astype('>f8').tofile(f)
    ant_wgt_ave.T.astype('>f8').tofile(f)
    f.close()

    if print_latlon==1:
        with open(fname_latlon,'wb') as f:
            pickle.dump([ntheta,nphi,theta,phi,wgt,lat_sub,lon_sub],f)


