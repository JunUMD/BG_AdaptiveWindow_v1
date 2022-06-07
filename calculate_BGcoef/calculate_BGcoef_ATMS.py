# -*- coding: utf-8 -*-
"""
Created on 9/21/2020

@author: Jun Zhou @UMD
"""
import matplotlib as mpl
#mpl.use('Agg')
import os
import numpy as np
import pickle
import sys
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from geopy.distance import geodesic
from math import ceil,floor,log,exp,sin,cos
import pp
from calculate_BGcoef.project2grid_single import project2grid_single
from calculate_BGcoef.determine_window_single import determine_window_single
from calculate_BGcoef.calculateUVS_single import calculateUVS_single
from calculate_BGcoef.determine_coef_single import determine_coef_single



class calculate_BGcoef_ATMS_class:
    def __init__(self,control,
                 output_dir,
                 delx,dely,
                 match_data_obj,
                 src_nedt,
                 window_dB,
                 noise_lev):
        print('Calculate BG coefficients ...')
        self.control=control
        self.src_nedt=src_nedt
        self.output_dir=output_dir
        self.delx=delx
        self.dely=dely
        self.xfovsiz=800
        self.yfovsiz=800
        self.winsize=match_data_obj.winsize
        self.halfwin=int((self.winsize-1)/2)
        isc_atm=match_data_obj.isc_atm
        bg=isc_atm-self.halfwin
        ed=isc_atm+self.halfwin+1        
        self.atm_chs=match_data_obj.select_chs
        self.atm_nch=match_data_obj.select_nch
        self.atm_nfov=match_data_obj.nfov
        self.atm_nscan=self.winsize
        self.atm_satzen=np.radians(match_data_obj.satzen[bg:ed,:])
        self.atm_satazi=np.radians(match_data_obj.satazi[bg:ed,:])
        self.atm_satran=match_data_obj.satran[bg:ed,:]
        self.atm_lat=np.radians(match_data_obj.lat[:,bg:ed,:])
        self.atm_lon=np.radians(match_data_obj.lon[:,bg:ed,:])
        self.atm_bw=match_data_obj.beamwidth 
        self.select_enlarge=match_data_obj.select_enlarge
        self.atm_fovedge_lat=np.degrees(match_data_obj.fovedge_lat[:,bg:ed,:,:])
        self.atm_fovedge_lon=np.degrees(match_data_obj.fovedge_lon[:,bg:ed,:,:])
        self.window_dB=window_dB
        self.noise_lev=noise_lev       

    def calculate_BGcoef(self):
        
        self.setupgrid() 

        self.project2grid_atm()        

        self.determine_window()
 
        self.calculate_UVS()

        self.determine_coef()


    def setupgrid(self):
        print('setupgrid ...')
        lonmin=min(np.min(self.atm_fovedge_lon),np.min(self.atm_fovedge_lon))-0.5
        lonmax=max(np.max(self.atm_fovedge_lon),np.max(self.atm_fovedge_lon))+0.5
        latmin=min(np.min(self.atm_fovedge_lat),np.min(self.atm_fovedge_lat))-0.5
        latmax=max(np.max(self.atm_fovedge_lat),np.max(self.atm_fovedge_lat))+0.5

        dx=geodesic((latmin,lonmin),(latmin,lonmax)).km
        dy=geodesic((latmin,lonmin),(latmax,lonmin)).km
        self.numx=ceil(dx/self.delx)+1
        self.numy=ceil(dy/self.dely)+1
        self.grid_origin_lat=latmin 
        self.grid_origin_lon=lonmin


    def project2grid_atm(self):
        print('Project to grid ...')

        for ich in [0,1]: 
            output_latlon_dir=os.path.join(self.output_dir,'project','atm_ch'+str(ich),'latlon')+'/'
            output_grid_dir=os.path.join(self.output_dir,'project','atm_ch'+str(ich),'grid')+'/'
            Path(output_latlon_dir).mkdir(parents=True,exist_ok=True)
            Path(output_grid_dir).mkdir(parents=True,exist_ok=True)
            if ich == 0:
                iscs=np.arange(self.atm_nscan)
            elif ich == 1:
                iscs=[self.halfwin] 
       
            for isc in iscs:
                tasks=[]
                inputs=[]
                for ifv in range(self.atm_nfov):
                    fname_latlon=output_latlon_dir+'isc'+str(isc)+'_ifv'+str(ifv)+'.pickle'
                    fname_grid=output_grid_dir+'isc'+str(isc)+'_ifv'+str(ifv)+'.bin'
                    if not (os.path.exists(fname_grid) and os.path.getsize(fname_grid) > 0 ):
                        tasks.append([ich,isc,ifv])
                        if ifv in [1,48,94]: 
                            print_latlon=1
                        else:
                            print_latlon=0
                        inputs.append([fname_latlon,fname_grid,
                                      self.atm_bw[ich],         # beamwidth radians 
                                      self.atm_lat[ich,isc,ifv],# fov center radians
                                      self.atm_lon[ich,isc,ifv],# fov center radians
                                      self.atm_satzen[isc,ifv], # fov center radians
                                      self.atm_satazi[isc,ifv], # fov center radians
                                      self.atm_satran[isc,ifv], # fov center radians
                                      self.delx,self.dely,      # 3km grid
                                      self.grid_origin_lat,self.grid_origin_lon, #lower left point of the grid
                                      self.numx,self.numy,
                                      self.xfovsiz,self.yfovsiz, #a fov occupy at most 800 grid
                                      self.select_enlarge[ich],
                                      print_latlon ])
                ntask=len(tasks)
                if ntask == 0: continue
                print('ch'+str(ich),'isc'+str(isc),'FOVs need to be processed: ',ntask)
                parallel=1
                if parallel == 1:
                    patch=24 
                    ncycle=floor(ntask/patch)
                    left=0
                    if ntask%patch != 0:
                        left=ntask-patch*ncycle
                        ncycle+=1
                    else:
                        left=patch
                    for icycle in range(ncycle):
                        if icycle == ncycle-1:
                            ibg=icycle*patch
                            ied=ibg+left
                        else:
                            ibg=icycle*patch
                            ied=(icycle+1)*patch
                        print(icycle,ibg,ied)                  

                        job_server = pp.Server(socket_timeout=86400)
                        jobs =  [(taskID,input1,job_server.submit(
                                 project2grid_single,   
                                 (input1,),     
                                 (),           
                                 ("numpy","pickle","math","geopy.distance",  
                                  "calculate_BGcoef.project2grid_fortranlib.ant_prj"))) 
                                 for taskID,input1 in enumerate(inputs[ibg:ied])]

                        job_server.wait()
                        job_server.print_stats()
                        job_server.destroy() 
                else:
                    for taskID,input1 in enumerate(inputs[0:1]):
                        project2grid_single(input1)

        #self.display_latlon_atm()
        #self.display_grid_atm()



    def determine_window(self):
        print('Determine window ...')

        self.output_win_dir=os.path.join(self.output_dir,
                                         'window_'+str(self.window_dB)+'dB')+'/'
        Path(self.output_win_dir).mkdir(parents=True,exist_ok=True)


        for ich in [0]: 

            save_dir=os.path.join(self.output_win_dir,'window','atmch'+str(ich))
            Path(save_dir).mkdir(parents=True,exist_ok=True)

            tasks=[]
            inputs=[]
            for ifr in range(self.atm_nfov):
                save_file=save_dir+'/tgt_ifr'+str(ifr)+'.bin'
                if not (os.path.exists(save_file) and os.path.getsize(save_file) > 0 ):             
                    tgt_grd_dir=os.path.join(self.output_dir,'project','atm_ch1','grid') 
                    atm_grd_dir=os.path.join(self.output_dir,'project','atm_ch'+str(ich),'grid')
                    tasks.append([ich,ifr])
                    inputs.append([save_file,tgt_grd_dir,self.halfwin,atm_grd_dir,ifr,self.xfovsiz,self.yfovsiz,
                                   self.atm_nscan,self.atm_nfov,self.window_dB])

            ntask=len(tasks)
            if ntask == 0: continue
            print('atm ch'+str(ich),'FORs need to be processed: ',ntask)
            parallel=1
            if parallel == 1:
                patch=24  
                ncycle=floor(ntask/patch)
                left=0
                if ntask%patch != 0:
                    left=ntask-patch*ncycle
                    ncycle+=1
                else:
                    left=patch
               
                for icycle in range(ncycle):
                    if icycle == ncycle-1:
                        ibg=icycle*patch
                        ied=ibg+left
                    else:
                        ibg=icycle*patch
                        ied=(icycle+1)*patch
                    print(icycle,ibg,ied) 


                    job_server = pp.Server(socket_timeout=86400)
                    jobs =  [(taskID,input1,job_server.submit(
                             determine_window_single, 
                             (input1,),        
                             (),              
                             ("numpy","pickle","math")))    
                             for taskID,input1 in enumerate(inputs[ibg:ied])]

                    job_server.wait()
                    job_server.print_stats()
                    job_server.destroy() 
            elif parallel==0:
                for taskID,input1 in enumerate(inputs[47:48]):
                    print(input1)
                    determine_window_single(input1)

        #self.display_window()



    def calculate_UVS(self):
        print('Calculate UVS ...')
        for ich in [0]: 

            save_dir=os.path.join(self.output_win_dir,'UVS','atmch'+str(ich))+'/'
            Path(save_dir).mkdir(parents=True,exist_ok=True)

            tasks=[]
            inputs=[]
            for ifr in range(self.atm_nfov):
                window_fname=os.path.join(self.output_win_dir,'window','atmch'+str(ich),
                                          'tgt_ifr'+str(ifr)+'.bin')
                save_file=save_dir+'/tgt_ifr'+str(ifr)+'.bin'
                if not (os.path.exists(save_file) and os.path.getsize(save_file) > 0 ):             
                    tgt_grd_dir=os.path.join(self.output_dir,'project','atm_ch1','grid')+'/'  
                    atm_grd_dir=os.path.join(self.output_dir,'project','atm_ch'+str(ich),'grid')+'/'

                    tasks.append([ich,ifr])
                    inputs.append([save_file,tgt_grd_dir,self.halfwin,atm_grd_dir,window_fname,ifr,
                                   self.xfovsiz,self.yfovsiz])

            ntask=len(tasks)
            if ntask == 0: continue
            print('ch'+str(ich),'ifr'+str(ifr),'FORs need to be processed: ',ntask)
  
            parallel=1
            if parallel == 1:
                patch=24 
                ncycle=floor(ntask/patch)
                left=0
                if ntask%patch != 0: 
                    left=ntask-patch*ncycle
                    ncycle+=1 
                else:
                    left=patch

                for icycle in range(ncycle):
                    if icycle == ncycle-1:
                        ibg=icycle*patch
                        ied=ibg+left
                    else:
                        ibg=icycle*patch
                        ied=(icycle+1)*patch
                    print(icycle,ibg,ied)

                    job_server = pp.Server(socket_timeout=86400)
                    jobs =  [(taskID,input1,job_server.submit(
                             calculateUVS_single,   
                             (input1,),       
                             (),          
                             ("calculate_BGcoef.calculateUVS_fortranlib.calculate_uvs",))) 
                             for taskID,input1 in enumerate(inputs[ibg:ied])]

                    job_server.wait()
                    job_server.print_stats()
                    job_server.destroy() 

            elif parallel==0:
                for taskID,input1 in enumerate(inputs[47:48]):
                    print(input1)
                    calculateUVS_single(input1)



    def determine_coef(self):
        print('Determine coefficients ...')
        for ich in [0]:  

            coef_save_dir=os.path.join(self.output_win_dir,
                                         'Coeff',
                                         'atmch'+str(ich))+'/'
            Path(coef_save_dir).mkdir(parents=True,exist_ok=True)

            win_dir=os.path.join(self.output_win_dir,'window','atmch'+str(ich))+'/'
            uvs_dir=os.path.join(self.output_win_dir,'UVS','atmch'+str(ich))+'/'

            tasks=[]
            inputs=[]
            for ifr  in range(self.atm_nfov):
                save_file=coef_save_dir+'/tgt_ifr'+str(ifr)+'.pickle'
                if not (os.path.exists(save_file) and os.path.getsize(save_file) > 0 ):             
                    tasks.append([ich,ifr])
                    win_file=win_dir+'/tgt_ifr'+str(ifr)+'.bin'
                    uvs_file=uvs_dir+'/tgt_ifr'+str(ifr)+'.bin'

                    inputs.append([self.src_nedt,ifr,self.xfovsiz,self.yfovsiz,
                                   win_file,uvs_file,self.noise_lev,save_file])

            ntask=len(tasks)
            if ntask == 0: continue
            print('atm ch'+str(ich),'FORs need to be processed: ',ntask)
            parallel=1
            if parallel == 1:
                patch=24
                ncycle=floor(ntask/patch)
                left=0
                if ntask%patch != 0:
                    left=ntask-patch*ncycle
                    ncycle+=1
                else:
                    left=patch

                for icycle in range(ncycle):
                    if icycle == ncycle-1:
                        ibg=icycle*patch
                        ied=ibg+left
                    else:
                        ibg=icycle*patch
                        ied=(icycle+1)*patch
                    print(icycle,ibg,ied)

                    job_server = pp.Server(socket_timeout=86400)
                    jobs =  [(taskID,input1,job_server.submit(
                             determine_coef_single,     
                             (input1,),      
                             (),             
                             ("calculate_BGcoef.determine_coef_pythonlib.determine_coef_class",))) 
                             for taskID,input1 in enumerate(inputs[ibg:ied])]

                    job_server.wait()
                    job_server.print_stats()
                    job_server.destroy()
            else:
                for taskID,input1 in enumerate(inputs[47:48]):
                    print(input1)
                    determine_coef_single(input1)
                
        self.display_coef()


               
    def display_latlon_atm(self):
        cm=plt.cm.get_cmap('rainbow')
        for ich in [0,1]:
            if ich==0:
                iscs=[0,int((self.atm_nscan-1)/2),self.atm_nscan-1]
            elif ich==1:
                iscs=[int((self.atm_nscan-1)/2)]
            for isc in iscs:
                minlon=np.min(self.atm_fovedge_lon[ich,isc,:,:])
                maxlon=np.max(self.atm_fovedge_lon[ich,isc,:,:])
                minlat=np.min(self.atm_fovedge_lat[ich,isc,:,:])
                maxlat=np.max(self.atm_fovedge_lat[ich,isc,:,:])

                Path(os.path.join(self.output_dir,'project','atm_ch'+str(ich),'latlon')).mkdir(parents=True,exist_ok=True) 
                gname=os.path.join(self.output_dir,'project','atm_ch'+str(ich),'latlon','isc'+str(isc)+'.png')
                 

                fig,ax=plt.subplots(1,1,figsize=(20, 10))
                m=Basemap(projection='cyl', 
                          llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,
                          resolution='c')
                m.drawmapboundary(fill_color='lightblue') 
                m.fillcontinents(color='wheat',lake_color='lightblue') 
                m.drawcoastlines()  
                parallels=np.arange(minlat,maxlat,2)   
                meridians=np.arange(minlon,maxlon,10)   
                m.drawparallels(parallels,
                                labels=[True,True,True,True])  
                m.drawmeridians(meridians,
                                labels=[True,True,True,True])
                for ifv in [1,48,94]:
                    fname=os.path.join(self.output_dir,'project','atm_ch'+str(ich),'latlon')+'/isc'+str(isc)+'_ifv'+str(ifv)+'.pickle'
                    with open(fname,'rb') as f:
                        [ntheta,nphi,theta,phi,wgt,lat_subgrid,lon_subgrid]=pickle.load(f) 

                    wgt[:,-1]=wgt[:,-2]
                    sc=m.scatter(np.degrees(lon_subgrid),np.degrees(lat_subgrid),c=wgt,
                                 cmap=cm,marker='.',s=0.1,alpha=1,zorder=10)
                    m.plot(self.atm_fovedge_lon[ich,isc,ifv,:],self.atm_fovedge_lat[ich,isc,ifv,:],color='k',zorder=11)
                    m.scatter(np.degrees(self.atm_lon[ich,isc,ifv]),
                              np.degrees(self.atm_lat[ich,isc,ifv]),color='k',marker='X',zorder=11)
                divider=make_axes_locatable(ax)
                cax=divider.append_axes('right',   
                                        size='5%',  
                                        pad=0.2)  
                cb=fig.colorbar(sc,cax=cax)
                plt.savefig(gname,dpi=300,bbox_inches='tight')
                plt.close()
                print(gname)
 

    def display_grid_atm(self):
        cm=plt.cm.get_cmap('rainbow')
        for ich in [0,1]:
            save_dir=os.path.join(self.output_dir,'project','atm_ch'+str(ich),'grid')
            Path(save_dir).mkdir(parents=True,exist_ok=True)

            if ich==0:
                iscs=[0,int((self.atm_nscan-1)/2),self.atm_nscan-1]
            elif ich==1:
                iscs=[int((self.atm_nscan-1)/2)]

            for isc in iscs:

                gname=os.path.join(save_dir,'isc'+str(isc)+'.png')
                fig,ax=plt.subplots(1,1,figsize=(20, 10))
                for ifv in [1,48,94]:
                    fname=os.path.join(self.output_dir,'project','atm_ch'+str(ich),'grid','isc'+str(isc)+'_ifv'+str(ifv)+'.bin')
                    f=open(fname,'rb')
                    [min_nx,max_nx,min_ny,max_ny,nx_foot,ny_foot]=np.fromfile(f,dtype='>i4',count=6)[:]
                    ant_wgt=np.fromfile(f,dtype='>f8',count=self.xfovsiz*self.yfovsiz).reshape(self.yfovsiz,self.xfovsiz).T
                    f.close()
                    x=[]
                    y=[]
                    z=[]
                    for idim in range(self.xfovsiz):
                        for jdim in range(self.yfovsiz):
                            if ant_wgt[idim,jdim]>0.0:
                                x.append(float(idim+min_nx)*self.delx)  
                                y.append(float(jdim+min_ny)*self.dely)
                                z.append(ant_wgt[idim,jdim])  
                    x_ct=float(nx_foot)*self.delx
                    y_ct=float(ny_foot)*self.dely
                    sc=ax.scatter(x,y,c=z,cmap=cm,marker='.',s=1,alpha=1,zorder=10)
                    #ax.scatter(x_ct,y_ct,color='k',marker='X',zorder=11)

                cb=fig.colorbar(sc)
                plt.savefig(gname,dpi=300,bbox_inches='tight')
                plt.close()
                print(gname)


    def display_window(self):
        cm=plt.cm.get_cmap('rainbow')
        for ich in [0]: 
            size=150
            for ifr in [0,47,95]: 

                f=open(os.path.join(self.output_win_dir,'window','atmch'+str(ich),'tgt_ifr'+str(ifr)+'.bin'),'rb')
                atm_windowsize=np.fromfile(f,dtype='>i4',count=1)[0]
                atm_windowindx=np.fromfile(f,dtype='>i4',count=atm_windowsize*2).reshape(2,atm_windowsize).T.tolist()
                idx_src=np.fromfile(f,dtype='>i2',count=1)[0]
                f.close()
                atm_flag=[]
                atm_lat_cen=[]
                atm_lon_cen=[]
                for [isc,ifv] in atm_windowindx:
                    atm_lat_cen.append(np.degrees(self.atm_lat[ich,isc,ifv]))
                    atm_lon_cen.append(np.degrees(self.atm_lon[ich,isc,ifv]))
                    if [isc,ifv]==atm_windowindx[idx_src-1]:
                        atm_flag.append(1)   # center
                    else:
                        atm_flag.append(0)   # others

                atm_lat_tgtedg=self.atm_fovedge_lat[1,self.halfwin,ifr,:]
                atm_lon_tgtedg=self.atm_fovedge_lon[1,self.halfwin,ifr,:]

                (minlonedg,maxlonedg,minlatedg,maxlatedg) = (np.min(atm_lon_tgtedg),np.max(atm_lon_tgtedg),
                                                             np.min(atm_lat_tgtedg),np.max(atm_lat_tgtedg))
                (minloncen,maxloncen,minlatcen,maxlatcen) = (np.min(atm_lon_cen),np.max(atm_lon_cen),
                                                             np.min(atm_lat_cen),np.max(atm_lat_cen))
                (minlon,maxlon,minlat,maxlat) = (min(minlonedg,minloncen)-0.5,max(maxlonedg,maxloncen)+0.5,
                                                 min(minlatedg,minlatcen)-0.5,max(maxlatedg,maxlatcen)+0.5)
                

                fig,ax=plt.subplots(1,1,figsize=(20, 20))
                m=Basemap(projection='cyl', 
                          llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,
                          resolution='c')
                m.drawmapboundary(fill_color='lightblue') 
                m.fillcontinents(color='wheat',lake_color='lightblue') 
                m.drawcoastlines()  
                parallels=np.arange(minlat,maxlat,2)  
                meridians=np.arange(minlon,maxlon,10)   
                m.drawparallels(parallels,
                                labels=[True,True,True,True])  
                m.drawmeridians(meridians,
                                labels=[True,True,True,True]) 

                sc_cir=m.plot(atm_lon_tgtedg,atm_lat_tgtedg,
                               color='b',linewidth=2,zorder=10)

                for flag,lat_cen,lon_cen in zip(atm_flag,atm_lat_cen,atm_lon_cen):
                    if flag==0:
                        m.scatter(lon_cen,lat_cen,c='k',marker='x',s=size,zorder=11)
                    elif flag == 1:
                        m.scatter(lon_cen,lat_cen,c='m',marker='X',s=250,zorder=11)

                gname=os.path.join(self.output_win_dir,'window','atmch'+str(ich),'tgt_ifr'+str(ifr)+'.png')
                plt.savefig(gname,dpi=300,bbox_inches='tight')
                plt.close()
                print(gname) 


    def display_coef(self):
        
        for ich in [0]: 
            for ifr in [1,47,94]: 
                path=os.path.join(self.output_win_dir,
                                  'Coeff',
                                  'atmch'+str(ich))+'/'
                coef_file=path+'tgt_ifr'+str(ifr)+'.pickle'
                
                with open(coef_file,'rb') as f:
                    [windowsize,windowindx,idx_src,alpha]=pickle.load(f)
                print(ifr,np.sum(alpha))

                idx=np.argmax(alpha)
                x=np.arange(len(alpha))-idx
                font_jun={'family':'Times New Roman','size':25}
                mpl.rc('font',**font_jun)
                fig=plt.figure(figsize=(10,5))
                ax=fig.add_subplot(111)
                ax.plot(x,alpha,'k-',linewidth=1)
                ax.scatter(x,alpha,c='k',s=20)
                fig.tight_layout()
                plt.savefig(path+'tgt_ifr'+str(ifr)+'.png',dpi=300,bbox_inches='tight')
                plt.close()   
