# -*- coding: utf-8 -*-
"""
Created on 9/21/2020

@author: Jun Zhou @UMD
"""
import matplotlib as mpl
#mpl.use('Agg')
from pathlib import Path
import os,glob,time
import pickle,sys
import h5py
import numpy as np
import statistics
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

class apply_BG_class:
    def __init__(self,case,datasource,input_dir,
                 local,orb,orbit_num,
                 src_ch,coef_dir,output_case_dir,nfov,
                 ta_vmin,ta_vmax,
                 dif_vmin,dif_vmax):
       
        print('Apply BG to case '+case+' '+datasource+' ...')
        
        self.datasource=datasource
        self.input_dir=input_dir
        [self.minlat,self.maxlat,self.minlon,self.maxlon]=[local['latmin'],local['latmax'],local['lonmin'],local['lonmax']]
        self.orb=orb,
        self.orbit_num=orbit_num
        self.select_chs=[src_ch]
        self.select_nch=len(self.select_chs)
        self.coef_dir=coef_dir
        self.output_case_dir=output_case_dir
        self.nfov=nfov

        self.select_bands=[]
        for ich in self.select_chs:
            if ich in [1]:
                iband=1
            elif ich in [2]:
                iband=2
            elif ich in np.arange(3,16):
                iband=3
            elif ich in [16]:
                iband=4
            elif ich in np.arange(17,23):
                iband=5
            self.select_bands.append(iband)
            
        self.ta_vmin=ta_vmin
        self.ta_vmax=ta_vmax
        self.dif_vmin=dif_vmin
        self.dif_vmax=dif_vmax

        self.draw_tgt,self.draw_src,self.draw_rmp,self.draw_dif=True,True,True,True


    def ingest(self):
   
        if self.datasource=='sdr':
            

            self.do_ingest()

  
        elif self.datasource=='simu':

            with open(os.path.join(self.input_dir,'Sim_Ta_src_tgt.pickle'),'rb') as f:
                [nscan,Ta_src,Ta_tgt,lat,lon,landfrac_src,landfrac_tgt]=pickle.load(f)

            self.nscan=nscan
            self.ta=np.zeros((self.select_nch,self.nscan,self.nfov)) 
            self.lat=np.zeros((self.select_nch,self.nscan,self.nfov)) 
            self.lon=np.zeros((self.select_nch,self.nscan,self.nfov)) 
            self.ta[0,:,:]=Ta_src
            self.lat[0,:,:]=lat
            self.lon[0,:,:]=lon
            self.ta_tgt=Ta_tgt  # lat lon are the same as src

            
    def prepcoef(self):

        self.alpha_all=[]
        self.windowsize_all=[]
        self.windowindx_all=[]
        self.idx_src_all=[]

        for ich in [0]: 
            self.alpha_all.append([])
            self.windowsize_all.append([])
            self.windowindx_all.append([])
            self.idx_src_all.append([])

            for ifv in range(self.nfov):

                with open(self.coef_dir+'/tgt_ifr'+str(ifv)+'.pickle','rb') as f:
                    [windowsize,windowindx,idx_src,alpha_bst]=pickle.load(f)


                self.alpha_all[ich].append(alpha_bst)               
                self.windowsize_all[ich].append(windowsize)               
                self.windowindx_all[ich].append(windowindx)               
                self.idx_src_all[ich].append(idx_src)               


 
    def apply(self):
        for ich in [0]:
            ta_rmp=np.full_like(self.ta[ich],np.nan)

            for ifr in range(self.nfov):
                windowsize=self.windowsize_all[ich][ifr]
                windowindx=self.windowindx_all[ich][ifr]
                idx_src=self.idx_src_all[ich][ifr]
                alpha=self.alpha_all[ich][ifr]
                isc_cen=windowindx[idx_src-1][0]
                ifr_cen=windowindx[idx_src-1][1]
                      
                for i in range(windowsize):
                    windowindx[i][0]=windowindx[i][0]-isc_cen
                    windowindx[i][1]=windowindx[i][1]-ifr_cen

                alpha=alpha/np.sum(alpha)
                        
                for isc in range(self.nscan):

                    temp=0.0
                    flag=0
                    for iwin in range(windowsize):
                        [i,j]=np.array(windowindx[iwin])+np.array([isc,ifr])
                        if i>=0 and i<self.nscan and j>=0 and j < self.nfov:
                            temp+=alpha[iwin]*self.ta[ich,i,j]
                        else:
                            flag=1
                            break
                    if flag==0: 
                        ta_rmp[isc,ifr]=temp
         

            save_file=self.output_case_dir+'ta_rmp.pickle'
            with open(save_file,'wb') as f:
                if self.datasource == 'sdr':
                    pickle.dump([self.lat[ich,:,:],self.lon[ich,:,:],self.ta[ich,:,:],ta_rmp],f)
                elif self.datasource == 'simu':
                    pickle.dump([self.lat[ich,:,:],self.lon[ich,:,:],self.ta[ich,:,:],ta_rmp,self.ta_tgt],f)



    def display(self):
        with open(self.output_case_dir+'ta_rmp.pickle','rb') as f:
            if self.datasource == 'sdr':
                [lat0,lon0,ta_src0,ta_rmp0]=pickle.load(f)
            elif self.datasource == 'simu':
                [lat0,lon0,ta_src0,ta_rmp0,ta_tgt0]=pickle.load(f)

        nsc=np.array(lat0).shape[0]        
        mask=np.zeros_like(lat0)
        
        for isc in range(nsc):
            if np.min(lat0[isc,1:94])>=self.minlat and np.max(lat0[isc,1:94])<=self.maxlat:
                mask[isc,1:94]=1
                
        idd=np.where(mask==1)
        lat=lat0[idd]
        lon=lon0[idd]
        ta_src=ta_src0[idd]
        ta_rmp=ta_rmp0[idd]

        [self.minlat_img,self.maxlat_img,self.minlon_img,self.maxlon_img]=[np.min(lat),np.max(lat),
                                                                           np.min(lon),np.max(lon)]

        if self.draw_src:
            print('draw Ta src ...')
            self.display_img(ta_src,
                             lat,
                             lon,
                             self.minlon_img,self.maxlon_img,self.minlat_img,self.maxlat_img,
                             self.output_case_dir+'TaOrg.png',
                             self.ta_vmin,self.ta_vmax)

        if self.draw_rmp:
            print('draw Ta rmp ...')
            save_fig=self.output_case_dir+'/TaRmp.png'
            self.display_img(ta_rmp,
                             lat,
                             lon,
                             self.minlon_img,self.maxlon_img,self.minlat_img,self.maxlat_img,
                             save_fig,self.ta_vmin,self.ta_vmax)

        if self.draw_dif:
            print('draw rmp minus src ...')
            save_fig=self.output_case_dir+'/Diff_Rmp_minus_Src.png'
            self.display_img(ta_rmp-ta_src,
                             lat,
                             lon,
                             self.minlon_img,self.maxlon_img,self.minlat_img,self.maxlat_img,
                             save_fig,vmin=self.dif_vmin,vmax=self.dif_vmax)

        if self.datasource == 'simu':

            ta_tgt=ta_tgt0[idd]

            if self.draw_tgt:
                print('draw Ta tgt ...')
                self.display_img(ta_tgt,
                                 lat,
                                 lon,
                                 self.minlon_img,self.maxlon_img,self.minlat_img,self.maxlat_img,
                                 self.output_case_dir+'TaTgt.png',
                                 self.ta_vmin,self.ta_vmax)

            if self.draw_dif:
                print('draw rmp minus tgt ...')
                save_fig=self.output_case_dir+'/Diff_Rmp_minus_Tgt.png'
                self.display_sct(ta_rmp-ta_tgt,
                                 lat,
                                 lon,
                                 self.minlon_img,self.maxlon_img,self.minlat_img,self.maxlat_img,
                                 save_fig,vmin=self.dif_vmin,vmax=self.dif_vmax)

            if self.draw_dif:
                print('draw src minus tgt ...')
                save_fig=self.output_case_dir+'/Diff_Src_minus_Tgt.png'
                self.display_img(ta_Src-ta_tgt,
                                 lat,
                                 lon,
                                 self.minlon_img,self.maxlon_img,self.minlat_img,self.maxlat_img,
                                 save_fig,vmin=self.dif_vmin,vmax=self.dif_vmax)


                 
    def do_ingest(self):
        self.create_filelist() 

        self.ta=[]
        self.lat=[]
        self.lon=[]
        self.satzen=[]
        self.satazi=[]
        self.satran=[]
        self.nscan=0
        
        for ele in self.select_chs:
            self.ta.append([])
            self.lat.append([])
            self.lon.append([])
            self.satzen.append([])
            self.satazi.append([])
            self.satran.append([])
       
        for ifile,filename in enumerate(self.sdr_filelist):
            [ta,lat,lon,satzen,satazi,satran,tim,vel,nscan]=self.read_granule(filename)
            for isc in range(nscan):
                if (vel[isc,-1]>0 and self.orb=='des') or (vel[isc,-1]<0 and self.orb=='asc'): continue
                self.nscan+=1
                for i,(ich,iband) in enumerate(zip(self.select_chs,self.select_bands)):
                    self.ta[i].append(ta[isc,:,ich-1])
                    self.lat[i].append(lat[isc,:,iband-1])
                    self.lon[i].append(lon[isc,:,iband-1])
                    self.satzen[i].append(satzen[isc,:])
                    self.satazi[i].append(satazi[isc,:])
                    self.satran[i].append(satran[isc,:])
        self.ta=np.array(self.ta)
        self.lat=np.array(self.lat)
        self.lon=np.array(self.lon)
        self.satzen=np.array(self.satzen)
        self.satazi=np.array(self.satazi)
        self.satran=np.array(self.satran)

    def read_granule(self,filename):
        f=h5py.File(self.input_dir+filename,'r')
        tem=np.array(f['/All_Data/ATMS-SDR_All/BrightnessTemperature']) 
        fac=np.array(f['/All_Data/ATMS-SDR_All/BrightnessTemperatureFactors']) 
        tim=np.array(f['/All_Data/ATMS-SDR_All/BeamTime']) 
        lat=np.array(f['/All_Data/ATMS-SDR-GEO_All/BeamLatitude']) 
        lon=np.array(f['/All_Data/ATMS-SDR-GEO_All/BeamLongitude']) 
        vel=np.array(f['/All_Data/ATMS-SDR-GEO_All/SCVelocity']) 
        satzen=np.array(f['/All_Data/ATMS-SDR-GEO_All/SatelliteZenithAngle']) 
        satazi=np.array(f['/All_Data/ATMS-SDR-GEO_All/SatelliteAzimuthAngle']) 
        satran=np.array(f['/All_Data/ATMS-SDR-GEO_All/SatelliteRange']) 
        f.close()
        nscan=tem.shape[0]
        ta=np.zeros_like(tem,dtype=float)
        for ich in range(22):
            for isc in range(nscan):
                for ifv in range(self.nfov):
                    ta[isc,ifv,ich]=float(tem[isc,ifv,ich])*fac[0]+fac[1]
        return [ta,lat,lon,satzen,satazi,satran,tim,vel,nscan]


    def create_filelist(self):
        if self.orbit_num=='non':
            fname_prefix=self.input_dir+'/GATMO-SATMS_*.h5'
        else:
            fname_prefix=self.input_dir+'/GATMO-SATMS_*_b'+str(self.orbit_num)+'*.h5'
        flist=glob.glob(fname_prefix)
        flist.sort(key=lambda x:x.split('_')[4])
        idx_str=os.path.basename(flist[0]).find('_d')
        uniq_key=[]
        uniq_filename=[]
        for fname in flist:
            filename=os.path.basename(fname)
            key=filename[0:idx_str+35]

            if key not in uniq_key:
                uniq_key.append(key)
                uniq_filename.append(filename)
            else:
                idx=uniq_key.index(key)
                uniq_filename[idx]=filename

        self.sdr_filelist=uniq_filename
        self.nfile=len(self.sdr_filelist)
        

    def display_sct(self,ta,lat,lon,lonmin,lonmax,latmin,latmax,save_file,vmin=-999,vmax=-999):
        font_jun={'size':20}
        mpl.rc('font',**font_jun)
        fig,ax=plt.subplots(1,1,figsize=(10, 10))
        m=Basemap(projection='cyl',resolution='c',
                  llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax)
        m.drawmapboundary(fill_color='lightblue')
        m.fillcontinents(color='wheat',lake_color='lightblue')
        parallels=np.arange(latmin,latmax+15,15)
        meridians=np.arange(lonmin,lonmax+30,30)
        m.drawparallels(parallels,
                        labels=[True,False,False,True]) 
        m.drawmeridians(meridians,
                        labels=[True,False,False,True]) 
        cm=plt.cm.get_cmap('jet')
        x=lon.ravel()
        y=lat.ravel()
        z=ta.ravel()
        idx=~np.isnan(z)
        x=x[idx]
        y=y[idx]
        z=z[idx]

        if vmin==-999: vmin=np.min(z)
        if vmax==-999: vmax=np.max(z)

        im=m.scatter(x,y,c=z,cmap=cm,marker='.',s=20,alpha=1,zorder=10,vmin=vmin,vmax=vmax)
        m.drawcoastlines(zorder=20)
        divider=make_axes_locatable(ax)
        cax=divider.append_axes('right',     
                                 size='5%', 
                                 pad=0.2)  
        cb=fig.colorbar(im,cax=cax)
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        plt.close()


    def display_img(self,ta,lat,lon,lonmin,lonmax,latmin,latmax,save_file,vmin=-999,vmax=-999):
        dd=0.1
        xg=np.arange(lonmin,lonmax+dd,dd)
        yg=np.arange(latmin,latmax+dd,dd)
        grid_x,grid_y=np.meshgrid(xg,yg)

        x=lon.ravel()
        y=lat.ravel()
        z=ta.ravel()
        idx=~np.isnan(z)
        x=x[idx]
        y=y[idx]
        z=z[idx]
        points=np.array([x,y]).T
        grid_z=griddata(points,z,(grid_x,grid_y),method='linear')

        if vmin==-999: vmin=np.min(grid_z)
        if vmax==-999: vmax=np.max(grid_z)


        font_jun={'size':15}
        mpl.rc('font',**font_jun)
        fig,ax=plt.subplots(1,1,figsize=(10, 10))
        m=Basemap(projection='cyl',resolution='f',
                  llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax)
        m.drawmapboundary(fill_color='lightblue')
        m.fillcontinents(color='wheat',lake_color='lightblue')
        parallels=np.arange(latmin,latmax+2,2)
        meridians=np.arange(lonmin,lonmax+5,5)
        m.drawparallels(parallels,
                        labels=[True,False,False,True]) 
        m.drawmeridians(meridians,
                        labels=[True,False,False,True]) 
        cm=plt.cm.get_cmap('rainbow')
        
        im=ax.imshow(grid_z,cmap=cm,origin='lower',
                     extent=[lonmin,lonmax,latmin,latmax],
                     vmin=vmin,vmax=vmax,zorder=10)
        m.drawcoastlines(zorder=20)
        divider=make_axes_locatable(ax)
        cax=divider.append_axes('right',    
                                 size='5%',   
                                 pad=0.2)    
        cb=fig.colorbar(im,cax=cax)

        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        plt.close()


