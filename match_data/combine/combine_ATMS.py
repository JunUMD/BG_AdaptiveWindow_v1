# -*- coding: utf-8 -*-
"""
Created on 10/09/2020

@author: Jun@UMD
"""
import matplotlib as mpl
#mpl.use('Agg')
import os
import pickle
from pathlib import Path
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from geopy.distance import geodesic

class combine_ATMS_ARTS_class:
    def __init__(self):
        self.isc_atm=0

class combine_ATMS_class:
    def __init__(self,output_dir,isc_atm,winsize,
                 ingest_ATMS_obj,
                 project_ATMS_obj):
        print('Combine ...')
        self.isc_atm=isc_atm
        self.winsize=winsize
        self.halfwin=int((self.winsize-1)/2)
        self.output_dir=output_dir

        self.atm_chs=ingest_ATMS_obj.select_chs
        self.atm_nch=ingest_ATMS_obj.select_nch
        self.atm_nscan=ingest_ATMS_obj.nscan
        self.atm_nfov=ingest_ATMS_obj.nfov
        self.atm_lat=ingest_ATMS_obj.lat  #[nch,nsc,nfv] degrees
        self.atm_lon=ingest_ATMS_obj.lon  #[nch,nsc,nfv]
        self.atm_tim=ingest_ATMS_obj.tim  #[nsc,nfv]  microsecond since 1/1/1958 IET

        self.atm_fovedge_lat=np.degrees(project_ATMS_obj.fovedge_lat) #[nch,nsc,nfv,nphi]
        self.atm_fovedge_lon=np.degrees(project_ATMS_obj.fovedge_lon) #[nch,nsc,nfv,nphi]

         

    def comebine_display(self):
        fname='combine.png'
        self.display_latlon(fname,
                            self.atm_fovedge_lon[0,self.isc_atm-self.halfwin:self.isc_atm+self.halfwin+1,:,:],
                            self.atm_fovedge_lat[0,self.isc_atm-self.halfwin:self.isc_atm+self.halfwin+1,:,:],
                            self.atm_fovedge_lon[1,self.isc_atm,:,:],
                            self.atm_fovedge_lat[1,self.isc_atm,:,:],
                            self.atm_lon[0,self.isc_atm-self.halfwin:self.isc_atm+self.halfwin+1,:],
                            self.atm_lat[0,self.isc_atm-self.halfwin:self.isc_atm+self.halfwin+1,:],
                            self.atm_lon[1,self.isc_atm,:],
                            self.atm_lat[1,self.isc_atm,:])
                                
                    
                         
    def display_latlon(self,fname,atmfovlon,atmfovlat,tgt_atmfovlon,tgt_atmfovlat,atmlon,atmlat,tgtlon,tgtlat):
        lonmin=np.min(atmfovlon)-0.5
        lonmax=np.max(atmfovlon)+0.5
        latmin=np.min(atmfovlat)-0.5
        latmax=np.max(atmfovlat)+0.5

        font_jun={'size':10}
        mpl.rc('font',**font_jun)
        fig,ax=plt.subplots(1,1,figsize=(10, 10))
        m=Basemap(projection='cyl',resolution='c',
                  llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax)
        m.drawmapboundary(fill_color='lightblue')
        m.fillcontinents(color='wheat',lake_color='lightblue')
        parallels=np.arange(latmin,latmax+5,5)
        meridians=np.arange(lonmin,lonmax+5,5)
        m.drawparallels(parallels,
                        labels=[True,False,False,True]) #,fmt=lat2str)
        m.drawmeridians(meridians,
                        labels=[True,False,False,True]) #,fmt=lon2str)    
        m.drawcoastlines(zorder=20)

        for isc in range(atmfovlon.shape[0]):
            for ifv in range(atmfovlon.shape[1]):
                m.plot(atmfovlon[isc,ifv,:],atmfovlat[isc,ifv,:],color='b',linewidth=1,alpha=1,zorder=10)
                m.scatter(atmlon[isc,ifv],atmlat[isc,ifv],marker='x',c='b',s=50,zorder=10)

        for ifv in range(tgt_atmfovlon.shape[0]):
            m.plot(tgt_atmfovlon[ifv,:],tgt_atmfovlat[ifv,:],color='r',linewidth=1,alpha=1,zorder=11)
            m.scatter(tgtlon[ifv],tgtlat[ifv],marker='x',c='r',s=20,zorder=11)


        plt.savefig(self.output_dir+fname,dpi=300,bbox_inches='tight')
        plt.close()

 
