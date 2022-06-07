# -*- coding: utf-8 -*-
"""
Created on 10/7/2020

@author: Jun@UMD
"""
import matplotlib as mpl
#mpl.use('Agg')
import os,sys
import glob
import numpy as np
import pickle
import h5py


class ingest_ATMS_class:
    def __init__(self,ymd,input_dir,output_dir,select_chs,orb,domain,select_orbnum):
        print('Ingest ...')
        self.geosdr_combine=1
        self.ymd=ymd
        self.input_dir=input_dir
        self.output_dir=output_dir
        self.select_chs=select_chs
        self.orb=orb
        self.select_orbnum=select_orbnum
        [self.latmin,self.latmax,self.lonmin,self.lonmax]=domain
        self.select_nch=len(select_chs)
        self.nfov=96
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


    def ingest(self):
        self.create_filelist()
        self.ta=[]
        self.lat=[]
        self.lon=[]
        self.satzen=[]
        self.satazi=[]
        self.satran=[]
        self.tim=[]
        self.nscan=0
 
        for ele in self.select_chs:
            self.ta.append([])
            self.lat.append([])
            self.lon.append([])
    
        idx_str=os.path.basename(self.sdr_filelist[0]).find('_b')+1
        for ifile,filename in enumerate(self.sdr_filelist):
            found=0
            if self.select_orbnum != '':
                if filename[idx_str:idx_str+6] != self.select_orbnum: continue
            [ta,lat,lon,satzen,satazi,satran,tim,vel,nscan]=self.read_granule(filename)
            for isc in range(nscan):
                if (vel[isc,-1]>0 and self.orb=='des') or (vel[isc,-1]<0 and self.orb=='asc'): continue
                if (lat[isc,48,0]>self.latmax or lat[isc,48,0]<self.latmin or
                    lon[isc,48,0]>self.lonmax or lon[isc,48,0]<self.lonmin): 
                    continue
                if found==0:
                    found=1
                    print(filename)
                self.nscan+=1
                self.satzen.append(satzen[isc,:])
                self.satazi.append(satazi[isc,:])
                self.satran.append(satran[isc,:])
                self.tim.append(tim[isc,:])
                for i,(ich,iband) in enumerate(zip(self.select_chs,self.select_bands)):
                    self.ta[i].append(ta[isc,:,ich-1])
                    self.lat[i].append(lat[isc,:,iband-1])
                    self.lon[i].append(lon[isc,:,iband-1])

        self.ta=np.array(self.ta) 
        self.lat=np.array(self.lat) 
        self.lon=np.array(self.lon) 
        self.satzen=np.array(self.satzen) 
        self.satazi=np.array(self.satazi) 
        self.satran=np.array(self.satran) 
        self.tim=np.array(self.tim)


    def read_granule(self,filename):
        f=h5py.File(self.input_dir+filename,'r')
        tem=np.array(f['/All_Data/ATMS-SDR_All/BrightnessTemperature']) #[12,96,22]
        fac=np.array(f['/All_Data/ATMS-SDR_All/BrightnessTemperatureFactors']) #[2]
        tim=np.array(f['/All_Data/ATMS-SDR_All/BeamTime']) #[12,96]
        lat=np.array(f['/All_Data/ATMS-SDR-GEO_All/BeamLatitude']) #[12,96,5]
        lon=np.array(f['/All_Data/ATMS-SDR-GEO_All/BeamLongitude']) #[12,96,5]
        vel=np.array(f['/All_Data/ATMS-SDR-GEO_All/SCVelocity']) #[12,3]
        satzen=np.array(f['/All_Data/ATMS-SDR-GEO_All/SatelliteZenithAngle']) #[12,96]
        satazi=np.array(f['/All_Data/ATMS-SDR-GEO_All/SatelliteAzimuthAngle']) #[12,96]
        satran=np.array(f['/All_Data/ATMS-SDR-GEO_All/SatelliteRange']) #[12,96]
        f.close()
        nscan=tem.shape[0]
        ta=np.zeros_like(tem,dtype=float)
        for ich in range(22):
            for isc in range(nscan):
                for ifv in range(self.nfov):
                    ta[isc,ifv,ich]=float(tem[isc,ifv,ich])*fac[0]+fac[1]
        return [ta,lat,lon,satzen,satazi,satran,tim,vel,nscan]


    def create_filelist(self):
        if self.geosdr_combine == 1:
            fname_prefix=self.input_dir+'/GATMO-SATMS_*.h5'
            flist=glob.glob(fname_prefix)
            flist.sort(key=lambda x:x.split('_')[4])
            idx_str=os.path.basename(flist[0]).find('_d')
            uniq_key=[]
            uniq_filename=[]
            for fname in flist:
                filename=os.path.basename(fname)
                key=filename[0:idx_str+35]

                if key.split('_')[2] != 'd'+self.ymd: 
                    continue

                if key not in uniq_key:
                    uniq_key.append(key)
                    uniq_filename.append(filename)
                else:
                    idx=uniq_key.index(key)
                    uniq_filename[idx]=filename

            self.sdr_filelist=uniq_filename
            self.nfile=len(self.sdr_filelist)


