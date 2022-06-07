# -*- coding: utf-8 -*-
"""
Created on 9/21/2020

@author: Jun@UMD
"""
import matplotlib as mpl
#mpl.use('Agg')
import os
import numpy as np
import pickle
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

from match_data.project.project_ATMS_lib.ant_prj import ant_prj_mod

class project_ATMS_ARTS_class:
    def __init__(self):
        self.beamwidth=0
        self.select_enlarge=0   


class project_ATMS_class:
    def __init__(self,output_dir,select_chs,orb,domain,
                 select_bws,select_enlarge,ingest_ATMS_obj):
        print('Project ...')
        self.output_dir=output_dir
        self.beamwidth=np.radians(select_bws)
        self.select_enlarge=select_enlarge
   
        self.nch=ingest_ATMS_obj.select_nch 
        self.chs=select_chs 
        self.nscan=ingest_ATMS_obj.nscan
        self.nfov=ingest_ATMS_obj.nfov
        # nscan,nfov
        self.sat_zen=np.radians(ingest_ATMS_obj.satzen)
        self.sat_azi=np.radians(ingest_ATMS_obj.satazi)
        self.sat_ran=ingest_ATMS_obj.satran

        # nch,nscan,nfov
        self.ta=ingest_ATMS_obj.ta
        self.lat=np.radians(ingest_ATMS_obj.lat)
        self.lon=np.radians(ingest_ATMS_obj.lon)


    def project(self):
        self.nphi=37 
        self.dphi=np.radians(360.0/(self.nphi-1))

        self.fovedge_lat=np.zeros((self.nch,self.nscan,self.nfov,self.nphi))
        self.fovedge_lon=np.zeros((self.nch,self.nscan,self.nfov,self.nphi))
       
        for ich in range(self.nch):
            for isc in range(self.nscan):
                for ifv in range(self.nfov):
                    ant_prj_mod.ant_prj(self.nphi,self.dphi,                  #radians
                                        self.beamwidth[ich]*self.select_enlarge[ich],      #radians
                                        self.lat[ich,isc,ifv],   #radians
                                        self.lon[ich,isc,ifv],   #radians
                                        self.sat_zen[isc,ifv],
                                        self.sat_azi[isc,ifv],   #radians
                                        self.sat_ran[isc,ifv])  # meters
                    self.fovedge_lat[ich,isc,ifv,:]=ant_prj_mod.lat
                    self.fovedge_lon[ich,isc,ifv,:]=ant_prj_mod.lon
