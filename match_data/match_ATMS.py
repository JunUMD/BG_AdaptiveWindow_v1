# -*- coding: utf-8 -*-
"""
Created on 9/21/2020

@author: Jun@UMD
"""
import matplotlib as mpl
#mpl.use('Agg')
from pathlib import Path
import sys,pickle,os,math
import numpy as np
from match_data.ingest.ingest_ATMS import ingest_ATMS_class
from match_data.project.project_ATMS_lib.project_ATMS import project_ATMS_class
from match_data.combine.combine_ATMS import combine_ATMS_class

class match_ATMS_class:
    def __init__(self,control,input_path,output_path,sensor,select_chs,select_bws,isc_atm,winsize):
        self.input_path=input_path
        self.output_path=output_path
        self.sensor=sensor
        self.select_chs=select_chs
        self.select_bws=select_bws
        self.winsize=winsize
        self.isc_atm=isc_atm
        self.control=control


    def match_data(self):
        if self.control==1:
            ymd='20190927'
            orb='both' 
            select_enlarge=[2.5,self.select_bws[0]*2.5/self.select_bws[1]] 
            domain=[-90,90,-180,180]
            select_orbnum='b09624' 
            input_dir=os.path.join(self.input_path,ymd,'sdr/')

            ingest_ATMS_obj=ingest_ATMS_class(ymd,input_dir,self.output_path,self.select_chs,orb,domain,select_orbnum)
            ingest_ATMS_obj.ingest()

            project_ATMS_obj=project_ATMS_class(self.output_path,self.select_chs,orb,domain,
                                                self.select_bws,select_enlarge,ingest_ATMS_obj)
            project_ATMS_obj.project()
       
            combine_obj=combine_ATMS_class(self.output_path,self.isc_atm,self.winsize,
                                           ingest_ATMS_obj,
                                           project_ATMS_obj)
            combine_obj.comebine_display()

            self.select_nch     = ingest_ATMS_obj.select_nch
            self.nfov           = ingest_ATMS_obj.nfov
            self.satzen         = ingest_ATMS_obj.satzen
            self.satazi         = ingest_ATMS_obj.satazi
            self.satran         = ingest_ATMS_obj.satran
            self.lat            = ingest_ATMS_obj.lat
            self.lon            = ingest_ATMS_obj.lon
            self.beamwidth      = project_ATMS_obj.beamwidth
            self.select_enlarge = project_ATMS_obj.select_enlarge
            self.fovedge_lat    = project_ATMS_obj.fovedge_lat
            self.fovedge_lon    = project_ATMS_obj.fovedge_lon
           
            with open(os.path.join(self.output_path,'match.pickle'),'wb') as f:
                pickle.dump(self,f)

 
        elif self.control==0:
            with open(os.path.join(self.output_path,'match.pickle'),'rb') as f:
                obj=pickle.load(f)
            self.winsize        = obj.winsize
            self.isc_atm        = obj.isc_atm
            self.select_chs     = obj.select_chs
            self.select_nch     = obj.select_nch
            self.nfov           = obj.nfov
            self.satzen         = obj.satzen
            self.satazi         = obj.satazi
            self.satran         = obj.satran
            self.lat            = obj.lat
            self.lon            = obj.lon
            self.beamwidth      = obj.beamwidth
            self.select_enlarge = obj.select_enlarge
            self.fovedge_lat    = obj.fovedge_lat
            self.fovedge_lon    = obj.fovedge_lon
 


    

