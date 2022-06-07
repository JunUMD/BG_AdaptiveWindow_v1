# -*- coding: utf-8 -*-
"""                 
Created on 12/25/2020
                    
@author: Jun Zhou @UMD             
"""
import matplotlib as mpl
#mpl.use('Agg')              
import os,pickle
from pathlib import Path
import numpy as np
import sys
from numpy.linalg import det,inv,pinv,svd,eig,cond
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from numpy import dot
from mpl_toolkits.basemap import Basemap


class determine_coef_class:
    def __init__(self,nedt,ifr,numxgrd,numygrd,
                 win_file,uvs_file,
                 noise_lev,save_file):
       
        self.nedt=nedt
        self.ifr=ifr
        self.numxgrd=numxgrd
        self.numygrd=numygrd
        self.win_file=win_file
        self.uvs_file=uvs_file
        self.noise_lev=noise_lev
        self.save_file=save_file

    def ingest(self):
        f=open(self.win_file,'rb')
        self.windowsize=np.fromfile(f,dtype='>i4',count=1)[0]
        self.windowindx=np.fromfile(f,dtype='>i4',count=self.windowsize*2).reshape(2,self.windowsize).T.tolist()
        self.idx_src=np.fromfile(f,dtype='>i2',count=1)[0]
        self.cen_indx=self.windowindx[self.idx_src-1]
        f.close()

        self.windowsize0=self.windowsize
        self.windowindx0=self.windowindx

        f=open(self.uvs_file,'rb')
        self.U=np.fromfile(f,dtype='>f8',count=self.windowsize)
        self.V=np.fromfile(f,dtype='>f8',count=self.windowsize)
        nn=self.windowsize**2
        self.S=np.fromfile(f,dtype='>f8',count=nn).reshape(self.windowsize,self.windowsize).T
        f.close()

       

    def process(self):
        w=0.001

        min_gama0=0
        max_gama0=np.radians(90)
        delta2=self.nedt**2

        step_gama=(max_gama0-min_gama0)/10
        gamas=np.arange(min_gama0,max_gama0+step_gama,step_gama)

        dif_bst=9999  
        icy=0  
        while True:
            icy+=1
            e2s=[]
            alphas=[]
            for igm,gama in enumerate(gamas):
                alpha=self.calculate(self.U,self.V,self.S,w,gama,delta2)
                alphas.append(alpha)
                e2s.append(np.sqrt(np.sum(np.array(alpha)**2)))

            id_bst=np.argmin(np.abs(np.array(e2s)-self.noise_lev))
            alpha_bst=alphas[id_bst]
            dif_bst=np.abs(e2s[id_bst]-self.noise_lev)
            if dif_bst<1e-5 or step_gama<1e-10: break

            min_gama=gamas[id_bst]-step_gama
            if min_gama<0: min_gama=0
            max_gama=gamas[id_bst]+step_gama
            if max_gama>max_gama0: 
                exit()
            step_gama=step_gama/10.0
            gamas=np.arange(min_gama,max_gama+step_gama,step_gama)

            
        with open(self.save_file,'wb') as f:
            pickle.dump([self.windowsize,self.windowindx,self.idx_src,alpha_bst],f)
           

    def calculate(self,U,V,S0,w,gama,delta2):
        n=S0.shape[0]
        S=S0*np.cos(gama)
        for i in range(n):
            S[i,i]+=w*np.sin(gama)*delta2

        alpha_val1 = inv(S)
        temparray1 = dot(U.T,alpha_val1)  
        val1       = dot(temparray1,V)    
        val2       = dot(temparray1,U)    
        lamda      = (1 - val1*np.cos(gama))
        alpha_val2 = lamda*U
        alpha_val3 = 1.0/val2
        alpha_val4 = np.cos(gama)*V + alpha_val3*alpha_val2 
        alpha      = dot(alpha_val1,alpha_val4)
        return alpha

