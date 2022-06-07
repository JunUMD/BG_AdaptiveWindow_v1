# -*- coding: utf-8 -*-
"""                 
Created on 12/25/2020
                    
@author: Jun@UMD             
"""
import matplotlib as mpl
#mpl.use('Agg')              
import os,pickle
from pathlib import Path
import numpy as np
from numpy.linalg import det,inv,pinv,svd,eig,cond
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from numpy import dot
import glob
from scipy.interpolate import griddata
import math

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from skimage.measure import EllipseModel

from Gaussian_fit import twoD_Gaussian_fit,twoD_Gaussian,\
                         oneD_Gaussian_fit,oneD_Gaussian

class display_G_class:
    def __init__(self,ifr,G_fv_dir,G_fv_figdir,delx,dely):
        self.ifr=ifr
        self.G_fv_dir=G_fv_dir
        self.G_fv_figdir=G_fv_figdir
        self.delx=delx
        self.dely=dely

    """ use this when the area of src and tgt are different"""
    def display_G(self,tab):
        with open(self.G_fv_dir+tab+'.pickle','rb') as f:
            [antwgt,antpos_nx,antpos_ny]=pickle.load(f)

        # record the idd of src, will be applied to tgt and sim---
        if tab=='Gsrc': 
            self.idd=np.where(antwgt>0)
            #print(self.idd)
            #self.nx_src=np.max(self.idd[0])
            #self.ny_src=np.max(self.idd[1])
            #print(self.nx_src,self.ny_src)
            #print(np.min(self.idd[0]),np.min(self.idd[1]))

        z=antwgt[self.idd].ravel()
        x=antpos_nx[self.idd].ravel()*self.delx
        y=antpos_ny[self.idd].ravel()*self.dely

        # normalize and centerize 
        idd=np.argmax(z)
        z=z/z[idd]
        x=x-x[idd]
        y=y-y[idd]

        # centerize more accurately ------------
        try:
            [popt,pcov]=twoD_Gaussian_fit(x,y,z)
        except:
            popt=np.full((3),np.nan)
            pcov=np.full((1,1),np.inf)
        if not math.isinf(pcov[0,0]):
            if popt[0]>0: z=z/popt[0]
            if popt[1]>0: x=x-popt[1]
            if popt[2]>0: y=y-popt[2]
        # --------------------------------------
        n_bm=0
        trd=0.001
        while n_bm < 25:
            idd=np.where(np.abs(z-0.5)<trd)
            n_bm=len(idd[0])
            trd += 0.001
        print('trd=',trd)
        print('n_bm=',n_bm)    
        x_bm=x[idd]
        y_bm=y[idd]
        z_bm=z[idd]
        bmwidth=np.nan
        dis_alongTr=np.nan
        dis_crossTr=np.nan

        if n_bm>5: # 5 independent parameters
            bmwidth=np.nan
            points=np.array([x_bm,y_bm]).T
            ell = EllipseModel()
            ell.estimate(points)
            xc, yc, semi_major, semi_minor, theta = ell.params
            dis_alongTr,dis_crossTr=semi_minor*2,semi_major*2
            if self.ifr==47: bmwidth=np.degrees(np.mean([semi_major, semi_minor])*2.0/824)

        # record the area of src, will be used for tgt and sim ---
        if tab=='Gsrc': self.area=[np.min(x),np.max(x),np.min(y),np.max(y)]

        cm=plt.cm.get_cmap('jet')
        font_jun={'size':30}
        mpl.rc('font',**font_jun)
        fig=plt.figure(figsize=(15,15))

        [xmin,xmax,ymin,ymax]=self.area
        dd=self.delx
        x_side=np.arange(xmin,xmax+dd,dd)
        y_side=np.arange(ymin,ymax+dd,dd)
        x_grd,y_grd=np.meshgrid(x_side,y_side)
        points=np.array([x,y]).T
        z_grd=griddata(points,z,(x_grd,y_grd),method='linear')
        ax=fig.add_subplot(111)
        im=ax.imshow(z_grd,cmap='jet',origin='lower',
                     extent=[xmin,xmax,ymin,ymax],
                     vmin=0, vmax=1)

        if n_bm>5:
            ax.scatter(x_bm,y_bm,s=20,marker='X',c='k')
            ell1 = Ellipse(xy = (xc,yc),
                           width = semi_major*2, height = semi_minor*2,
                           angle = np.degrees(theta),
                           color='k', fill=False, linewidth=0.5, linestyle='--', alpha=1)
            ax.add_patch(ell1)

        divider=make_axes_locatable(ax)
        cax=divider.append_axes('right',     #position
                                size='5%',   #width
                                pad=0.2)     #seperation between graph and colorbar
        cb=fig.colorbar(im,cax=cax)
        info='Beam width = '+str(bmwidth)+'\n'+\
             'along tr = '+str(dis_alongTr)+'\n'+\
             'cross tr = '+str(dis_crossTr)
        ax.set_xlabel(info)

        plt.savefig(self.G_fv_figdir+tab+'.png',dpi=300,bbox_inches='tight')
        plt.close()
        print(self.G_fv_figdir+tab+'.png')



    """ use this when the area of tgt is the same as that of src"""
    def display_G_0(self,tab):
        with open(self.G_fv_dir+tab+'.pickle','rb') as f:
            [antwgt,antpos_nx,antpos_ny]=pickle.load(f)

        # record the idd of src, will be applied to tgt and sim---
        # some glitch at break, set this restraint to avoid showing them in figure
        # 0.02 ~ -17dB
        if tab=='Gsrc_wgtave': self.idd=np.where(antwgt>0.02)

        z=antwgt[self.idd].ravel()
        x=antpos_nx[self.idd].ravel()*self.delx
        y=antpos_ny[self.idd].ravel()*self.dely

        # normalize and centerize 
        idd=np.argmax(z)
        z=z/z[idd]
        x=x-x[idd]
        y=y-y[idd]

        # centerize more accurately ------------
        try:
            [popt,pcov]=twoD_Gaussian_fit(x,y,z)
        except:
            popt=np.full((3),np.nan)
            pcov=np.full((1,1),np.inf)
        if not math.isinf(pcov[0,0]):
            fit_info='Gaussian fitting succeed!'
            if popt[0]>0: z=z/popt[0]
            if popt[1]>0: x=x-popt[1]
            if popt[2]>0: y=y-popt[2]
        else:
            fit_info='Gaussian fitting failed!'
        # --------------------------------------

        idd=np.where(np.abs(z-0.5)<0.01) #0.01
        x_bm=x[idd]
        y_bm=y[idd]
        z_bm=z[idd]
        n_bm=len(x_bm)
        bmwidth=np.nan
        dis_alongTr=np.nan
        dis_crossTr=np.nan

        if n_bm>5: # 5 independent parameters
            bmwidth=np.nan
            points=np.array([x_bm,y_bm]).T
            ell = EllipseModel()
            ell.estimate(points)
            xc, yc, semi_major, semi_minor, theta = ell.params
            dis_alongTr,dis_crossTr=semi_minor*2,semi_major*2
            if self.ifr==47: 
                ll=np.mean([semi_major, semi_minor])*2.0
                if tab=='Gsrc_wgtave': self.alt=ll/np.radians(5.2) # use src to estimate the satellite altitude
                bmwidth=np.degrees(ll/self.alt)

        # record the area of src, will be used for tgt and sim ---
        if tab=='Gsrc_wgtave': self.area=[np.min(x),np.max(x),np.min(y),np.max(y)]

        cm=plt.cm.get_cmap('jet')
        font_jun={'size':30}
        mpl.rc('font',**font_jun)
        fig=plt.figure(figsize=(15,15))

        [xmin,xmax,ymin,ymax]=self.area
        dd=self.delx
        x_side=np.arange(xmin,xmax+dd,dd)
        y_side=np.arange(ymin,ymax+dd,dd)
        x_grd,y_grd=np.meshgrid(x_side,y_side)
        points=np.array([x,y]).T
        z_grd=griddata(points,z,(x_grd,y_grd),method='linear')
        ax=fig.add_subplot(111)
        im=ax.imshow(z_grd,cmap='jet',origin='lower',
                     extent=[xmin,xmax,ymin,ymax],
                     vmin=0, vmax=1)

        if n_bm>5:
            ax.scatter(x_bm,y_bm,s=20,marker='X',c='k')
            ell1 = Ellipse(xy = (xc,yc),
                           width = semi_major*2, height = semi_minor*2,
                           angle = np.degrees(theta),
                           color='k', fill=False, linewidth=0.5, linestyle='--', alpha=1)
            ax.add_patch(ell1)

        divider=make_axes_locatable(ax)
        cax=divider.append_axes('right',     #position
                                size='5%',   #width
                                pad=0.2)     #seperation between graph and colorbar
        cb=fig.colorbar(im,cax=cax)
        info='Beam width = '+str(bmwidth)+'\n'+\
             'along tr = '+str(dis_alongTr)+'\n'+\
             'cross tr = '+str(dis_crossTr)+'\n'+\
             fit_info
        ax.set_xlabel(info)

        plt.savefig(self.G_fv_figdir+tab+'.png',dpi=300,bbox_inches='tight')
        plt.close()
        print(self.G_fv_figdir+tab+'.png')

