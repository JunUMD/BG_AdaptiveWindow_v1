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
from mpl_toolkits.basemap import Basemap
import glob
from scipy.interpolate import griddata
import math

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from calculate_BGcoef.adjust_window_pythonlib.Gaussian_fit import twoD_Gaussian_fit,twoD_Gaussian,\
                                                                  oneD_Gaussian_fit,oneD_Gaussian
from calculate_BGcoef.adjust_window_pythonlib.FitCircle import FitCircle_class
from calculate_BGcoef.adjust_window_pythonlib.module_fit_ellipse import Ellipse_type

class display_adjust_class:
    def __init__(self,delx,dely,ifr,winadj_subdir,winadj_subfigdir,
                 tgt_edg_lat,tgt_edg_lon,
                 tgt_cen_lat,tgt_cen_lon,
                 src_cen_lat,src_cen_lon):
        self.delx=delx
        self.dely=dely
        self.ifr=ifr
        self.winadj_subdir=winadj_subdir
        self.winadj_subfigdir=winadj_subfigdir
        # degree
        self.tgt_edg_lat=tgt_edg_lat  #[nphi]
        self.tgt_edg_lon=tgt_edg_lon
        self.tgt_cen_lat=tgt_cen_lat
        self.tgt_cen_lon=tgt_cen_lon
        self.src_cen_lat=src_cen_lat  #[nsc,nfv]
        self.src_cen_lon=src_cen_lon

        self.fit_bm='cir'  # or ell

        fname=os.path.join(self.winadj_subdir,'winsize_Q0_e2_gm.pickle')
        with open(fname,'rb') as f:
            [winsize_allwin,Q0_allwin,e2_allwin,gm_allwin,
             alpha_allwin,sumalpha_allwin,cond_S_allwin]=pickle.load(f)
        self.winsize_all=winsize_allwin  # including gm == 0

    def display_winsize_Q0_e2_gama(self):
        fname=os.path.join(self.winadj_subdir,'winsize_Q0_e2_gm.pickle')
        with open(fname,'rb') as f:
            [winsize_allwin,Q0_allwin,e2_allwin,gm_allwin,
             alpha_allwin,sumalpha_allwin,cond_S_allwin]=pickle.load(f)
        winsize_allwin,Q0_allwin,e2_allwin,gm_allwin,\
        alpha_allwin,sumalpha_allwin,cond_S_allwin = np.array(winsize_allwin),\
                                       np.array(Q0_allwin),\
                                       np.array(e2_allwin),\
                                       np.array(gm_allwin), \
                                       np.array(alpha_allwin), \
                                       np.array(sumalpha_allwin),\
                                       np.array(cond_S_allwin)
        self.winsize_allwin=winsize_allwin  # including gm == 0
        idd=np.where(gm_allwin>0)
        self.winsize_all=winsize_allwin #[idd]
        self.e2_all=e2_allwin #[idd]
        self.gm_all=gm_allwin #[idd]
        self.Q0_all=Q0_allwin #[idd]
        self.alpha_all=alpha_allwin
        self.sumalpha_all=sumalpha_allwin
        self.cond_S_all=cond_S_allwin

        #self.e2_all=np.log10(self.e2_all)
        #self.Q0_all=np.log10(self.Q0_all)

        font_jun={'family':'Times New Roman','size':30}
        mpl.rc('font',**font_jun)
        fig=plt.figure(figsize=(50,8))
        ax1=fig.add_subplot(142)
        ax2=ax1.twinx()
        ax1.plot(self.winsize_all,self.Q0_all,'b-',linewidth=2,marker='o',markersize=10)
        ax2.plot(self.winsize_all,self.e2_all,'r-',linewidth=2,marker='o',markersize=10)
        ax1.set_xlabel('win size')
        ax1.set_ylabel('Q0 of the best point',color='b')
        ax2.set_ylabel('e2 of the best point',color='r')

        ax=fig.add_subplot(143)
        ax.plot(self.winsize_all,self.sumalpha_all,'k-',linewidth=2,marker='o',markersize=10)
        ax.set_xlabel('win size')
        ax.set_ylabel('sumalpha bst')

        ax=fig.add_subplot(141)
        ax.plot(self.winsize_all,np.log10(self.cond_S_all),'g-',linewidth=2,marker='o',markersize=10)
        ax.set_xlabel('win size')
        ax.set_ylabel('Cond num of S')

        ax=fig.add_subplot(144)
        ax.plot(self.winsize_all,self.gm_all,'m-',linewidth=2,marker='o',markersize=10)
        ax.set_xlabel('win size')
        ax.set_ylabel('gm bst')
        plt.savefig(self.winadj_subfigdir+'winsize_Q0_e2_gm.png',dpi=300,bbox_inches='tight')
        plt.close()

    def display_winsize_bmwidth(self):
        fname=os.path.join(self.winadj_subdir,'Estimate_bmwidth_dis.pickle')
        with open(fname,'rb') as f:
            [bmwidth,dis_alongTr,dis_crossTr]=pickle.load(f)

        x0=self.winsize_all

        bmwidth,dis_alongTr,dis_crossTr=np.array(bmwidth),np.array(dis_alongTr),np.array(dis_crossTr)

        font_jun={'family':'Times New Roman','size':40}
        mpl.rc('font',**font_jun)

        if self.ifr==47:
            idd=np.where(np.logical_and(~np.isnan(bmwidth),self.gm_all>0))
            y1=bmwidth[idd]
            y2=self.e2_all[idd]
            x=self.winsize_allwin[idd]            
            fig=plt.figure(figsize=(20,40))
            ax1=fig.add_subplot(211)
            ax2=ax1.twinx()
            ax1.plot(x,y1,'b-',linewidth=2,marker='o',markersize=10)
            ax2.plot(x,y2,'r-',linewidth=2,marker='o',markersize=10)
            ax1.set_xlabel('winsize')
            ax1.set_ylabel('beam width')
            ax2.set_ylabel('e2')
            ax=fig.add_subplot(212)
        else:
            fig=plt.figure(figsize=(20,20))
            ax=fig.add_subplot(111)
        idd=np.where(np.logical_and(~np.isnan(dis_alongTr),self.gm_all>0))
        y=dis_alongTr[idd]
        x=self.winsize_all[idd]
        ax.plot(x,y,'b-',linewidth=2,marker='o',markersize=10,label='along-track')
        idd=np.where(np.logical_and(~np.isnan(dis_crossTr),self.gm_all>0))
        y=dis_crossTr[idd]
        x=self.winsize_all[idd]
        ax.plot(x,y,'g-',linewidth=2,marker='o',markersize=10,label='cross-track')
        ax.legend()
        ax.set_xlabel('winsize')
        ax.set_ylabel('along cross track dis')
        plt.savefig(self.winadj_subfigdir+'Beamwidth_dis_e2.png',dpi=300,bbox_inches='tight')
        plt.close()    

    def display_window(self):
        for i,winsiz in enumerate(self.winsize_all):
            if i==0:
                fname=self.winadj_subdir+'Winsize_'+str(winsiz)+'_full.pickle'
                tab='Window_win'+str(winsiz)+'_full.png'
            else:
                fname=self.winadj_subdir+'Winsize_'+str(winsiz)+'.pickle'
                tab='Window_win'+str(winsiz)+'.png'

            print(i,os.path.basename(fname))
            with open(fname,'rb') as f:
                [alpha_bst,Q0_bst,e2_bst,gm_bst,
                 windowsize,windowindx,idx_src,S,cond_num]=pickle.load(f)

            if i==0: 
                self.windowsize0=windowsize
                self.windowindx0=windowindx

            fsave=self.winadj_subfigdir+tab
            self.display_img(alpha_bst,Q0_bst,e2_bst,gm_bst,S,cond_num,
                             windowsize,windowindx,idx_src,fsave)

    def display_img(self,alpha_bst,Q0_bst,e2_bst,gm_bst,S,cond_num,
                    windowsize,windowindx,idx_src,fsave):

        winsize=windowsize

        font_jun={'size':50}
        mpl.rc('font',**font_jun)
        cm=plt.cm.get_cmap('jet')
        fig=plt.figure(figsize=(40,40))
      
        ax=fig.add_subplot(221)
        im=ax.imshow(S,cmap='jet',origin='lower',
                     extent=[0,windowsize,0,windowsize])
        divider=make_axes_locatable(ax)
        cax=divider.append_axes('right',     #position
                                size='5%',   #width
                                pad=0.2)     #seperation between graph and colorbar
        cb=fig.colorbar(im,cax=cax)
        ax.set_xlabel('log10(Cond) = '+str(np.log10(cond_num)))

        ax=fig.add_subplot(222)
        lon_cen=self.src_cen_lon[windowindx[idx_src-1][0],windowindx[idx_src-1][1]]
        lat_cen=self.src_cen_lat[windowindx[idx_src-1][0],windowindx[idx_src-1][1]]
        lon0=[]
        lat0=[]
        lon=[]
        lat=[]
        for [isc,ifv] in self.windowindx0:
            lon0.append(self.src_cen_lon[isc,ifv])
            lat0.append(self.src_cen_lat[isc,ifv])
            if [isc,ifv] in windowindx:
                lon.append(self.src_cen_lon[isc,ifv])
                lat.append(self.src_cen_lat[isc,ifv])
        minlon=min(np.min(self.tgt_edg_lon),np.min(lon0))-0.1
        maxlon=max(np.max(self.tgt_edg_lon),np.max(lon0))+0.1
        minlat=min(np.min(self.tgt_edg_lat),np.min(lat0))-0.1
        maxlat=max(np.max(self.tgt_edg_lat),np.max(lat0))+0.1
        cm=plt.cm.get_cmap('rainbow')
        m=Basemap(projection='cyl', #map projection,cyl is default
                  llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,
                  resolution='c') #c: crude; l:low; i:intermediate, h:high, f:full; sometimes can only >=h
        m.drawmapboundary(fill_color='lightblue') #the whole map is filled with the specified color
        m.fillcontinents(color='wheat',lake_color='lightblue') #the continents are filled
        m.drawcoastlines()  #draw coastline, from gshhs
        sc=m.plot(self.tgt_edg_lon,self.tgt_edg_lat,color='b',linewidth=2,zorder=10)
        m.scatter(lon0,lat0,c='k',marker='X',s=200,zorder=11)
        m.scatter(lon,lat,c='orangered',marker='X',s=200,zorder=12)
        m.scatter(lon_cen,lat_cen,c='m',marker='X',s=400,zorder=13)

        ax=fig.add_subplot(223)
        x=np.arange(winsize)
        ax.plot(x,alpha_bst,'-m',linewidth=2,marker='o',markersize=15,zorder=12,label='igm=best')
        ax.legend()
        ax.set_xlabel('window index \n sum(Coef)='+str(np.sum(alpha_bst))+'\n gama bst = '+str(gm_bst))
        ax.set_ylabel('alpha')

        plt.savefig(fsave,dpi=300,bbox_inches='tight')
        plt.close()

    def display_Gsim(self):
        bmwidth=[]
        dis_alongTr=[]
        dis_crossTr=[]
        for i,winsiz in enumerate(self.winsize_all):
            #if winsiz > 120: continue

            if i==0:         
                fname=self.winadj_subdir+'Gsim_win'+str(winsiz)+'_full.pickle'
                tab='Gsim_win'+str(winsiz)+'_full'
            else:
                fname=self.winadj_subdir+'Gsim_win'+str(winsiz)+'.pickle'
                tab='Gsim_win'+str(winsiz)
            
            print(i,os.path.basename(fname)) 

            if os.path.isfile(fname):
                print(tab)
                param=self.display_G(tab)
                bmwidth.append(param[0])
                dis_alongTr.append(param[1])
                dis_crossTr.append(param[2])

        fname=os.path.join(self.winadj_subdir,'Estimate_bmwidth_dis.pickle')
        with open(fname,'wb') as f:
            pickle.dump([bmwidth,dis_alongTr,dis_crossTr],f)


    def display_G(self,tab):
        print(tab) 
        with open(self.winadj_subdir+tab+'.pickle','rb') as f:
            [antwgt,antpos_nx,antpos_ny]=pickle.load(f)

        idd=np.where(antwgt>0)
        z=antwgt[idd].ravel()
        x=antpos_nx[idd].ravel()*self.delx
        y=antpos_ny[idd].ravel()*self.dely
       
        # normalize and centerize 
        idd=np.argmax(z)
        z=z/z[idd]
        x=x-x[idd]
        y=y-y[idd]

        # centerize more accurately
        [popt,pcov]=twoD_Gaussian_fit(x,y,z)
        if not math.isinf(pcov[0,0]):
            if popt[0]>0: z=z/popt[0]
            if popt[1]>0: x=x-popt[1]
            if popt[2]>0: y=y-popt[2]

        idd=np.where(np.abs(z-0.5)<0.02)
        x_bm=x[idd]
        y_bm=y[idd]
        z_bm=z[idd]
        n_bm=len(x_bm)
        bmwidth=np.nan
        dis_alongTr=np.nan
        dis_crossTr=np.nan

        if n_bm>2:
           if self.fit_bm=='cir': 
               # for circle fitting --------
               fitcircle_obj=FitCircle_class(x_bm,y_bm)
               fitcircle_obj.fit()
               x_cr,y_cr=fitcircle_obj.x_fit2,fitcircle_obj.y_fit2
               dis_alongTr,dis_crossTr=fitcircle_obj.R_2*2.0,fitcircle_obj.R_2*2.0
               if self.ifr==47: bmwidth=np.degrees(fitcircle_obj.R_2*2.0/824)
           else:
               # for ellipse fitting ---------
               ellipse_obj=Ellipse_type()
               center, param_axis, param_orien=ellipse_obj.fit_ellipse([x_bm,y_bm,z_bm])
               dis_alongTr,dis_crossTr=param_axis[:]*2  #???               
               if self.ifr==47: bmwidth=np.degrees(np.mean(param_axis)*2.0/824)

        if tab=='Gsrc': self.area=[np.min(x),np.max(x),np.min(y),np.max(y)]       

        cm=plt.cm.get_cmap('jet')
        font_jun={'size':30}
        mpl.rc('font',**font_jun)
        fig=plt.figure(figsize=(15,15))
        #ax=fig.add_subplot(111)
        #ax.scatter(x,y,c=z,cmap='jet',s=10,marker='o')
        #if n_bm>2:
        #    ax.scatter(x_bm,y_bm,s=20,marker='X',c='k')
        #    if self.fit_bm=='cir':
        #        ax.plot(x_cr,y_cr,'k-')
        #    else:
        #        width, height = param_axis
        #        ell1 = Ellipse(xy = center, width = width, height = height, angle = param_orien,
        #                   color='k', fill=False, linewidth=0.5, linestyle='--', alpha=1)
        #        ax.add_patch(ll1)
        #ax.set_xlim(self.area[0],self.area[1])
        #ax.set_ylim(self.area[2],self.area[3])

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
        if n_bm>2:
            ax.scatter(x_bm,y_bm,s=20,marker='X',c='k')
            if self.fit_bm=='cir':
                ax.plot(x_cr,y_cr,'k-')
            else:
                width, height = param_axis
                ell1 = Ellipse(xy = center, width = width, height = height, angle = param_orien,
                           color='k', fill=False, linewidth=2, linestyle='--', alpha=1)
                ax.add_patch(ell1)
        divider=make_axes_locatable(ax)
        cax=divider.append_axes('right',     #position
                                size='5%',   #width
                                pad=0.2)     #seperation between graph and colorbar
        cb=fig.colorbar(im,cax=cax)
        if self.ifr==47:
            info='Beam width = '+str(bmwidth)+'\n'+\
                 'along tr = '+str(dis_alongTr)+'\n'+\
                 'cross tr = '+str(dis_crossTr)
        else:
            info='along tr = '+str(dis_alongTr)+'\n'+\
                 'cross tr = '+str(dis_crossTr)

        ax.set_xlabel(info)


        plt.savefig(self.winadj_subfigdir+tab+'.png',dpi=300,bbox_inches='tight')
        plt.close()
       
        print('Estimated Beamwidth = ',bmwidth)        
        print('          Along-trc = ',dis_alongTr)        
        print('          Cross-trc = ',dis_crossTr)        
        return bmwidth,dis_alongTr,dis_crossTr



