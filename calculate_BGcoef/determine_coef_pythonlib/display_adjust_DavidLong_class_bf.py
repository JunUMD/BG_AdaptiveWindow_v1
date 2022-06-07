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
#from skimage.measure import EllipseModel

from calculate_BGcoef.adjust_window_pythonlib.Gaussian_fit import twoD_Gaussian_fit,twoD_Gaussian,\
                                                                  oneD_Gaussian_fit,oneD_Gaussian
#from calculate_BGcoef.adjust_window_pythonlib.FitCircle import FitCircle_class
#from calculate_BGcoef.adjust_window_pythonlib.module_fit_ellipse import Ellipse_type


class display_adjust_DavidLong_class:
    def __init__(self,delx,dely):
        self.delx=delx
        self.dely=dely


    def display_compare_winsize(self):

        task='Remap_ATMS_5p2_3p3'

        if task=='Remap_ATMS_5p2_3p3':
            wins=['square3','square5','square7','square11','in0','out0p01']
            dBs=list(np.arange(-1,-13,-2))+[-999]

            input_dir='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_Ch1_5p2_Ch1_3p3/BGcoef/'
            com_dir='/WindowAdjust_noiselv2.5_stringSearch/atmch0/'
  
            nwin=0
            labels=[]
            winsiz_all=[]
            for win in wins:
                for dB in dBs:
                    if win!='out0p01' and dB!=-999: continue
                    nwin+=1
                    winsiz_all.append([])
                    labels.append(win+str(dB))
                    for fov in range(96):
                        fname=input_dir+'window_'+win+com_dir+'fov'+str(fov)+'/gm_Q0e2_'+str(dB).replace('-','neg')+'dB.pickle' 
                        print(fname)
                        with open(fname,'rb') as f:
                            [gamas,alphas,sum_alphas,Q0s,e2s,
                             windowsize,windowindx,idx_src,S,cond_num]=pickle.load(f)
                        winsiz_all[nwin-1].append(windowsize)

            save_dir='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_Ch1_5p2_Ch1_3p3/BGcoef/comp_winsize/'
            Path(save_dir).mkdir(parents=True,exist_ok=True) 


            colors=['k-','k--','k-.','grey',
                'r-','g-','b-','yellow',
                'm-','darkviolet','deepskyblue','darkorange']
      
            linewidths=[2,2,2,2,2,2,2,2,2,6,4,2] 
             
            font_jun={'family':'Times New Roman','size':30}
            mpl.rc('font',**font_jun)
            fig=plt.figure(figsize=(15,10))
            ax=fig.add_subplot(111)
            x=np.arange(96)
            for i,(color,linewidth,label) in enumerate(zip(colors,linewidths,labels)):
                print(label,color)
                ax.plot(x,winsiz_all[i],color,linewidth=linewidth,label=label)

            ax.set_xlabel('fov')
            ax.set_ylabel('winsize')
            ax.legend() 
           
            plt.savefig(save_dir+'Compare_winsize.png',dpi=300,bbox_inches='tight')
            print(save_dir+'Compare_winsize.png')
            plt.close() 


    def initialize_fov(self,ifr,winadj_subdir,winadj_subfigdir,
                       tgt_edg_lat,tgt_edg_lon,
                       tgt_cen_lat,tgt_cen_lon,
                       src_cen_lat,src_cen_lon,
                       windowsize0,windowindx0):
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
     
        self.windowsize0=windowsize0
        self.windowindx0=windowindx0

    def display_parameters(self,noise_lv):
        with open(self.winadj_subdir+'dBs_winsize_Q0_e2_gm.pickle','rb') as f:
            [dBs,winsize_allwin,Q0_allwin,e2_allwin,gm_allwin,
             alpha_allwin,sumalpha_allwin,cond_S_allwin]=pickle.load(f)

        font_jun={'family':'Times New Roman','size':40}
        mpl.rc('font',**font_jun) 
        fig=plt.figure(figsize=(50,20))
        ax=fig.add_subplot(231)
        ax.plot(dBs,np.log10(cond_S_allwin),'-m',linewidth=2,marker='o',markersize=10,label='cond S')
        ax.legend()
        ax.set_xlabel('dB')
        ax.set_ylabel('Cond S')

        ax=fig.add_subplot(232)
        ax.plot(dBs,Q0_allwin,'-b',linewidth=2,marker='o',markersize=10,label='Q0')
        ax.legend()
        ax.set_xlabel('dB')
        ax.set_ylabel('Q0')

        ax=fig.add_subplot(233)
        ax.plot(dBs,e2_allwin,'-r',linewidth=2,marker='o',markersize=10,label='e2')
        ax.legend()
        ax.set_ylim(noise_lv-0.001,noise_lv+0.001)
        ax.set_xlabel('dB')
        ax.set_ylabel('e2')

        ax=fig.add_subplot(234)
        ax.plot(dBs,gm_allwin,'-k',linewidth=2,marker='o',markersize=10,label='gama')
        ax.legend()
        ax.set_xlabel('dB')
        ax.set_ylabel('gama')

        ax=fig.add_subplot(235)
        ax.plot(dBs,sumalpha_allwin,'-g',linewidth=2,marker='o',markersize=10,label='sum alpha')
        ax.legend()
        ax.set_ylim(1-0.001,1+0.001)
        ax.set_xlabel('dB')
        ax.set_ylabel('sum alpha')

        plt.savefig(self.winadj_subfigdir+'Parameter_allwin.png',dpi=300,bbox_inches='tight')
        plt.close() 
        print(self.winadj_subfigdir+'Parameter_allwin.png')
           

    def display_window(self,dB):
        dB_str=str(dB).replace('-','neg')+'dB'
        with open(self.winadj_subdir+'Window_'+dB_str+'.pickle','rb') as f:
            [windowsize,windowindx,idx_src,S,cond_num]=pickle.load(f)

        print(self.winadj_subdir+'Window_'+dB_str+'.pickle')      

        winsize=windowsize

        font_jun={'size':50}
        mpl.rc('font',**font_jun)
        cm=plt.cm.get_cmap('jet')
        fig=plt.figure(figsize=(20,30))

        ax=fig.add_subplot(211)
        im=ax.imshow(S,cmap='jet',origin='lower',
                     extent=[0,windowsize,0,windowsize])
        divider=make_axes_locatable(ax)
        cax=divider.append_axes('right',     #position
                                size='5%',   #width
                                pad=0.2)     #seperation between graph and colorbar
        cb=fig.colorbar(im,cax=cax)
        ax.set_xlabel('log10(Cond) = '+str(np.log10(cond_num)))

        ax=fig.add_subplot(212)
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

        plt.savefig(self.winadj_subfigdir+'Window_fov'+str(self.ifr)+'_win'+str(winsize)+'_'+dB_str+'.png',dpi=300,bbox_inches='tight')
        plt.close()
        print(self.winadj_subfigdir+'Window_fov'+str(self.ifr)+'_win'+str(winsize)+'_'+dB_str+'.png')

    def display_alpha(self,dB):
        dB_str=str(dB).replace('-','neg')+'dB'
        with open(self.winadj_subdir+'Alpha_'+dB_str+'.pickle','rb') as f:
            [alpha]=pickle.load(f)

        print(self.winadj_subdir+'Alpha_'+dB_str+'.pickle')

        font_jun={'family':'Times New Roman','size':40}
        mpl.rc('font',**font_jun) 
        x=np.arange(len(alpha))
        fig=plt.figure(figsize=(20,10))
        ax=fig.add_subplot(111)
        ax.plot(x,alpha,'-m',linewidth=2,marker='o',markersize=15,zorder=12)
        ax.set_xlabel('window index')
        ax.set_ylabel('alpha')
        plt.savefig(self.winadj_subfigdir+'Alpha_fov'+str(self.ifr)+'_'+dB_str+'.png',dpi=300,bbox_inches='tight')
        plt.close()

    def display_G(self,tab):
        with open(self.winadj_subdir+tab+'.pickle','rb') as f:
            [antwgt,antpos_nx,antpos_ny]=pickle.load(f)
 
        # record the idd of src, will be applied to tgt and sim---
        if tab=='Gsrc': self.idd=np.where(antwgt>0)

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

        idd=np.where(np.abs(z-0.5)<0.01)
        x_bm=x[idd]
        y_bm=y[idd]
        z_bm=z[idd]
        n_bm=len(x_bm)
        bmwidth=np.nan
        dis_alongTr=np.nan
        dis_crossTr=np.nan

        if n_bm>5: # 5 independent parameters
            bmwidth=np.nan
            ## for circle fitting --------
            #fitcircle_obj=FitCircle_class(x_bm,y_bm)
            #fitcircle_obj.fit()
            #x_cr,y_cr=fitcircle_obj.x_fit2,fitcircle_obj.y_fit2
            #dis_alongTr,dis_crossTr=fitcircle_obj.R_2*2.0,fitcircle_obj.R_2*2.0
            #if self.ifr==47: bmwidth=np.degrees(fitcircle_obj.R_2*2.0/824)
            # for Ellipse fitting -------
            #with open('/mars/jzhou/XY.pickle','wb') as f:
            #    pickle.dump([x_bm,y_bm],f)           

 
            #ellipse_obj=Ellipse_type()
            #center, param_axis, param_orien=ellipse_obj.fit_ellipse([x_bm,y_bm,z_bm])
            #dis_alongTr,dis_crossTr=np.array(param_axis[:]) 

            points=np.array([x_bm,y_bm]).T
            #ell = EllipseModel()
            #ell.estimate(points)
            #xc, yc, semi_major, semi_minor, theta = ell.params
            #dis_alongTr,dis_crossTr=semi_minor*2,semi_major*2
            #if self.ifr==47: bmwidth=np.degrees(np.mean([semi_major, semi_minor])*2.0/824)


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
            #ax.plot(x_cr,y_cr,'k-')
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

        plt.savefig(self.winadj_subfigdir+tab+'.png',dpi=300,bbox_inches='tight')
        plt.close()
        print(self.winadj_subfigdir+tab+'.png')



    # only for string search-----------
    def display_compare_Q0e2(self):
        fov='fov'+str(self.ifr) 

        task='Remap_ATMS_5p2_3p3'

        if task=='Remap_ATMS_5p2_3p3':
            wins=['square3','out0p01']
            dBs=[-999,-9]
            #dBs=list(np.arange(-1,-13,-2))+[-999]

            input_dir='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_ARTSprj_GeoCrr_Ch1_5p2_Ch1_3p3/BGcoef/'
            com_dir='/WindowAdjust_noiselv999_stringSearch/atmch0/'+fov+'/'

            save_dir=input_dir+'/comp_Q0e2/'
            Path(save_dir).mkdir(parents=True,exist_ok=True) 

            nwin=0
            labels=[]
            fnames=[]
            for win,dB in zip(wins,dBs):
                nwin+=1
                fname=input_dir+'window_'+win+com_dir+'/gm_Q0e2_'+str(dB).replace('-','neg')+'dB.pickle'

                if (os.path.exists(fname) and os.path.getsize(fname) > 0 ):
                    fnames.append(fname)
                    labels.append(win+str(dB))
                    print(fname)

            colors=['k-','r-']
            #colors=['k-','k--','k-.',
            #    'r-','g-','lime','mediumspringgreen',
            #    'aquamarine','turquoise','teal','aqua','lightskyblue',
            #    'deepskyblue','dodgerblue','steelblue']
  

        elif task=='Remap_ATMS_5p2_3p3_less':
            win5='out0p01'
            input_dir='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_Ch1_5p2_Ch1_3p3/BGcoef/'
            com_dir='/WindowAdjust_noiselv2.5_stringSearch/atmch0/'+fov+'/'

            save_dir='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_Ch1_5p2_Ch1_3p3/BGcoef/comp_Q0e2/'
            Path(save_dir).mkdir(parents=True,exist_ok=True) 

            fnames=[
                input_dir+'window_'+win5+com_dir+'gm_Q0e2_neg1dB.pickle',
                input_dir+'window_'+win5+com_dir+'gm_Q0e2_neg3dB.pickle',
                input_dir+'window_'+win5+com_dir+'gm_Q0e2_neg5dB.pickle',
                input_dir+'window_'+win5+com_dir+'gm_Q0e2_neg7dB.pickle',
                input_dir+'window_'+win5+com_dir+'gm_Q0e2_neg9dB.pickle',
                input_dir+'window_'+win5+com_dir+'gm_Q0e2_neg11dB.pickle',
                input_dir+'window_'+win5+com_dir+'gm_Q0e2_neg999dB.pickle']

            labels=[
                win5+'neg1dB',win5+'neg3dB',win5+'neg5dB',win5+'neg7dB',
                win5+'neg9dB',win5+'neg11dB',win5+'neg999dB']

            colors=['k',
                'r','g','lime','mediumspringgreen',
                'aquamarine','turquoise','teal','aqua','lightskyblue',
                'deepskyblue','dodgerblue','steelblue']

        elif task=='Remap_ATMS_5p2_1p1':
            win1='square3'
            win2='out0p01'
            input_dir='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_5p2_1p1/BGcoef/'
            com_dir='/WindowAdjust_noiselv999_stringSearch/atmch0/'+fov+'/'

            save_dir='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_5p2_1p1/BGcoef/comp_Q0e2/'
            Path(save_dir).mkdir(parents=True,exist_ok=True)

            fnames=[input_dir+'window_'+win2+com_dir+'gm_Q0e2_neg999dB.pickle',
                    input_dir+'window_'+win2+com_dir+'gm_Q0e2_neg999dB.pickle'] 
            labels=[win1,win2]
            colors=['k-','r-']

        elif task=='Remap_ATMS_5p2_1p1_fov':
            input_dir='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_5p2_1p1/BGcoef/window_out0p01/WindowAdjust_noiselv999_stringSearch/atmch0/fov'
            fovs=[self.ifr,95-self.ifr]  #[47]+list(np.arange(88,96)) 
            save_dir='/mars/jzhou/data/BG_AdaptiveWindow_output/Remap_ATMS_5p2_1p1/BGcoef/comp_Q0e2_rightFovs/'
            Path(save_dir).mkdir(parents=True,exist_ok=True)
            
            colors=['k-',
                    'r-','g-','lime','mediumspringgreen',
                    'aquamarine','turquoise','teal','aqua','lightskyblue',
                    'deepskyblue','dodgerblue','steelblue']
            labels=['fov'+str(fov+1) for fov in fovs]
            fnames=[]
            for fov in fovs:
                fnames.append(input_dir+str(fov)+'/gm_Q0e2_neg999dB.pickle')

        font_jun={'family':'Times New Roman','size':30}
        mpl.rc('font',**font_jun)
        fig=plt.figure(figsize=(15,10))
        ax=fig.add_subplot(111)
        for i,(fname,color,label) in enumerate(zip(fnames,colors,labels)):
            print(label,color)
            with open(fname,'rb') as f:
                [gamas,alphas,sum_alphas,Q0s,e2s,
                windowsize,windowindx,idx_src,S,cond_num]=pickle.load(f)
            if i==0:
                id_bst=np.argmin(np.abs(np.array(e2s)-5))
                print(e2s[id_bst],gamas[id_bst],np.degrees(gamas[id_bst]))
          
            
            ax.plot(Q0s,e2s,color,linewidth=2,label=label)
            #ax.scatter(Q0s[id_bst],e2s[id_bst],c='grey',marker='X',s=100) 

        ax.set_xlabel('Q0')
        ax.set_ylabel('e2')
        ax.legend() 
        if task=='Remap_ATMS_5p2_1p1_fov':
            plt.savefig(save_dir+'Compare_Q0e2_noiselv5_fov'+str(self.ifr)+'_'+str(95-self.ifr)+'.png',dpi=300,bbox_inches='tight')
            print(save_dir+'Compare_Q0e2.png')
            plt.close()
        else:
            plt.savefig(save_dir+'Compare_Q0e2_'+fov+'.png',dpi=300,bbox_inches='tight')
            print(save_dir+'Compare_Q0e2_'+fov+'.png')
            plt.close()


    def display_compare_alpha(self):
        fov='fov'+str(self.ifr)
        dB=-9
        input_nogeo='/mars/jzhou/data/BG_AdaptiveWindow_output/'+\
                     'Remap_ATMS_ARTSprj_noGeoCrr_Ch1_5p2_Ch1_3p3/'+\
                     'BGcoef/window_out0p01/WindowAdjust_noiselv5_cycleSearch/'+\
                     'atmch0/'+fov+'/'
        input_geo='/mars/jzhou/data/BG_AdaptiveWindow_output/'+\
                     'Remap_ATMS_ARTSprj_GeoCrr_Ch1_5p2_Ch1_3p3/'+\
                     'BGcoef/window_out0p01/WindowAdjust_noiselv5_cycleSearch/'+\
                     'atmch0/'+fov+'/'
        save_dir='/mars/jzhou/data/BG_AdaptiveWindow_output/compare_alpha/'        
        Path(save_dir).mkdir(parents=True,exist_ok=True)

        dB_str=str(dB).replace('-','neg')+'dB'
        with open(input_nogeo+'Alpha_'+dB_str+'.pickle','rb') as f:
            [alpha_nogeo]=pickle.load(f)
        with open(input_geo+'Alpha_'+dB_str+'.pickle','rb') as f:
            [alpha_geo]=pickle.load(f)

        font_jun={'family':'Times New Roman','size':40}
        mpl.rc('font',**font_jun) 
        fig=plt.figure(figsize=(20,10))
        ax=fig.add_subplot(111)
        x=np.arange(len(alpha_nogeo))
        nx_str1='nogeo nx='+str(len(x))
        #ax.plot(x,alpha_nogeo,'-k',linewidth=2,marker='o',markersize=15,zorder=12,label='NoGeoCrr')
        ax.plot(x,alpha_nogeo,'-k',linewidth=2,zorder=12,label='NoGeoCrr')
        x=np.arange(len(alpha_geo))
        nx_str2='geo nx='+str(len(x))
        #ax.plot(x,alpha_geo,'-r',linewidth=2,marker='o',markersize=15,zorder=12,label='GeoCrr')
        ax.plot(x,alpha_geo,'-r',linewidth=1,zorder=12,label='GeoCrr')
        ax.set_xlabel('window index \n'+nx_str1+'\n'+nx_str2)
        ax.set_ylabel('alpha')
        ax.legend()
        plt.savefig(save_dir+'Alpha_fov'+str(self.ifr)+'_'+dB_str+'.png',dpi=300,bbox_inches='tight')
        plt.close()
       
    
