# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 12:19:08 2021

@author: Leo.Ji
"""

import numpy as np
#import pandas as pd
import os
import pickle

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


class Ellipse_type():
    def __init__(self, ):
        self.f_decimal = np.float32    # numerical accuracy
                

        return
    
    
    def fit_ellipse(self, data):
            
        points = np.vstack((data[0], data[1])).astype(self.f_decimal)
        
        res = self.get_ellipse(points)
        center, param_axis, param_orien = self.check_ellipse(res)    
        self.test_plot_Ellipse( (center, param_axis, param_orien), data )

        return center, param_axis, param_orien
    
    
    def get_ellipse(self, data):
        # https://stackoverflow.com/questions/17952171/not-sure-how-to-fit-data-with-a-gaussian-python
        center = data.mean(axis=1)
#        if True in (center>180): center[center>180] -= 360.
        
#        sigma = np.cov(data, rowvar=0)
        sigma = np.cov(data)
        
        # compute eigenvalues and associated eigenvectors
        vals, vecs = np.linalg.eigh(sigma)
        
        # eigenvalues give length of ellipse along each eigenvector
        major, minor = 4 * np.sqrt(vals)  #### this is the whole lengths of axises, not the half lengths
        
        # compute "tilt" of ellipse using first eigenvector
        x, y = vecs[:, 0]
        orientation = np.arctan2(y, x)
        
        return [center, [major, minor], orientation]
    
    
    def check_ellipse(self, ellipses):
        
#        param_axises, param_oriens = np.array(ellipses[1]).T.astype(self.f_decimal), np.array(ellipses[2]).T.astype(self.f_decimal)
        if ellipses[1][0] < ellipses[1][1]:
            ellipses[1][0], ellipses[1][1] = ellipses[1][1].astype(self.f_decimal), ellipses[1][0].astype(self.f_decimal)   # major axis should be larger than minor axis
            ellipses[2] -= np.sign(ellipses[2])*0.5*np.pi    # at this time, orientation should be turn 0.5*np.pi
                        
        ellipses[2] = self.constraint_orientation(ellipses[2])   # constraint orientation in [-0.5*np.pi, 0.5*np.pi]
        return ellipses  # centers, param_axises, param_oriens
    

    def constraint_orientation(self, oriens):
        if np.abs(oriens) >= 2*np.pi:  oriens %= 2*np.pi    # [-2*np.pi , 2*np.pi ]
        if np.abs(oriens) >= 1.5*np.pi: oriens -= np.sign(oriens)*2*np.pi     # [-np.pi, np.pi]
        if np.abs(oriens) >= 0.5*np.pi: oriens -= np.sign(oriens)*np.pi     # [-np.pi, np.pi]
        
        return np.degrees(oriens).astype(self.f_decimal)


    def test_plot_Ellipse(self, param_ellipse, points):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
#        ell1 = Ellipse(xy = param_ellipse[0], width = max(param_ellipse[1]), height = max(min(param_ellipse[1]),min(0.2,max(param_ellipse[1]))), 
#                       angle = param_ellipse[2], facecolor= 'yellow', alpha=1)
        width, height = param_ellipse[1]
        
        ell1 = Ellipse(xy = param_ellipse[0], width = width, height = height, angle = param_ellipse[2],
                       color='blue', fill=False, linewidth=0.5, linestyle='--', alpha=1)
        ax.add_patch(ell1)
        ax.scatter(points[0], points[1], marker='.', c=points[2], s=20, alpha=1, zorder=2)
        
        plt.axis('scaled')        
        plt.axis('equal')   #changes limits of x or y axis so that equal increments of x and y have the same length        
        plt.show()
        plt.savefig('/mars/jzhou/test_ellipse.png',dpi=300,bbox_inches='tight')
        return

    
    def test_scatter(self, data):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(data[0], data[1], marker='.', c=data[2], s=20, alpha=1, zorder=2)
        plt.show()
        
        return
