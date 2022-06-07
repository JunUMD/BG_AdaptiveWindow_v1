# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:23:32 2020

@author: jzhou128
"""

from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt

class FitCircle_class:

    def __init__(self,x,y):
        self.x=x
        self.y=y
        
    def calc_R(self,xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((self.x-xc)**2 + (self.y-yc)**2)

    def f_2(self,c):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        Ri = self.calc_R(*c)
        return Ri - Ri.mean()
    
    def fit(self):
        x_m = np.mean(self.x)
        y_m = np.mean(self.y)
        center_estimate = x_m, y_m
        center_2, ier = optimize.leastsq(self.f_2, center_estimate)

        self.xc_2, self.yc_2 = center_2
        Ri_2       = self.calc_R(*center_2)
        self.R_2        = Ri_2.mean()
        residu_2   = sum((Ri_2 - self.R_2)**2)
        
        theta_fit = np.linspace(-np.pi, np.pi, 360)
        self.x_fit2 = self.xc_2 + self.R_2 * np.cos(theta_fit)
        self.y_fit2 = self.yc_2 + self.R_2 * np.sin(theta_fit)
        
        return self.xc_2, self.yc_2, self.R_2,self.x_fit2,self.y_fit2
    
    def display(self):
        
        fig=plt.figure(figsize=(10,10))
        ax=fig.add_subplot(111)
        ax.scatter(self.x,self.y,c='k',s=20,marker='X')
        ax.scatter(self.x_fit2,self.y_fit2,c='r',s=5)
        plt.show()
        plt.close()


if __name__=="__main__":      
    x=[-1,0,0,1]
    y=[0,-1,1,0]
    fitcircle_obj=FitCircle_class(x,y)
    fitcircle_obj.fit()
    fitcircle_obj.display()
