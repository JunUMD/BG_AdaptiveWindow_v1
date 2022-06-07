from scipy import optimize
import numpy as np

def twoD_Gaussian(data, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x,y)=data
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)+ c*((y-yo)**2)))
    return g.ravel()

#def twoD_Gaussian(data, amplitude, xo, yo, sigma_x, sigma_y, offset):
#    (x,y)=data
#    xo = float(xo)
#    yo = float(yo)
#    g=offset+amplitude*np.exp(-((x-xo)**2/(2*sigma_x**2)+(y-yo)**2/(2*sigma_y**2)))    
#    return g.ravel()


def oneD_Gaussian(x, amplitude, mean, stddev):
    return amplitude*np.exp(-(x - mean)**2/(2*stddev**2))

def twoD_Gaussian_fit(x,y,G):

    initial_guess = (1,0,0,1,1,np.radians(8.7),0)
    #initial_guess = (1,np.mean(x),np.mean(y),1,1,0)

    popt, pcov = optimize.curve_fit(twoD_Gaussian,(x,y),G,p0=initial_guess,maxfev=5000)
    G_fit=twoD_Gaussian((x,y),*popt)

    return [popt,pcov]

def oneD_Gaussian_fit(x,G):
    initial_guess = (1,np.mean(x),1)

    popt, pcov = optimize.curve_fit(oneD_Gaussian,x,G,p0=initial_guess)

    return [popt]

    



