# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import matplotlib.pyplot as plt


dset=np.loadtxt('0.1 mM D6-D10 980 ms IT 1st Experiment 0072_1311057U1.txt',skiprows=8);
centre = 550;
width = 250;

def window_hanning(x_data,y_data,centre,width):
    y = np.hanning(0.05*len(x_data))
    cen=y.argmax()
    x = np.arange(-cen,cen+1,1)
    ybuf = np.ones(width)
    yybuf = np.concatenate((y[0:cen],ybuf,y[cen:-1]))
    xxbuf=np.arange(-(cen+width/2),(cen+width/2),1 )
    yy = np.interp(dset[:,0],xxbuf+centre,yybuf) 
    red_ydata = yy*y_data    
    
    plt.figure()
    plt.plot(x_data,y_data,'.')
    plt.plot(x_data,yy)
    plt.plot(x_data,red_ydata,'.')
    plt.show()
    return red_ydata


def window_kaiser(x_data,y_data,centre,width,beta=2.5):    
    y = np.kaiser(width,beta)
    y=np.insert(y,0,0); y[-1]=0;
    padindex=1
    x = np.arange(-len(y)/2,(len(y)/2)+padindex,1)
    ybuf=np.pad(y, (0, padindex), 'constant')
    yy = np.interp(x_data,x+centre,ybuf) 
    red_ydata = yy*y_data
    
    plt.figure()
    #plt.plot(y)
    plt.plot(x_data,y_data,'.')
    plt.plot(x_data,yy)
    plt.plot(x_data,red_ydata,'.')
    plt.show()
    return red_ydata 
    
#ydone = window_kaiser(dset[:,0],dset[:,4],centre,width,3)

#ydone=window_hanning(dset[:,0],dset[:,4],centre,width)