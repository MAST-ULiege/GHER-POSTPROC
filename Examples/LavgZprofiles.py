import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

#from mpl_toolkits.basemap import Basemap  
#from multiprocessing import Pool
#import gsw                                                                                                                                                                                                         
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import sys
import os

import G3D_class

ztab=-1*np.concatenate([np.arange(0,10,2), np.arange(10,40,5),np.arange(50,120,10),np.arange(120,320,50)])
firstyear=1980
lastyear =1981
tempstep=4
colormap = plt.cm.hsv





for yy in range(firstyear,lastyear):
    G1     = G3D_class.G3D('../Out_CART/'+str(yy)+'.nc')
    maskDS= (G1.bat<120) | (G1.bat.mask)

    a1,zp  = G1.avgprofile('T1age', ztab=ztab,maskin=maskDS.squeeze())

    fig=plt.figure(figsize=(15,15))
    ax=plt.subplot(1, 1, 1)
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 1, a1.shape[0]/tempstep)])
    labels = []
 
    for d in range(0,a1.shape[0], tempstep):
        plt.plot(a1[d],zp)
        labels.append(G1.dates[d].strftime("%d %B")) 
 
    plt.legend(labels, ncol=6, loc='upper center', 
               bbox_to_anchor=[0.5, 1.1], 
               columnspacing=1.0, labelspacing=0.0,
               handletextpad=0.0, handlelength=1.5,
               fancybox=True, shadow=True)
    plt.ylabel('depth - [m]')
    plt.ylim(-200,0)
    plt.grid(True)
    plt.xlabel('age - [d]')
    plt.xlim(0,1000)
    fig.savefig(G1.figoutputdir+'feather'+str(yy)+'.png')
