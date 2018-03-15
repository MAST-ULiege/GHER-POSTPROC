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

ztab   = -1*np.concatenate([np.arange(0,10,2), np.arange(10,40,5),np.arange(50,120,10),np.arange(120,320,50)])

firstyear=1
lastyear =1

yy=1

G1     = G3D_class.G3D('../Out_CARTseq/r'+str(yy)+'.nc')
maskDS = (G1.bat<120) | (G1.bat.mask)
maskDS3D=np.tile(maskDS,(31,1,1))

G1.gload('T1age')
G1.testz()

months = [g.month for g in G1.dates]
year   = [g.month for g in G1.year]

for yy in [1982]:
    for mm in range(1,13):
        tindexes = np.where((np.array(months)==4 )& (np.array(year)==1982))[0]
        locbotage=ma.mean(G1.T1age[tindexes,11],axis=0)

        fig=plt.figure(figsize=(6, 6))
        ax=plt.subplot(1, 1, 1)                                                                                                              
        plt.contourf(G1.lon,G1.lat,locbotage,cmap='RdPu')
