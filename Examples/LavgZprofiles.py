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
import calendar
import G3D_class

ztab=-1*np.concatenate([np.arange(0,200,5)])
colormap = plt.cm.hsv
firstyear=1
lastyear =5

# Compile the complite Age profile time series
for yy in range(firstyear,lastyear+1):
    G1     = G3D_class.G3D('../../data/postfpCARTseq/r'+str(yy)+'.nc')
    maskDS = (G1.bat<120) | (G1.bat.mask)
    # a1 has two dimension : ztab and original time.
    a1,zp  = G1.avgprofile('T1age', ztab=ztab,maskin=maskDS.squeeze())
    if yy==firstyear:
        a1a=a1
        datesa=G1.dates
    else:
        a1a=np.append(a1a,a1,0)
        datesa=datesa+G1.dates
    if not(yy==lastyear):
        del G1


# Loop over years, 1 plot per year

minyear = np.min([ i.year for i in datesa])
maxyear = np.max([ i.year for i in datesa])

for yy in range(minyear,maxyear):
    # Select the good part of datesa and a1a
    indx     = [i for i,x in enumerate(datesa) if x.year==yy]
    locdates = [datesa[i] for i in indx]
    loca     = a1a[indx]

    fig=plt.figure(figsize=(15,15))
    ax=plt.subplot(1, 1, 1)
    
    labels=[]

    for mm in range(1,13):
        mindx   = [i for i,x in enumerate(locdates) if x.month==mm]
        locprof = loca[mindx].mean(axis=0)
        plt.plot(locprof,zp)
        labels.append( calendar.month_abbr[mm] ) 
 
    plt.legend(labels, loc='right', 
               fancybox=True, shadow=True)
    plt.title(str(yy)+' monthly age profiles')

    plt.ylabel('depth - [m]')
    plt.ylim(-200,0)
    plt.grid(True)
    plt.xlabel('age - [d]')
    plt.xlim(0,1200)
    fig.savefig(G1.figoutputdir+'Bfeather'+str(yy)+'.png')
