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
lastyear =4
UpperVoldDepth=150


for yy in range(firstyear,lastyear+1):
    G1     = G3D_class.G3D('../Out_CARTseq/r'+str(yy)+'.nc')
    maskDS = (G1.bat<120) | (G1.bat.mask)
    maskDS3D=np.tile(maskDS,(31,1,1))

    G1.gload('T1age')
    G1.testz()

    NAgeClasses = 100
    AgeClasses  = np.linspace(0,500,NAgeClasses )
    AgeClasses  = [0,5,25,50,100,150,300,500]
    AgeVolumes  = np.zeros([len(G1.dates),len(AgeClasses)])
    AgeVolumes2  = np.zeros([len(G1.dates),len(AgeClasses)])

    Vol = G1.dx*G1.dy*G1.dz*1e-9
    UpperVol = ma.masked_where( maskDS3D| (G1.z<-UpperVoldDepth),Vol)
    
    daysincebeg=np.zeros(len(G1.dates))

    if yy==firstyear:
        datebeg=G1.dates[0]

    for t in range(len(G1.dates)):
# make a vector with the volume of water For each age class
        localagevector = ma.masked_where(maskDS3D|(G1.z<-UpperVoldDepth),G1.T1age[t])
    
        for ageClassindex in range(len(AgeClasses)-1):
            bi  =  ma.masked_where( (localagevector<AgeClasses[ageClassindex+1]), Vol)
            bi2 = ma.masked_where( (localagevector>=AgeClasses[ageClassindex+1]), Vol)

            AgeVolumes[t,ageClassindex]=bi.sum()
            AgeVolumes2[t,ageClassindex]=bi2.sum()

        daysincebeg[t]=(G1.dates[t]-datebeg).days

    if yy==firstyear:
        AVa=AgeVolumes
        AVa2=AgeVolumes2
#        datesa=daysincebeg
        datesa=G1.dates
    else:
        AVa=np.append(AVa,AgeVolumes,0)
        AVa2=np.append(AVa2,AgeVolumes2,0)
#        datesa=np.append(datesa,daysincebeg,0)
        datesa=datesa+G1.dates

locator = mdates.AutoDateLocator()
formator = mdates.AutoDateFormatter(locator)

AVa=AVa/UpperVol.sum()*100
AVa2=AVa2/UpperVol.sum()*100

####################                                                                                                                                                                                               
# 1st figure :  
####################                                       



                                                                                                                                                        
fig=plt.figure(figsize=(15, 15))
ax=plt.subplot(3, 1, 1)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.plot(datesa,AVa)
plt.legend(labels=["> "+str(d)+" d" for d in AgeClasses[1:]], ncol=4, loc='upper center',
           bbox_to_anchor=[0.5, 1.1],
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)

plt.title('Volumes of Waters older than .. - [% of volume]')
plt.ylabel('% of Upper Volume - [%]')
plt.xlabel('Time- [d]')
plt.grid(True)

ax=plt.subplot(3, 1, 2)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.plot(datesa,AVa2)
plt.legend(labels=["< "+str(d)+" d" for d in AgeClasses[1:]], ncol=4, loc='upper center',
           bbox_to_anchor=[0.5, 1.1],
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)

plt.title('Volumes of Waters younger than .. - [% of volume]')
plt.ylabel('% of Upper Volume - [%]')
plt.xlabel('Time- [d]')
plt.grid(True)

ax=plt.subplot(3, 1, 3)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.plot(datesa,AVa+AVa2)
plt.legend(labels=[">< "+str(d)+" d" for d in AgeClasses[1:]], ncol=4, loc='upper center',
           bbox_to_anchor=[0.5, 1.1],
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)

plt.title('Volumes of Waters older than .. - [% of volume]')
plt.ylabel('% of Upper Volume - [%]')
plt.xlabel('Time- [d]')
plt.grid(True)
fig.savefig(G1.figoutputdir+'AgeVolumesPercCheck.png')
