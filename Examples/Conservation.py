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

for yy in range(firstyear,lastyear+1):
    G1     = G3D_class.G3D('../Out_CARTseq/r'+str(yy)+'.nc')
    maskDS = (G1.bat<120) | (G1.bat.mask)

    G1.gload('T1age')
    G1.gstore('T1age')
    G1.testz()

    NAgeClasses = 100
    AgeClasses  = np.linspace(0,1000,NAgeClasses )
    AgeVolumes  = np.zeros([len(G1.dates),NAgeClasses])

    Vol = G1.dx*G1.dy*G1.dz*1e-9

    daysincebeg=np.zeros(len(G1.dates))

    if yy==firstyear:
        datebeg=G1.dates[0]

    for t in range(len(G1.dates)):
# make a vector with the volume of water For each age class
        localagevector = G1.T1age[t]
    
        for ageClassindex in range(len(AgeClasses)-1):
            bi = ma.masked_where( (localagevector<AgeClasses[ageClassindex]) | (localagevector>=AgeClasses[ageClassindex+1]), Vol)
            AgeVolumes[t,ageClassindex]=bi.sum()

        daysincebeg[t]=(G1.dates[t]-datebeg).days

    if yy==firstyear:
        AVa=AgeVolumes
        datesa=daysincebeg
    else:
        AVa=np.append(AVa,AgeVolumes,0)
        datesa=np.append(datesa,daysincebeg,0)

locator = mdates.AutoDateLocator()
formator = mdates.AutoDateFormatter(locator)

AVa=AVa/Vol.sum()*100

####################                                                                                                                                                                                               
# 1st figure :  
####################                                                                                                                                                                                               
fig=plt.figure(figsize=(15, 15))

ax=plt.subplot(1, 1, 1)
#ax.xaxis_date()
#ax.xaxis.set_major_locator(locator)
#ax.xaxis.set_major_formatter(formator)
#plt.contourf(datesa, AgeClasses, AVa.T,levels=np.linspace(0,10,100),cmap='GnBu')
plt.contourf(datesa, AgeClasses, AVa.T,levels=np.linspace(0,1.5,100),cmap='gist_ncar_r')
plt.colorbar()
plt.plot([0.0, datesa.max()], [0.0, datesa.max()], 'r-', lw=2) 
plt.title('Volumes for age of Waters - [% of volume]')
plt.ylabel('Age - [d]')
plt.xlabel('Time- [d]')
plt.grid(True)
fig.savefig(G1.figoutputdir+'AgeVolumes.png')
