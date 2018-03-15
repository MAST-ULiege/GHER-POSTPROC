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

firstyear=1
lastyear =4


for yy in range(firstyear,lastyear+1):
    G1     = G3D_class.G3D('../Out_CARTseq/r'+str(yy)+'.nc')
    maskDS= (G1.bat<120) | (G1.bat.mask)

    tt,zp  = G1.avgprofile('TEM'  , ztab=ztab,maskin=maskDS.squeeze())
    a1,zp  = G1.avgprofile('T1age', ztab=ztab,maskin=maskDS.squeeze())
    ss,zp  = G1.avgprofile('SAL', ztab=ztab,maskin=maskDS.squeeze())
    if yy==firstyear:
        tta=tt
        a1a=a1
        ssa=ss
        datesa=G1.dates
    else:
        tta=ma.append(tta,tt,0)
        a1a=ma.append(a1a,a1,0)
        ssa=ma.append(ssa,ss,0)
        datesa=datesa+G1.dates


locator = mdates.AutoDateLocator()
formator = mdates.AutoDateFormatter(locator)

####################                                                                                                                                                                                               
# 1st figure :  
####################                                                                                                                                                                                               
fig=plt.figure(figsize=(15, 8))
Tlevels=np.linspace(6,26,20)

ax=plt.subplot(1, 1, 1)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.contourf(datesa,zp,a1a.T,levels=[0,5,10,20,50,80,120,140,160,180,200,220,240,260,300,350,400,450,500,800],cmap="RdPu")
plt.colorbar()
plt.contour(datesa,zp,tta.T,levels=[8.35,15],linewidth=2)
plt.contour(datesa,zp,ssa.T,levels=[20,20.2],color='k')
#plt.contour(datesa,zp,tta.T,levels=Tlevels,linestyles='--')
#plt.contour(datesa,zp,a1a.T,levels=np.linspace(0,500,10))
#plt.colorbar()
plt.title('Temp for H>120')
plt.ylabel('depth - [m]')
plt.ylim(-300,0)
plt.grid(True)
fig.savefig(G1.figoutputdir+'TEM_AGE_avgZ_SEQ.png')



####################                                                                                 
# 2nd figure : mean age at different depths
#################### 

fig=plt.figure(figsize=(15, 8))
DepthsforPlot=[30,50,70,100]
ax=plt.subplot(1, 1, 1)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.plot(datesa,a1a[:, [np.where(ztab==(-l))[0][0] for l in DepthsforPlot ]   ])
plt.title('Mean Age at various depth')
plt.ylabel('Age - [d]')
plt.grid(True)
plt.legend(labels=[str(d) for d in DepthsforPlot], ncol=4, loc='upper center',
           bbox_to_anchor=[0.5, 1.1],
           columnspacing=1.0, labelspacing=0.0,
           handletextpad=0.0, handlelength=1.5,
           fancybox=True, shadow=True)

fig.savefig(G1.figoutputdir+'AgeAtDepth.png')
