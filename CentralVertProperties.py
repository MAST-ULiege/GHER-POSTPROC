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

yearstart=1980
yearend=2000

for years in range(yearstart,yearend):
    G     = G3D_class.G3D('../Out_002/'+str(years)+'.nc')
    tt,z  = G.profileatxy('TEM',120,50)
    ss,z  = G.profileatxy('SAL',120,50)
    dd,z  = G.profileatxy('DEN',120,50)
    dates = [dt.datetime(1858,11,17)+dt.timedelta(days=int(t)) for t in G.time]
    if (years==yearstart):
        ttall    = tt
        ssall    = ss
        ddall    = dd
        datesall = dates
    else:
        ttall    = ma.concatenate([ttall,tt])
        ssall   = ma.concatenate([ssall,ss])
        ddall    = ma.concatenate([ddall,dd])
        datesall = datesall+dates

lonlab,latlab=G.getlonlat(120,50)


locator = mdates.AutoDateLocator()
formator = mdates.AutoDateFormatter(locator)

####################                                                                                                                                                                                               
# 1st figure :  
####################                                                                                                                                                                                               

fig=plt.figure(figsize=(15, 8))

ax=plt.subplot(1, 1, 1)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.contourf(datesall,z,ttall.T,25)
plt.colorbar()
plt.contour(datesall,z,ttall.T,levels=[8.35],color='black')
plt.title('Temp @ Lon : %s E; Lat : % sN'%(lonlab,latlab))
plt.ylabel('depth - [m]')
plt.ylim([-500,0])
fig.savefig(G.figoutputdir+'TEMatXY.png')

####################                                                                                                                                                                                               
# 1st figure :
####################                                                                                                                                                                                               

fig=plt.figure(figsize=(15, 8))

ax=plt.subplot(1, 1, 1)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.contourf(datesall,z,ssall.T,25)
plt.title('Salinity @ Lon : %s E; Lat : % sN'%(lonlab,latlab))
plt.ylabel('depth - [m]')
plt.ylim([-500,0])
plt.colorbar()
fig.savefig(G.figoutputdir+'SALatXY.png')

####################                                                                                                                                                                                               
# 3st figure :
####################                                                                                                                                                                                               

fig=plt.figure(figsize=(15, 8))

ax=plt.subplot(1, 1, 1)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.contourf(datesall,z,ddall.T,25)
plt.title('Density @ Lon : %s E; Lat : % sN'%(lonlab,latlab))
plt.ylabel('depth - [m]')
plt.ylim([-500,0])
plt.colorbar()
fig.savefig(G.figoutputdir+'DENatXY.png')





