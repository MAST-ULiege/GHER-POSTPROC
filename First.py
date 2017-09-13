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

G=G3D_class.G3D('../data/CART1CLIP/1980.nc')
t1,zforplot1=G.avgprofile('T1age',ztab=-np.array([0,2,5,10,20,50]))
t2,zforplot2=G.avgprofile('T1age')

dates = [dt.datetime(1858,11,17)+dt.timedelta(days=int(t)) for t in G.time]
fig=plt.figure(figsize=(10, 8))

####################                                                                                                                                                                                                
# 1st figure : Age profile                                                                                                                                                                                        
####################                                                                                                                                                                                               

ax=plt.subplot(2, 2, 1)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
plt.pcolor(dates,zforplot1,t1.T)
plt.title('T1age')
plt.ylabel('depth - [m]')
plt.clim([0,len(G.time)])
plt.ylim([-500,0])
plt.colorbar()

ax=plt.subplot(2, 2, 2)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
plt.pcolor(dates,zforplot1,(1-t1.T/len(G.time)))
plt.title('Ventilation Rate')
plt.ylabel('depth - [m]')
plt.clim([0,1])
plt.ylim([-500,0])
plt.colorbar()

ax=plt.subplot(2, 2, 3)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
plt.pcolor(dates,zforplot2,t2.T)
plt.title('T1age')
plt.ylabel('depth - [m]')
plt.clim([0,len(G.time)])
plt.ylim([-500,0])
plt.colorbar()

ax=plt.subplot(2, 2, 4)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
plt.pcolor(dates,zforplot2,(1-t2.T/len(G.time)))
plt.title('Ventilation Rate')
plt.ylabel('depth - [m]')
plt.ylim([-500,0])
plt.clim([0,1])
plt.colorbar()

fig.savefig('../MeanAgeProfile.png')
