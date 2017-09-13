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
t1,zforplot1=G.avgprofile('TEM',ztab=-np.array([0,10,20,30,50,100,200,300,400,500,700,1000]))
t2,zforplot2=G.avgprofile('TEM')

dates = [dt.datetime(1858,11,17)+dt.timedelta(days=int(t)) for t in G.time]
fig=plt.figure(figsize=(10, 8))

####################                                                                                                                                                                                                
# 1st figure : Age profile                                                                                                                                                                                        
####################                                                                                                                                                                                               

ax=plt.subplot(1, 2, 1)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
plt.pcolor(dates,zforplot1,t1.T)
plt.title('Temp')
plt.ylabel('depth - [m]')
plt.clim([7,12])
plt.ylim([-500,0])
plt.colorbar()

ax=plt.subplot(1, 2, 2)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
plt.pcolor(dates,zforplot2,t2.T)
plt.title('Temp')
plt.ylabel('depth - [m]')
plt.clim([7,12])
plt.ylim([-500,0])
plt.colorbar()


fig.savefig('../MeanTemProfile.png')
