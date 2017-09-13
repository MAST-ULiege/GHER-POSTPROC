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

G=G3D_class.G3D('../Out_001/1980.nc')

tt,z = G.profileatxy('TEM',120,50)

print tt.shape

lonlab,latlab=G.getlonlat(120,50)

dates = [dt.datetime(1858,11,17)+dt.timedelta(days=int(t)) for t in G.time]
fig=plt.figure(figsize=(10, 8))

####################                                                                                                                                                                                               
# 1st figure :  
####################                                                                                                                                                                                               

print z.shape
print tt.shape

ax=plt.subplot(1, 1, 1)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
plt.pcolor(dates,z,tt.T)
plt.title('Temp @ Lon : %s E; Lat : % sN'%(lonlab,latlab))
plt.ylabel('depth - [m]')
plt.clim([7,12])
plt.ylim([-500,0])
plt.colorbar()


fig.savefig(G.figoutputdir+'MeanTemProfile.png')
