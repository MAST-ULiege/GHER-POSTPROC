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

maskDS= (G.bat<120) | (G.bat.mask)  # Mask should be True where masked
maskSH= (G.bat>120) | (G.bat.mask)  # Mask should be True where masked

#tDS,zforplotDS=G.avgprofile('T1age',maskin=maskDS.squeeze())
tDS=G.avgprofileSIGMA('T1age',maskin=maskDS.squeeze())


#tSH,zforplotSH=G.avgprofile('T1age',maskin=maskSH.squeeze())
tSH=G.avgprofileSIGMA('T1age',maskin=maskSH.squeeze())

####################   
# 1st figure : Age profile                                                                                                                                                                                        
####################                                                                                                                                                                                               


fig=plt.figure(figsize=(10, 8))

ax=plt.subplot(2, 2, 1)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
#plt.pcolor(G.dates,zforplotSH,tSH.T)
plt.pcolor(tSH.T)
plt.title('T1age SH')
plt.ylabel('depth - [m]')
plt.clim([0,len(G.time)])
plt.ylim([-500,0])
plt.colorbar()

ax=plt.subplot(2, 2, 2)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
#plt.pcolor(G.dates,zforplotSH,(1-tSH.T/len(G.time)))
plt.pcolor((1-tSH.T/len(G.time)))
plt.title('Ventilation Rate SH')
plt.ylabel('depth - [m]')
plt.clim([0,1])
plt.ylim([-500,0])
plt.colorbar()


ax=plt.subplot(2, 2, 3)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
#plt.pcolor(G.dates,zforplotDS,tDS.T)
plt.pcolor(tDS.T)
plt.title('T1age DS')
plt.ylabel('depth - [m]')
plt.clim([0,len(G.time)])
plt.ylim([-500,0])
plt.colorbar()

ax=plt.subplot(2, 2, 4)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
#plt.pcolor(G.dates,zforplotDS,(1-tDS.T/len(G.time)))
plt.pcolor((1-tDS.T/len(G.time)))
plt.title('Ventilation Rate DS')
plt.ylabel('depth - [m]')
plt.clim([0,1])
plt.ylim([-500,0])
plt.colorbar()


fig.savefig(G.figoutputdir+'MeanAgeProfile.png')
