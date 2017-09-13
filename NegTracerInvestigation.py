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

# local function definition for this very specific use : 
def instance_negT1(self):
    print('Im in *** ')
    self.gload('T1C')
    print('Zorro is masked, and also T1c %s'%(ma.isMaskedArray(self.T1C)))
    self.negT1=ma.masked_where(self.T1C>=0.0,self.T1C)
    print('all is fine')

#adding the local method to the G3D class
G3D_class.G3D.instance_negT1=instance_negT1

# Let us start now : Instantiate the class

G=G3D_class.G3D('../Out/1980.nc')

#####################################
# 1st figure : Spatially integrated.  
#####################################                                                                                                                                                                                               
intN= G.integratespatial('negT1')
fig=plt.figure(figsize=(15, 8))
ax=plt.subplot(1,1, 1)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
plt.plot(G.dates,intN)
plt.title('integrated negative T1C')
plt.ylabel('')
fig.savefig(G.figoutputdir+'IntegratedNegativeTracer.png')

####################                                                                                                                                                                                               
# 2nd figure : horizontal avg sigma                                                                                                   
####################                                                                                                                                                                                               

vertsig=G.avgprofileSIGMA('negT1')
fig=plt.figure(figsize=(15, 8))
ax=plt.subplot(1,1, 1)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
plt.pcolor(G.dates,np.arange(1,32),vertsig.T)
plt.title('integrated negative T1C')
plt.ylabel('')
plt.colorbar()
fig.savefig(G.figoutputdir+'NTavgsigma.png')

#################
# 3rd vertical integraton 

G.VNEG=G.vertint('negT1')
fig=plt.figure(figsize=(15, 8))
ax=plt.subplot(1,1, 1)
plt.pcolor(G.lon,G.lat,G.VNEG[-1])
plt.title('integrated negative T1C')
plt.colorbar()
fig.savefig(G.figoutputdir+'NTvertint.png')

