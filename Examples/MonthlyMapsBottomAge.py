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
import cmocean 
import G3D_class
import calendar

firstyear=1
lastyear =5

# Compile the complite Age profile time series
for yy in range(firstyear,lastyear+1):
    G1     = G3D_class.G3D('../../data/postfpCARTseq/r'+str(yy)+'.nc')
    maskDS = (G1.bat>120) | (G1.bat.mask)
    # a1 has two dimension : ztab and original time.
    G1.gload('T1age')
    G1.T1age=G1.maskvar(G1.T1age,maskDS.squeeze())
    
    botmaxage = G1.T1age[:,11]
    
    if yy==firstyear:
        bmaa=botmaxage
        datesa=G1.dates
    else:
        bmaa=ma.append(bmaa,botmaxage,0)
        datesa=datesa+G1.dates
    if not(yy==lastyear):
        del G1
    else:
        G1.testz()
        
minyear = np.min([ i.year for i in datesa])
maxyear = np.max([ i.year for i in datesa])
        
for yy in range(minyear,maxyear):
    # Select the good part of datesa and bmaa
    indx     = [i for i,x in enumerate(datesa) if x.year==yy]
    locdates = [datesa[i] for i in indx]
    loca     = bmaa[indx]
    
    # Maps of max bottom age
    fig = plt.figure(figsize=(30,15))
    ax  = plt.subplot(1, 2, 1)
    mdd = plt.contourf(G1.lon[1:150],G1.lat,loca.max(axis=0)[:,1:150],cmap=cmocean.cm.deep, levels=np.linspace(0,500,20))
    plt.colorbar(mdd)
    plt.title('Max Age of Bottom Waters : ' + str(yy))
    
    # Maps of when max age occurs (days of the year) 
    idxma = loca.argmax(0)
    dayl  = [ i.timetuple().tm_yday for i in locdates ]
    dayd  = dict(enumerate(dayl))
    
    b=ma.copy(idxma)
    for old, new in dayd.items():
        b[idxma == old] = new

    ax=plt.subplot(1, 2, 2)
    plt.pcolor(G1.lon[1:150],G1.lat,ma.masked_where(loca[0].mask,b)[:,1:150],cmap=cmocean.cm.phase, vmin=0, vmax=366)
    clb=plt.colorbar(ticks=np.linspace(15,350,12).round())
    clb.ax.set_yticklabels(calendar.month_abbr[1:13]) 
    plt.title('Occurence')
    fig.savefig(G1.figoutputdir+'BottomMaxAge'+str(yy)+'.png')
    
    
    



'''
ztab   = -1*np.concatenate([np.arange(0,10,2), np.arange(10,40,5),np.arange(50,120,10),np.arange(120,320,50)])

firstyear=1
lastyear =1

yy=1

G1     = G3D_class.G3D('../Out_CARTseq/r'+str(yy)+'.nc')
maskDS = (G1.bat<120) | (G1.bat.mask)
maskDS3D=np.tile(maskDS,(31,1,1))

G1.gload('T1age')
G1.testz()

months = [g.month for g in G1.dates]
year   = [g.month for g in G1.year]

for yy in [1982]:
    for mm in range(1,13):
        tindexes = np.where((np.array(months)==4 )& (np.array(year)==1982))[0]
        locbotage=ma.mean(G1.T1age[tindexes,11],axis=0)

        fig=plt.figure(figsize=(6, 6))
        ax=plt.subplot(1, 1, 1)                                                                                                              
        plt.contourf(G1.lon,G1.lat,locbotage,cmap='RdPu')
'''
