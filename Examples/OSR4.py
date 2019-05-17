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

firstyear=1992
lastyear =1996

# Compile the complite Age profile time series
for yy in range(firstyear,lastyear+1):
    print(yy)
    G1     = G3D_class.G3D('/home/ulg/mast/lvandenb/bsmfc/GHER_5km_120/Out/RAN_V4/'+str(yy)+'.nc', 'local_GHER_V4.yml')
    maskDS = (G1.bat>120) | (G1.bat.mask)

    G1.testvar('O2bottom')
#    G1.testvar('pO2sat', k='bottom')
    del G1.DOX

    G1.O2bottom=G1.maskvar(G1.O2bottom,maskDS.squeeze())
    
    if yy==firstyear:
        G=G1
    else:
        G.dates     = ma.append(G.dates,    G1.dates   ,0)
        G.O2bottom  = ma.append(G.O2bottom, G1.O2bottom,0)
#        G.pO2satkbottom  = ma.append(G.pO2satkbottom, G1.pO2satkbottom,0)
    del G1


G.makeclim('O2bottom')
G.mapMonthlyClim('O2bottom',figsuffix='SHELF',cmapname='oxy', subdomain="NWS", Clim=[0,300])

#G.makeclim('pO2satkbottom')
#G.mapMonthlyClim('pO2satkbottom',figsuffix='SHELF', subdomain="NWS",cmapname="balance",Clim=[0,200])    

#####################
# Fig 1 Time series

G.HYPOXbottom=ma.where(condition=G.O2bottom<60,x=1.0,y=0.0)

G.H=G.HYPOXbottom.sum(axis=(1,2,3))

G.H=G.H*G.dx*G.dy/1e6
G.plotseries(varname='H',title = 'Area for [O2]<60 uM - [km2]')
   

#################################
# Fig 2
# One line per year, are VS DoY
    
fig = plt.figure(figsize=(8*1.618,8))
#locator = mdates.AutoDateLocator()
    
coll = cmocean.cm.ice_r((float(maxyear)- np.array(range(minyear,maxyear))  )/(float(maxyear)-float(minyear)))

for yi,yy in enumerate(range(minyear,maxyear)):
    # Select the good part of datesa and bmaa
    indx     = [i for i,x in enumerate(G.dates) if x.year==yy]
    locdates = [G.dates[i] for i in indx]
    plt.plot(G.H[indx]/1e3, label=yy, color=coll[yi])
    plt.ylabel('Hypoxic Area - [10^3 km2]')
    plt.xlabel('Day of the year - [d]')
    plt.xlim([0,366])

    plt.legend()

## TODO - ADD CLIMATOLOGY 

fig.savefig(G.figoutputdir+'Harea_peryear.png')

################################################
## COMPUTE integrals (Capet et al. 2013)

minyear = np.min([ i.year for i in G.dates])
maxyear = np.max([ i.year for i in G.dates])

maxApyear = np.zeros(maxyear-minyear+1)
Hperyear  = np.zeros(maxyear-minyear+1)
Hnormyear  = np.zeros(maxyear-minyear+1)
Dperyear  = np.zeros(maxyear-minyear+1)

for yy in range(minyear,maxyear+1):
    indx     = [i for i,x in enumerate(G.dates) if x.year==yy]
    Hperyear[yy-minyear]  = np.sum(G.H[indx],0)  # This is in km2*d !! assuming dt is 1 day !! 
    maxApyear[yy-minyear] = np.max(G.H[indx],0)

Dperyear     = Hperyear/maxApyear
Hnormperyear = Hperyear/Dperyear.mean()

###################
# FIGURE 3 SCATTER

plt.figure(figsize=(8*1.618,8))
fig,ax=plt.subplots()
## Make contour plot of the H-index relationship                                                                                                                                                           
# .. grid H and D from min and max +-std                                                                                                                                                                   
# compute corresponding Hindex                                                                                                                                                                             
# plt.contourf                                            

xlist = np.linspace(40, 130,50)
ylist = np.linspace(1, 18, 50)
X, Y = np.meshgrid(xlist, ylist)

Z= Y*X/Dperyear.mean()
#Z=np.divide(X,Y)

cs=plt.contourf(X,Y,Z, cmap= cmocean.cm.thermal)

bbox = dict(boxstyle="round", fc="0.8")

def bboxfy(y):
    bbox = dict(boxstyle="round", fc=cmocean.cm.ice_r((float(maxyear)-float(y))/(float(maxyear)-float(minyear))))
    return(bbox)

#ax.scatter(Hperyear, Dperyear, c=range(minyear,maxyear+1),  edgecolors='none', s=200, cmap=cmocean.cm.ice)
for i, yr in enumerate(   range(minyear,maxyear+1)  ):
    ax.annotate(str(yr)[2:4], (Dperyear[i], Hnormperyear[i]/1e3), bbox=bboxfy(yr), color='white', horizontalalignment='center', verticalalignment='center')

plt.xlabel('Hypoxia Duration - [d]')
plt.ylabel('Max. Hypoxic Area - [10^3 km2]')
cbaxes = fig.add_axes([0.25, 0.94, 0.5, 0.02]) 
clb=plt.colorbar(cs,  cax=cbaxes, orientation="horizontal")#, position='top') 
clb.ax.set_title('Hypox-Index - [10^3 km2]')

fig.savefig(G.figoutputdir+'H_scatter.png')


####################################################
import matplotlib as mpl

mpl.style.use('ggplot')
plt.rcParams['axes.facecolor']='w'
plt.rcParams['grid.color']='lightgrey'#   b0b0b0    # grid color                                                                                                                                            
plt.rcParams['grid.linestyle']='--'         # solid                                                                                                                                                         
plt.rcParams['axes.edgecolor']='black'

fig=plt.figure(figsize=(8*1.618,8))
ax=plt.subplot(1, 1, 1)
#ax.xaxis_date()
#locator  = mdates.AutoDateLocator()
#ax.xaxis.set_major_formatter(locator)

annual_dates= [dt.datetime(y,6,1) for y in range(minyear,maxyear+1)]

cs=plt.plot(annual_dates, Hnormperyear/1e3, color = 'lightseagreen',linewidth = 2, marker='o',markersize = 10)
plt.title('Annual BH index')
fig.savefig(G.figoutputdir+'BH_Annual.png')

fit1,cov        = np.polyfit( range(minyear,maxyear)  , Hnormperyear/1e3,1, cov=True)
trend1 = fit1[0]
stderroronslope = np.sqrt(np.diag(cov))[0]

print( 'Annual Trend is %s +/- %s'%(trend1, stderroronslope ) 

'''
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
