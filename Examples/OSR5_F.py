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
import N3D_class


## ## ## ## ## ###  ## ## ## ## ## ###  ## ## ## ## ## ### 
## Arguments Management ##
import  getopt,sys
 
arglist = sys.argv[1:]
unixOptions = "y:"
gnuOptions = ["yml="]
try:
    arguments, values = getopt.getopt(arglist, unixOptions, gnuOptions)
except getopt.error as err:
    # output error, and return with an error code                                                                                                                                                                  
    print (str(err))
    sys.exit(2)

ymlfile   = None
for a,v  in arguments:
        if a in ("-y",'--yml'):
            ymlfile = v

if ymlfile is None:
    print( 'Dude, I need some -y or --yml argument for the .yml file')
    sys.exit()
## ## ## ## ## ##### ## ## ## ## ###  ## ## ## ## ## ###         

mlist = N3D_class.FullLoad(ymlfile)

for mm in mlist:
    G1  = N3D_class.N3D(mm,ymlfile)
    maskDS = (G1.bat>60) | (G1.bat.mask)

    G1.testvar('O2bottom')
    del G1.DOX

    G1.O2bottom=G1.maskvar(G1.O2bottom,maskDS.squeeze())
    
    if mm==mlist[0]:
        G=G1
    else:
        G.dates     = ma.append(G.dates,    G1.dates   ,0)
        G.O2bottom  = ma.append(G.O2bottom, G1.O2bottom,0)
    del G1

G.timeclean()

G.makeclim('O2bottom')
G.mapMonthlyClim('O2bottom',figsuffix='SHELF',cmapname='oxy', subdomain="NWS", Clim=[0,300])

#####################
# Fig 1 Time series

# HYPOXbottom is ones where bottom O2 is below threshold, zeros elsewhere.
G.HYPOXbottom=ma.where(condition=G.O2bottom<60,x=1.0,y=0.0)
# Spatial integration accounts for grid cell's dimension. 
G.H=G.integratespatial('HYPOXbottom').squeeze()/1e6
G.plotseries('H',title = 'Area for [O2]<60 uM - [km2]')

#################################
# Fig 2
# One line per year, are VS DoY
    
fig = plt.figure(figsize=(8*1.618,8))
#locator = mdates.AutoDateLocator()

minyear = np.min([ i.year for i in G.dates])                                                                                                                                                                       
maxyear = np.max([ i.year for i in G.dates]) 

coll = cmocean.cm.ice_r((float(maxyear)- np.array(range(minyear,maxyear+1))  )/(float(maxyear)-float(minyear)))

for yi,yy in enumerate(range(minyear,maxyear+1)):
    # Select the good part of datesa and bmaa
    indx     = [i for i,x in enumerate(G.dates) if x.year==yy]
    locdates = [G.dates[i] for i in indx]
    locyday  = [ d.timetuple().tm_yday for d in locdates]
    plt.plot( locyday , G.H[indx]/1e3, label=yy, color=coll[yi])
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

# timestepindays = float((G.dates[1]-G.dates[0]).seconds)/86400 # Asumming Constant time steps
timestepindays = float((G.dates[1]-G.dates[0]).days)+float((G.dates[1]-G.dates[0]).seconds)/86400

print('time step is %s days'%timestepindays)

for yy in range(minyear,maxyear+1):
    indx     = [i for i,x in enumerate(G.dates) if x.year==yy]
    Hperyear[yy-minyear]  = np.sum(G.H[indx],0)*timestepindays 
    maxApyear[yy-minyear] = np.max(G.H[indx],0)

Dperyear     = Hperyear/maxApyear
Hnormperyear = Hperyear/Dperyear.mean()


print(Dperyear)
print(Hnormperyear)
###################
# FIGURE 3 SCATTER

plt.figure(figsize=(8*1.618,8))
fig,ax=plt.subplots()
## Make contour plot of the H-index relationship                                                                                                                                                           
# .. grid H and D from min and max +-std                                                                                                                                                                   
# compute corresponding Hindex                                                                                                                                                                             
# plt.contourf                                            

xlist = np.linspace(10, 200,50)
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
    ax.annotate(str(yr)[2:4], (Dperyear[i], maxApyear[i]/1e3), bbox=bboxfy(yr), color='white', horizontalalignment='center', verticalalignment='center')

plt.xlabel('Hypoxia Duration - [d]')
plt.ylabel('Max. Hypoxic Area - [10^3 km2]')
cbaxes = fig.add_axes([0.25, 0.94, 0.5, 0.02]) 
clb=plt.colorbar(cs,  cax=cbaxes, orientation="horizontal")#, position='top')

clb.ax.set_title('Hypox-Index - [10^3 km2]')
fig.savefig(G.figoutputdir+'H_scatter.png')

###################################################
import matplotlib as mpl

mpl.style.use('ggplot')
plt.rcParams['axes.facecolor']='w'
plt.rcParams['grid.color']='lightgrey'
 # b0b0b0  
 # grid color
plt.rcParams['grid.linestyle']='--'  
 # solid
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

'''
fit1,cov = np.polyfit( range(minyear,maxyear)  , Hnormperyear/1e3,1, cov=True)
trend1 = fit1[0]
stderroronslope = np.sqrt(np.diag(cov))[0]
print( 'Annual Trend is %s +/- %s'%(trend1, stderroronslope )

'''

tfile = open(G.figoutputdir+'data.txt','w') 
tfile.write( 'year ; BH-index ; Duration; max Area \n')
tfile.write( '%s \n' %range(minyear,maxyear+1))
tfile.write( '%s \n' %(Hnormperyear/1e3))
tfile.write( '%s \n' %Dperyear)
tfile.write( '%s \n' %maxApyear)
tfile.close() 
