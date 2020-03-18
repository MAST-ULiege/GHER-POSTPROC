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

#mlist=mlist[]

vlist = ['NOS','SIO','PHO','DOX','NPP']

for mm in mlist:
    Nb  = N3D_class.N3D(mm,ymlfile)

    for v in vlist: 
#        exec('Nb.'+v+'ksurface=Nb.vertmean(v, zinf=10)')
        Nb.testvar(v,k='surface')

    if mm==mlist[0]:
        N=Nb
    else:
        N.dates        = ma.append(N.dates       , Nb.dates,0)
        N.time         = ma.append(N.time        , Nb.time,0)
        for v in vlist:
            exec('N.'+v+'ksurface = ma.append(N.'+v+'ksurface       , Nb.'+v+'ksurface,0)') 
    del Nb
 
#for v in vlist:
#    exec('N.'+v+'=N.'+v+'*86400')
#    N.makeclim(v+'ksurface')

N.timeclean()

climdic={ 'DOX':[200,400],
          'NOS':[0,5],
          'SIO':[0,25],
          'PHO':[0,1], 
          'NPP':[-1e-5,1e-5]}

for v in vlist:
#    N.mapMonthlyClim(v+'ksurface',Clim=climdic[v],figsuffix='_10m')#, subdomain="NWS"), Clim=climdic[v], cmapname='balance',extend='both')
    N.SeasonStrip(v+'ksurface', Clim=climdic[v])#, subdomain="NWS", Clim=climdic[v], cmapname='balance',extend='both')



climdic={ 'DOX':[200,400],
          'NOS':[0,20],
          'SIO':[0,50],
          'PHO':[0,5],
          'NPP':[-1e-4,1e-4]}

for v in vlist:
#    N.mapMonthlyClim(v+'ksurface',Clim=climdic[v],figsuffix='_10m')#, subdomain="NWS"), Clim=climdic[v], cmapname='balance',extend='both')
    N.SeasonStrip(v+'ksurface', figsuffix='NWS', subdomain="NWS", Clim=climdic[v])

