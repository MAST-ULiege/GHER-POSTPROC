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

#mlist= mlist[:-10]

vlist = ['botfluxDOX','botfluxNOS','botfluxNHS','botfluxDIC','fCSED','sCSED','NCrSED','fSSED','sSSED']

for mm in mlist:
    Nb  = N3D_class.N3D(mm.replace('ptrc','btrc'),ymlfile)

    for v in vlist: 
        Nb.testvar(v)

    if mm==mlist[0]:
        N=Nb
    else:
        N.dates        = ma.append(N.dates       , Nb.dates,0)
        N.time         = ma.append(N.time        , Nb.time,0)
        for v in vlist:
            exec('N.'+v+' = ma.append(N.'+v+'       , Nb.'+v+',0)') 
    del Nb
 
N.timeclean()

for v in vlist:
    if (v in ['botfluxDOX','botfluxNOS','botfluxNHS','botfluxDIC']):
            exec('N.'+v+'=N.'+v+'*86400')
    N.makeclim(v)

climdic={ 'botfluxDOX':[-80,80],
          'botfluxNOS':[-5,5],
          'botfluxNHS':[-10,10],
          'botfluxDIC':[-50,50],
          'fCSED':None,
          'sCSED':None,
          'NCrSED':None,
          'fSSED':None,
          'sSSED':None}


for v in vlist:
    N.mapMonthlyClim(v,figsuffix='SHELF', subdomain="NWS", Clim=climdic[v], cmapname='balance',extend='both')
    N.SeasonStrip(v, figsuffix='SHELF', subdomain="NWS", Clim=climdic[v], cmapname='balance',extend='both')

