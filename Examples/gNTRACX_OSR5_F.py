import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

#from mpl_toolkits.basemap import Basemap 
#from multiprocessing import Pool                                                                        
#import gsw                                                                                                                                                                                             \
                                                                                                                                                                                                                   
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

#mlist=mlist[:10]

for mm in mlist:
    Ni  = N3D_class.N3D(mm,ymlfile)

    NmaskDS= (Ni.bat<120 ) & ~(Ni.bat.mask)  # Mask should be True where masked
    Ni.viNOS = Ni.vertint('NOS')
    Ni.apNOS = Ni.avgprofileSIGMA(varname='NOS',maskin=NmaskDS)
    Ni.viNHS = Ni.vertint('NHS')
    Ni.apNHS = Ni.avgprofileSIGMA(varname='NHS',maskin=NmaskDS)
    Ni.viDEN = Ni.vertint('DENITRIFICATION')
    Ni.apDEN = Ni.avgprofileSIGMA(varname='DENITRIFICATION',maskin=NmaskDS)
    Ni.viAMA = Ni.vertint('AMANOX')
    Ni.apAMA = Ni.avgprofileSIGMA(varname='AMANOX',maskin=NmaskDS)
#    Ni.viOXN = Ni.vertint('OXIDATIONBYNOS')
#    Ni.apOXN = Ni.avgprofileSIGMA(varname='OXIDATIONBYNOS',maskin=NmaskDS)

    if mm==mlist[0]:
        N=Ni
    else:
        N.dates      = ma.append(N.dates   , Ni.dates,0)
        N.time       = ma.append(N.time    , Ni.time,0)
        N.apNOS      = ma.append(N.apNOS    , Ni.apNOS,0)
        N.apNHS      = ma.append(N.apNHS    , Ni.apNHS,0)
        N.apDEN      = ma.append(N.apDEN    , Ni.apDEN,0)
        N.apAMA      = ma.append(N.apAMA    , Ni.apAMA,0)
 #       N.apOXN      = ma.append(N.apOXN    , Ni.apOXN,0)
        N.viNOS      = ma.append(N.viNOS     , Ni.viNOS,0)
        N.viNHS      = ma.append(N.viNHS     , Ni.viNHS,0)
        N.viAMA      = ma.append(N.viAMA     , Ni.viAMA,0)
        N.viDEN      = ma.append(N.viDEN     , Ni.viDEN,0)
 #       N.viOXN      = ma.append(N.viOXN     , Ni.viOXN,0)
        

    del Ni

N.timeclean()

for v in ['NOS','NHS','DEN','AMA']:#,'OXN']:
    N.plotprofile('ap'+v,z=-N.z[0,:,0,0])
#    N.makeclim('vi'+v)
#    N.mapMonthlyClim('vi'+v,figsuffix='SHELF',cmapname='haline', subdomain="NWS")
#    N.mapMonthlyClim('vi'+v,figsuffix='WHOLE',cmapname='haline')
#    N.mapStrip('vi'+v,figsuffix='_D',daysbetween=30,diff=True)

for v in ['NOS','NHS','DEN','AMA']:
    N.SeasonStrip('vi'+v)


