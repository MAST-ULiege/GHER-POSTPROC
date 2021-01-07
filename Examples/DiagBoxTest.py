import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

#from mpl_toolkits.basemap import Basemap 
#from multiprocessing import Pool                                                                        
#import gsw                                                                                                                                                                                             \
                                                                                                                                                                                                                   
import matplotlib
#matplotlib.use('pdf')
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
'''
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
'''
ymlfile='./local_NEMO_diags.yml'
## ## ## ## ## ##### ## ## ## ## ###  ## ## ## ## ## ###         


mlist = N3D_class.FullLoad(YAML_FILE =ymlfile, begdate=dt.datetime(2014,1,1), enddate=dt.datetime(2014,12,31))
flist2D = ['AirSeaOxygenFlux_Source','AirSeaOxygenFlux_Sink','Photosynthesis2D','PhytoNitrateReduction2D','Resp_Phyto2D','Resp_Zoo2D','Resp_Gel2D','Resp_Bac2D','Oxidation2D','BenthicOxygenFlux']

for v in flist2D:
    for mm in mlist[:]:
        Ni  = N3D_class.N3D(mm,ymlfile)
        NmaskDS= ~(Ni.bat.mask)  # Mask should be True where masked                                                                                                                                 
        Ni.testvar(v)
        
        if mm==mlist[0]:
            N=Ni
            reg = N.loadregionmask('./regions.nc','kopelevich')
        else:
            N.dates      = ma.append(N.dates   , Ni.dates,0)
            N.time       = ma.append(N.time    , Ni.time,0)
            setattr(N,v,ma.append(getattr(N,v), getattr(Ni,v),0))

    del Ni

    #m5 = (reg != 5.0) & ~ (G1.bat.mask)

    N.timeclean(begdate=dt.datetime(2014,1,1), enddate=dt.datetime(2014,12,31))

#reg=N.bat.copy() 
#reg[N.bat<=120]=1 
#reg[N.bat>120]=2

#reg.mask=N.bat.mask
reg=N.loadregionmask('regio_from_O2diags.nc','regio')

Bmean = N.FluxBars(flist2D,reg, figsuffix='means',type='mean',factor=86400, unit='mmol/m²/d') 
Bint  = N.FluxBars(flist2D,reg, figsuffix='ints' ,type='int' ,factor=1e-12, unit='Gmol') 


np.savetxt('Integrated.txt',Bint, header=' '.join(flist2D))
np.savetxt('Mean.txt',Bint, header=' '.join(flist2D))