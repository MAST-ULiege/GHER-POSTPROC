# We only import librairies needed for plotting
# Other librairies are imported in the class definition file, G3D_class.py,
# which contains all process and variables function definition.

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import numpy.ma as ma
import N3D_class
import G3D_class

# We instantiate an object of the class G3D, just by giving the path to the netcdf file to work with
# Up to now I'm working with 4D netcdf files containing several variables. 
# Outputs from different files can be merged easily, as can be seen in other examples


for mm in range(1,13):
    Ni  = N3D_class.N3D('BS_1d_20100101_20171231_ptrc_T_2010'+format(mm,'02')+'-2010'+format(mm,'02')+'.nc','local_NEMO_004.yml')

    Ni.testvar('O2bottom')
    NmaskDS= (Ni.bat<120 ) & ~(Ni.bat.mask)  # Mask should be True where masked

    Ni.apO2=Ni.avgprofileSIGMA(varname='DOX',maskin=NmaskDS)
        
    if mm==1:
        N=Ni
    else:
        N.dates     = ma.append(N.dates   , Ni.dates,0)
        N.time      = ma.append(N.time    , Ni.time,0)
        N.O2bottom  = ma.append(N.O2bottom, Ni.O2bottom,0)
        N.apO2      = ma.append(N.apO2    , Ni.apO2,0)

    del Ni

N.makeclim('O2bottom')
N.mapMonthlyClim('O2bottom',figsuffix='SHELF',cmapname='oxy', subdomain="NWS", Clim=[0,300])
N.mapMonthlyClim('O2bottom',figsuffix='WHOLE',cmapname='oxy', Clim=[0,30])
N.mapMonthlyClim('O2bottom',figsuffix='WHOLEb',cmapname='oxy', Clim=[0,3])

N.plotprofile('apO2',z=-N.z[0,:,0,0],cmapname='oxy',Clim=[0,300])
N.plotprofile('apO2',z=-N.z[0,:,0,0],cmapname='oxy',zlim=[-200,0])
N.plotprofile('apO2',z=-N.z[0,:,0,0],cmapname='oxy',Clim=[0,3],zlim=[-2200,-1000],figout='apO2b')
