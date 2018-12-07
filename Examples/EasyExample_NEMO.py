# We only import librairies needed for plotting
# Other librairies are imported in the class definition file, G3D_class.py,
# which contains all process and variables function definition.

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt

import N3D_class
import G3D_class

# We instantiate an object of the class G3D, just by giving the path to the netcdf file to work with
# Up to now I'm working with 4D netcdf files containing several variables. 
# Outputs from different files can be merged easily, as can be seen in other examples

N  = N3D_class.N3D('OutNEMO/BS_1d_19900101_20071231_grid_T_199401-199401.nc','local_NEMO.yml')
G  = G3D_class.G3D('OutGHER/r5.nc', 'local_GHER.yml')

# All loaded variables are attributes of the G3D instance.
# For example the variable "bat" is defined directly when the object is instantiated. 
# Other are loaded only when needed.
# Variables are python Masked_array, so they have an attribute mask which is an arreay of booleans 

# Here we want to define a mask based on bathymetry 
NmaskDS= (N.bat<50 ) & ~(N.bat.mask)  # Mask should be True where masked
NmaskSH= (N.bat>=50) & ~(N.bat.mask)  # Mask should be True where masked                                                                                                                                         

GmaskSH= (G.bat>=50) & ~(G.bat.mask)
GmaskDS= (G.bat<50)  & ~(G.bat.mask)
# All processing functions are called as function of the G3D instance. 
# Variable name is given as an argument. Some functions allows more argument.

# This would give the basin averaged time serie of salinity
#T1 = N.avgspatial('SAL')

# The avgspatial function enables an optional mask argument
# Note also , that we can use a variable name that is not defined in the netcdf file. 
# In this case the toolbox will automatically look for the function "instance_SSS"

NsssDS=N.avgspatial('SSS',NmaskDS)
NsssSH=N.avgspatial('SSS',NmaskSH)

GsssDS=G.avgspatial('SSS',GmaskDS)
GsssSH=G.avgspatial('SSS',GmaskSH)

# The following is general python plotting .. 
# the "dates" attribute is also loaded automatically 

####################                                                                                                                                                                                               
# 1st figure :  
####################                                                                                                                                                                                               

locator = mdates.AutoDateLocator()
formator = mdates.AutoDateFormatter(locator)

fig=plt.figure(figsize=(15, 8))

ax=plt.subplot(1, 1, 1)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.plot(N.dates,NsssSH, label = 'Average surface salinity on the shelf - NEMO', color='r',linestyle='-')
plt.plot(G.dates,GsssSH, label = 'Average surface salinity on the shelf - GHER', color='r',linestyle='--')
plt.plot(N.dates,NsssDS, label = 'Average surface salinity in the open sea - NEMO', color='b',linestyle='-')
plt.plot(G.dates,GsssDS, label = 'Average surface salinity in the open sea - GHER', color='b',linestyle='--')
plt.title('Sea Surface Salinity')
plt.ylabel('Practical Salinity - []')
plt.legend(loc='best')
fig.savefig(N.figoutputdir+'Simple.png')

