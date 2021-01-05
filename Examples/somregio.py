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
'''
ymlfile='./local_NEMO_diags.yml'
## ## ## ## ## ##### ## ## ## ## ###  ## ## ## ## ## ###

mlist = N3D_class.FullLoad(YAML_FILE =ymlfile, begdate=dt.datetime(2014,1,1), enddate=dt.datetime(2014,12,31), filetype='grid_T')

vlist = ['vosaline','votemper']
vlist2D = ['SST','SSS','sossheig']

# 1. make clims
for v in vlist2D:
    for mm in mlist[:]:
        Ni  = N3D_class.N3D(mm,ymlfile)
        NmaskDS=  ~(Ni.bat.mask) #& (Ni.bat.mask<120) # Mask should be True where masked
        Ni.testvar(v)
        #setattr(Ni,v,Ni.maskvar(getattr(Ni,v),NmaskDS))
        #setattr(Ni,'vi'+v,Ni.vertmean(v))
        #setattr(Ni,'vi'+v,getattr(Ni,v))
    
        if mm==mlist[0]:
            N=Ni
        else:
            N.dates      = ma.append(N.dates   , Ni.dates,0)
            N.time       = ma.append(N.time    , Ni.time,0)
            setattr(N,v,ma.append(getattr(N,v), getattr(Ni,v),0))
        del Ni

    N.timeclean(begdate=dt.datetime(2014,1,1), enddate=dt.datetime(2017,12,31))

    N.makeclim(v)
    N.mapMonthlyClim(v,cmapname='haline')
    if v==vlist2D[0]:
        Na=N
    else:
        setattr(Na,'clim_'+v, getattr(N,'clim_'+v))

# 2. SOM 1
substep = 100
som_shape = (10,1)
som2_shape = (5,1)



vfa= [getattr(Na,'clim_'+v).compressed()[0:-1:substep] for v in vlist2D ]

som1data = np.stack(vfa, axis=1)

som1means= np.mean(som1data, axis=0)
som1stds = np.std(som1data, axis=0)
som1datan = (som1data - som1means) / som1stds

print(som1data.shape)

from minisom import MiniSom    
som = MiniSom(som_shape[0], som_shape[1], len(vlist2D), sigma=0.5, learning_rate=0.5) # initialization of 6x6 SOM
som.train(som1datan, 10000) # trains the SOM with 100 iterations

# each neuron represents a cluster
winner_coordinates = np.array([som.winner(x) for x in som1datan]).T
# with np.ravel_multi_index we convert the bidimensional
# coordinates to a monodimensional index
cluster_index = np.ravel_multi_index(winner_coordinates, som_shape)

'''
import matplotlib.pyplot as plt
# plotting the clusters using the first 2 dimentions of the data
for c in np.unique(cluster_index):
    plt.scatter(som1data[cluster_index == c, 0],
                som1data[cluster_index == c, 1], label='cluster='+str(c), alpha=.7)

# plotting centroids
for centroid in som.get_weights():
    plt.scatter(centroid[:, 0]*som1stds[0]+som1means[0], centroid[:, 1]*som1stds[1]+som1means[1], marker='x', 
                s=80, linewidths=35, color='k', label='centroid')
plt.legend()
'''

# 3. SOM 2
# Here we first rextract the outcome of the first SOM step for visualization.
# We now have a lon*lat*time matrix of indexes. 
vfc = [getattr(Na,'clim_'+v).compressed() for v in vlist2D ]
som1datac = np.stack(vfc, axis=1)
som1datac = (som1datac- som1means) / som1stds 

winner_coordinates = np.array([som.winner(x) for x in som1datac]).T
cluster_index = np.ravel_multi_index(winner_coordinates, som_shape)

bb=getattr(Na,'clim_'+vlist2D[0]).copy()
bb[~bb.mask]=cluster_index

'''
pcolor(bb[0].squeeze())
som2data = 
'''

Na.somstep1 = bb
Na.gstore('somstep1')

# Prepare data for the second SOM step.
# Not sure it still needed to normalize the data in this case. 

vft= [bb[t].compressed() for t in range(len(N.climdates)) ]
som2data = np.stack(vft, axis=1)

#som2means= np.mean(som1data, axis=0)
#som2stds = np.std(som1data, axis=0)
#som1datan = (som1data - som1means) / som1stds

som2 = MiniSom(som2_shape[0], som2_shape[1], som2data.shape[1], sigma=0.5, learning_rate=0.5) # initialization of 6x6 SOM
som2.train(som2data, 10000) # trains the SOM with 100 iterations

winner_coordinates = np.array([som2.winner(x) for x in som2data]).T
cluster_index = np.ravel_multi_index(winner_coordinates, som2_shape)

bb=getattr(Na,'clim_'+vlist2D[0])[0].copy()
bb[~bb.mask]=cluster_index


Na.regio=bb[:,None,:,:]
Na.gstore('regio')

#N.saveregionmask(bb,'regio.nc') 