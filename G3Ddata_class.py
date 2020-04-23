import numpy as np

import numpy.ma as ma
from netCDF4 import Dataset
import gsw
import yaml
import os.path
import datetime as dt
import cmocean
import calendar
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from   scipy import interpolate 
import scipy.io as sio
from scipy.interpolate import RegularGridInterpolator

#from mpl_toolkits.basemap import Basemap

class G3Ddata(object): 
    '''This is a class for model-data comparison.''' 
######################################################################

    def __init__(self,infile, variable, vardataname, mode='mat'):

        self.infile = infile
        self.mode   = mode

        if self.mode =='mat':
            print(' *** Loading Data ****  \n')
            if os.path.isfile(self.infile):
                print(self.infile + ' -> OK')
                self.found=True
            else:
                print(self.infile + ' can not be found')
                self.found=False
                return

            test = sio.loadmat(infile)

            self.lon      = test['lon'].squeeze()
            self.lat       = test['lat'].squeeze()
        # HERE FOR LATER : use dict for data and model variable. Save relavant value in the object using model variable (or equivalent to instance_...) 
            self.obs      = test[vardataname].squeeze()
            self.time     = test['time'].squeeze()
            self.depth     = test['depth'].squeeze()
            self.variable = variable
            self.model    = np.empty_like(self.obs)

            self.dates = [dt.datetime(1858,11,17)+dt.timedelta(seconds=int(t*86400)) for t in self.time]

        if self.mode =='nc':
            with Dataset(infile,'r') as inf:
                for name, variable in inf.variables.items():
                    self.__setattr__(name,inf.variables[name][:]) 
            self.dates = [dt.datetime(1858,11,17)+dt.timedelta(seconds=int(t*86400)) for t in self.time]

    def subset(self, lons=None, lats=None, datess=None, depths=None):
        if lons is not None:  
            idx= (self.lon>lons[0]) & (self.lon<lons[1])
            self.filter(idx)
        if lats is not None:
            idx= (self.lat>lats[0]) & (self.lat<lats[1])
            self.filter(idx)
        if datess is not None:
            idx= [ (i>datess[0]) & (i<datess[1]) for i in self.dates]
            self.filter(idx)
        if depths is not None: 
           idx= (self.depth>depths[0]) & (self.depth<depths[1])
           self.filter(idx)

    def filter(self, idx):
        otl = len(self.dates)
        bb  = [ a for a in self.__dict__  if isinstance(self.__getattribute__(a),np.ndarray) ]
        bbb = [ b for b in bb if self.__getattribute__(b).shape[0]==otl]

        for b in bbb:
            print('Filtering '+b)
            self.__setattr__(b, self.__getattribute__(b)[idx])

        self.dates = [ self.dates[i] for i,r in enumerate(idx) if r]
            
    def addmodel(self,G,vavar=None):
        if vavar is None:
            vavar=self.variable
        # G should be a model object of class G3D or N3D
        # Add a model column to the data matrix    
        dayfromO   = [ (d-G.dates[0]).days for d in G.dates]
        self.dayfromO = [ (d-G.dates[0]).days for d in self.dates]
        interpolator = RegularGridInterpolator((dayfromO, G.z[0,:,0,0], G.lat,G.lon), G.__getattribute__(vavar), bounds_error=False, fill_value=np.nan)
        points=[  [ self.dayfromO[i], self.depth[i], self.lat[i], self.lon[i]] for i in range(len(self.dates))]
        self.model = interpolator(points)        

    def addmodelbat(self,G):
        interpbat = RegularGridInterpolator((G.lat,G.lon), G.bat.squeeze(), bounds_error=False, fill_value=np.nan)
        points = [ [self.lat[i], self.lon[i]] for i in range(len(self.dates))]
        self.modelbat = interpbat(points) 

    def plot_time(self, binwidth='months', Clim=None, title=None, figout=None, figoutputdir='./'):
        locator = mdates.AutoDateLocator()
        formator = mdates.AutoDateFormatter(locator)

        fig=plt.figure(figsize=(10, 4))
        ax=plt.subplot(1, 1, 1)
        ax.xaxis_date()
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formator)
        plt.plot(self.dates, self.obs, 'r.',label='Obs.'  , alpha=0.8)
        plt.plot(self.dates, self.model,'b.', label='Mod.', alpha=0.8)
        if Clim is None :
            Clim = [ min(np.nanmin(self.obs), np.nanmin(self.model)), max(np.nanmax(self.obs), np.nanmax(self.model))]
            print(Clim)
        plt.ylim(Clim)
        plt.legend()
        if title is None:
            title=self.variable
        plt.title(title)
        if figout is None:
            figout='Obs_'+self.variable
        fig.savefig(figoutputdir+'/'+figout+'.png')
        plt.close()



    def ncsave(self,fname):
        with Dataset(fname,'w') as outf:
            # create record dimension, open
            # All values of the good length using recod dimension .. 
            # Deal with the rest later. 
            outf.createDimension('record',None)
            otl = len(self.dates)
            bb  = [ a for a in self.__dict__  if isinstance(self.__getattribute__(a),np.ndarray) ]
            bbb = [ b for b in bb if self.__getattribute__(b).shape[0]==otl]

            for b in bbb:
                print('Writing '+b)
                outf.createVariable(b, np.float32, 'record')
                outf.variables[b][:] = self.__getattribute__(b)


    def ncload(self,fname):
        with Dataset(fname,'r') as inf:
            for name, dimension in inf.dimensions.items():
                if self.verbose: print ('\n XXX '+name+'  '+str(dimension))
                if ((ndim==2) and (name in [self.latdimname, self.londimname])):
                    diagf.createDimension(name, 1)
                else:
                    diagf.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

            for name, variable in inf.variables.items():
 
                if name  in [self.depthdimname, self.latdimname, self.londimname]: 
                    if ((ndim==2) and (name in [self.latdimname, self.londimname])):
                        x = diagf.createVariable(name, variable.datatype, variable.dimensions)
                        exec('diagf.variables[name][:] = self.'+name)
                    else:
                        x = diagf.createVariable(name, variable.datatype, variable.dimensions)
                        diagf.variables[name][:] = inf.variables[name][:]

            # create record dimension, open
            # All values of the good length using recod dimension .. 
            # Deal with the rest later. 
            outf.createDimension('record',None)
            otl = len(self.dates)
            bb  = [ a for a in self.__dict__  if isinstance(self.__getattribute__(a),np.ndarray) ]
            bbb = [ b for b in bb if self.__getattribute__(b).shape[0]==otl]

            for b in bbb:
                print('Writing '+b)
                outf.createVariable(b, np.float32, 'record')
                outf.variables[b][:] = self.__getattribute__(b)
