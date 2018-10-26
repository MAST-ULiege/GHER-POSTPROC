import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import gsw
import yaml
import os.path
import datetime as dt
import cmocean
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

class G3D(object): 
    '''This is a class for model outputs exploration. It is based on GHER3D model outputs, but child classes are available for other models.'''

######################################################################

    def __init__(self,infile):

        self.infile=infile

        #The diag filename can be used to store computed diagnostic
        self.diagfile= self.infile[:-3]+".diag.nc"

        print(' *******  \n')
        if os.path.isfile(self.infile):
            print(self.infile + ' -> OK')
            self.found=True
        else:
            print(self.infile + ' can not be found')
            self.found=False
            return
        try:
            YAML_FILE = 'local.yml'
            print("\nLaunching with YAML file: %s" % YAML_FILE)
            # Read yaml configuration file
            with open(YAML_FILE, 'r') as stream:
                config = yaml.load(stream)
        except Exception:
            print("".join(("\n A file called local.yml should be present","'\n")))
        

        try:
            self.model    = config['MODEL'] 
        except :
            self.model    = 'GHER'

    # this should be read somewhere
        self.batfile      = config['BATFILE']
        self.dx           = config['DX']
        self.dy           = config['DY']

    # to get intervals in the double sigma context of GHER
        self.sigI         = config['SIGI']
        self.sigII        = config['SIGII']
        self.hlim         = config['HLIM']
        self.verbose      = config['VERBOSE']
        self.figoutputdir = config['PLOTDIR']
        self.resultdir    = config['RESULTDIR']

        self.instance_bat()
        self.testtime()

######################################################################
# VARIABLE : BAT

    def instance_bat(self):
        if os.path.isfile(self.batfile):
            print(self.batfile + ' -> OK')
        else:
            print(self.batfile + ' can not be found')

        with Dataset(self.batfile,'r') as nc:
            self.bat = nc.variables['bat'][:]
            
######################################################################
# VARIABLE : Z

    def instance_z(self):
        # For now let's build a constant z,
        # dynamic z considering ETA can be done later
        
        if self.model is 'NEMO':
            self.z   = nc.variables['deptht'][:]

        with Dataset(self.infile,'r') as nc:
            self.z   = nc.variables['depth'][1] # Let it be 3D for now           
            self.lon = nc.variables['longitude'][:]
            self.lat = nc.variables['latitude'][:]

######################################################################
# VARIABLE : ZI
            
    def instance_zi(self):
        # For now let's build a constant dz,
        # dynamic dz considering ETA can be done later
        self.zi=ma.zeros(self.z.shape+np.array([1,0,0]))
        for k in xrange(0,12):
            self.zi[k]= - ma.masked_where( self.bat<=self.hlim,  (( self.bat - self.hlim ) * (1-self.sigII[k]) ) + self.hlim )
        for k in xrange(12,self.zi.shape[0]):
            self.zi[k]= - np.minimum(self.bat,self.hlim) * (1-self.sigI[k-12+1])
            
######################################################################
# VARIABLE : DZ

    def instance_dz(self):
        # For now let's build a constant z,
        # dynamic z considering ETA can be done later
        self.dz = ma.abs(self.zi[1:]-self.zi[:-1]) 

######################################################################
# UTILITY : LOAD
        
#    def gload(self,varname):
#        with Dataset(self.infile,'r') as nc:
#            try:
#                exec('self.'+varname+ '= nc.variables[''varname''][:]')
#                print( 'Just loaded %s'%(varname))
#            except: 
#                print( '%s not found in %s'%(varname,self.infile))
#                print( '-> try to compute')
#             #   try:
#                exec('self.instance_'+varname+'()')
#               # except:
#               #     print('self.instance_'+varname+' is not defined.')

######################################################################
# UTILITY : LOAD LOCAL 
        
    def gload(self,varname,ti=None,k=None,i=None,j=None):
        try:
            # looking for the variable in the original file
            with Dataset(self.infile,'r') as nc:
                if (ti is None) and (k is None) and (i is None) and (j is None):
                    exec('self.'+varname+ '= nc.variables[varname][:]')
                    print( 'Just loaded '+ (varname) + ' for full')
                elif (ti is None) and (k is None) and (i is not None) and (j is not None):
                    print( 'Loading '+ (varname) + ' for i:'+str(i)+' and j:'+str(j))
                    exec('self.'+varname+'i'+str(i)+'j'+str(j)+'= nc.variables[varname][:,:,j,i]')
                    print( 'Just loaded '+ (varname) + ' for i:'+str(i)+' and j:'+str(j))
                elif (ti is None) and (k is not None) and (i is None) and (j is None):
                    exec('self.'+varname+'k'+str(k)+'= nc.variables[varname][:,k]')
                    print( 'Just loaded '+ (varname) +' for ki:'+str(k))
                else:
                    print(' Stange case encountered in  G3D_class.py : def gload ')
                    
        except: 
            try:
                # looking for the variable in the diag file (created by function gstore)
                with Dataset(self.diagfile,'r') as nc:
                    if (ti is None) and (k is None) and (i is None) and (j is None):
                        exec('self.'+varname+ '= nc.variables[varname][:]')
                        print( 'Just loaded '+ (varname) + ' for full')
                    elif (ti is None) and (k is None) and (i is not None) and (j is not None):
                        print( 'Loading '+ (varname) + ' for i:'+str(i)+' and j:'+str(j))
                        exec('self.'+varname+'i'+str(i)+'j'+str(j)+'= nc.variables[varname][:,:,j,i]')
                        print( 'Just loaded '+ (varname) + ' for i:'+str(i)+' and j:'+str(j))
                    elif (ti is None) and (k is not None) and (i is None) and (j is None):
                        exec('self.'+varname+'k'+str(k)+'= nc.variables[varname][:,k]')
                        print( 'Just loaded '+ (varname) +' for ki:'+str(k))
                    else:
                        print(' Stange case encountered in  G3D_class.py : def gload ')


            except Exception as e: 
                # looking for an instance function allowing to define this variable
                print(e)
                print( '\n %s not found in %s'%(varname,self.infile))
                print( '-> Calling')
                #try:
                print('     self.instance_'+varname+'(i='+str(i)+',j='+str(j)+',k='+str(k)+')')
                exec('self.instance_'+varname+'(i=i,j=j,k=k)')
                #except:
                #    print('self.instance_'+varname+' is not defined.')
                    
######################################################################
# UTILITY : Save a map of regions (as integers)
#
    def saveregionmask(self,regions,outname):
        regionvarame='region'
        
        with Dataset(self.infile,'r') as inf, Dataset(self.resultdir+outname+'.nc','w') as regf:
            # copy attributes
            for name in inf.ncattrs():
                regf.setncattr(name, inf.getncattr(name))

            #copy dimensions
            for name, dimension in inf.dimensions.items():
                if name  in ['longitude','latitude']: 
                    regf.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

            for name, variable in inf.variables.items():
                # take out the variable you don't want
                if name  in ['longitude','latitude']: 
                    x = regf.createVariable(name, variable.datatype, variable.dimensions)
                    regf.variables[name][:] = inf.variables[name][:]

            regf.createVariable(regionvarame, np.float32,('latitude', 'longitude'),zlib=True)
            regf.variables[regionvarame][:]=regions
            
            
######################################################################
# UTILITY : Load a map of regions (as integers)

    def loadregionmask(self,regionFileName,regionvarname='region'):
        with Dataset(regionFileName,'r') as nc:
            regionmap = nc.variables[regionvarname][:]
        return(regionmap)
                
######################################################################
# UTILITY : STORE
#
# create a copy of the file "in.nc" as "in.ext.nc" copy dimension and attributes and add the requested value
         
    def gstore(self,varname, ztab=None):
        
        
        try:
            with Dataset(self.diagfile,'r') as nc:
                lon= nc.variables['time'][:]
        except:
            # -> the diag file does not exist.
            # We should create one with same dimension and attributes
            with Dataset(self.infile,'r') as inf, Dataset(self.diagfile,'w') as diagf:
                # copy attributes
                for name in inf.ncattrs():
                    diagf.setncattr(name, inf.getncattr(name))

                #copy dimensions
                for name, dimension in inf.dimensions.items():
                    diagf.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

                for name, variable in inf.variables.items():
                    # take out the variable you don't want
                    if name  in ['longitudes','lattitude','time','depth','level']: 
                        x = diagf.createVariable(name, variable.datatype, variable.dimensions)
                        diagf.variables[name][:] = inf.variables[name][:]

        exec('ndim=len(self.'+varname+'.squeeze().shape)')
        print('\n Storing now '+varname+' ('+ str(ndim)+' dimensions) on '+self.infile)
        exec('print(self.'+varname+'.shape)')

        with Dataset(self.diagfile,'a') as nc:
            print('ndim : '+ str(ndim) )
            try:
                if ndim==4:      # assuming here : time, level, lat,lon
                    nc.createVariable(varname, np.float32, ('time',     'level', 'latitude', 'longitude'),zlib=True)
                elif ndim == 3:  # assuming here : time, lat,lon 
                    nc.createVariable(varname, np.float32, ('time', 'singleton', 'latitude', 'longitude'),zlib=True)
                elif ((ndim == 2) and ~(ztab is None)) :  # assuming here : time, ztab
                    try:
                        nc.variables['ztab'][:]
                    except:
                        nc.createDimension('ztab', len(ztab))
                        nc.createVariable('ztab', np.float32, 'ztab')
                        nc.variables['ztab'][:] = ztab
                    nc.createVariable(varname, np.float32, ('time', 'ztab'),zlib=True)
                elif ndim == 1:
                    nc.createVariable(varname, np.float32, ('time'),zlib=True)
            except:
                print ('Seems that '+varname+' already exists on '+self.infile+'. \n I overwrites')

            
            exec('nc.variables[varname][:]=self.'+varname)
                
######################################################################
# PROCESS : FULL INTEGRATION

    def integratespatial (self,varname):
    # TODO allows for a optional mask, that could be 2D, 3D or 4D. 
    #    if dz.shape!=field.shape:
    #        print("Wrong Shapes")
    #        return
    #    if len(field.shape)==3:
    #        print("3D field")
    #    
        self.testz()
        self.testvar(varname)
        self.testtime()

        integrated=ma.empty(len(self.time))
        exec('loc=self.'+varname)
        print('dz: %s  and field: %s'%( len(self.dz.shape),len(loc.shape)))
        if (len(self.dz.shape)==3)and(len(loc.shape)==4):
            print("4D")
            for t in xrange(len(self.time)):
                bi=loc[t]*self.dz*self.dy*self.dx
                integrated[t] = ma.sum(bi)
        # return should be 1D (time)    
        return(integrated)
        
######################################################################
# PROCESS : TIME AVERAGE
    def avgtime (self, varname, maskin=None):
        self.testvar(varname)
        avg=self.varname.mean(axis=0)
        return(avg)
        
######################################################################
# PROCESS : FULL AVERAGE

    def avgspatial (self,varname,maskin=None):
        # in the sense of volumetric mean
        self.testz()
        self.testvar(varname)        
                        
        avg=ma.empty(len(self.time))
        exec('loc=self.'+varname+'.copy()')
        print('In avgptial \n dz: %s  \n Field: %s'%( len(self.dz.shape),len(loc.shape)))

        # DIMENSIONAL CASE 
        # 2D VARIABLE
        print(loc.shape)
        if ( loc.shape[1]==1 ):
            print('2D variable')
            # MASKING
            if (maskin is not None):
                print('Found Mask with ' + str(maskin.sum()) + ' masked points.')
                if len(maskin.shape)==3:
                    print('2D mask .. OK')
                    for t in xrange(len(self.time)):
                        loc[t]=ma.masked_where(maskin,loc[t])

            # AVERAGING
            for t in xrange(len(self.time)):
                bi=loc[t]*self.dy*self.dx
                vol=ma.masked_where(bi.mask,self.dy*self.dx*np.ones(bi.shape))
                avg[t] = ma.sum(bi)/ma.sum(vol)

        # 3D VARIABLE
        elif (loc.shape[1]>1):
            if (maskin is not None):
                print('Masking of 3D vars not implemented yet .. Complete G3D_class.py') 
            print("3D variable")
            for t in xrange(len(self.time)):
                bi=loc[t]*self.dz*self.dy*self.dx
                vol=ma.masked_where(bi.mask,self.dz*self.dy*self.dx)
                avg[t] = ma.sum(bi)/ma.sum(vol)
        else:
            print('this dimensional case is not condidered yet')
                
        # return should be 1D (time)    
        return(avg)
        
######################################################################
# UTILITY : test z        
    def testz(self):
        try:
            self.dz
        except :
            print('dz not found -> loading')
            self.instance_z()
            self.instance_zi()
            self.instance_dz()
            
######################################################################
# UTILITY : test variable
            
    def testvar(self,varname,doload=True):
        try:
            exec('self.'+varname)
            isthere=True
        except:
            print('%s not found '%(varname))
            isthere=False
            if doload:
                print ('Loading %s'%(varname) )
                self.gload(varname)
                isthere=True

        return(isthere)
                
            
######################################################################
# UTILITY : test time
            
    def testtime(self):
        try:
            self.time
        except:
            print('%s not found -> loading'%('time'))
            self.gload('time')
            self.dates = [dt.datetime(1858,11,17)+dt.timedelta(days=int(t)) for t in self.time]

######################################################################
# UTILITY : Apply a mask on a variable, handling dimension cases. 

    def maskvar(self,loc, maskin):
    
        if (len(loc.shape)==4) and (len(maskin.squeeze().shape)==2):
            print('Masking '+ str(len(loc.shape))+ 'D variable with '+str(len(maskin.shape))+'D mask')
            for t in xrange(loc.shape[0]):
                for k in  xrange(loc.shape[1]):
                    loc[t,k,:,:]=ma.masked_where(maskin.squeeze(),loc[t,k,:,:]) # ma.expand_dims(ma.expand_dims(maskin,0),0)
        else:
            print('NOT Masking '+ str(len(loc.shape))+ 'D variable with '+str(len(maskin.shape))+'D mask')
                
        return loc

        
######################################################################
# PROCESS : Horizontally-averaged profile in z coordinate, taking sigma coordinate in consideration

    def avgprofile(self,varname,
                   ztab=-1*np.concatenate([np.arange(0,10,2), np.arange(10,40,5),np.arange(50,120,10),np.arange(120,300,50),np.arange(300,1000,200)]),
                   maskin=None
                   ):
        # The idea is to get an average profile as a function of time
        # return is 2D : [z,time]
        # sigma-space makes it a bit complicate
        #ztab=-1*np.concatenate([np.arange(0,10,2), np.arange(10,40,5),np.arange(50,120,10),np.arange(120,300,50),np.arange(300,1000,200)])
        #ztab=-1*np.arange(0,100,5)
        
        self.testz()
        self.testvar(varname)
        self.testtime()
        
        avg = ma.empty((len(self.time),ztab.shape[0]-1))
        
        exec('loc=self.'+varname)

        if maskin is not None:
            loc = self.maskvar(loc,maskin)
       
        if (len(self.dz.shape)==3)and(len(loc.shape)==4):            
            gridZU = self.zi[1:]
            gridZD = self.zi[:-1]
            for k in xrange(ztab.shape[0]-1):
                print('%s / %s'%(k+1,ztab.shape[0]-1))
                dzloc= ma.maximum(ma.zeros(self.dz.shape), np.minimum(gridZU, ztab[k])-np.maximum(gridZD, ztab[k+1]))
                for t in xrange(len(self.time)):
                    vol=ma.masked_where(loc[t].mask,dzloc*self.dy*self.dx)
                    bi=loc[t]*vol
                    #avg[t,k]= ma.sum(bi)/ma.sum(vol)
                    avg[t,k]= ma.masked_where( (float(bi.count())/float(bi.size))<0.1,ma.sum(bi)/ma.sum(vol))
                    if (t==10):
                        print('k %s : from %s to %s' %(k,ztab[k],ztab[k+1]))
                        print('total volume considered for this layer :  %s km^3' %(ma.sum(vol)/1e9))
                        print('mean age :  %s' %(avg[t,k]))
                    
        zforplot=(ztab[:-1]+ztab[1:])/2
        return avg, zforplot
        
######################################################################
# PROCESS : Horizontally-averaged profile in DENSITY coordinate, taking sigma coordinate in consideration

    def avgprofile_den(self,varname,
                   ztab=-1*np.concatenate([np.arange(1010,1018,.2)]),
                   maskin=None
                   ):
        # The idea is to get an average profile as a function of time
        # return is 2D : [den,time]
        # sigma-space makes it a bit complicate
        
        self.testz()
        self.testvar(varname)
        self.testvar('DEN')
        self.testtime()
        
        avg = ma.empty((len(self.time),ztab.shape[0]-1))
        
        exec('loc=self.'+varname)

        if maskin is not None:
            loc = self.maskvar(loc,maskin)
       
        if (len(self.dz.shape)==3)and(len(loc.shape)==4):            
            gridZU = self.zi[1:]
            gridZD = self.zi[:-1]
            for k in xrange(ztab.shape[0]-1):
                print('%s / %s'%(k+1,ztab.shape[0]-1))
                dzloc= ma.maximum(ma.zeros(self.dz.shape), np.minimum(gridZU, ztab[k])-np.maximum(gridZD, ztab[k+1]))
                for t in xrange(len(self.time)):
                    vol=ma.masked_where(loc[t].mask,dzloc*self.dy*self.dx)
                    bi=loc[t]*vol
                    avg[t,k]= ma.sum(bi)/ma.sum(vol)
                    if (t==10):
                        print('k %s : from %s to %s' %(k,ztab[k],ztab[k+1]))
                        print('total volume considered for this layer :  %s km^3' %(ma.sum(vol)/1e9))
                        print('mean age :  %s' %(avg[t,k]))
                    
        zforplot=(ztab[:-1]+ztab[1:])/2
        return avg, zforplot
######################################################################                                                                                        
# PROCESS : Horizontally-averaged profile in sigma coordinate

    def avgprofileSIGMA(self,varname, maskin=None):
        self.testvar(varname)
        self.testtime()
        exec('loc=self.'+varname)
        
        if maskin is not None:
            loc = self.maskvar(loc,maskin)

        avg=ma.average(loc,(2,3))
        return avg



######################################################################                                                                                                                                              
# PROCESS : Horizontally integrated profile in sigma coordinate 


    def horintprofileSIGMA(self,varname, maskin=None):
        self.testvar(varname)
        self.testtime()
        exec('loc=self.'+varname)

        if maskin is not None:
            loc = self.maskvar(loc,maskin)

        hint=ma.empty( (loc.shape[0],loc.shape[1]) )
    
        for t in xrange(len(self.time)):
            hint[t]=ma.sum(loc[t]*self.dz*self.dx*self.dy,(1,2))

        return hint

######################################################################
# UTILITY : test spatial coordinates

    def test_coord(self,c1,c2):
        if (type(c1) is int) & (type(c2) is int):
           #considered as indexes                                                                                                                                   
            i=c1
            j=c2
        else:
           #considered coordinates   
            X=c1-self.lon
            i=abs(X).argmin()
            Y=c2-self.lat
            j=abs(Y).argmin()
        return i,j
        
######################################################################    
# UTILITY : provide lon lat coordinates for indexes
    
    def getlonlat(self,c1,c2):
        i,j=self.test_coord(c1,c2)
        lon=self.lon[i]
        lat=self.lat[j]
        return lon,lat

######################################################################
#
# PROCESS : Depth profile at given coordinate
#          
#  inputs : 
#          * varname 
#          * indexes : integer will be interpretated as grid indexes, floats as coordinates (lon,lat) 

    def profileatxy(self,varname,c1,c2):

        self.testz()
        isvar=self.testvar(varname,False)
        self.testtime()
        i,j=self.test_coord(c1,c2)

        zforplot=self.z[:,j,i]
        if isvar:
            # for some reason the 4d variable is already loaded in memory
            exec('loc=self.'+varname)
            lloc=loc[:,:,j,i]
        else:
            # It's not loaded, but we'll load only what we need
            self.gload(varname,i=i,j=j)
            exec('lloc=self.'+varname+'i'+str(i)+'j'+str(j))

        return lloc,zforplot
        
        
######################################################################
#
# PROCESS : Value at given coordinate
#          
#  inputs : 
#          * varname 
#          * indexes : integer will be interpretated as grid indexes, floats as coordinates (lon,lat) 

    def valatxyz(self,varname,c1,c2,c3):

        self.testz()
        isvar=self.testvar(varname,False)
        self.testtime()
        i,j=self.test_coord(c1,c2)

        zloc=self.z[:,j,i]
        if c3=='bottom':
            k=11 # 11 is the python index for the bottom layer in the first sigma layer !! very case specific XX  
        else: 
            k=(np.abs(zloc-c3)).argmin()
        if isvar:
            # for some reason the 4d variable is already loaded in memory
            exec('lloc=self.'+varname)
            lloc=loc[:,:,j,i]
        else:
            # It's not loaded, but we'll load only what we need
            self.gload(varname,i=i,j=j)
            exec('lloc=self.'+varname+'i'+str(i)+'j'+str(j))
        if lloc.shape[1]==1:
            k=0
        lloc=lloc[:,k]

        return lloc

##########################################################################
# PROCESS : Vertical Integration

    def vertint(self,varname,zinf=-10000,zsup=2):

        self.testz()
        self.testvar(varname)
        self.testtime()

        exec('loc=self.'+varname)
        vint=ma.empty( (loc.shape[0],loc.shape[2],loc.shape[3]) )

        for t in xrange(len(self.time)):
            vint[t]=ma.sum(loc[t]*self.dz,0)
        
        return vint

#########################################################################                                                                                                                                          
# PROCESS : Vertical Mean
  
    def vertmean(self,varname,zinf=-10000,zsup=2):

        self.testz()
        self.testvar(varname)
        self.testtime()

        exec('loc=self.'+varname)
        vmean=ma.empty( (loc.shape[0],loc.shape[2],loc.shape[3]) )

        for t in xrange(len(self.time)):
            vmean[t] = ma.sum ( loc[t]*self.dz , 0)
            vol      = ma.sum (        self.dz , 0)
            vmean[t] = vmean[t]/vol

        return vmean
        
#########################################################################                                                                                                                                          
# PROCESS : Vertical Mean - Isopycnals
  
    def vertmean_den(self,varname,rhoinf=1000,rhosup=1026):

        self.testz()
        self.testvar(varname)
        self.testvar('DEN')
        self.testtime()

        exec('loc=self.'+varname)
        vmean=ma.empty( (loc.shape[0],1,loc.shape[2],loc.shape[3]) )

        for t in xrange(len(self.time)):
            mask3D   = (self.DEN[t]<=rhoinf)|(self.DEN[t]>rhosup)
            loct=ma.masked_where(mask3D,loc[t])
            vmeanloc = ma.sum ( loct*self.dz , 0)
            vol      = ma.sum ( ma.masked_where(mask3D,self.dz), 0)
            vmean[t] = vmeanloc/vol

        return vmean        
        
############################################################################
# PROCESS : Build Climatology
#
# From a multi-year G object, compile a climatology for the required variable. 
#
# A climatology object has 365 daily step. 
# For each day, a field is the average of all input field within [d-climspan: d+climspan], with climspan=7 by default.

    def makeclim(self,varname,climspan=7):
        import calendar
        self.testz()
        # G might be reconstruct for several input file. In that case however, the var instance should exist, as gload won't work.
        self.testvar(varname) 
        self.testtime()
        
        exec('loc=self.'+varname)
        
        clim=ma.empty([366,loc.shape[1],loc.shape[2],loc.shape[3]])
        climdates=[i.replace(year=2000) for i in self.dates]
        num_days = [calendar.monthrange(2000, month)[1] for month in range(1,12+1)]
        clim.dates=[]

        for month in range(1,12+1):
            for day in range(1, num_days[month-1]+1):
                clim.dates.append(dt.datetime(2000, month, day))
                
        for ti in range(clim.shape[0]):
            tds=np.array([ abs( (clim.dates[ti]-i).days ) for i in climdates ])
            # for early Jan days or late december
            tds[tds>300]=367-tds[tds>300]
            
            idx=[ i for i,delt in enumerate(tds) if (delt<climspan)]
            ww =[climspan-i for i in tds if (i<climspan) ]
            
            clim[ti]=ma.average(loc[idx],axis=0,weights=ww)
            
        exec('self.clim_'+varname+'=clim')
        self.climdates=clim.dates
        return clim
        
############################################################################
# PLOTS : Plot Map 

    def mapClimTimeMean(self, varname):
        exec('loc=self.clim_'+varname)
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        
        m = Basemap(llcrnrlat=self.lat[0],\
                    urcrnrlat=self.lat[-1],\
                    llcrnrlon=self.lon[0],\
                    urcrnrlon=self.lon[-1],\
                    resolution='l')
                    
        m.drawcoastlines()
        llon,llat=np.meshgrid(G.lon,G.lat)
        x, y = m(llon,llat)        
                
        m.contourf(x,y,G.clim_T1age[:,11].mean(axis=0), np.linspace(0,600,10),cmap=cmocean.cm.deep, extend="both")
        plt.title()
        m.colorbar()#ticks=np.linspace(Tmin,Tmax,num=11)) matplotlib.rcParams.update({'font.size': 18})
        fig.savefig(figout)

        
############################################################################
# PLOTS : Plot Map 

    def map2D(self, varname,cmapname="deep", region=None, figout='Map.png', title=None, Clim=None, scatmat=None):
        self.testz()
        exec('loc=self.'+varname)
        exec('cmap=cmocean.cm.'+cmapname)
        
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_axes([0.1,0.1,0.8,0.8])

        if region==None:
            i1=0; i2=-1; j1=0; j2=-1
        elif region == 'Shelf':
            i1=0; i2=150; j1=35; j2=-1            
            
        if (title==None): title=varname
            
        if Clim==None : Clim=[loc.min(),loc.max()]
        m = Basemap(llcrnrlat=self.lat[j1],\
                    urcrnrlat=self.lat[j2],\
                    llcrnrlon=self.lon[i1],\
                    urcrnrlon=self.lon[i2],\
                    resolution='l')
                    
        m.drawcoastlines()
        parallels = np.arange(40.,48.,1.)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
        meridians = np.arange(20.,43.,1.)
        m.drawmeridians(meridians,labels=[1,0,0,0],fontsize=10)

        llon,llat=np.meshgrid(self.lon[i1:i2],self.lat[j1:j2])
        x, y = m(llon,llat)        
        
       # print('x.shape %s'%(x.shape))
       # print('y.shape %s'%(y.shape))
       # print('llon.shape %s'%(llon.shape))
       # print('llat.shape %s'%(llon.shape))
                
        cs=m.contourf(x,y,loc.squeeze()[j1:j2,i1:i2], np.linspace(Clim[0],Clim[1],100),cmap=cmap, extend="both")
        cbar = m.colorbar(cs,location='bottom',pad="5%",ticks=np.linspace(Clim[0],Clim[1],10))  #          matplotlib.rcParams.update({'font.size': 18})
        cbar.set_label(title)
        if (scatmat!=None):
            m.scatter(scatmat[:,1],scatmat[:,0],50,scatmat[:,2],edgecolors='black',cmap=cmap,vmin=Clim[0],vmax=Clim[1])
        m.contour(x,y,self.bat[0,j1:j2,i1:i2],levels=[40,80,120], colors='k',linestyles='dashed')
        
        #plt.title(title)
        #m.colorbar()
        fig.savefig(figout)

############################################################################
# VARIABLE : Density

    def instance_DEN(self,i=None,j=None,k=None):

        self.testz()
        tlat=np.tile(self.lat, (self.z.shape[0],self.z.shape[2],1))
        tlat=np.transpose(tlat,(0,2,1))
        tlon=np.tile(self.lon, (self.z.shape[0], self.z.shape[1],1))

        if (i is None) and (j is None) and (k is None):
            p  = gsw.p_from_z(self.z,tlat)
            self.testvar('SAL')
            self.testvar('TEM')
            self.DEN = ma.empty(self.SAL.shape)
            for t in xrange(self.SAL.shape[0]):
                SA = gsw.SA_from_SP(self.SAL[t],p,tlon,tlat)
                self.DEN[t]=gsw.rho(SA,self.TEM[t],p)
            
            if False:
                self.gstore('DEN')

            del self.SAL
            del self.TEM

        elif (i is not None) and (j is not None) and (k is None):
            # need to be computed only for one profile
            p  = gsw.p_from_z(self.z[:,j,i],tlat[:,j,i])
            print("p.shape")
            print(p.shape)
#            self.DEN = ma.empty(self.SAL.shape)
#            exec('self.testvar(SALi'+str(i)+'j'+str(j)+')')
#            exec('self.testvar(TEMi'+str(i)+'j'+str(j)+')')
            self.gload('SAL',i=i,j=j)
            self.gload('TEM',i=i,j=j)
            exec('SALloc=self.SALi'+str(i)+'j'+str(j))
            exec('TEMloc=self.TEMi'+str(i)+'j'+str(j))
            for t in xrange(len(self.time)):
                if t==0:
                    print("SALloc.shape")
                    print(SALloc[t].shape)
                    print("TEMloc.shape")
                    print(TEMloc[t].shape)
                SA = gsw.SA_from_SP(SALloc[t],p,self.lon[i],self.lat[j])
                DENloc=gsw.rho(SA,TEMloc[t],p)
                exec('self.DENi'+str(i)+'j'+str(j)+'=DENloc')

            if False:
                self.gstore('DEN')

############################################################################
# UTILITY : unload to free some memory


############################################################################
# VARIABLE SSS

    def instance_SSS(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.gload('SAL',k=30) 
            self.SSS=self.SALk30
            del self.SALk30
            self.SSS=ma.expand_dims(self.SSS,1)
            
        else:
            print(' Case not ready, please code .. ')

#        self.SSSinfo.units='p.s.u.'
#        self.SSSinfo.longname='Surface Salinity'
#        self.SSSinfo.dims='2D'

############################################################################
# VARIABLE SST

    def instance_SST(self, i=None,j=None, k=None):
        if (i is None) and (j is None) :
            self.gload('TEM',k=30) 
            self.SST=self.TEMk30
            del self.TEMk30
            self.SST=ma.expand_dims(self.SST,1)
            
        else:
            print(' Case not ready, please code .. ')

#        self.SSSinfo.units='C'
#        self.SSSinfo.longname='Surface Temperature'
#        self.SSSinfo.dims='2D'



############################################################################                                                                                                                                       # VARIABLE CCCn

    def instance_CCCn(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.testvar('DEN')
            self.testvar('TEM')
            self.CC   = ma.masked_where( (self.DEN<1014.5)|(self.TEM > 8.35), self.DEN * (self.TEM-8.35) * 3985 / 1e6 )
            self.CCCn = self.vertint('CC')
            self.CCCn = ma.expand_dims(self.CCCn,1)
        else:
            print(' Case not ready, please code .. ')


############################################################################
# VARIABLE PEA : Potential Energy Anomaly

    def instance_PEA(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.testvar('TEM')
            self.testvar('SAL')
            self.testvar('DEN')

            self.AVRDEN = np.tile(self.vertmean('DEN'),(1,31,1,1))

            for t in xrange(len(self.time)):
                self.PEAv[t] = 9.81 *self.z*(self.AVRDEN[t]-self.DEN[t])

            self.PEA = self.vertint(PEAv)
            self.PEA= ma.expand_dims(self.PEA,1)

        else:
            print(' Case not ready, please code .. ')


############################################################################                                                                                                             
# VARIABLE : Mixed layer depth  
#
# NOT READY !!! 
#
#        def instance_MLD(t):
#            #  Might probably be ipmroved a lot.                                                                                                                                                                  #             #  Currently, there are two 1D interpolation by i,j
#            Dloc=DEN[t]
#            zloc=z[t]
#            MLDloc=Dloc[30].copy()*0 # for init                                                                                                             
#            for i in xrange(Dloc.shape[1]):
#                for j in xrange(Dloc.shape[2]):
#                    if ma.is_masked(MLDloc[i,j]):
#                        continue
#                    f = interp1d(zloc[:,i,j], Dloc[:,i,j])
#                    d3=f(-3)
#                    f = interp1d(Dloc[:,i,j], zloc[:,i,j])
#            MLDloc[i,j]=f(d3+deltasig)#
#
#            return MLDloc


############################################################################   
# VARIABLE tracer age ... should be made generic for tracer number

    def instance_T1age(self, i=None, j=None, k=None):

        if (i is None) and (j is None) and (k is None):
            self.testvar('T1A')
            self.testvar('T1C')
            self.T1age=ma.masked_where(self.T1C<1e-2, self.T1A/self.T1C)
        else:
            self.gload('T1A',i=i,j=j)
            self.gload('T1C',i=i,j=j)
            exec('T1Aloc=self.T1Ai'+str(i)+'j'+str(j))
            exec('T1Cloc=self.T1Ci'+str(i)+'j'+str(j))
            exec('self.T1agei'+str(i)+'j'+str(j)+'=ma.masked_where(T1Cloc<1e-2, T1Aloc/ T1Cloc)')
            
############################################################################   
# VARIABLE Current Magnitude - [m/s]

    def instance_CUR(self, i=None, j=None, k=None):

        if (i is None) and (j is None) and (k is None):
            self.testvar('U')
            self.testvar('V')
            self.CUR = np.sqrt( self.U*self.U + self.V*self.V) 
        else:
            self.gload('U',i=i,j=j)
            self.gload('V',i=i,j=j)
            exec('Uloc=self.Ui'+str(i)+'j'+str(j))
            exec('Vloc=self.Vi'+str(i)+'j'+str(j))
            exec('self.CURi'+str(i)+'j'+str(j)+'=np.sqrt( Uloc*Uloc + Vloc*Vloc)')            

############################################################################   
# VARIABLE Oxygen Solubility - [mmol/m3]
# http://www.helcom.fi/documents/action%20areas/monitoring%20and%20assessment/manuals%20and%20guidelines/manual%20for%20marine%20monitoring%20in%20the%20combine%20programme%20of%20helcom_partb_annexb8_appendix3.pdf

    def instance_O2sat(self, i=None, j=None, k=None):

        I= -135.90205
        J= 1.575701e5
        K= -6.642308e7
        L= 1.243800e10
        M= -8.621949e11
        N= 0.017674
        P= -10.754
        Q= 2140.7
        
        if (i is None) and (j is None) and (k is None):
            self.testvar('TEM')
            self.testvar('SAL')
            TK = self.TEM+273.15
            self.O2sat = np.exp(I + J/TK + K/TK**2 + L/TK**3 + M/TK**4 - self.SAL*(N + P/TK + Q/TK**2))
        else:
            self.gload('TEM',i=i,j=j)
            self.gload('SAL',i=i,j=j)
            exec('TK=self.TEMi'+str(i)+'j'+str(j)+' + 273.15')
            exec('S=self.SALi'+str(i)+'j'+str(j))
            loc = np.exp(I + J/TK + K/TK**2 + L/TK**3 + M/TK**4 - S*(N + P/TK + Q/TK**2))
            exec('self.O2sati'+str(i)+'j'+str(j)+'=loc')            

############################################################################   
# VARIABLE Oxygen Saturation Rate - [%]
# http://www.helcom.fi/documents/action%20areas/monitoring%20and%20assessment/manuals%20and%20guidelines/manual%20for%20marine%20monitoring%20in%20the%20combine%20programme%20of%20helcom_partb_annexb8_appendix3.pdf

    def instance_pO2sat(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.testvar('DOX')
            self.testvar('O2sat')
            self.pO2sat = self.DOX/self.O2sat*100
        else:
            self.gload('DOX',i=i,j=j)
            self.gload('O2sat',i=i,j=j)
            exec('self.pO2sati'+str(i)+'j'+str(j)+'=self.DOXi'+str(i)+'j'+str(j)+'/self.O2sati'+str(i)+'j'+str(j)+'*100')
            
#############################################################################
# VARIABLE Consecutive days below given concentration - [d]
# I can only code the very specific case I need right now...

    def instance_H60(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.testvar('DOX')
        else:
            self.gload('DOX',i=i,j=j)            
            exec('DOXloc=self.DOXi'+str(i)+'j'+str(j))
            botDOX=DOXloc[:,11]
            loc=np.zeros_like(botDOX)
            for di,d in enumerate(botDOX):
                if d>60:
                    loc[i]=0
                else:
                    if di==0  :
                        loc[di]=1
                    elif loc[di-1]==0:
                        loc[di]=1
                    else:
                        loc[di]=loc[di-1]+1
            exec('self.H60i'+str(i)+'j'+str(j)+'=loc[:,None]') 
            print('done')           

#############################################################################
# VARIABLE Consecutive days below given concentration - [d]
# I can only code the very specific case I need right now...

    def instance_H120(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.testvar('DOXk11')
            botDOX=self.DOXk11
            loc=np.zeros_like(botDOX)
#            for i in botDOX.shape[2]:
#                for j in botDOX.shape[]:
            for di,d in enumerate(botDOX):
            # loc[di] = ma.where (d>120,0, if .... ) 
                for i in range(d.shape[0]):
                    for j in range(d.shape[1]):
                        if d[i,j]>120:
                            loc[di,i,j]=0 # 29/08/2018 corrected from loc[i] AFTER provision of forcings for habitat. Hope it doesn't affect..
                        else:
                            if di==0  :
                                loc[di,i,j]=1
                            elif loc[di-1,i,j]==0:
                                loc[di,i,j]=1
                            else:
                                loc[di,i,j]=loc[di-1,i,j]+1
            self.H120=loc[:,None,:,:]
            
        else:
            self.gload('DOX',i=i,j=j)            
            exec('DOXloc=self.DOXi'+str(i)+'j'+str(j))
            botDOX=DOXloc[:,11]
            loc=np.zeros_like(botDOX)
            for di,d in enumerate(botDOX):
                if d>120:
                    loc[i]=0
                else:
                    if di==0  :
                        loc[di]=1
                    elif loc[di-1]==0:
                        loc[di]=1
                    else:
                        loc[di]=loc[di-1]+1
            exec('self.H120i'+str(i)+'j'+str(j)+'=loc[:,None]') 
            print('done')           
