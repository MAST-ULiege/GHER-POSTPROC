import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import gsw
import yaml
import os.path
import datetime as dt
import cmocean
import calendar
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#from mpl_toolkits.basemap import Basemap

class G3D(object): 
    '''This is a class for model outputs exploration. It is based on GHER3D model outputs, but child classes are available for other models.'''

######################################################################

    def __init__(self,infile, YAML_FILE = 'local.yml'):

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
#            YAML_FILE = 'local.yml'
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

        self.ksurface     = len(self.sigII)+len(self.sigI)-2-1
        self.kbottom      = ma.where(self.bat>self.hlim,0,len(self.sigII)-1)

    # instantiate the dictionnary with model parameters
#        try: 
        self.initparams(config['PARAMFILE'])

 #       except:
  #      print('I did not found any model paramter file')

######################################################################
# UTILITY : Instantiate a dictionary with model parameter values. 

    def initparams(self,paramfile):
            # first the param file (default values)
        d={}
        fp = open(paramfile, 'r') 
        line='a'
        while line:
            line = fp.readline()
#            print(line)
            if '::' in line.split('!')[0]:
                paramline=line.split('!')[0].split('::')[1]
                p,v=paramline.split('=')
                try: # deals with all numeric values
                    daytosecond=86400.0
                    exec('v='+v)
                    d[p.strip()]=float(v)
                except: # probably a useless boolean or string, so I don't go further at the moment
                    p,v=paramline.split('=')
                    d[p.strip()]=v
        self.paramd=d

######################################################################
# VARIABLE : BAT

    def instance_bat(self):
        if os.path.isfile(self.batfile):
            print(self.batfile + ' -> OK')
        else:
            print(self.batfile + ' can not be found')

        with Dataset(self.batfile,'r') as nc:
            self.bat = nc.variables['bat'][:]
        self.bat=self.bat[None,:,:,:]
            
######################################################################
# VARIABLE : Z

    def instance_z(self):
        # For now let's build a constant z,
        # dynamic z considering ETA can be done later
        with Dataset(self.infile,'r') as nc:
            self.z   = nc.variables['depth'][1][None,:,:,:] # Make it 4D, even if we don't use elevation for now           
            self.lon = nc.variables['longitude'][:]
            self.lat = nc.variables['latitude'][:]

######################################################################
# VARIABLE : ZI
            
    def instance_zi(self):
        # For now let's build a constant dz,
        # dynamic dz considering ETA can be done later
        self.zi=ma.zeros(self.z.shape+np.array([0,1,0,0]))
        for k in xrange(0,12):
            self.zi[0,k]= - ma.masked_where( self.bat<=self.hlim,  (( self.bat - self.hlim ) * (1-self.sigII[k]) ) + self.hlim )
        for k in xrange(12,self.zi.shape[0]):
            self.zi[0,k]= - np.minimum(self.bat,self.hlim) * (1-self.sigI[k-12+1])
            
######################################################################
# VARIABLE : DZ

    def instance_dz(self):
        # For now let's build a constant z,
        # dynamic z considering ETA can be done later
        self.dz = ma.abs(self.zi[0,1:]-self.zi[0,:-1]) 

######################################################################
# UTILITY : LOAD LOCAL 
# The idea is that ouput always 4 dimension [time,depth,lat,lon], even if some are singleton.
# TODO simplify the code to avoid copying lines uselessly
# TODO 

    def gload(self,varname,k=None,i=None,j=None):
        try:
            # looking for the variable in the original file
            with Dataset(self.infile,'r') as nc:
                if  (k is None) and (i is None) and (j is None):
                    l=nc.variables[varname][:]
                    if (len(l.shape)==3):
                        print(varname + ' is 2+1D')
                        exec('self.'+varname+ '= nc.variables[varname][:,None,:,:]')
                    else:
                        exec('self.'+varname+ '= nc.variables[varname][:]')

                    print( 'Just loaded '+ (varname) + ' for full')
                elif (k is None) and (i is not None) and (j is not None):
                    print( 'Loading '+ (varname) + ' for i:'+str(i)+' and j:'+str(j))
#                    try:
                    l=nc.variables[varname][:]
                    if (len(l.shape)==4):
                        print(varname + 'varname is 3+1D') 
                        exec('self.'+varname+'i'+str(i)+'j'+str(j)+'= l[:,:,j,i]')#nc.variables[varname][:,:,j,i]')
                    elif (len(l.shape)==3):
                        print(varname + ' is 2+1D')
                        exec('self.'+varname+'i'+str(i)+'j'+str(j)+'= l[:,None,j,i]')#nc.variables[varname][:,:,j,i]')
                    else:
                        print('cannot understand diomension of %s'%varname)

                        

                    print( 'Just loaded '+ (varname) + ' for i:'+str(i)+' and j:'+str(j))
                elif (k is not None) and (i is None) and (j is None):
                    if (k=="surface"):
                        exec('self.'+varname+'ksurface= nc.variables[varname][:,self.ksurface]')
                        exec('self.'+varname+'ksurface=self.'+varname+'ksurface[:,None,:,:]')
                    elif (k=='bottom'):
#                        exec('self.'+varname+'kbottom=ma.empty([len(self.time),1,self.bat.shape[2], self.bat.shape[3]])')
                        exec('self.'+varname+ '= nc.variables[varname][:]')
                        exec('self.'+varname+'kbottom=ma.empty_like(self.'+varname+'[:,self.ksurface])[:,None,:,:]')
                        for i in range(self.bat.shape[2]):
                            for j in range(self.bat.shape[3]):
                                if (not ma.is_masked(self.kbottom[0,0,i,j])):
                                    exec('self.'+varname+'kbottom[:,0,i,j]= self.'+varname+'[:,self.kbottom[0,0,i,j],i,j]')
                                    
                    else:                        
                        exec('self.'+varname+'k'+str(k)+'= nc.variables[varname][:,k]')
                        exec('self.'+varname+'k'+str(k)+'=self.'+varname+'k'+str(k)+'[:,None,:,:]')
                    print( 'Just loaded '+ (varname) +' for k:'+str(k))
                else:
                    print(' Strange case encountered in  G3D_class.py : def gload (i:%,j:%;k:%)'%(i,j,k))
                    
        except: 
            try:
                # looking for the variable in the diag file (created by function gstore)
                with Dataset(self.diagfile,'r') as nc:
                    if (k is None) and (i is None) and (j is None):
                        exec('self.'+varname+ '= nc.variables[varname][:]')
                        print( 'Just loaded '+ (varname) + ' for full')
                    elif (k is None) and (i is not None) and (j is not None):
                        print( 'Loading '+ (varname) + ' for i:'+str(i)+' and j:'+str(j))
                        exec('self.'+varname+'i'+str(i)+'j'+str(j)+'= nc.variables[varname][:,:,j,i]')
                        print( 'Just loaded '+ (varname) + ' for i:'+str(i)+' and j:'+str(j))
                    elif (k is not None) and (i is None) and (j is None):
                        exec('self.'+varname+'k'+str(k)+'= nc.variables[varname][:,k]')
                        print( 'Just loaded '+ (varname) +' for ki:'+str(k))
                    else:
                        print(' Stange case encountered in  G3D_class.py : def gload (i:%,j:%;k:%)'%(i,j,k))


            except Exception as e: 
                # looking for an instance function allowing to define this variable
                print(e)
                print( '\n %s not found in %s'%(varname,self.infile))
                print( '-> Calling')
                #try:
                print('     self.instance_'+varname+'(i='+str(i)+',j='+str(j)+',k='+str(k)+')')
                exec('self.instance_'+varname+'(i=i,j=j,k=k)')

#                exec('loc=self.+varname+')



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
        exec('avg=self.'+varname+'.mean(axis=0)')
        return(avg)
        
######################################################################
# PROCESS : FULL AVERAGE

    def avgspatial (self,varname,maskin=None):
        # in the sense of volumetric mean
        self.testz()
        self.testvar(varname)        
                        
        avg=ma.empty(len(self.dates))
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
                if len(maskin.squeeze().shape)==2:
                    print('static 2D mask .. OK')
                    for t in xrange(len(self.time)):
                        loc[t]=ma.masked_where(maskin.squeeze()[None,:,:],loc[t])
                elif ((len(maskin.squeeze().shape)==3) and (maskin.shape[0]==self.time.shape)):
                    print('dynamic 2D mask .. OK')
                    for t in xrange(len(self.time)):
                        loc[t]=ma.masked_where(maskin.squeeze()[t],loc[t])
                else:
                    print('Mask dimension not understood: %s'%(maskin.shape) )

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
            
    def testvar(self,varname,doload=True,i=None,j=None, k=None):
        try:
            if (i is None) and (j is None) and (k is None):
                exec('self.'+varname)
            elif (k is not None):
                exec('self.'+varname+'k'+str(k))
            elif (i is not None) and (j is not None):
                exec('self.'+varname+'i'+str(i)+'j'+str(j))
            isthere=True
        except:
            print('%s not found '%(varname))
            isthere=False
            if doload:
                print ('Loading %s'%(varname) )
                self.gload(varname,i=i,j=j,k=k)
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
            gridZU = self.zi[0,1:]
            gridZD = self.zi[0,:-1]
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
            gridZU = self.zi[0,1:]
            gridZD = self.zi[0,:-1]
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

        zforplot=self.z[0,:,j,i]
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
        isvar=self.testvar(varname,doload=False)
        self.testtime()
        i,j=self.test_coord(c1,c2)

        zloc=self.z[0,:,j,i]
        if c3=='bottom':
            k=self.kbottom
        else: 
            k=(np.abs(zloc-c3)).argmin()
        if isvar:
            # for some reason the 4d variable is already loaded in memory
            exec('lloc=self.'+varname)
            lloc=lloc[:,:,j,i]
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
        vint=ma.empty( (loc.shape[0],1,loc.shape[2],loc.shape[3]) )

        for t in xrange(len(self.time)):
            vint[t]=ma.sum(loc[t]*self.dz[0],0)
        
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

    def mapMonthlyClim(self, varname,title=None,cmapname='haline',Clim=None,figsuffix='', batlines=True, subdomain=None):

        exec('loc=self.clim_'+varname)
        loclon=self.lon
        loclat=self.lat
        locbat=self.bat[0,0]
        if (subdomain is not None):
            if (subdomain=='NWS'):
                limlonll,limlatll = self.test_coord(28.0,43.0)
                limlonur,limlatur = self.test_coord(33.5,47.0)
            elif (subdomain=='BOSP'):
                limlonll,limlatll = self.test_coord(27.7,41.0)
                limlonur,limlatur = self.test_coord(30.8,41.7)
            else:
                print('subdomain unknwon')
            print('coord for subdomain : l1:%s, l2:%s, L1:%s, L2:%s'%(limlatll,limlatur,limlonll,limlonur))
            loc    = loc[:,:,limlatll:limlatur,limlonll:limlonur]
            loclat = self.lat[limlatll:limlatur]
            loclon = self.lon[limlonll:limlonur]
            locbat = locbat[limlatll:limlatur,limlonll:limlonur]
        if (loc.shape[1]>1):
            print('!! use mapMonthlyClim for 2D clims only !!')
            print('!! proceeding now for the surface layer ''')
            loc=loc[:,self.ksurface][:,None,:,:]
        if (title==None): title=varname
        if Clim==None : Clim=[loc.min(),loc.max()]
        exec('cmap=cmocean.cm.'+cmapname)

        fig, aaxes = plt.subplots(4,3,figsize=(10, 12))

        parallels = np.arange(np.floor(min(loclat)),np.ceil(max(loclat)),1.)
        meridians = np.arange(np.floor(min(loclon)),np.ceil(max(loclon)),1.)
        llon,llat = np.meshgrid(loclon,loclat)

        for monthi in range(1,12+1):
            indx = [i for i,x in enumerate(self.climdates) if (x.month==monthi)]
            
            try: # if BaseMap is installed and OK
                m    = Basemap(llcrnrlat=loclat[0],urcrnrlat=loclat[-1],llcrnrlon=loclon[0],urcrnrlon=loclon[-1],\
                                   resolution='i',ax=aaxes[int(np.ceil((monthi-1)/3)),(monthi-1)%3 ])
                xx, yy = m(llon,llat)
                m.drawcoastlines()
                m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
                m.drawmeridians(meridians,labels=[1,0,0,0],fontsize=10)
                cs = m.contourf(xx,yy,loc[indx,0].mean(axis=0),cmap=cmocean.cm.deep, extend="both")
                if batlines: 
                    m.contour(xx,yy,G.bat[0,:,:140],levels=[40,80,120], colors='k',linestyles='dashed')
            except:
                cs = aaxes[int(np.ceil((monthi-1)/3)),(monthi-1)%3 ].contourf(loclon, loclat,loc[indx,0].mean(axis=0),\
                                                                                  levels= np.linspace(Clim[0],Clim[1],20),\
                                                                                  cmap=cmap, extend="both")
                if batlines:
                    aaxes[int(np.ceil((monthi-1)/3)),(monthi-1)%3 ].contour(loclon, loclat,locbat,levels=[40,80,120], colors='k',linestyles='dashed')

            aaxes[int(np.ceil((monthi-1)/3)),(monthi-1)%3 ].set_title(calendar.month_name[monthi])

        fig.subplots_adjust(hspace=0.1,wspace=0.1, bottom=0.1, right=0.95, left=0.05, top=0.95)
        cbar_ax = fig.add_axes([0.1, 0.04, 0.8, 0.03])
        cbar    = fig.colorbar(cs,ticks=np.linspace(Clim[0],Clim[1],10),cax=cbar_ax, orientation="horizontal")
        cbar.set_label(varname)

        fig.savefig(self.figoutputdir+'MonthlyClim_'+varname+figsuffix+'.png')
        fig.close()

############################################################################
# PLOTS : Plot Map 

    def mapStrip(self, varname,title=None,cmapname='haline',Clim=None,figsuffix='', batlines=True, subdomain=None, daysbetween=31,extend="max", diff=False, difflag='init'):
        
        exec('loc=self.'+varname+'.copy()')

        loclon=self.lon
        loclat=self.lat
        locbat=self.bat[0,0]
        if (subdomain is not None):
            if (subdomain=='NWS'):
                limlonll,limlatll = self.test_coord(28.0,43.0)
                limlonur,limlatur = self.test_coord(33.5,47.0)
            elif (subdomain=='BOSP'):
                limlonll,limlatll = self.test_coord(27.7,41.0)
                limlonur,limlatur = self.test_coord(30.9,41.7)
            else:
                print('subdomain unknwon')
            print('coord for subdomain : l1:%s, l2:%s, L1:%s, L2:%s'%(limlatll,limlatur,limlonll,limlonur))
            loc    = loc[:,:,limlatll:limlatur,limlonll:limlonur]
            loclat = self.lat[limlatll:limlatur]
            loclon = self.lon[limlonll:limlonur]
            locbat = locbat[limlatll:limlatur,limlonll:limlonur]

        if (loc.shape[1]>1):
            print('!! use mapStrip for 2D-stacks only !!')
            print('!! proceeding now for the surface layer .. ''')
            loc=loc[:,self.ksurface][:,None,:,:]
        if (title==None): title=varname
        if Clim==None : Clim=[loc.min(),loc.max()]
        exec('cmap=cmocean.cm.'+cmapname)
        # usefull to check specific evolution of a variable
        if diff:
            locinit=loc[0].copy()
            for t in range(loc.shape[0]):
                if difflag=='init':
                    loc[t]=loc[t]-locinit
                else:
                    locorig=loc.copy()
                    try: 
                        loc[t]=loc[t]-locorig[max(t-difflag,0)]
                    except TypeError:
                        Print( 'difflag should be ''init'' or an integer (for now timesteps)')
            Clim=[-max(abs(loc.min()),abs(loc.max())), max(abs(loc.min()),abs(loc.max()))]
            exec('cmap=cmocean.cm.'+'balance')

        # computing number of sub-plots
        nframe=int(np.floor(len(self.dates)/daysbetween)) # TODO, should instead allow for different output time-steps, but still consider days and makes average when neeeded
        cols = int(np.ceil(np.sqrt(nframe)))
        rows = int(1 + (nframe - 1)//np.ceil(np.sqrt(nframe)))
 
        fig, aaxes = plt.subplots(nrows=rows, ncols=cols, figsize=(10, 12), squeeze=False)
 
        parallels = np.arange(np.floor(min(loclat)),np.ceil(max(loclat)),1.)
        meridians = np.arange(np.floor(min(loclon)),np.ceil(max(loclon)),1.)
        llon,llat = np.meshgrid(loclon,loclat)

        for fi in range(1,nframe+1):
            indx = range ((fi-1)*daysbetween,fi*daysbetween)
            spi=int(np.ceil((fi-1)/cols))
            spj=(fi-1)%cols
         #   print('ind %s, ind %s, spi %s, spj %s'%(fi,indx,spi,spj))
            try: # if BaseMap is installed and OK
                m    = Basemap(llcrnrlat=loclat[0],urcrnrlat=loclat[-1],llcrnrlon=loclon[0],urcrnrlon=loclon[-1],\
                                   resolution='i',ax=aaxes[int(np.ceil((fi-1)/cols)),(fi-1)%cols ])
                xx, yy = m(llon,llat)
                m.drawcoastlines()
                m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
                m.drawmeridians(meridians,labels=[1,0,0,0],fontsize=10)
                cs = m.contourf(xx,yy,loc[indx,0].mean(axis=0),cmap=cmocean.cm.deep, extend=extend)
                if batlines: 
                    m.contour(xx,yy,self.bat[0,:,:140],levels=[40,80,120], colors='k',linestyles='dashed')
               # if not (spi==0):
                    
            except:
                cs = aaxes[spi,spj].contourf(loclon, loclat,loc[indx,0].mean(axis=0),\
                                                                                  levels= np.linspace(Clim[0],Clim[1],20),\
                                                                                  cmap=cmap, extend=extend)
                if batlines:
                    aaxes[spi,spj ].contour(loclon, loclat,locbat,levels=[40,80,120, 500, 1000], colors='k',linestyles='dashed')

            aaxes[spi,spj].set_title(self.dates[indx[0]].strftime('%d%b')+'-'+self.dates[indx[-1]].strftime('%d%b'))
            if not (spj==0):
                 aaxes[spi,spj].get_yaxis().set_visible(False)
            if not (spi==(cols-1)):
                 aaxes[spi,spj].get_xaxis().set_visible(False)
                 
        fig.subplots_adjust(hspace=0.1,wspace=0.1, bottom=0.1, right=0.95, left=0.05, top=0.95)
        cbar_ax = fig.add_axes([0.1, 0.04, 0.8, 0.03])
        cbar    = fig.colorbar(cs,ticks=np.linspace(Clim[0],Clim[1],10),cax=cbar_ax, orientation="horizontal")
        cbar.set_label(varname)

        fig.savefig(self.figoutputdir+'Strip_'+varname+figsuffix+'.png')
        fig.close()
        
############################################################################
# PLOTS : Plot Map 

    def map2D(self, varname,cmapname="haline", subdomain=None, figout=None, title=None, Clim=None, scatmat=None,extend='max', batlines=True):
        self.testz()
        exec('loc=self.'+varname)
        loclon=self.lon
        loclat=self.lat
        locbat=self.bat[0,0]
        if (subdomain=="NWS"):
            limlon,limlat = self.test_coord(33.5,43.0)
            loc    = loc[:,:,limlat:,:limlon]
            loclat = self.lat[limlat:]
            loclon = self.lon[:limlon]
            locbat = locbat[limlat:,:limlon]
        exec('cmap=cmocean.cm.'+cmapname)
        if figout==None:
            figout='MAP_'+varname
        if (title==None): title=varname
        if Clim==None : Clim=[loc.min(),loc.max()]

        parallels = np.arange(40.,48.,1.) # TODO: generalize
        meridians = np.arange(20.,43.,1.)
        llon,llat=np.meshgrid(loclon,loclat)

        fig = plt.figure(figsize=(10,8))
        ax = fig.add_axes([0.1,0.1,0.8,0.8])

        try: # if BaseMap is installed and OK
            m = Basemap(llcrnrlat=loclat[0],\
                            urcrnrlat=loclat[-1],\
                            llcrnrlon=loclon[0],\
                            urcrnrlon=loclon[-1],\
                            resolution='l')
            m.drawcoastlines()
            m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
            m.drawmeridians(meridians,labels=[1,0,0,0],fontsize=10)
            x, y = m(llon,llat)        
            cs=m.contourf(x,y,loc.squeeze(), np.linspace(Clim[0],Clim[1],100),cmap=cmap, extend=extend)
            if batlines:
                m.contour(x,y,self.bat[0,j1:j2,i1:i2],levels=[40,80,120], colors='k',linestyles='dashed')
        except:
            cs = ax.contourf(loclon, loclat,loc.squeeze(),levels= np.linspace(Clim[0],Clim[1],20),cmap=cmap, extend=extend)
            if batlines:
                ax.contour(loclon, loclat,locbat,levels=[40,80,120,500,1000], colors='k',linestyles='dashed')
        cbar_ax = fig.add_axes([0.1, 0.04, 0.8, 0.03])
#        cbar    = fig.colorbar(cs,ticks=np.linspace(Clim[0],Clim[1],10),cax=cbar_ax, orientation="horizontal")
        cbar    = fig.colorbar(cs,ticks=np.linspace(Clim[0],Clim[1],11),cax=cbar_ax, orientation="horizontal") 
        cbar.set_label(title)
        if (scatmat!=None):
            try:
                m.scatter(scatmat[:,1],scatmat[:,0],50,scatmat[:,2],edgecolors='black',cmap=cmap,vmin=Clim[0],vmax=Clim[1])
            except:
                ax.scatter(scatmat[:,1],scatmat[:,0],50,scatmat[:,2],edgecolors='black',cmap=cmap,vmin=Clim[0],vmax=Clim[1])
        
        #plt.title(title)
        #m.colorbar()
        fig.savefig(figout)
        fig.savefig(self.figoutputdir+figout+'.png')
        fig.close()

############################################################################
# PLOTS : Plot Vertical Profile

    def plotprofile (self, varname, z=None,cmapname="haline", figout=None, title=None, Clim=None, zlim=None):
        if figout==None:
            figout=varname
        exec('loc=self.'+varname)
        if z.all()==None:
            z=range(loc.shape[1])
        if zlim!=None:
            zimin=abs(z-zlim[0]).argmin()
            zimax=abs(z-zlim[1]).argmin()
            loc=loc[:,min(zimin,zimax):max(zimin,zimax)]
            z=z[min(zimin,zimax):max(zimin,zimax)]
        if Clim==None: 
            Clim=[loc.min(),loc.max()]
        exec('cmap=cmocean.cm.'+cmapname)

        print(loc.shape)
#        print(z.shape)

        locator = mdates.AutoDateLocator()
        formator = mdates.AutoDateFormatter(locator)
        
        fig=plt.figure(figsize=(15, 8))
        ax=plt.subplot(1, 1, 1)
        ax.xaxis_date()
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formator)
        cs=plt.contourf(self.dates, z, loc.transpose(), np.linspace(Clim[0],Clim[1],20), extend="max", cmap=cmap )
        cbar_ax = fig.add_axes([0.1, 0.04, 0.8, 0.03])
        cbar    = fig.colorbar(cs,ticks=np.linspace(Clim[0],Clim[1],11),cax=cbar_ax, orientation="horizontal")
        fig.savefig(self.figoutputdir+figout+'.png')
        fig.close()

############################################################################
# PLOTS : Plot Time Series

    def plotseries (self, varname, figout=None, title=None, Clim=None):
        if figout==None:
            figout=varname
        exec('loc=self.'+varname)
        if Clim==None: 
            Clim=[loc.min(),loc.max()]

        locator = mdates.AutoDateLocator()
        formator = mdates.AutoDateFormatter(locator)
        
        fig=plt.figure(figsize=(10, 4))
        ax=plt.subplot(1, 1, 1)
        ax.xaxis_date()
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formator)
        cs=plt.plot(self.dates, loc)
        plt.title(title)
        fig.savefig(self.figoutputdir+'TimeSeries_'+figout+'.png')
        fig.close()

############################################################################
# VARIABLE : Density

    def instance_DEN(self,i=None,j=None,k=None):

        self.testz()
        tlat=np.tile(self.lat, (self.z.shape[1],self.z.shape[3],1))
        tlat=np.transpose(tlat,(0,2,1))
        tlon=np.tile(self.lon, (self.z.shape[1], self.z.shape[2],1))

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
            p  = gsw.p_from_z(self.z[0,:,j,i],tlat[:,j,i])
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
# VARIABLE Sea Bottom Salinity

    def instance_SBS(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.gload('SAL',k='bottom')
            self.SBS=self.SALkbottom
            del self.SALkbottom
        else:
            print(' Case not ready, please code .. ')

############################################################################
# VARIABLE Sea Bottom Temperature
                                                                                                                                                 
    def instance_SBT(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.gload('TEM',k='bottom')
            self.SBT=self.TEMkbottom
            del self.TEMkbottom
        else:
            print(' Case not ready, please code .. ')


############################################################################                                                                                                      
# VARIABLE Sea Bottom Oxygen
    def instance_O2bottom(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.gload('DOX',k='bottom')
            self.O2bottom=self.DOXkbottom
            del self.DOXkbottom
        else:
            print(' Case not ready, please code .. ')

############################################################################
# VARIABLE SSS

    def instance_SSS(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.gload('SAL',k='surface') 
            self.SSS=self.SALksurface
            del self.SALksurface            
        else:
            print(' Case not ready, please code .. ')

############################################################################
# VARIABLE SST

    def instance_SST(self, i=None,j=None, k=None):
        if (i is None) and (j is None) :
            self.gload('TEM',k='surface') 
            self.SST=self.TEMksurface
            del self.TEMksurface            
        else:
            print(' Case not ready, please code .. ')

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
    def instance_MLD(self, i=None,j=None,k=None):
         #  Might probably be ipmroved a lot.                                                                                                                                                                  #            #  Currently, there are two 1D interpolation by i,j
        self.testvar('DEN',i=i,j=j)
        self.testvar('DENat3',i=i,j=j)

        self.MLD = ma.ones(self.DEN[:,self.ksurface])[:,None,:,:]*3

        for t in range(len(self.dates)):        
            for i in xrange(Dloc.shape[1]):
                for j in xrange(Dloc.shape[2]):
                    if ma.is_masked(MLDloc[i,j]):
                        continue
                    f = interp1d(zloc[:,i,j], Dloc[:,i,j])
                    d3=f(-3)
                    f = interp1d(Dloc[:,i,j], zloc[:,i,j])

        MLDloc[i,j]=f(d3+deltasig)

        return MLDloc

###########################################################################

    def instance_DENat3(self, i=None,j=None,k=None):
        self.testvar('DEN',i=i,j=j)
        self.DENat3=self.DEN[:,self.ksurface].copy()[:,None,:,:]
        for t in range(len(self.dates)):
            for i in range(self.bat.shape[2]):
                for j in range(self.bat.shape[3]):
                    if ma.is_masked(self.DENat3[t,0,i,j]):
                        continue
                    self.DENat3[t,0,i,j]=np.interp(3,self.z[0,:,i,j],self.DEN[t,:,i,j])



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

############################################################################
# VARIABLE : totN, Total Nitrogen
    def instance_totN(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.totN = self.sumvar({'Nphyto':1,'Nzoo':1,'Ngel':1,'Nbac':1,'Ndis':1,'Norg':1},i,j,k)
            '''
            self.testvar('Nphyto')
            self.testvar('Nzoo')
            self.testvar('Ngel')
            self.testvar('Nbac')
            self.testvar('Ndis')
            self.testvar('Norg')
            self.totN=self.Nphyto+self.Nzoo+self.Ngel+self.Nbac+self.Ndis+self.Norg
            '''
        else:
            print('NEED TO BE COMPLETED : instance_totN')

############################################################################
# VARIABLE : Nphyto, Nitrogen in phyto form

    def instance_Nphyto(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.Nphyto = self.sumvar({'NDI':1,'NEM':1,'NFL':1},i,j,k)
        else:
            print('NEED TO BE COMPLETED : instance_Ndis')
            '''
            self.testvar('NDI')
            self.Nphyto = self.NDI
            if self.sparemem: del self.NDI
            self.testvar('NEM')
            self.Nphyto = self.Nphyto+self.NEM
            if self.sparemem: del self.NEM
            self.testvar('NFL')
            self.Nphyto = self.Nphyto+self.NFL
            if self.sparemem: del self.NFL
        else:
            print('NEED TO BE COMPLETED : instance_Nphyto')
            '''
############################################################################ 
# VARIABLE : Nzoo, Nitrogen in zoo form

    def instance_Nzoo(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.Nzoo = self.sumvar({'MES':self.paramd['NCrMesoZoo'] ,'MIC':self.paramd['NCrMicroZoo']},i,j,k)
            '''
            self.testvar('MES')
            self.Nzoo = self.MES*self.paramd['NCrMesoZoo']
            if self.sparemem: del self.MES
            self.testvar('MIC')
            self.Nzoo = self.Nzoo+self.MIC*self.paramd['NCrMicroZoo']
            if self.sparemem: del self.MIC
            '''
        else:
            print('NEED TO BE COMPLETED : instance_Nzoo')

############################################################################ 
# VARIABLE : Ngel, Nitrogen in gel form

    def instance_Ngel(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.Ngel = self.sumvar({'GEL':self.paramd['NCrGelatinous'] ,'NOC':self.paramd['NCrNoctiluca']},i,j,k) 
            '''
            self.testvar('NOC')
            self.Ngel = self.NOC*self.paramd['NCrNoctiluca']
            if self.sparemem: del self.NOC
            self.testvar('GEL')
            self.Ngel = self.NGel+self.NOC*self.paramd['NCrNoctiluca']
            if self.sparemem: del self.GEL
            '''
        else:
            print('NEED TO BE COMPLETED : instance_Ngel')

############################################################################
# VARIABLE : Ndis, Nitrogen in dissolved inorganic form
    '''
    def instance_Ndis(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.testvar('NOS')
            self.testvar('NHS')
            self.Ndis = self.NHS+self.NOS
#            del self.NOS, self.NHS
        else:
            print('NEED TO BE COMPLETED : instance_Ndis')
    '''

###########################################################################

    def instance_Ndis(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.Ndis = self.sumvar({'NOS':1,'NHS':1},i,j,k)
        else:
            print('NEED TO BE COMPLETED : instance_Ndis')

###########################################################################

    def sumvar(self, sumdic, i=None, j=None, k=None):
        out=ma.zeros(self.time.shape+self.z.shape[1:])
        for vvar, factor in sumdic.items():
            self.testvar(vvar)
            exec('out = ma.add(out, self.'+vvar+'*factor)')
            if self.sparemem: exec('del self.'+vvar)
        return(out)

############################################################################
# VARIABLE : Norg, Nitrogen in detritic organic form

    def instance_Norg(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.Norg = self.sumvar({'DNL':1,'DNS':1,'PON':1},i,j,k)
            '''
            self.testvar('DNL')
            self.testvar('DNS')
            self.testvar('PON')
            self.Norg = self.DNS+self.DNL+self.PON
            del self.DNS, self.DNL, self.PON
            '''
        else:
            print('NEED TO BE COMPLETED : instance_Norg')
            
############################################################################ 
# VARIABLE : Nbac, Nitrogen in bacteria form

    def instance_Nbac(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.testvar('BAC')
            self.Nbac = self.BAC*self.paramd['NCrBac']
            del self.BAC
        else:
            print('NEED TO BE COMPLETED : instance_Ndis')

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
            self.testvar('DOXk11') # TODO : replace with O2bottom
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

############################################################################ 
# VARIABLE : Ns, Nitrogen in dissolved inorganic form                     

    def instance_N2loss(self, i=None, j=None, k=None):
        if (i is None) and (j is None) and (k is None):
            self.testvar('AMANOX')
            self.testvar('DENITRIFICATION')
            self.testvar('OXIDATIONBYNITRATE')
            self.N2loss = self.AMANOX+self.DENITRIFICATION+self.OXIDATIONBYNITRATE
        else:
            print('NEED TO BE COMPLETED : instance_N2loss')


