import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import gsw
import yaml
import os.path
import datetime as dt

import G3D_class

class N3D(G3D_class.G3D): 
    '''This is a class for NEMO model outputs exploration. It is based on G3D class but overides specific methods (for z, lon, lat definitions)'''

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
        self.verbose      = config['VERBOSE']
        self.figoutputdir = config['PLOTDIR']

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
            mbathy = nc.variables['mbathy'][:]     # 2D
            print('mbathy shape') 
            print(mbathy.shape)
            self.z  = nc.variables['gdept'][:]     # 3D
            print('gdept shape')
            print(self.z.shape)

            self.bat=ma.empty((1,1,mbathy.shape[1],mbathy.shape[2]))
            for i in range(mbathy.shape[2]):
                for j in range(mbathy.shape[1]):
                    self.bat[0,0,j,i]=self.z[0,mbathy[0,j,i],j,i]

            self.dx= nc.variables['e1t'][:]
            self.dy= nc.variables['e2t'][:]
            
            self.landmask= nc.variables['tmask'][:] # 3D 
######################################################################
# VARIABLE : Z

    def instance_z(self):
        # For now let's build a constant z,
        # dynamic z considering ETA can be done later

        with Dataset(self.batfile,'r') as nc:
            self.z   = nc.variables['gdept'][:] # Let it be 3D for now           
            self.lon = nc.variables['glamt'][:]
            self.lat = nc.variables['gphit'][:]
            self.dz  = nc.variables['e3t'][:]
######################################################################
# VARIABLE : ZI
            
    def instance_zi(self):
        # For now let's build a constant dz,
        with Dataset(self.batfile,'r') as nc:
            self.zi  = nc.variables['gdepw'][:]
      
######################################################################
# UTILITY : test z        
    def testz(self):
        try:
            self.dz
        except :
            print('dz not found -> loading')
            self.instance_z()
            
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
                # NEMO NEEDS EXTRA MASKING
                print('remasking ' +varname)
                exec('self.'+varname+'=ma.masked_array(self.'+varname+')')
                for t in range(len(self.time)):
                    exec('self.'+varname+'[t]=ma.masked_where(self.landmask[0,:self.'+varname+'.shape[1]], self.'+varname+'[t])')


        return(isthere)
                
            
######################################################################
# UTILITY : test time
            
    def testtime(self):
        try:
            self.time
        except:
            print('%s not found -> loading'%('time'))
            self.gload('time_counter')
            self.time=self.time_counter
            del self.time_counter
            with Dataset(self.infile,'r') as nc:
                t0 = nc.variables['time_counter'].time_origin 
            self.dates = [dt.datetime.strptime(t0,' %Y-%b-%d %H:%M:%S')+dt.timedelta(days=int(t)) for t in self.time]

######################################################################
# UTILITY : 
            
    def instance_SAL(self,i=None,j=None, k=None):
#        self.gload('vosaline',i,j,k)
        self.testvar('vosaline')
        print(type(self.vosaline))
        self.SAL=self.vosaline
        if k is not None:
            exec('self.SALk'+str(k)+'=self.SAL[:,k]')

        del self.vosaline

    def instance_TEM(self):
        self.gload('votemper')
        self.TEM=self.votamper
        del self.votemper


    def instance_SSS(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.gload('SAL',k=0)
            self.SSS=self.SALk0
            del self.SALk0
            self.SSS=ma.expand_dims(self.SSS,1)

        else:
            print(' Case not ready, please code .. ')
