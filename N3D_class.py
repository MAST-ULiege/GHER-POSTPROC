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
    def __init__(self,infile,YAML_FILE='local.yml'):
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

        self.ksurface=0


        with Dataset(self.batfile,'r') as nc:
            self.kbottom = nc.variables['mbathy'][:]  # 3D
            self.kbottom = self.kbottom[:,None,:,:]   # 4D
            self.z = nc.variables['gdept_0'][:]       # 3D # Wrong for now, we don't consider free surface !! 
            print('gdept_0 shape')
            print(self.z.shape)

            self.bat=ma.empty_like(self.kbottom[:])
            for i in range(self.bat.shape[3]):
                for j in range(self.bat.shape[2]):
                    self.bat[0,0,j,i]=self.z[0,self.kbottom[0,0,j,i],j,i]

            self.bat= ma.masked_where(self.kbottom==0,self.bat)

            self.dx= nc.variables['e1t'][:]
            self.dy= nc.variables['e2t'][:]
            
            self.landmask= nc.variables['tmask'][:]   # 3D 

        self.kbottom = ma.where(self.kbottom!=0, self.kbottom-1,self.kbottom)
######################################################################
# VARIABLE : Z
    def instance_z(self):
        # Dynamic z considering ETA

        with Dataset(self.batfile,'r') as nc:
            self.z   = nc.variables['gdept_0'][:]     # 3D # Wrong for now, we don't consider free surface !!           
            self.lon = nc.variables['glamt'][:][0,0,:]
            self.lat = nc.variables['gphit'][:][0,:,0]
            self.dz  = nc.variables['e3t_0'][:]
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
#####################################################################
# UTILITY : test variable
    def testvar(self,varname, doload=True,i=None,j=None, k=None):
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
                if (k is not None):
                    varname=varname+'k'+str(k)
                isthere=True
       #         if (i is None) and (j is None) and (k is None):
                    # NEMO NEEDS EXTRA MASKING
                print('remasking ' +varname)
                exec('self.'+varname+'=ma.masked_array(self.'+varname+',mask=False)')
                # This assumes that all variables comes with 4 dimension, eventually with only one level
                exec('nt=self.'+varname+'.shape[0]')
                for t in range(nt):
                    exec('self.'+varname+'[t]=ma.masked_where(self.landmask[0,:self.'+varname+'.shape[1]]==0, self.'+varname+'[t])')
                #elif (k is not None):
                #    print('remasking ' +varname+str(k))
#                else:
#                    print('remasking for N vars is not yet implemented for cases where i or j are not None')
        return(isthere)
######################################################################
# UTILITY : test time
    def testtime(self):
        try:
            self.time
        except:
            print('%s not found -> loading'%('time'))
            self.gload('time_counter') # In NEMO time_counter is the number of seconds since referene time 
            self.time=self.time_counter
            del self.time_counter
            with Dataset(self.infile,'r') as nc:
                t0 = nc.variables['time_counter'].time_origin 
            print('Time Origin : %s'%(t0))    
            self.dates = [dt.datetime.strptime(t0,'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=int(t)) for t in self.time]
######################################################################
# UTILITY : 
    def instance_SAL(self,i=None,j=None, k=None):
        self.testvar('vosaline',i=i,j=j, k=k)
        if k is not None:
            exec('self.SALk'+str(k)+'=self.vosalinek'+str(k))
            exec('del self.vosalinek'+str(k))
        else: 
            self.SAL=self.vosaline
            del self.vosaline

    def instance_TEM(self,i=None,j=None,k=None):
        self.testvar('votemper',i=i,j=j, k=k)
        if k is not None:
            exec('self.TEMk'+str(k)+'=self.votemperk'+str(k))
            exec('del self.votemperk'+str(k))
        else:
            self.TEM=self.votemper
            del self.votemper





