# -*- coding: utf-8 -*-
import io
import numpy as np
from netCDF4 import num2date, date2num
from datetime import datetime, timedelta
import pandas as pd
from ..DropSizeDistribution import DropSizeDistribution
from . import common

def read_OTT_MtView(filename):
    '''
    Takes a filename pointing to a parsivel raw file and returns
    a drop size distribution object.
    Usage:
    dsd = read_parsivel(filename)
    Returns:
    DropSizeDistrometer object

    '''
    reader = OTTreaderRaw(filename)
    dsd = DropSizeDistribution(reader)
    return dsd


class OTTreaderRaw(object):
    """
    ParsivelReader class takes takes a filename as it's only argument(for now).
    This should be a parsivel raw datafile(output from the parsivel).

    """
    def __init__(self, filename):
        self.filename = filename
        self.rain_rate = []
        self.Z = []
        self.num_particles = []
        self._base_time = []
        self.nd = []
        self.vd = []
        self.raw = []
        self.code = []
        self.time = []
        self.dates=[]
        self.ndt = []
        self.pcm = np.reshape(self.pcm_matrix, (32, 32))
        self._read_file_MtView()
        self._get_epoch_time()
        self._prep_data()        
        self.bin_edges = np.hstack((0, self.diameter['data'] + np.array(self.spread['data'])/2))
        self.bin_edges = common.var_to_dict(
            'bin_edges',
            self.bin_edges,
            'mm', 'Bin Edges')
        self.filtered_raw_matrix=[]
        self.ndf=[]
        self._apply_pcm_matrix()
      




    def _read_file_MtView(self):
        """  Read the Parsivel Data file and store it in internal structure.
        Returns: None

        """
        read_file=pd.read_csv(self.filename,low_memory=False,skip_blank_lines=True,parse_dates=True) #Read monthly merged CSV file
        read_file=read_file[(read_file.iloc[:,25]>10)] #Filter for minimum number of particles (20 particles for 1 min)
        read_file=read_file[(read_file.iloc[:,26]>0.05)] #Filter for Rainfall rate (1-min time step >0.1mm/hr)
        read_file.index=pd.to_datetime(read_file.ix[:,0])
        self.ndt=read_file.index
        self.file=read_file
        self.rain_rate=list(read_file.iloc[:,26]) #Extract Raw rain rate from OTT
        self.Z=list(read_file.iloc[:,24]) # Extract Reflectivity (Z)
        self.num_particles=list(read_file.iloc[:,25]) #Extact No of particles
        nd=read_file.iloc[:,29:61] # Extract raw Nd
        #nd=nd[~nd.isnull()]  #Remove rows with NaN
        nd[(nd>=-9.999)&(nd<=0)]=0
        nd_filter=np.power(10,nd)
        self.nd=np.asarray(nd_filter)
        self.vd=np.asarray(read_file.iloc[:,61:93]) #Extract raw vd
        self.raw=np.asarray(read_file.iloc[:,93:1117]) #Extract raw 32 x 32 matrix
      
        


    def _read_file_bom(self):
        """  Read the Parsivel Data file and store it in internal structure.
        Returns: None
        Quality check of the data and apply filters
        QC1
        QC2
        QC3
        """
        
        r_file=pd.read_csv(self.filename,low_memory=False,skip_blank_lines=True,parse_dates=True)
        f_file=r_file[(r_file.ix[:,6]>20) & (r_file.ix[:,1]>=0)]#Filter for minimum no of particles 20 for 1 minute
        self.file=f_file
        f_file.index=pd.to_datetime(f_file.ix[:,0])
        self.rain_rate=list(f_file.ix[:,1]) #Extract Raw rain rate from OTT
        self.Z=list(f_file.ix[:,3]) # Extract Reflectivity (Z)
        self.num_particles=list(f_file.ix[:,6]) #Extact No of particles
        raw_nd=f_file.ix[:,9:41] # N(D) is in log10 units
        raw_nd[(raw_nd>=-9.999) &(raw_nd<=0)]=0
        nd=np.power(10,raw_nd)
        self.nd=np.asarray(nd)
        self.vd=np.asarray(f_file.ix[:,41:73]) #Extract raw vd
        self.raw=np.asarray(f_file.ix[:,73:1097]) #Extract 1024 raw vector- 32 x 32 matrix




    def _apply_pcm_matrix(self):
        """ Apply Data Quality matrix from Ali Tokay
            Based on Jaffrain et.al, 2011
            
            Returns: None

        """
       
        dt=30 # deltaT= 60 seconds        
        self.filtered_raw_matrix = np.ndarray(shape=(len(self.raw),32, 32), dtype=float) 
        self.ndf=np.zeros(shape=(len(self.raw),32),dtype=float)

        for i in range(len(self.raw)):
            self.filtered_raw_matrix[i] = np.multiply(self.pcm, np.reshape(self.raw[i], (32, 32)))
            
        sa=180*(30-0.5*self.diameter['data'])*10**-6  # Effective sampling area =180mm x (30mm -0.5Di) Jaffrain et.al, 2011
        
        for i in range(len(self.raw)):
            for dia in range(2,22):
                Nacc=0
                Ntemp=0
                for vel in range(0,32):
                    f_raw_mat=self.filtered_raw_matrix[i].T
                    Ntemp=f_raw_mat[dia,vel]/(self.velocity['data'][vel]*self.spread['data'][dia]*dt*sa[dia])
                    Nacc=Nacc+Ntemp
                if Nacc<=0:
                    self.ndf[i,dia]=1
                else:
                    self.ndf[i,dia]=Nacc


    def _prep_data(self):
        self.fields = {}

        self.fields['rain_rate'] = common.var_to_dict(
            'Rain rate', np.ma.array(self.rain_rate), 'mm/h', 'Rain rate')
        self.fields['reflectivity'] = common.var_to_dict(
            'Reflectivity', np.ma.masked_equal(self.Z,-9.999), 'dBZ',
            'Equivalent reflectivity factor')
        self.fields['Nd'] = common.var_to_dict(
            'Nd', np.ma.masked_equal(self.nd, np.power(10,-9.999)), 'm^-3 mm^-1',
            'Liquid water particle concentration')
        self.fields['Nd']['data'].set_fill_value(0)

        self.fields['num_particles'] = common.var_to_dict(
            'Number of Particles', np.ma.array(self.num_particles),
            '', 'Number of particles')
        self.fields['terminal_velocity'] = common.var_to_dict(
            'Terminal Fall Velocity', self.vd,  # np.ndarray(self.vd),
            'm/s', 'Terminal fall velocity for each bin')
        self.time = {'data': np.array(self.dates), 'units': common.EPOCH_UNITS,
                  'title': 'Time', 'long_name': 'time'}


    
    def _get_epoch_time(self):
        """
        Convert datetime to number using netCDF4.date2num
        """
        dts=[]
        for dt in range(len(self.file)):
            x=date2num(self.file.index[dt],units=common.EPOCH_UNITS)
            dts.append(x)
        self.dates=self.ndt
        
 
    
    
    
    
    

    diameter = common.var_to_dict(
        'diameter',
        np.array(
            [0.062, 0.187, 0.312, 0.437, 0.562, 0.687, 0.812, 0.937, 1.062, 1.187, 1.375, 1.625, 1.875, 2.125,
             2.375, 2.75, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17, 19, 21.5, 24.5]),
        'mm', 'Particle diameter of bins')

    
    
    spread = common.var_to_dict(
        'spread',
        [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
         0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0],
        'mm', 'Bin size spread of bins')


    velocity = common.var_to_dict(
        'velocity',
        np.array(
            [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.1, 1.3, 1.5, 1.7, 1.9,
             2.2, 2.6, 3, 3.4, 3.8, 4.4, 5.2, 6.0, 6.8, 7.6, 8.8, 10.4, 12.0, 13.6, 15.2,
             17.6, 20.8]),
        'm s^-1', 'Terminal fall velocity for each bin')

    v_spread = [.1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .2, .2, .2, .2, .2, .4,
                .4, .4, .4, .4, .8, .8, .8, .8, .8, 1.6, 1.6, 1.6, 1.6, 1.6, 3.2, 3.2]


    pcm_matrix =(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)