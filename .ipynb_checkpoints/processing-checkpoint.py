import os
from datetime import datetime

import netCDF4
import cftime
import numpy as np

import common

def brisbane_disdro_to_radar_moments(file_list, scatterer):
    
    DBZ_list = []
    ZDR_list = []
    KDP_list = []
    ATTEN_list = []
    TIME_list = []
    RAIN_list = []

    #parse text files
    data_line_len = 50
    dsd_dictlist = []
    
    for infile in file_list:
        print('processing', infile)
        with open(infile) as f:
            #init
            bin_name_list = []
            bin_value_list = []
            bin_mean_drop = []
            #parse date from filename
            fn = os.path.basename(infile)
            year = '20' + fn[1:3]
            doy = fn[3:6]
            #parse each line
            for line in f:
                #skip empty lines
                if len(line) == 1:
                    continue
                elif len(line) < data_line_len:
                    #parse rainfall
                    rainfall = float(line)
                    dsd_dict = {'dt_start':line_dt_start, 'dt_end':line_dt_end, 'rainfall':rainfall,
                                'bin_mean_drop':bin_mean_drop, 'bin_name_list':bin_name_list, 'bin_value_list':bin_value_list}
                    dsd_dictlist.append(dsd_dict)
                    #init
                    bin_name_list = []
                    bin_value_list = []
                    bin_mean_drop = []
                else:
                    #parse line
                    line_list = line.split()
                    #parse date
                    line_dt_start = datetime.strptime(year + '-' + doy + ' ' + ':'.join(line_list[0:3]),'%Y-%j %H:%M:%S')
                    line_dt_end = datetime.strptime(year + '-' + doy + ' ' + ':'.join(line_list[3:6]),'%Y-%j %H:%M:%S')
                    #parse bin
                    bin_mean_drop.append((float(line_list[6])+float(line_list[7]))/2)
                    bin_name_list.append('-'.join(line_list[6:8]) + 'mm')
                    bin_value_list.append(float(line_list[8]))
                
    #scattering calcs
    print('running scattering calc')
    for item in dsd_dictlist:
        nd = np.array(item['bin_value_list'])
        mean_diam_drop_class = np.array(item['bin_mean_drop'])  

        if np.sum(nd) == 0:
            continue

        #calc radar moments
        dbz, zdr, kdp, atten_spec = common.scatter_off_2dvd_packed(mean_diam_drop_class, nd, scatterer)
        DBZ_list.append(dbz)
        ZDR_list.append(zdr)
        KDP_list.append(kdp)
        ATTEN_list.append(atten_spec)
        TIME_list.append(item['dt_start'])
        RAIN_list.append(item['rainfall'])
        
    dbz_array = np.array(DBZ_list)
    zdr_array = np.array(ZDR_list)
    kdp_array = np.array(KDP_list)
    att_array = np.array(ATTEN_list)
    rain_array = np.array(RAIN_list)
    time_array = np.array(TIME_list)
    
    return dbz_array, zdr_array, kdp_array, att_array, time_array, rain_array


def darwin_disdro_to_radar_moments(file_list, scatterer):
    
    #init lists
    DBZ_list = []
    ZDR_list = []
    KDP_list = []
    ATTEN_list = []
    RAIN_list = []
    TIME_list = []
    
    #read DSD data
    for infile in file_list:
        print('processing', infile)
        with netCDF4.Dataset(infile, 'r') as ncid:
            time = cftime.num2pydate(ncid['time'][:], ncid['time'].units)
            mean_diam_drop_class = ncid['mean_diam_drop_class'][:]
            num_drop = ncid['num_drop'][:]        
            ndensity = ncid['nd'][:]
            liq_water = ncid['liq_water'][:]
            Z = ncid['Z'][:]    
            nclambda = ncid['lambda'][:]
            n_0 = ncid['n_0'][:]
            rain = ncid['rain_rate'][:]
        
    #fir each sample, use number density
    print('running scattering calcs')
    cnt = 0
    for nd in ndensity:
        if np.sum(nd) == 0:
            continue
        cnt += 1
        
        #calc radar moments
        dbz, zdr, kdp, atten_spec = common.scatter_off_2dvd_packed(mean_diam_drop_class, nd, scatterer)    
        DBZ_list.append(dbz)
        ZDR_list.append(zdr)
        KDP_list.append(kdp)
        ATTEN_list.append(atten_spec)
        TIME_list.append(time[cnt].date())
        RAIN_list.append(rain[cnt])
        
    dbz_array = np.array(DBZ_list)
    zdr_array = np.array(ZDR_list)
    kdp_array = np.array(KDP_list)
    att_array = np.array(ATTEN_list)
    rain_array = np.array(RAIN_list)
    time_array = np.array(TIME_list)
    
    return dbz_array, zdr_array, kdp_array, att_array, time_array, rain_array


def broadmeadows_disdro_to_radar_moments(file_list, scatterer):
    
    #init
    empty_value = -9.999
    mean_diam_drop_class = np.array([0.062, 0.187, 0.312, 0.437, 0.562, 0.687, 0.812, 
                                      0.937, 1.062, 1.187, 1.375, 1.625, 1.875, 2.125,
                                      2.375, 2.75, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 
                                      7.5, 8.5, 9.5, 11, 13, 15, 17, 19, 21.5, 24.5])
    
    #init lists
    DBZ_list = []
    ZDR_list = []
    KDP_list = []
    ATTEN_list = []
    RAIN_list = []
    TIME_list = []
    
    #parse all text files
    data_line_len = 50
    dsd_dictlist = []
    
    for infile in file_list:
        print('processing', infile)
        #init lists
        dt_list = []
        nd_list = []
        rr_list = []
        #find total number of times
        num_lines = sum(1 for line in open(infile,'r'))
        #read file
        with open(infile) as f:
            #for every line
            for count, line in enumerate(f):
                if count>0: # skip first line
                    #split line
                    data_read = line.split(',')
                    #only process data string lines
                    if len(data_read) >= 1000:
                        #append to lists
                        dt_list.append(datetime.strptime(data_read[0],'%Y-%m-%d %H:%M:%S'))
                        nd_list.append(list(map(float,data_read[11:43]))) #Field N(d) 90
                        rr_list.append(float(data_read[1]))
                    else:
                        #catch short lists
                        print('error', data_read)
                        break
    
    print('running scattering calcs')
    for i, nd in enumerate(nd_list):
        #convert list to array
        nd_array = np.array(nd)
        #replace empty value with zeros
        nd_array[nd_array==empty_value] = 0
        #skip when no data is present
        if np.sum(nd_array) == 0:
            continue

        #calc radar moments
        normal_nd_array = nd_array/np.gradient(mean_diam_drop_class)
        dbz, zdr, kdp, atten_spec = common.scatter_off_2dvd_packed(mean_diam_drop_class, normal_nd_array, scatterer)

        #append outputs
        DBZ_list.append(dbz)
        ZDR_list.append(zdr)
        KDP_list.append(kdp)
        ATTEN_list.append(atten_spec)
        RAIN_list.append(rr_list[i])
        TIME_list.append(dt_list[i])

    dbz_array = np.array(DBZ_list)
    zdr_array = np.array(ZDR_list)
    kdp_array = np.array(KDP_list)
    att_array = np.array(ATTEN_list)
    rain_array = np.array(RAIN_list)
    time_array = np.array(TIME_list)
    
    return dbz_array, zdr_array, kdp_array, att_array, time_array, rain_array

def mtview_disdro_to_radar_moments(file_list, scatterer):

    empty_value = -9.999
    mean_diam_drop_class = np.array([0.062, 0.187, 0.312, 0.437, 0.562, 0.687, 0.812, 
                                  0.937, 1.062, 1.187, 1.375, 1.625, 1.875, 2.125,
                                  2.375, 2.75, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 
                                  7.5, 8.5, 9.5, 11, 13, 15, 17, 19, 21.5, 24.5])
    
    
    #init lists
    DBZ_list = []
    ZDR_list = []
    KDP_list = []
    ATTEN_list = []
    RAIN_list = []
    TIME_list = []
    
    #parse all text files
    data_line_len = 50
    dsd_dictlist = []

    for infile in file_list:
        print('processing', infile)
        #init lists
        dt_list = []
        nd_list = []
        rr_list = []
        #find total number of times
        num_lines = sum(1 for line in open(infile,'r'))
        #read file
        with open(infile) as f:
            #for every line
            for count, line in enumerate(f):
                if count>0: # skip first line
                    #split line
                    data_read = line.split(',')
                    #only process data string lines
                    if len(data_read) >= 1000:
                        #append to lists
                        try:
                            dt_list.append(datetime.strptime(data_read[0],'%Y-%m-%d %H:%M:%S'))
                        except:
                            dt_list.append(datetime.strptime(data_read[0],'%d/%m/%Y %H:%M'))
                        rr_list.append(float(data_read[26]))
                        nd_list.append(list(map(float,data_read[29:61]))) #Field N(d) 90
                    else:
                        #catch short lists
                        print('error', data_read)
                        break
    
    print('running scatteirng calcs')
    for i, nd in enumerate(nd_list):
        #convert list to array
        nd_array = np.array(nd)
        #replace empty value with zeros
        nd_array[nd_array==empty_value] = 0
        #skip when no data is present
        if np.sum(nd_array) == 0:
            continue

        #calc radar moments
        normal_nd_array = nd_array/np.gradient(mean_diam_drop_class)
        dbz, zdr, kdp, atten_spec = common.scatter_off_2dvd_packed(mean_diam_drop_class, normal_nd_array, scatterer)

        #append outputs
        DBZ_list.append(dbz)
        ZDR_list.append(zdr)
        KDP_list.append(kdp)
        ATTEN_list.append(atten_spec)
        RAIN_list.append(rr_list[i])
        TIME_list.append(dt_list[i])

    dbz_array = np.array(DBZ_list)
    zdr_array = np.array(ZDR_list)
    kdp_array = np.array(KDP_list)
    att_array = np.array(ATTEN_list)
    rain_array = np.array(RAIN_list)
    time_array = np.array(TIME_list)
    
    return dbz_array, zdr_array, kdp_array, att_array, time_array, rain_array