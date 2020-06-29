""" 
This script reads OTT disdrometer data from BoM. BoM had 3 OTT/ 3 Thies installed at Broadmeadows starting from 2014/07 till the mid of 2017. 
Data frequency : 1 min
Duration from : 2014/07 till middle of 2017
Field are extracted based on the Manual provided by BoM
Separate file per Year
All necessary field including N(d) and V(d) derived from the raw data

"""

import os,glob
import pandas as pd
import datetime


def list_flatten(l, a=None):
    if a is None:
        a = []
    for i in l:
        if isinstance(i, list):
            list_flatten(i, a)
        else:
            a.append(i)
    return a


folder=r'\\ad.monash.edu\home\User063\jaya0001\Desktop\Jaya_work\Data_analysis\disdrometer\OTT1'
os.chdir(folder)
for year in [2014,2015,2016,2017]:
    file_list=glob.glob('*'+str(year)+'*.txt') #group by each year
    all_raw_disdrometer=pd.DataFrame()
    all_data_disdrometer=pd.DataFrame()
    nd_disdrometer=pd.DataFrame()
    vd_disdrometer=pd.DataFrame()
    raw_disdrometer=pd.DataFrame()
    
    for i in range(0,len(file_list)):
        print file_list[i]
        all_date=[]
        rain_intensity=[]
        rain_acc=[]
        radar_reflectivity=[]
        MOR_vis=[]
        data_interval=[]
        no_of_particles=[]
        st_name=[]
        st_no=[]
        nd=[]
        vd=[]
        raw_data=[]
        dt_disdrometer=pd.DataFrame()
        with open(file_list[i],'rt') as f:
            for count,line in enumerate(f):
                if count>1: # To skip first 2 rows of data
                    try:
                       dt_field=line.split('-')
                       yr=dt_field[1]
                       mon=dt_field[2]
                       day=dt_field[3].split(' ')[0]
                       time=dt_field[3].split(' ')[1].split(':')
                       hr=time[0]
                       mins=time[1]
                       sec=time[2][0:2]
                       dt=datetime.datetime(int(yr),int(mon),int(day),int(hr),int(mins),int(sec))
                    except:
                       data_read=line.split(';')
                       if (len(data_read)>=1104):
                           all_date.append(dt)
                           rain_intensity.append(float(data_read[0]))   # Rain Intensity field 1
                           rain_acc.append(float(data_read[1])) # Rain amount accumulated field 2
                           radar_reflectivity.append(float(data_read[5])) # Radar reflectivity field 7
                           MOR_vis.append(float(data_read[6])) #MOR visibility in precipitation field 8
                           data_interval.append(float(data_read[7])) # Sample interval field 9
                           no_of_particles.append(float(data_read[8])) # No of detected particles field 11
                           st_name.append(data_read[12]) #Station Name field 22
                           st_no.append(data_read[13]) #Station Number field 23
                           
                           nd.append(map(float,data_read[16:48])) #Field N(d) 90
                           vd.append(map(float,data_read[48:80])) #Field v(d) 91
                           raw_data.append(map(float,data_read[80:1104])) #Raw data field 93 (32 x 32)
            dt_disdrometer['Rain Intensity']=rain_intensity
            dt_disdrometer['Rain Accumulated']=rain_acc
            dt_disdrometer['Radar Reflectivity']=radar_reflectivity
            dt_disdrometer['MOR Visibility']=MOR_vis
            dt_disdrometer['Inteval']=data_interval
            dt_disdrometer['No of Partciles']=no_of_particles
            dt_disdrometer['Station Name']=st_name
            dt_disdrometer['Station No']=st_no
            #dt_disdrometer['Nd']=nd
            #dt_disdrometer['Vd']=vd
            dt_disdrometer.index=all_date
            dt_disdrometer.index=dt_disdrometer.index.map(lambda x:x.replace (second=0))
            ndisd=pd.DataFrame(nd,index=dt_disdrometer.index)#columns={'Nd(1)','Nd(2)','Nd(3)','Nd(4)','Nd(5)','Nd(6)','Nd(7)','Nd(8)','Nd(9)','Nd(10)','Nd(11)','Nd(12)','Nd(13)','Nd(14)','Nd(15)','Nd(16)','Nd(17)','Nd(18)','Nd(19)','Nd(20)','Nd(21)','Nd(22)','Nd(23)','Nd(24)','Nd(25)','Nd(26)','Nd(27)','Nd(28)','Nd(29)','Nd(30)','Nd(31)','Nd(32)'})
            vdisd=pd.DataFrame(vd,index=dt_disdrometer.index)#columns={'Vd(1)','Vd(2)','Vd(3)','Vd(4)','Vd(5)','Vd(6)','Vd(7)','Vd(8)','Vd(9)','Vd(10)','Vd(11)','Vd(12)','Vd(13)','Vd(14)','Vd(15)','Vd(16)','Vd(17)','Vd(18)','Vd(19)','Vd(20)','Vd(21)','Vd(22)','Vd(23)','Vd(24)','Vd(25)','Vd(26)','Vd(27)','Vd(28)','Vd(29)','Vd(30)','Vd(31)','Vd(32)'})
            #raw_disdrometer['Rain Intensity']=rain_intensity
            #raw_disdrometer['Station Name']=st_name               
            #raw_disdrometer.index=all_date
            #raw_disdrometer['Raw data']=raw_data
            #raw_data=pd.DataFrame(raw_data,index=dt_disdrometer.index) #Raw data field
        
        all_data_disdrometer=all_data_disdrometer.append(dt_disdrometer)
        nd_disdrometer=nd_disdrometer.append(ndisd)
        vd_disdrometer=vd_disdrometer.append(vdisd)
        #raw_disdrometer=raw_disdrometer.append(raw_data)
    
    #all_data=pd.concat([all_data_disdrometer,nd_disdrometer,vd_disdrometer],axis=1)
    all_data=pd.concat([all_data_disdrometer,nd_disdrometer,vd_disdrometer],axis=1)
    all_data=all_data[~all_data.index.duplicated(keep='first')]
    my_data=all_data.reindex(pd.date_range(all_data.index[0],all_data.index[-1],freq='1T'),fill_value='NaN')
    my_data.to_csv('OTT1_'+str(year)+'.csv',sep=',') #will be faster if we don't use the network drive to write the file
    del my_data
               















  
    
    
    

#f_read=pd.read_table(file_list[0],header=None,index_col=None,skiprows=3,sep=';')

'''
for f in [file_list[1]]:
    print f
    f_read=pd.read_table(f,header=None,index_col=None)
    ff_read=f_read.ix[:,1]
    dt_frame=pd.DataFrame()
    all_data=pd.DataFrame()

    
    i=0    
    while i <10:
        raw_data=[]
        try:
            dt_read=ff_read.ix[i,:].split('|')[0][:-4]
            dt=datetime.datetime.strptime(dt_read,'%Y-%m-%d %H:%M:%S')
            print dt
            data_read=ff_read.ix[i,:].split('|')[1].split(';')
            rain_intensity=data_read[0]
            field_nd=data_read[16:48]
            field_vd=data_read[48:80]
            dt_frame['Nd']=field_nd
            dt_frame['Vd']=field_vd
            dt_frame['Date']=dt
            dt_frame['Intensity']=rain_intensity
            all_data=all_data.append(dt_frame)
            
            data_read1=ff_read.ix[i+1,:].split('|')[1].split(';')
            data_read2=ff_read.ix[i+2,:].split('|')[1].split(';')
            data_read3=ff_read.ix[i+3,:].split('|')[1].split(';')
            data_read4=ff_read.ix[i+4,:].split('|')[1].split(';')
            raw_data.append([data_read[80:-1],data_read1[0:],data_read2[0:-1],data_read3[0:-1],data_read4[0:-2]])
            #print len(raw_data)
            i=i+5
        except:
            i=i+1
            print i
    all_files=all_files.append(all_data)

x=list_flatten(raw_data)
a=np.array(x).reshape(32,32)
#all_files.to_csv('data.csv')   

 '''           




    








