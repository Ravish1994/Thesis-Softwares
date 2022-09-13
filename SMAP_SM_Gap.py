import h5py
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

'''----------------------------------------Visualizing SMAP Soil Moisture Data of AM and PM-------------------------------------------'''

'''
   Input : path of all h5py files 
           a: day number in the year 2020 
           
   Output : Visualizing plots of Soil Moisture Data of AM and PM

''' 

def SMAP_SM_Gap(path,a):
    Data = h5py.File(path,'r')
    dat1 = Data['Soil_Moisture_Retrieval_Data_AM']    
    dat2 = Data['Soil_Moisture_Retrieval_Data_PM']
    
    # AM Soil Moisture data
    lat1 = dat1['latitude'][:]
    dff1 = pd.DataFrame(lat1)
    lat1 = np.array(dff1.replace(to_replace=-9999,value=np.nan))
    lat1 = lat1[97:131,676:723]
    
    DF1 = pd.DataFrame(lat1.flatten())
    DF1.columns = ['lat']
    
    lon1 = dat1['longitude'][:]
    dff2 = pd.DataFrame(lon1)
    lon1 = np.array(dff2.replace(to_replace=-9999,value=np.nan))
    lon1 = lon1[97:131,676:723]
    lon1 = lon1.flatten()
    DF1['lon'] = lon1

    DATA1 = dat1['soil_moisture'][:] # Unit cm^3/cm^3
    df3 = pd.DataFrame(DATA1)
    Data_SM1 = np.array(df3.replace(to_replace=-9999,value=np.nan))
    Data_SM_1 = Data_SM1[97:131,676:723]
    Data_SM1 = Data_SM_1.flatten()
    DF1['SM_AM']   = Data_SM1
    DF1 = DF1.dropna()
    
    lat1 = DF1['lat']
    lon1 = DF1['lon']
    Data_SM1 = DF1['SM_AM']
    lat_1,lon_1 = np.meshgrid(lat1,lon1)    
    
    # PM Soil Moisture data
    lat2 = dat2['latitude_pm'][:]
    dff2 = pd.DataFrame(lat2)
    lat2 = np.array(dff2.replace(to_replace=-9999,value=np.nan))
    lat2 = lat2[97:131,676:723]
    
    DF2 = pd.DataFrame(lat2.flatten())
    DF2.columns = ['lat']
    
    lon2 = dat2['longitude_pm'][:]
    dff2 = pd.DataFrame(lon2)
    lon2 = np.array(dff2.replace(to_replace=-9999,value=np.nan))
    lon2 = lon2[97:131,676:723]
    lon2 = lon2.flatten()
    DF2['lon'] = lon2

    DATA2 = dat2['soil_moisture_pm'][:] # Unit cm^3/cm^3
    df3 = pd.DataFrame(DATA2)
    Data_SM2 = np.array(df3.replace(to_replace=-9999,value=np.nan))
    Data_SM_2 = Data_SM2[97:131,676:723]
    Data_SM2 = Data_SM_2.flatten()
    DF2['SM_PM']   = Data_SM2
    DF2 = DF2.dropna()
    
    lat2 = DF2['lat']
    lon2 = DF2['lon']
    Data_SM2 = DF2['SM_PM']
    lat_2,lon_2 = np.meshgrid(lat2,lon2)    
    
    if ((len(Data_SM2)>0) and (len(Data_SM1)>0)):
        Data_SM1 = griddata((lat1,lon1),Data_SM1,(lat_1,lon_1), method = 'nearest')
        Data_SM2 = griddata((lat2,lon2),Data_SM2,(lat_2,lon_2), method = 'nearest')
    
        plt.figure(figsize=(20,5))
        cmap = 'gist_ncar'
        plt.subplot(1,2,1)
        plt.contourf(lon_1,lat_1,Data_SM1,cmap=cmap)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        cbar = plt.colorbar(shrink=1)
        plt.title(f'{a}: Contour Plot of Soil Moisture AM')
        
        plt.subplot(1,2,2)
        plt.contourf(lon_2,lat_2,Data_SM2,cmap=cmap)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.title(f'{a}: Contour Plot of Soil Moisture PM')
        cbar = plt.colorbar(shrink=1)
        
        
'''
   Input : m: starting day of the year
           n: ending day of the year
   
   Output : Visualizing plots of Soil Moisture Data of AM and PM

''' 

def SMAP_SM_Gaps_2020(m,n):
    for i in range(m,n+1):
        if i<=31:
            d = i
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_01\SMAP_L3_SM_P_2020010{d}_R17000_001.h5'
                a = f'0{d} Jan 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_01\SMAP_L3_SM_P_202001{d}_R17000_001.h5'
                a = f'{d} Jan 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
        elif i<=60:
            d = i-31
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_02\SMAP_L3_SM_P_2020020{d}_R17000_001.h5'
                a = f'0{d} Feb 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_02\SMAP_L3_SM_P_202002{d}_R17000_001.h5'
                a = f'{d} Feb 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
        elif i<=91:
            d = i-60
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_03\SMAP_L3_SM_P_2020030{d}_R17000_001.h5'
                a = f'0{d} Mar 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_03\SMAP_L3_SM_P_202003{d}_R17000_001.h5'
                a = f'{d} Mar 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
        elif i<=121:
            d = i-91
            if i <= 99:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_2020040{d}_R17000_001.h5'
                a = f'0{d} April 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:                
                if d<=9:
                    path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_2020040{d}_R17000_001.h5'
                    a = f'0{d} April 2020 Soil Moisture'
                    SMAP_SM_Gap(path2,a)
                else:
                    path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_202004{d}_R17000_001.h5'
                    a = f'{d} April 2020 Soil Moisture'
                    SMAP_SM_Gap(path2,a)
        elif i<=152:
            d = i-121
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_05\SMAP_L3_SM_P_2020050{d}_R17000_001.h5'
                a = f'0{d} May 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_05\SMAP_L3_SM_P_202005{d}_R17000_001.h5'
                a = f'{d} May 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
        elif i<=182:            
            d = i-152
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_06\SMAP_L3_SM_P_2020060{d}_R17000_001.h5'
                a = f'0{d} June 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_06\SMAP_L3_SM_P_202006{d}_R17000_001.h5'
                a = f'{d} June 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
        elif i<=213:            
            d = i-182
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_07\SMAP_L3_SM_P_2020070{d}_R17000_001.h5'
                a = f'0{d} July 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_07\SMAP_L3_SM_P_202007{d}_R17000_001.h5'
                a = f'{d} July 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
        elif i<=244:
            d = i-213
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_08\SMAP_L3_SM_P_2020080{d}_R17000_001.h5'
                a = f'0{d} August 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_08\SMAP_L3_SM_P_202008{d}_R17000_001.h5'
                a = f'{d} August 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
        elif i<=274:
            d = i-244
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_09\SMAP_L3_SM_P_2020090{d}_R17000_001.h5'
                a = f'0{d} Sept 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_09\SMAP_L3_SM_P_202009{d}_R17000_001.h5'
                a = f'{d} Sept 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
        elif i<=305:
            d = i-274
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_10\SMAP_L3_SM_P_2020100{d}_R17000_001.h5'
                a = f'0{d} Oct 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_10\SMAP_L3_SM_P_202010{d}_R17000_001.h5'
                a = f'{d} Oct 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
        elif i<=335:
            d = i-305
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_11\SMAP_L3_SM_P_2020110{d}_R17000_001.h5'
                a = f'0{d} Nov 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_11\SMAP_L3_SM_P_202011{d}_R17000_001.h5'
                a = f'{d} Nov 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
        elif i<=366:
            d = i-335
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_12\SMAP_L3_SM_P_2020120{d}_R17000_001.h5'
                a = f'0{d} Dec 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_12\SMAP_L3_SM_P_202012{d}_R17000_001.h5'
                a = f'{d} Dec 2020 Soil Moisture'
                SMAP_SM_Gap(path2,a)              