'''----------------------------------------------- Important Libraries-----------------------------------------------------------------'''

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from scipy.interpolate import griddata

'''------------------------------------------------Creating String for the Month-------------------------------------------------'''
'''
   Input  : Month Number
   Output : String of Month number
   
   Example, 
   Input  = 1, 11, 10
   Output = '01', '11', '10'

'''

def Month_No(i):
    if i <=9:
        Month_No = '0'+str(i)
    else:
        Month_No = str(i)
    return Month_No  

'''
   Input :  Month Number of the year 2020, Week Number in a Particular Month
   Output:  Generating weekly paths according to the week number

'''

def Weekly_path_Generating(Month_No,Week_No):
    str1 = f'D:\EG\Project Data\ESA_CCI_DATA\Month_{Month_No}\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020{Month_No}'
    str2 = '000000-TCDR-v202012.0.0.nc'
    
    if Week_No == 1:
        path1  = str1 + '01' + str2
        path2  = str1 + '02' + str2
        path3  = str1 + '03' + str2
        path4  = str1 + '04' + str2
        path5  = str1 + '05' + str2
        path6  = str1 + '06' + str2
        path7  = str1 + '07' + str2
        
    elif Week_No == 2:
        path1  = str1 + '08' + str2
        path2  = str1 + '09' + str2
        path3  = str1 + '10' + str2
        path4  = str1 + '11' + str2
        path5  = str1 + '12' + str2
        path6  = str1 + '13' + str2
        path7  = str1 + '14' + str2

    elif Week_No == 3:
        path1  = str1 + '15' + str2
        path2  = str1 + '16' + str2
        path3  = str1 + '17' + str2
        path4  = str1 + '18' + str2
        path5  = str1 + '19' + str2
        path6  = str1 + '20' + str2
        path7  = str1 + '21' + str2

    elif Week_No == 4:
        path1  = str1 + '22' + str2
        path2  = str1 + '23' + str2
        path3  = str1 + '24' + str2
        path4  = str1 + '25' + str2
        path5  = str1 + '26' + str2
        path6  = str1 + '27' + str2
        path7  = str1 + '28' + str2
        
    return path1,path2,path3,path4,path5,path6,path7


'''
    Input : Grids of the two days Soil Moisture 
    Ouput : Stacked values of Soil Moisture on ESA CCI grid points

'''
def Fill_Gap(SM1,SM2):
    m,n = SM1.shape
    for i in range(m):
        for j in range(n):
            if SM1[i][j]==0:
                SM1[i][j] = SM2[i][j]
    return SM1

def SM_Daily_ESACCI(path):
    D = Dataset(path)
    Data1 = D.variables['sm'][0]
    df1 = pd.DataFrame(Data1)
    SM_Data = df1.replace(to_replace=-9999,value=0)
    SM = np.array(SM_Data)
    SM = SM[260:270,1025:1035]
    return SM

'''
    Input : Path of the daily SMAP Soil Moisture data 
    Ouput : Daily Soil Moisture 

'''
def SM_Daily_ESACCI_lat_lon(path):
    D = Dataset(path)
    lat = D.variables['lat']
    lon = D.variables['lon']
    ESACCI_Lat = np.array(lat)
    ESACCI_Lon = np.array(lon)
    lat = ESACCI_Lat[260:270]
    lon = ESACCI_Lon[1025:1035]
    return lat,lon

'''-------------------------------------------------------Reading Weekly data------------------------------------------------'''

def Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7):
    SM_1 = SM_Daily_ESACCI(path1)
    SM_2 = SM_Daily_ESACCI(path2)
    SM_3 = SM_Daily_ESACCI(path3)
    SM_4 = SM_Daily_ESACCI(path4)
    SM_5 = SM_Daily_ESACCI(path5)
    SM_6 = SM_Daily_ESACCI(path6)
    SM_7 = SM_Daily_ESACCI(path7)
    return SM_1,SM_2,SM_3,SM_4,SM_5,SM_6,SM_7

'''------------------------------------------------Visualizing Daily Soil Moisture Data ----------------------------------------------'''

def Daily_SM(Month_No,Week_No):
    path1,path2,path3,path4,path5,path6,path7 = Weekly_path_Generating(Month_No,Week_No)
    SM1,SM2,SM3,SM4,SM5,SM6,SM7 = Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7)    
    path = 'D:\EG\Project Data\ESA_CCI_DATA\Month_01\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-20200101000000-TCDR-v202012.0.0.nc'
    lat,lon = SM_Daily_ESACCI_lat_lon(path)         
    SM_act = [SM1,SM2,SM3,SM4,SM5,SM6,SM7]     
    for i in range(1,8):   
        SMa = SM_act[i-1]
        plt.figure(figsize=(10,5))
        cmap = 'gist_ncar'
        plt.contourf(lon,lat,SMa/100,cmap = cmap)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        cbar = plt.colorbar(shrink=1)
        cbar.set_label(f'Daily Soil Moisture')
        plt.title(f'''Contour plot of Soil Moisture
        Month No.{Month_No} Week No{Week_No} Day-{i}''')

'''----------------------------------------------Visualizing weekly Average Soil Moisture Data-------------------------------------'''       
        
def Weekly_avg_SM(Month_No,Week_No):
    path1,path2,path3,path4,path5,path6,path7 = Weekly_path_Generating(Month_No,Week_No)
    SM1,SM2,SM3,SM4,SM5,SM6,SM7 = Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7)    
    path = 'D:\EG\Project Data\ESA_CCI_DATA\Month_01\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-20200101000000-TCDR-v202012.0.0.nc'
    lat,lon = SM_Daily_ESACCI_lat_lon(path) 
    
    SM1 = pd.DataFrame(SM1)
    SM2 = pd.DataFrame(SM2)
    SM3 = pd.DataFrame(SM3)
    SM4 = pd.DataFrame(SM4)
    SM5 = pd.DataFrame(SM5)
    SM6 = pd.DataFrame(SM6)
    SM7 = pd.DataFrame(SM7)
    
    SM1 = SM1.replace(to_replace=np.nan,value=0)
    SM2 = SM2.replace(to_replace=np.nan,value=0)
    SM3 = SM3.replace(to_replace=np.nan,value=0)
    SM4 = SM4.replace(to_replace=np.nan,value=0)
    SM5 = SM5.replace(to_replace=np.nan,value=0)
    SM6 = SM6.replace(to_replace=np.nan,value=0)
    SM7 = SM7.replace(to_replace=np.nan,value=0)
    
    SM_avg = (SM1+SM2+SM3+SM4+SM5+SM6+SM7)/700
    plt.figure(figsize=(10,5))
    cmap = 'gist_ncar'
    plt.contourf(lon,lat,SM_avg,cmap = cmap)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    cbar = plt.colorbar(shrink=1)
    cbar.set_label(f'Weekly Avg Soil Moisture')
    plt.title(f'''Contour weekly average plot of Soil Moisture
    Month No.{Month_No} Week No{Week_No}''')
 
'''---------------------------------------Visualizing weekly stacked Soil Moisture Data---------------------------------------------''' 

def Weekly_stacking_SM(Month_No,Week_No):
    path1,path2,path3,path4,path5,path6,path7 = Weekly_path_Generating(Month_No,Week_No)
    SM1,SM2,SM3,SM4,SM5,SM6,SM7 = Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7)    
    path = 'D:\EG\Project Data\ESA_CCI_DATA\Month_01\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-20200101000000-TCDR-v202012.0.0.nc'
    lat,lon = SM_Daily_ESACCI_lat_lon(path) 
    
    SM1 = pd.DataFrame(SM1)
    SM2 = pd.DataFrame(SM2)
    SM3 = pd.DataFrame(SM3)
    SM4 = pd.DataFrame(SM4)
    SM5 = pd.DataFrame(SM5)
    SM6 = pd.DataFrame(SM6)
    SM7 = pd.DataFrame(SM7)
    
    SM1 = SM1.replace(to_replace=np.nan,value=0)
    SM2 = SM2.replace(to_replace=np.nan,value=0)
    SM3 = SM3.replace(to_replace=np.nan,value=0)
    SM4 = SM4.replace(to_replace=np.nan,value=0)
    SM5 = SM5.replace(to_replace=np.nan,value=0)
    SM6 = SM6.replace(to_replace=np.nan,value=0)
    SM7 = SM7.replace(to_replace=np.nan,value=0)
    
    SM2 = Fill_Gap(SM1,SM2)
    SM3 = Fill_Gap(SM2,SM3)
    SM4 = Fill_Gap(SM3,SM4)
    SM5 = Fill_Gap(SM4,SM5)
    SM6 = Fill_Gap(SM5,SM6)
    SM7 = Fill_Gap(SM6,SM7)
    
    plt.figure(figsize=(10,5))
    cmap = 'gist_ncar'
    plt.contourf(lon,lat,SM7,cmap = cmap) 
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    cbar = plt.colorbar(shrink=1)
    cbar.set_label(f'Weekly stacked Soil Moisture')
    plt.title(f'''Contour weekly stacked plot of Soil Moisture
    Month No.{Month_No} Week No{Week_No}''')