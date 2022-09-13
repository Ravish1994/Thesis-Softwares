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
    
    str1 = f'D:\EG\Project Data\SMAP_DATA\Month_{Month_No}\SMAP_L3_SM_P_2020{Month_No}'
    str2 = '_R17000_001.h5'
    
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
    Ouput : Stacked values of Soil Moisture on SMAP grid points

'''

def Fill_Gap(S1,S2):
    m,n = S1.shape
    for i in range(m):
        for j in range(n):
            if S1[i][j]==0:
                S1[i][j] = S2[i][j]
    return S1

'''
    Input : Path of the daily SMAP Soil Moisture data 
    Ouput : Stacked values of AM and PM Soil Moisture 

'''
def SM_Daily_SMAP(path):
    Data = h5py.File(path,'r')
    d1 = Data['Soil_Moisture_Retrieval_Data_AM']    
    d2 = Data['Soil_Moisture_Retrieval_Data_PM']
    
    D1 = d1['soil_moisture'][:] # Unit cm^3/cm^3
    df1 = pd.DataFrame(D1)
    SM1 = np.array(df1.replace(to_replace=-9999,value=0))
    
    lat1 = d1['latitude'][:]
    dff1 = pd.DataFrame(lat1)
    lat1 = np.array(dff1.replace(to_replace=-9999,value=0))
    
    lon1 = d1['longitude'][:]
    dff2 = pd.DataFrame(lon1)
    lon1 = np.array(dff2.replace(to_replace=-9999,value=0))
    
    D2 = d2['soil_moisture_pm'][:] # Unit cm^3/cm^3
    df2 = pd.DataFrame(D2)
    SM2 = np.array(df2.replace(to_replace=-9999,value=0))
    
    lat2 = d2['latitude_pm'][:]
    dff2 = pd.DataFrame(lat2)
    lat2 = np.array(dff2.replace(to_replace=-9999,value=0))
    
    lon2 = d2['longitude_pm'][:]
    dff2 = pd.DataFrame(lon2)
    lon2 = np.array(dff2.replace(to_replace=-9999,value=0))
    
    SM = Fill_Gap(SM1,SM2)
    SM = SM[97:131,676:723]
    lat = Fill_Gap(lat1,lat2)
    lat = lat[97:131,676:723]
    lon = Fill_Gap(lon1,lon2)
    lon = lon[97:131,676:723]
    return lat,lon,SM


'''-------------------------------------------------- Masking the Data in our Catchment----------------------------------------'''

def Remove_Null(lat,lon,SM):
    DF = pd.DataFrame(lat)
    DF.columns = ['lat']
    DF['lon'] = lon
    DF['SM'] = SM
    
    mask = ((DF['lat']>=21) & (DF['lat']<=31) & (DF['lon']>=73) & (DF['lon']<=89) & (DF['SM']>=0))
    DF = DF[mask]
    
    DF  = DF.dropna()
    lat = DF['lat']
    lon = DF['lon']
    SM  = DF['SM']
    return lat,lon,SM



'''-------------------------------------------------------Reading Weekly data------------------------------------------------'''


def Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7):
    lat1,lon1,SM_1 = SM_Daily_SMAP(path1)
    lat2,lon2,SM_2 = SM_Daily_SMAP(path2)
    lat3,lon3,SM_3 = SM_Daily_SMAP(path3)
    lat4,lon4,SM_4 = SM_Daily_SMAP(path4)
    lat5,lon5,SM_5 = SM_Daily_SMAP(path5)
    lat6,lon6,SM_6 = SM_Daily_SMAP(path6)
    lat7,lon7,SM_7 = SM_Daily_SMAP(path7)
    return lat1,lat2,lat3,lat4,lat5,lat6,lat7,lon1,lon2,lon3,lon4,lon5,lon6,lon7,SM_1,SM_2,SM_3,SM_4,SM_5,SM_6,SM_7




'''----------------------------------Visualizing Daily Soil Moisture Data after stacking AM and PM----------------------------------'''



def Daily_SM(Month_No,Week_No):
    path1,path2,path3,path4,path5,path6,path7 = Weekly_path_Generating(Month_No,Week_No)
    lat1,lat2,lat3,lat4,lat5,lat6,lat7,lon1,lon2,lon3,lon4,lon5,lon6,lon7,SM_1,SM_2,SM_3,SM_4,SM_5,SM_6,SM_7 = Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7)
    lat1,lon1,SM1= lat1.flatten(),lon1.flatten(),SM_1.flatten()
    lat2,lon2,SM2= lat2.flatten(),lon2.flatten(),SM_2.flatten()
    lat3,lon3,SM3= lat3.flatten(),lon3.flatten(),SM_3.flatten()
    lat4,lon4,SM4= lat4.flatten(),lon4.flatten(),SM_4.flatten()
    lat5,lon5,SM5= lat5.flatten(),lon5.flatten(),SM_5.flatten()
    lat6,lon6,SM6= lat6.flatten(),lon6.flatten(),SM_6.flatten()
    lat7,lon7,SM7= lat7.flatten(),lon7.flatten(),SM_7.flatten()
    
    lat1,lon1,SM1 = Remove_Null(lat1,lon1,SM1)
    lat2,lon2,SM2 = Remove_Null(lat2,lon2,SM2)
    lat3,lon3,SM3 = Remove_Null(lat3,lon3,SM3)
    lat4,lon4,SM4 = Remove_Null(lat4,lon4,SM4)
    lat5,lon5,SM5 = Remove_Null(lat5,lon5,SM5)
    lat6,lon6,SM6 = Remove_Null(lat6,lon6,SM6)
    lat7,lon7,SM7 = Remove_Null(lat7,lon7,SM7)
    
    Lat1,Lon1 = np.meshgrid(lat1,lon1)
    Lat2,Lon2 = np.meshgrid(lat2,lon2)
    Lat3,Lon3 = np.meshgrid(lat3,lon3)
    Lat4,Lon4 = np.meshgrid(lat4,lon4)
    Lat5,Lon5 = np.meshgrid(lat5,lon5)
    Lat6,Lon6 = np.meshgrid(lat6,lon6)
    Lat7,Lon7 = np.meshgrid(lat7,lon7)
    
    SM1 = griddata((lat1,lon1),SM1,(Lat1,Lon1),method='cubic')
    SM2 = griddata((lat2,lon2),SM2,(Lat2,Lon2),method='cubic')
    SM3 = griddata((lat3,lon3),SM3,(Lat3,Lon3),method='cubic')
    SM4 = griddata((lat4,lon4),SM4,(Lat4,Lon4),method='cubic')
    SM5 = griddata((lat5,lon5),SM5,(Lat5,Lon5),method='cubic')
    SM6 = griddata((lat6,lon6),SM6,(Lat6,Lon6),method='cubic')
    SM7 = griddata((lat7,lon7),SM7,(Lat7,Lon7),method='cubic') 
    
    SM_t = [SM1,SM2,SM3,SM4,SM5,SM6,SM7] 
    Lat_t = [Lat1,Lat2,Lat3,Lat4,Lat5,Lat6,Lat7]
    Lon_t = [Lon1,Lon2,Lon3,Lon4,Lon5,Lon6,Lon7]
    SM_act = [SM_1,SM_2,SM_3,SM_4,SM_5,SM_6,SM_7]
    
    for i in range(1,8):
        SM    = SM_t[i-1]
        m1,n1 = SM.shape
        
        ## Masking the array without disturbing its shape and size
        for i1 in range(m1):
            for j1 in range(n1):
                a = SM[i1][j1]
                if (a>0) and (a<1):
                    SM[i1][j1] = SM[i1][j1]
                else:
                    SM[i1][j1] = np.nan
                    
        Lat = Lat_t[i-1]
        Lon = Lon_t[i-1]    
        SMa = SM_act[i-1]
        m,n = SMa.shape
        
        ## Masking the array without disturbing its shape and size
        for i2 in range(m):
            for j2 in range(n):
                a1 = SMa[i2][j2]
                if (a1>0) and (a1<1):
                    SMa[i2][j2] = SMa[i2][j2]
                else:
                    SMa[i2][j2] = np.nan

                
        plt.figure(figsize=(5,5))
        cmap = 'gist_ncar'
        plt.figure(figsize=(10,5))
        plt.contourf(Lon,Lat,SM,cmap = cmap)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        cbar = plt.colorbar(shrink=1)
        cbar.set_label(f'Daily Soil Moisture')
        plt.title(f'''Contour plot of Soil Moisture Month No.{Month_No} Week No{Week_No} Day-{i}''')
        
        
        
'''----------------------------------------------Visualizing weekly Average Soil Moisture Data-------------------------------------''' 


def Weekly_Avg_SM(Month_No,Week_No):
    path1,path2,path3,path4,path5,path6,path7 = Weekly_path_Generating(Month_No,Week_No)
    lat1,lat2,lat3,lat4,lat5,lat6,lat7,lon1,lon2,lon3,lon4,lon5,lon6,lon7,SM_1,SM_2,SM_3,SM_4,SM_5,SM_6,SM_7 = Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7)
    
    lon2 = Fill_Gap(lon1,lon2)
    lon3 = Fill_Gap(lon2,lon3)
    lon4 = Fill_Gap(lon3,lon4)
    lon5 = Fill_Gap(lon4,lon5)
    lon6 = Fill_Gap(lon5,lon6)
    lon7 = Fill_Gap(lon6,lon7) 

    lat2 = Fill_Gap(lat1,lat2)
    lat3 = Fill_Gap(lat2,lat3)
    lat4 = Fill_Gap(lat3,lat4)
    lat5 = Fill_Gap(lat4,lat5)
    lat6 = Fill_Gap(lat5,lat6)
    lat7 = Fill_Gap(lat6,lat7)
       
    SM_avg = (SM_1+SM_2+SM_3+SM_4+SM_5+SM_6+SM_7)/7 
    SM_avg1 = pd.DataFrame(SM_avg)
    SM_Avgf = np.array(SM_avg1.replace(to_replace=0,value=np.nan))
    
    lat7,lon7,SM_avgf = lat7.flatten(),lon7.flatten(),SM_Avgf.flatten()
    lat7,lon7,SM_avgf = Remove_Null(lat7,lon7,SM_avgf)
    Lat7,Lon7 = np.meshgrid(lat7,lon7)
    
    SM_avgf = griddata((lat7,lon7),SM_avgf,(Lat7,Lon7),method='nearest')
    SM_avgf = np.array(pd.DataFrame(SM_avgf).replace(to_replace=0,value=np.nan))
    
    m,n = SM_avgf.shape
    ## Masking the array without disturbing its shape and size
    for i2 in range(m):
        for j2 in range(n):
            a1 = SM_avgf[i2][j2]
            if (a1>0) and (a1<1):
                SM_avgf[i2][j2] = SM_avgf[i2][j2]
            else:
                SM_avgf[i2][j2] = np.nan 
    
    cmap = 'gist_ncar'  
    plt.figure(figsize=(10,5))
    plt.contourf(Lon7,Lat7,SM_avgf,cmap = cmap)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    cbar = plt.colorbar(shrink=1)
    cbar.set_label(f'Weekly Avg Soil Moisture')
    plt.title(f'''Average Soil Moisture Month No.{Month_No} Week No{Week_No}''')
    

'''---------------------------------------Visualizing weekly stacked Soil Moisture Data---------------------------------------------''' 

def Weekly_Stack_SM(Month_No,Week_No):
    path1,path2,path3,path4,path5,path6,path7 = Weekly_path_Generating(Month_No,Week_No)
    lat1,lat2,lat3,lat4,lat5,lat6,lat7,lon1,lon2,lon3,lon4,lon5,lon6,lon7,SM_1,SM_2,SM_3,SM_4,SM_5,SM_6,SM_7 = Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7)
    
    lon2 = Fill_Gap(lon1,lon2)
    lon3 = Fill_Gap(lon2,lon3)
    lon4 = Fill_Gap(lon3,lon4)
    lon5 = Fill_Gap(lon4,lon5)
    lon6 = Fill_Gap(lon5,lon6)
    lon7 = Fill_Gap(lon6,lon7) 
    
    lat2 = Fill_Gap(lat1,lat2)
    lat3 = Fill_Gap(lat2,lat3)
    lat4 = Fill_Gap(lat3,lat4)
    lat5 = Fill_Gap(lat4,lat5)
    lat6 = Fill_Gap(lat5,lat6)
    lat7 = Fill_Gap(lat6,lat7)
    
    SM_2 = Fill_Gap(SM_1,SM_2)
    SM_3 = Fill_Gap(SM_2,SM_3)
    SM_4 = Fill_Gap(SM_3,SM_4)
    SM_5 = Fill_Gap(SM_4,SM_5)
    SM_6 = Fill_Gap(SM_5,SM_6)
    SM_7 = Fill_Gap(SM_6,SM_7)
     
    SM_stk1 = pd.DataFrame(SM_7)
    SM_Stkf = np.array(SM_stk1.replace(to_replace=0,value=np.nan))
    
    lat7,lon7,SM_stkf = lat7.flatten(),lon7.flatten(),SM_Stkf.flatten()
    lat7,lon7,SM_stkf = Remove_Null(lat7,lon7,SM_stkf)
    Lat7,Lon7 = np.meshgrid(lat7,lon7)
    
    SM_stkf = griddata((lat7,lon7),SM_stkf,(Lat7,Lon7),method='nearest')
    SM_stkf = np.array(pd.DataFrame(SM_stkf).replace(to_replace=0,value=np.nan))
    
    m,n = SM_stkf.shape
    ## Masking the array without disturbing its shape and size
    for i2 in range(m):
        for j2 in range(n):
            a1 = SM_stkf[i2][j2]
            if (a1>0) and (a1<1):
                SM_stkf[i2][j2] = SM_stkf[i2][j2]
            else:
                SM_stkf[i2][j2] = np.nan 
  
    cmap = 'gist_ncar'
    plt.figure(figsize=(10,5))
    plt.contourf(Lon7,Lat7,SM_stkf,cmap = cmap)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    cbar = plt.colorbar(shrink=1)
    cbar.set_label(f'Weekly stacked Soil Moisture')
    plt.title(f'''Weekly stacked Soil Moisture of Month No.{Month_No} Week No{Week_No}''')
    