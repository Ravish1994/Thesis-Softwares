from netCDF4 import Dataset
import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import h5py

'''----------------------------------Taking out single SMAP Soil Moisture Data within a particular grid cell-------------------------'''

'''
   Input : path: Path of the SMAP Daily Raw files
           lat and lon: Centroid of latitude and longitude of the SMAP single grid cell
   
   Output : Single SMAP soil moisture within that SMAP grid cell

'''

def SMAP_Data(path,lat,lon):
    
    # Importing the SMAP hdf files 
    SMAP_data = h5py.File(path,'r')  
    
    # PM data
    df1     = SMAP_data['Soil_Moisture_Retrieval_Data_PM']
    SM_lat1 = np.array(pd.DataFrame(df1['latitude_pm']))
    SM_lon1 = np.array(pd.DataFrame(df1['longitude_pm']))
    SM1     = np.array(pd.DataFrame(df1['soil_moisture_pm']))
    S1      = pd.DataFrame(SM_lat1.reshape(406*964,1))
    S2      = pd.DataFrame(SM_lon1.reshape(406*964,1))
    SM_1    = pd.DataFrame(SM1.reshape(406*964,1))
    
    # AM data
    df2     = SMAP_data['Soil_Moisture_Retrieval_Data_AM']
    SM_lat2 = np.array(pd.DataFrame(df2['latitude']))
    SM_lon2 = np.array(pd.DataFrame(df2['longitude']))        
    SM2     = np.array(pd.DataFrame(df2['soil_moisture']))
    T1      = pd.DataFrame(SM_lat2.reshape(406*964,1))
    T2      = pd.DataFrame(SM_lon2.reshape(406*964,1))  
    SM_2    = pd.DataFrame(SM2.reshape(406*964,1))

    # Concatenating AM and PM data  
    S_1 = pd.concat([S1,T1])
    S_2 = pd.concat([S2,T2])
    SM  = pd.concat([SM_1,SM_2]) 
    
    # Creating dataframe of SMAP_Lat, SMAP_Lon and SMAP_SM
    S_1.columns = ['Lat']
    S_1['lon']  = S_2
    S_1['SM']   = SM  
    
    # Masking and taking out SMAP soil moisture of that particular cell
    mask = ((S_1['Lat']>(lat-0.18)) & (S_1['Lat']<(lat+0.18)) & (S_1['lon']>(lon-0.18)) & (S_1['lon']<(lon+0.18)) & (S_1['SM']>0))
    Df = S_1[mask]
    SM2 = np.array(Df['SM'])
    if len(SM2)>0:
        SM1 = np.mean(SM2)
    else:
        SM1 = np.nan
    return SM1*100


'''-----------------------------Taking out average SMAP Soil Moisture Data within its 0.5 x 0.5 degree grid cell-------------------------'''

'''
   Input : path: Path of the SMAP Daily Raw files
           lat and lon: Centroid of latitude and longitude of the 0.5 x 0.5 degree grid cell
   
   Output : Saving CSV file for a single pixel averaged SMAP soil moisture values for 366 days of year 2020 within that grid cell

'''
def SMAP_SM_availability(m,n,lat,lon):
    SMAP_Data1 = []
    Day_No = []
    for i in range(m,n+1):
        Day_No1 = i
        Day_No.append(Day_No1)
        if i<=31:
            d = i
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_01\SMAP_L3_SM_P_2020010{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_01\SMAP_L3_SM_P_202001{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
        elif i<=60:
            d = i-31
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_02\SMAP_L3_SM_P_2020020{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_02\SMAP_L3_SM_P_202002{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
        elif i<=91:
            d = i-60
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_03\SMAP_L3_SM_P_2020030{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_03\SMAP_L3_SM_P_202003{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
        elif i<=121:
            d = i-91
            if i <= 99:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_2020040{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                if d<=9:
                    path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_2020040{d}_R17000_001.h5'
                    SM = SMAP_Data(path2,lat,lon)
                    SMAP_Data1.append(SM)
                    
                else:
                    path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_202004{d}_R17000_001.h5'
                    SM = SMAP_Data(path2,lat,lon)
                    SMAP_Data1.append(SM)
                    
        elif i<=152:
            d = i-121
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_05\SMAP_L3_SM_P_2020050{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_05\SMAP_L3_SM_P_202005{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
        elif i<=182:
            d = i-152
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_06\SMAP_L3_SM_P_2020060{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_06\SMAP_L3_SM_P_202006{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
        elif i<=213:
            d = i-182
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_07\SMAP_L3_SM_P_2020070{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_07\SMAP_L3_SM_P_202007{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
        elif i<=244:
            d = i-213
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_08\SMAP_L3_SM_P_2020080{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_08\SMAP_L3_SM_P_202008{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
        elif i<=274:
            d = i-244
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_09\SMAP_L3_SM_P_2020090{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_09\SMAP_L3_SM_P_202009{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
        elif i<=305:
            d = i-274
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_10\SMAP_L3_SM_P_2020100{d}_R17000_001.h5'
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_10\SMAP_L3_SM_P_202010{d}_R17000_001.h5'                
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
        elif i<=335:
            d = i-305
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_11\SMAP_L3_SM_P_2020110{d}_R17000_001.h5'                
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_11\SMAP_L3_SM_P_202011{d}_R17000_001.h5'                
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
        elif i<=366:
            d = i-335
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_12\SMAP_L3_SM_P_2020120{d}_R17000_001.h5'               
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_12\SMAP_L3_SM_P_202012{d}_R17000_001.h5'                
                SM = SMAP_Data(path2,lat,lon)
                SMAP_Data1.append(SM) 
    
    Day_Nodf = pd.DataFrame(Day_No)
    Day_Nodf.columns = ['Day_No']
    Day_Nodf['SMAP_SM'] = pd.DataFrame(SMAP_Data1)
    Day_Nodf.to_csv(f'SMAP_SM_Variations_{lat}_{lon}.csv',index = False)
  

'''-------------------------Taking out average ESA CCI Soil Moisture Data within its 0.5 x 0.5 degree grid cell-------------------------'''

'''
   Input : path: Path of the ESA CCI Daily Raw files
           lat and lon: Centroid of latitude and longitude of the 0.5 x 0.5 degree grid cell
   
   Output : Saving CSV file for a single pixel averaged ESA CCI soil moisture values for 366 days of year 2020 within that grid cell

'''

'''
Subfunction

'''

def SM_Data(path,latc,lonc):
    D = Dataset(path)
    Data1 = D.variables['sm'][0]
    df1 = pd.DataFrame(Data1)
    SM_Data = df1.replace(to_replace=-9999,value=np.nan)
    df = np.array(SM_Data)
    df_SM = df[235:275,1012:1078]
    
    lat = D.variables['lat']
    lon = D.variables['lon']
    ESACCI_Lat = np.array(lat)
    ESACCI_Lon = np.array(lon)
    X = ESACCI_Lat[235:275]
    Y = ESACCI_Lon[1012:1078]
    X1 = pd.DataFrame(X)
    Y1 = pd.DataFrame(Y)
    X2,Y2 = np.meshgrid(X1,Y1)
    Y2_ = np.transpose(Y2)
    X2_ = np.transpose(X2)
    lat = X2_.flatten()
    lon = Y2_.flatten()
    SM = (np.array(np.transpose(df_SM))).flatten()
    
    Df = pd.DataFrame(lat)
    Df.columns = ['lat']
    Df['lon']  = lon
    Df['SM']   = SM/100
    
    mask = ((Df['lat']>(latc-0.5)) & (Df['lat']<(latc+0.5)) & (Df['lon']>(lonc-0.5)) & (Df['lon']<(lonc+0.5)))
    Df = Df[mask]
    SM = Df['SM']
    if len(SM)>0:
        SM1 = (np.sum(SM))/len(SM)
    else:
        SM1 = np.nan
    return SM1


'''
Main function

'''

def SM_ESA_CCI(m,n,lat,lon):
    SM = []
    Day = []
    for i in range(m,n+1):
        Day.append(float(i))
        if i <= 31:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_01'
            d = i
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020010{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
               
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202001{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
        elif i <= 60:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_02'
            d = i-31
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020020{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202002{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
              
        elif i <= 91:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_03'
            d = i-60
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020030{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202003{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
        elif i <= 121:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_04'
            d = i-91
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020040{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202004{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
              
        elif i <= 152:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_05'
            d = i-121
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020050{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202005{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
        elif i <= 182:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_06'
            d = i-152
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020060{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202006{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
        elif i <= 213:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_07'
            d = i-182
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020070{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202007{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
        elif i <= 244:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_08'
            d = i-213
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020080{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202008{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
              
        elif i <= 274:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_09'
            d = i-244
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020090{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202009{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)

        elif i <= 305:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_10'
            d = i-274
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020100{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202010{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)

        elif i <= 335:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_11'
            d = i-305
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020110{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202011{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)

        elif i <= 366:
            Month = 'D:\EG\Project Data\ESA_CCI_DATA\Month_12'
            d = i-335
            if d<=9:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2020120{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
                
            else:
                filename = f'\C3S-SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-202012{d}000000-TCDR-v202012.0.0.nc'
                path = Month+filename
                Df = SM_Data(path,lat,lon)
                SM.append(Df)
    SM1 = pd.DataFrame(Day)
    SM1.columns = ['Day']
    SM1['SM'] = SM
    SM1.to_csv(f'ESACCI_SM_{lat}_{lon}.csv',index = False)
    
'''
   Input : lat and lon: Centroid of latitude and longitude of the 0.5 x 0.5 degree grid cell
           
   Output : Visualizing Correlation with plots 

'''
    
def ESA_SMAP_SM_Correlations(lat,lon):
    ESA_SM = pd.read_csv(f'D:\EG\Project Data\CYGNSS_Data_in_0p36Dg\ESACCI_RF_SM\ESACCI_SM_{lat}_{lon}.csv')
    SM_ESA = ESA_SM['SM']
    SMAP_SM = pd.read_csv(f'D:\EG\Project Data\CYGNSS_Data_in_0p36Dg\SMAP_RF_SM\SMAP_SM_Variations_{lat}_{lon}.csv')
    SM_SMAP = SMAP_SM['SMAP_SM'] 
    Df = pd.DataFrame(SM_ESA)
    Df.columns = ['ESA_SM']
    Df['SMAP_SM'] = SM_SMAP
    D = np.array(Df.corr())
    Day = SMAP_SM['Day_No']
    
    act   = pd.DataFrame(Df['SMAP_SM']).replace(to_replace=np.nan,value=0)
    pred  = pd.DataFrame(Df['ESA_SM']).replace(to_replace=np.nan,value=0)
    rmse  = mean_squared_error(act/100, pred/20, squared=False)
    
    plt.figure(figsize=(30,10))
    plt.scatter(Day,SM_SMAP/100,label='SMAP Soil Moisture by Modelling Brightness Temperature using Different Model')
    plt.scatter(Day,SM_ESA,label='ESA CCI Soil Moisture by Modelling Brightness Temperature using Different Model')
    plt.xlabel('Day_Number')
    plt.ylabel('Soil_Moisture')
    plt.legend()
    plt.title(f'Latitude: {lat} Longitude: {lon} Correlation: {np.round(((D[0][1])*100),2)}%  RMSE:{round(rmse,4)}')
    plt.rc('legend',fontsize=20)    