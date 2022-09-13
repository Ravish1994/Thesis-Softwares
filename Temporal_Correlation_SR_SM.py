from netCDF4 import Dataset
import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import h5py
import seaborn as sns
from sklearn.metrics import mean_squared_error
import statistics

'''--------------------------------------------Function to Determine Surface Reflectivity-------------------------------------------'''

def Surface_Reflectivity2(Data):
    P1  =  Data['gps_ant_gain_db_i'][:]           # Gain of the transmitting antenna:(in dBi(Decibel Isotropic))    
    P2  =  Data['sp_rx_gain'][:]                  # Gain of the recieving antenna:(in dBi)    
    P3  =  Data['tx_to_sp_range'][:]              # Distance(m), specular point to the transmitter    
    P4  =  Data['rx_to_sp_range'][:]              # Distance(m), specular point to the reciever    
    P5  =  Data['gps_tx_power_db_w'][:]           # GPS transmitting power(RHCP Power in dB)           
    DDM_peak =  Data['peak of power_analog'][:]   # 17 Delay × 11 Doppler(in Watt) bins corresp. to 50 km^2    
    P_r = DDM_peak*90
    T_rl = (P_r*(4*np.pi*(P3+P4))**2)/(P5*P1*P2*(0.19)**2)  # Effective Surface Reflectivity in dB
    return T_rl


'''
   Input : path: Path of the CYGNSS Daily Raw files
           lat and lon: Centroid of latitude and longitude of the SMAP single grid cell
   
   Output : Averaged Surface Reflectivity and its Standard Deviation within that SMAP grid cell

'''
def SR_CYGNSS(path,lat,lon):
    
    ## Reading raw data file of CYGNSS
    Cygnss_csv1 = pd.read_csv(path)
    
    ## Masking the Raw data files by removing outliers
    # DDM_SNR        > 2
    # SP_Rx_Gain     > 0
    # SP_inclination < 40
    
    mask = ((Cygnss_csv1['ddm_snr']>2) & (Cygnss_csv1['sp_rx_gain']>0) & (Cygnss_csv1['sp_inc_angle']<40))
    Cygnss_csv = Cygnss_csv1[mask]
    Data = Cygnss_csv
    
    # Calculating Effective Surface Reflectivity
    SR2 = Surface_Reflectivity2(Data)                   
    DF2 = pd.DataFrame(Cygnss_csv['sp_lat'])
    DF2.columns = ['sp_lat']
    DF2['sp_lon'] = Cygnss_csv['sp_lon']
    DF2['SR_eff'] = SR2
    
    # Masking SR values between 0 and 20 to remove outliers
    mask1 = ((DF2['SR_eff']<20) & (DF2['SR_eff']>0))    
    DF3 = DF2[mask1]
    S_1 = DF3
    
    # Averaging surface reflectivity values within SMAP grid cell of 36 Km with lat,lon as its centroid
    mask = ((S_1['sp_lat']>(lat-0.18)) & (S_1['sp_lat']<(lat+0.18)) & (S_1['sp_lon']>(lon-0.18)) & (S_1['sp_lon']<(lon+0.18)))
    S_1  = S_1[mask]
    SR1   = S_1['SR_eff']
    if len(SR1)>1:                             ## Ignoring the cells having less than 2 data points
        SR  = np.mean(np.array(SR1))
        STD = np.std(np.array(SR1))
    else:
        SR  = np.nan
        STD = np.nan
    return SR,STD


'''----------------------------------Taking out single SMAP Soil Moisture Data within that particular grid cell-------------------------'''

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


'''-------------------------Temporal Correlation between CYGNSS surface reflectivity with SMAP soil moisture---------------------------'''

'''
   Input : m: starting day of the year
           n: ending day of the year
           lat and lon: Centroid of latitude and longitude of the SMAP single grid cell
           
   Output : Saving CSV file for a single pixel containing Day_no, Standard deviation of SR in a particular cell with day and SM

'''

def CYGNSS_SMAP_Data_availability(m,n,lat,lon):
    CYGNSS_Data = []
    CYGNSS_STD = []
    SMAP_Data1 = []
    Day_No = []
    for i in range(m,n+1):
        Day_No1 = i
        Day_No.append(Day_No1)
        if i<=31:
            d = i
            if d <= 9:                
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_00{i}.csv'
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_01\SMAP_L3_SM_P_2020010{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM    = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_01\SMAP_L3_SM_P_202001{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
        elif i<=60:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
            d = i-31
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_02\SMAP_L3_SM_P_2020020{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_02\SMAP_L3_SM_P_202002{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
        elif i<=91:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
            d = i-60
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_03\SMAP_L3_SM_P_2020030{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_03\SMAP_L3_SM_P_202003{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
        elif i<=121:
            d = i-91
            if i <= 99:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_2020040{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
                if d<=9:
                    path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_2020040{d}_R17000_001.h5'
                    
                    SR,SD = SR_CYGNSS(path1,lat,lon)
                    SM = SMAP_Data(path2,lat,lon)
                    CYGNSS_Data.append(SR)
                    CYGNSS_STD.append(SD)
                    SMAP_Data1.append(SM)
                    
                else:
                    path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_202004{d}_R17000_001.h5'
                    
                    SR,SD = SR_CYGNSS(path1,lat,lon)
                    SM = SMAP_Data(path2,lat,lon)
                    CYGNSS_Data.append(SR)
                    CYGNSS_STD.append(SD)
                    SMAP_Data1.append(SM)
                    
        elif i<=152:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-121
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_05\SMAP_L3_SM_P_2020050{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_05\SMAP_L3_SM_P_202005{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
        elif i<=182:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-152
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_06\SMAP_L3_SM_P_2020060{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_06\SMAP_L3_SM_P_202006{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
        elif i<=213:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-182
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_07\SMAP_L3_SM_P_2020070{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_07\SMAP_L3_SM_P_202007{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
        elif i<=244:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-213
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_08\SMAP_L3_SM_P_2020080{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_08\SMAP_L3_SM_P_202008{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
        elif i<=274:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-244
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_09\SMAP_L3_SM_P_2020090{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_09\SMAP_L3_SM_P_202009{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
        elif i<=305:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-274
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_10\SMAP_L3_SM_P_2020100{d}_R17000_001.h5'
                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_10\SMAP_L3_SM_P_202010{d}_R17000_001.h5'                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
        elif i<=335:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-305
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_11\SMAP_L3_SM_P_2020110{d}_R17000_001.h5'                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_11\SMAP_L3_SM_P_202011{d}_R17000_001.h5'                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
        elif i<=366:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-335
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_12\SMAP_L3_SM_P_2020120{d}_R17000_001.h5'               
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM    = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_12\SMAP_L3_SM_P_202012{d}_R17000_001.h5'                
                SR,SD = SR_CYGNSS(path1,lat,lon)
                SM    = SMAP_Data(path2,lat,lon)
                CYGNSS_Data.append(SR)
                CYGNSS_STD.append(SD)
                SMAP_Data1.append(SM) 
    
    Day_Nodf              = pd.DataFrame(Day_No)
    Day_Nodf.columns      = ['Day_No']
    Day_Nodf['Cygnss_SR'] = pd.DataFrame(CYGNSS_Data)
    Day_Nodf['Cygnss_SD'] = pd.DataFrame(CYGNSS_STD)
    Day_Nodf['SMAP_SM']   = pd.DataFrame(SMAP_Data1)
    Day_Nodf.to_csv(f'Annual_Variations_{lat}_{lon}.csv',index = False)


'''
   Input : path of all csv files saved previously
           n: numbers of files 
           
   Output : Visualizing Correlation with band plots of standard deviation

''' 
    
def Plotting_Variations(path,n):
    lat = [21,22,23,24,25,26,27,28,29,30,31]
    lon = [73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89]
    Lat = []
    Lon = []
    for i in range(len(lat)):
        lat1 = lat[i]
        for j in range(len(lon)):
            lon1 = lon[j]
            Lat.append(lat1)
            Lon.append(lon1)
    Lon = (np.array(Lon)).flatten()
    Lat = (np.array(Lat)).flatten()
    
    Path = []
    for i in range(len(Lat)):
        lat = str(Lat[i])
        lon = str(Lon[i])
        File_Name = lat+lon
        Path1 = path + '\Annual_Variations_' + lat + '_' + lon + '.csv'
        Path.append(Path1)
        
    Path = Path[0:n]
    for Path in Path:
        Path_1 = Path
        Lat_Lon = Path_1[-9:-4]
        lat = (Lat_Lon[0:2])
        lon = (Lat_Lon[3:])
        
        df          = pd.read_csv(Path_1)
        Day_No      = df['Day_No']
        CYGNSS_Data = df['Cygnss_SR']
        SMAP_Data1  = df['SMAP_SM']
        std1        = df['Cygnss_SD']
        
        df1            = pd.DataFrame(CYGNSS_Data)/20
        df1.columns    = ['CYGNSS_SR']
        SMAP_Data1     = df['SMAP_SM']/100
        df1['SMAP_SM'] = SMAP_Data1
        df1            = df1.dropna()
        CR             = (np.array(df1.corr()))*100
        CR1            = np.round(CR[0][1],2)
        
        act   = pd.DataFrame(SMAP_Data1).replace(to_replace=np.nan,value=0)
        pred  = pd.DataFrame(CYGNSS_Data).replace(to_replace=np.nan,value=0)
        rmse  = mean_squared_error(act/100, pred/20, squared=False)
        

        plt.figure(figsize=(40,15))
        plt.scatter(Day_No,SMAP_Data1,label = 'SMAP Soil Moisture')
        plt.plot(Day_No,CYGNSS_Data/20,label    = 'Normalized CYGNSS Derived Surface Reflectivity')
        plt.fill_between(Day_No,(CYGNSS_Data-std1)/20,(CYGNSS_Data+std1)/20)
        plt.ylabel('''Surface reflectiviy in dB''',fontsize=20)
        plt.xlabel('Days of Year 2020',fontsize=20)
        plt.title(f'Latitude:{lat} Longitude:{lon} Correlation:{CR1}% RMSE:{round(rmse,4)}',fontsize=20)
        plt.legend(fontsize=20) 
        
        sns.lmplot(x='CYGNSS_SR',y='SMAP_SM',data=df1,aspect=4,height=6)
        plt.xlabel('Normalized CYGNSS Derived Surface Reflectivity')
        plt.ylabel('SMAP Soil Moisture in cm3/cm3')
        plt.title(f'''Fitting a Regression Line using lmplot''')