'''----------------------------------------------- Important Libraries-----------------------------------------------------------------'''

from netCDF4 import Dataset
import numpy as np
import pandas as pd
from pandas import DataFrame
from scipy.interpolate import griddata
import h5py
import matplotlib.pyplot as plt
import seaborn as sns


'''------------------------------------------Importing CYGNSS CSV files & Removing outliers -------------------------------------------
Removing outliers:
                                       1. Removing Observations with SNR value < 2 dB
                                       2. Removing Observations with a Receiver Antenna gain value < 0 dB
                                       3. Removing Observations with an Incidence angle > 40 degree
                                       
---------------------------------& Calculating Effective Surface Reflectivity for CYGNSS Co-ordinates---------------------------

    
         ((4π)**2*Pr*(R_r+R_t )**2)
Pr_eff = -----------------------------
              (λ**2*P_t*G_t*G_r )

R_r    = Range between specular point to the reciever(meter)
R_t    = Range between specular point to the transmitter(meter)
P_t    = GPS transmitting power(dB)
G_t    = Gain of the transmitting antenna(dB)
G_r    = Gain of the receiving antenna(dB)
λ      = 0.19 m( L1 frequency)

1. Taking DDM_peak as Pr : 
Input:  Path of the CYGNSS CSV file
Output: Effective Surface Reflectivity for CYGNSS Co-ordinates

'''
def Surface_Reflectivity2(Data):
    P1  =  Data['gps_ant_gain_db_i'][:]           # Gain of the transmitting antenna:(in dBi(Decibel Isotropic))    
    P2  =  Data['sp_rx_gain'][:]                  # Gain of the recieving antenna:(in dBi)    
    P3  =  Data['tx_to_sp_range'][:]              # Distance(m), specular point to the transmitter    
    P4  =  Data['rx_to_sp_range'][:]              # Distance(m), specular point to the reciever    
    P5  =  Data['gps_tx_power_db_w'][:]           # GPS transmitting power(RHCP Power in dB)           
    DDM_peak =  Data['peak of power_analog'][:]   # 17 Delay × 11 Doppler(in Watt) bins corresp. to 50 km^2    
    P_r  = DDM_peak*90
    T_rl = (P_r*(4*np.pi*(P3+P4))**2)/(P5*P1*P2*(0.19)**2)  # Effective Surface Reflectivity in dB
    return T_rl

'''------------------------------------------Importing CYGNSS CSV files & Removing outliers -------------------------------------------
Removing outliers:
                                       1. Removing Observations with SNR value < 2 dB
                                       2. Removing Observations with a Receiver Antenna gain value < 0 dB
                                       3. Removing Observations with an Incidence angle > 40 degree'''
                                       

def SR_CYGNSS(path):
    Cygnss_csv1 = pd.read_csv(path)
    mask = ((Cygnss_csv1['ddm_snr']>2) & (Cygnss_csv1['sp_rx_gain']>0) & (Cygnss_csv1['sp_inc_angle']<40))
    Cygnss_csv = Cygnss_csv1[mask]
    Data = Cygnss_csv
    SR2 = Surface_Reflectivity2(Data)             # Calculating Effective Surface Reflectivity
    DF2 = pd.DataFrame(Cygnss_csv['sp_lat'])
    DF2.columns = ['sp_lat']
    DF2['sp_lon'] = Cygnss_csv['sp_lon']
    DF2['SR_eff'] = SR2
    mask1 = ((DF2['SR_eff']<20) & (DF2['SR_eff']>0))
    DF3 = DF2[mask1]
    return DF3

'''-------------------------------------1. Importing the SMAP hdf files & 
                                        2. Concatinating both AM and PM data for individual day
                                        3. Masking the Concatenated SMAP data  
                                        (i)   Latitude_SMAP               
                                        (ii)  Longitude_SMAP               
                                        (iii) Soil Moisture_SMAP     
                                           within ganga catchment for SMAP co-ordinates within AOI------------------------------'''

'''
Input  : Path of the SMAP h5py files
Output : Masked SMAP soil moisture data within ganga catchment with SMAP co-ordinates

'''

def SMAP_Data(path):
    SMAP_data = h5py.File(path,'r')                 # Importing the SMAP hdf files  
    # PM data
    df1 = SMAP_data['Soil_Moisture_Retrieval_Data_PM']
    SM_lat1 = np.array(pd.DataFrame(df1['latitude_pm']))
    SM_lon1 = np.array(pd.DataFrame(df1['longitude_pm']))
    SM1 = np.array(pd.DataFrame(df1['soil_moisture_pm']))
    S1 = pd.DataFrame(SM_lat1.reshape(406*964,1))
    S2 = pd.DataFrame(SM_lon1.reshape(406*964,1))
    SM_1 = pd.DataFrame(SM1.reshape(406*964,1))
    
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
    
    # Masking the Concatenated SMAP data within ganga catchment for SMAP co-ordinates within AOI    
    mask = ((S_1['Lat']>21) & (S_1['Lat']<31) & (S_1['lon']>73) & (S_1['lon']<89) & (S_1['SM']>0))
    Df   = S_1[mask]
    return Df


def Corr(path1,path2,a,b):      
    
    DF2           = SR_CYGNSS(path1)    # CYGNSS data 
    CYGNSS_lat    = DF2['sp_lat']
    CYGNSS_lon    = DF2['sp_lon']
    CYGNSS_SR     = DF2['SR_eff'] 
    
    Df            = SMAP_Data(path2)    # SMAP data
    SMAP_Lat      = Df['Lat']
    SMAP_Lon      = Df['lon']
    SMAP_SM       = Df['SM'] 
    
    ## Interpolating CYGNSS surface reflectivity along SMAP coordinates
    SM_SMP1       = griddata((CYGNSS_lat,CYGNSS_lon),CYGNSS_SR,(SMAP_Lat,SMAP_Lon), method = a)
    
    ## Dataframe with SMAP lat and lon and SR and SM
    SmN           = pd.DataFrame(SM_SMP1)
    SmN.columns   = ['SR']
    SmN['SM']     = np.array(SMAP_SM)
    SmN['sp_lat'] = np.array(Df['Lat'])
    SmN['sp_lon'] = np.array(Df['lon'])
    
    ## Masking SR in between 0 and 20 dB
    mask1 = (SmN['SR']>0)
    SmN   = SmN[mask1]    
    mask2 = (SmN['SR']<20)
    SmN   = SmN[mask2]
    
    ## Applying change detection model on SR to retrieve Soil Moisture
    Min       = np.min(SmN['SR'])     ## SR corresponding to the driest condition
    Max       = np.max(SmN['SR'])     ## SR corresponding to the wettest condition
      
    SmN['SR'] = (SmN['SR']-Min)/(Max-Min)         ## Finding SM
    SmN       = SmN.dropna()
    mask3     = (SmN['SR']>0) & (SmN['SR']<0.7)   ## Taking only upto 0.7 Soil Moisture as more than this may be not due to soil surface
    SmN       = SmN[mask3]
    
    ## Finding correlation between CYGNSS SM and SMAP SM
    Cor1         = pd.DataFrame(SmN['SR'])
    Cor1.columns = ['SR_Derived_SM']
    Cor1['SM']   = SmN['SM']
    Corr         = np.array(Cor1.corr())
    CR           = np.round(Corr[0][1]*100,2)  
    
    
    ## Visualizing CYGNSS and SMAP SM
    Lat = np.array(SmN['sp_lat'])
    Lon = np.array(SmN['sp_lon'])
    lat,lon = np.meshgrid(Lat,Lon)
    SR  = np.array(SmN['SR'])
    SM  = np.array(SmN['SM'])
    SRSM = griddata((Lat,Lon),SR,(lat,lon), method='nearest')
    
    m1,n1 = SRSM.shape
        
    ## Masking the array without disturbing its shape and size
    for i1 in range(m1):
        for j1 in range(n1):
            a2 = SRSM[i1][j1]
            if (a2>0) and (a2<1):
                SRSM[i1][j1] = SRSM[i1][j1]
            else:
                SRSM[i1][j1] = np.nan
                    
    plt.figure(figsize=(20,5))
    plt.subplot(1,2,1)
    cmap = 'gist_ncar'
    plt.contourf(lon,lat,SRSM,cmap = cmap)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    cbar = plt.colorbar(shrink=1)
    cbar.set_label(f'Daily Soil Moisture')
    plt.title(f'''CYGNSS Derived Soil Moisture''')
    
    SMAPSM = griddata((Lat,Lon),SM,(lat,lon), method='nearest')
    plt.subplot(1,2,2)
    plt.contourf(lon,lat,SMAPSM,cmap = cmap)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    cbar = plt.colorbar(shrink=1)
    cbar.set_label(f'Daily Soil Moisture')
    plt.title(f'''SMAP Soil Moisture: Day:{b}, Correlation:{CR}%''')
    
    Points = []
    for i in range(len(SMAPSM.flatten())):
        Points.append(i)
        
    plt.figure(figsize=(40,10))
    plt.plot(Points[:2000],SMAPSM.flatten()[:2000],label = 'SMAP Soil Moisture')
    plt.plot(Points[:2000],SRSM.flatten()[:2000],label    = 'CYGNSS Derived Soil Moisture')
    plt.ylabel('''Surface reflectiviy and Soil Moisture''',fontsize=20)
    plt.xlabel('Points',fontsize=20)
    plt.legend(fontsize=20) 
    
    sns.lmplot(x='SR',y='SM',data=SmN,aspect=4,height=6)
    plt.xlabel('CYGNSS Derived Soil Moisture')
    plt.ylabel('SMAP Soil Moisture in cm3/cm3')
    plt.title(f'''Fitting a Regression Line using lmplot''')
    

'''------------------------------------------Function for Correlating all the 2020 data------------------------------------------

Input  : m-n(Period of time for which we are correlating the data)
         a = method of interpolation

Output : Correlation between the CYGNSS Soil Moisture derived using Change detection model and the Soil Moisture data of SMAP

'''    
    
def CYGNSS_SR_SMAP_SM_Corr(m,n,a):
    for i in range(m,n+1):
        if i<=31:
            d = i
            if d <= 9:                
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_00{i}.csv'
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_01\SMAP_L3_SM_P_2020010{d}_R17000_001.h5'
                b = f'0{i} Jan 2020'
                Corr(path1,path2,a,b)
            else:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_01\SMAP_L3_SM_P_202001{d}_R17000_001.h5'
                b = f'{i} Jan 2020'
                Corr(path1,path2,a,b)
        elif i<=60:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
            d = i-31
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_02\SMAP_L3_SM_P_2020020{d}_R17000_001.h5'
                b = f'0{d} Feb 2020'
                Corr(path1,path2,a,b)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_02\SMAP_L3_SM_P_202002{d}_R17000_001.h5'
                b = f'{d} Feb 2020'
                Corr(path1,path2,a,b)
        elif i<=91:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
            d = i-60
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_03\SMAP_L3_SM_P_2020030{d}_R17000_001.h5'
                b = f'0{d} March 2020'
                Corr(path1,path2,a,b)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_03\SMAP_L3_SM_P_202003{d}_R17000_001.h5'
                b = f'{d} March 2020'
                Corr(path1,path2,a,b)
        elif i<=121:
            d = i-91
            if i <= 99:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_2020040{d}_R17000_001.h5'
                b = f'0{d} April 2020'
                Corr(path1,path2,a,b)
            else:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
                if d<=9:
                    path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_2020040{d}_R17000_001.h5'
                    b = f'0{d} April 2020'
                    Corr(path1,path2,a,b)
                else:
                    path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_202004{d}_R17000_001.h5'
                    b = f'{d} April 2020'
                    Corr(path1,path2,a,b)
        elif i<=152:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-121
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_05\SMAP_L3_SM_P_2020050{d}_R17000_001.h5'
                b = f'0{d} May 2020'
                Corr(path1,path2,a,b)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_05\SMAP_L3_SM_P_202005{d}_R17000_001.h5'
                b = f'{d} May 2020'
                Corr(path1,path2,a,b)
        elif i<=182:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-152
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_06\SMAP_L3_SM_P_2020060{d}_R17000_001.h5'
                b = f'0{d} June 2020'
                Corr(path1,path2,a,b)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_06\SMAP_L3_SM_P_202006{d}_R17000_001.h5'
                b = f'{d} June 2020'
                Corr(path1,path2,a,b)
        elif i<=213:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-182
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_07\SMAP_L3_SM_P_2020070{d}_R17000_001.h5'
                b = f'0{d} July 2020'
                Corr(path1,path2,a,b)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_07\SMAP_L3_SM_P_202007{d}_R17000_001.h5'
                b = f'{d} July 2020'
                Corr(path1,path2,a,b)
        elif i<=244:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-213
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_08\SMAP_L3_SM_P_2020080{d}_R17000_001.h5'
                b = f'0{d} August 2020'
                Corr(path1,path2,a,b)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_08\SMAP_L3_SM_P_202008{d}_R17000_001.h5'
                b = f'{d} August 2020'
                Corr(path1,path2,a,b)
        elif i<=274:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-244
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_09\SMAP_L3_SM_P_2020090{d}_R17000_001.h5'
                b = f'0{d} Sept 2020'
                Corr(path1,path2,a,b)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_09\SMAP_L3_SM_P_202009{d}_R17000_001.h5'
                b = f'{d} Sept 2020'
                Corr(path1,path2,a,b)
        elif i<=305:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-274
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_10\SMAP_L3_SM_P_2020100{d}_R17000_001.h5'
                b = f'0{d} Oct 2020'
                Corr(path1,path2,a,b)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_10\SMAP_L3_SM_P_202010{d}_R17000_001.h5'
                b = f'{d} Oct 2020'
                Corr(path1,path2,a,b)
        elif i<=335:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-305
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_11\SMAP_L3_SM_P_2020110{d}_R17000_001.h5'
                b = f'0{d} Nov 2020'
                Corr(path1,path2,a,b)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_11\SMAP_L3_SM_P_202011{d}_R17000_001.h5'
                b = f'{d} Nov 2020'
                Corr(path1,path2,a,b)
        elif i<=366:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-335
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_12\SMAP_L3_SM_P_2020120{d}_R17000_001.h5'
                b = f'0{d} Dec 2020'
                Corr(path1,path2,a,b)
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_12\SMAP_L3_SM_P_202012{d}_R17000_001.h5'
                b = f'{d} Dec 2020'
                Corr(path1,path2,a,b)                              
                