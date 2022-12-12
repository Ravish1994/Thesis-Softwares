'''
1. Importing CYGNSS CSV files

2. Removing outliers from CYGNSS data

3. Calculating Effective Surface Reflectivity for CYGNSS Co-ordinates

4. Importing the SMAP hdf files

5. Masking the SMAP data  
   (i)   Latitude_SMAP               
   (ii)  Longitude_SMAP               
   (iii) Soil Moisture_SMAP
   for both AM and PM data and Concatinating all the data for individual day within ganga catchment for SMAP co-ordinates within AOI
   
6. Interpolating CYGNSS Effective Surface Reflectivity along SMAP co-ordinates

7. Correlating the Surface reflectivity of CYGNSS interpolated along SMAP co-ordinates and the Soil Moisture data of SMAP

'''
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from pandas import DataFrame

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

1. Taking DDM_SNR as Pr : 
Input:  Path of the CYGNSS CSV file
Output: Effective Surface Reflectivity for CYGNSS Co-ordinates

'''

def Surface_Reflectivity1(Data):
    v1    =  Data['gps_tx_power_db_w']                     # P_rt 
    v2    =  Data['gps_ant_gain_db_i']                     # G_t              
    P3    =  v1*v2                                         # P_rt*G_t     
    P1    =  Data['sp_rx_gain'][:]               # G_r    
    P2    =  Data['tx_to_sp_range'][:]           # R_ts    
    P3    =  Data['rx_to_sp_range'][:]           # R_sr    
    P_r   =  Data['ddm_snr'][:]                  # DDM_SNR        
    T_rl1 =  (P_r*(4*np.pi*(P2+P3))**2)/(P1*P3*(0.19)**2)  # Effective Surface Reflectivity in dB
    return T_rl1

'''
2. P_r    = The peak value of the analog DDM in dB
Input:  Path of the CYGNSS CSV file
Output: Effective Surface Reflectivity for CYGNSS Co-ordinates

'''
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
import h5py
def SMAP_Data(path):
    SMAP_data = h5py.File(path,'r')                 # Importing the SMAP hdf files
    
    # Reshaping data
    
    # PM data
    df1 = SMAP_data['Soil_Moisture_Retrieval_Data_PM']
    SM_lat1 = np.array(pd.DataFrame(df1['latitude_pm']))
    SM_lon1 = np.array(pd.DataFrame(df1['longitude_pm']))
    SM1 = np.array(pd.DataFrame(df1['soil_moisture_pm']))
    S1 = pd.DataFrame(SM_lat1.reshape(406*964,1))
    S2 = pd.DataFrame(SM_lon1.reshape(406*964,1))
    SM_1 = pd.DataFrame(SM1.reshape(406*964,1))
    
    # AM data
    df2 = SMAP_data['Soil_Moisture_Retrieval_Data_AM']
    SM_lat2 = np.array(pd.DataFrame(df2['latitude']))
    SM_lon2 = np.array(pd.DataFrame(df2['longitude']))        
    SM2 = np.array(pd.DataFrame(df2['soil_moisture']))
    T1 = pd.DataFrame(SM_lat2.reshape(406*964,1))
    T2 = pd.DataFrame(SM_lon2.reshape(406*964,1))  
    SM_2 = pd.DataFrame(SM2.reshape(406*964,1))

    # Concatenating AM and PM data  
    S_1 = pd.concat([S1,T1])
    S_2 = pd.concat([S2,T2])
    SM = pd.concat([SM_1,SM_2]) 
    
    # Creating dataframe of SMAP_Lat, SMAP_Lon and SMAP_SM
    S_1.columns = ['Lat']
    S_1['lon'] = S_2
    S_1['SM'] = SM  
    
    # Masking the Concatenated SMAP data within ganga catchment for SMAP co-ordinates within AOI
    
    mask = ((S_1['Lat']>21) & (S_1['Lat']<31) & (S_1['lon']>73) & (S_1['lon']<89) & (S_1['SM']>0))
    Df = S_1[mask]
    return Df

'''------------------------Interpolating CYGNSS Effective Surface Reflectivity along SMAP co-ordinates-------------------------'''
'''
Input  : Path of the CYGNSS CSV & SMAP h5py files
Output : Correlation between the Surface reflectivity of CYGNSS interpolated along SMAP co-ordinates and the Soil Moisture data of SMAP

'''
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

def Corr(path1,path2,a):      #**********************************Main correlating function*****************************************
    DF2 = SR_CYGNSS(path1)
    Df = SMAP_Data(path2)
    SMAP_Lat = Df['Lat']
    SMAP_Lon = Df['lon']
    SMAP_SM  = Df['SM']
    CYGNSS_lat = DF2['sp_lat']
    CYGNSS_lon = DF2['sp_lon']
    CYGNSS_SR = DF2['SR_eff']  
    SM_SMP1 = griddata((CYGNSS_lat,CYGNSS_lon),CYGNSS_SR,(SMAP_Lat,SMAP_Lon), method = a)
    SmN = pd.DataFrame(SM_SMP1)
    SmN.columns = ['SR']
    SmN['SM'] = np.array(SMAP_SM)
    S1M = SmN.dropna()
    Corr = np.array(S1M.corr())
    CR = np.round(Corr[0][1],2)
    
#     SmN['SMAP_Lat'] = np.array(SMAP_Lat)
#     X = SmN['SMAP_Lat']
#     Y = SmN['SM']
#     Z = SmN['SR']

#     plt.figure(figsize=(10,8))
#     plt.scatter(X,Y,label='SMAP Soil Moisture')
#     plt.scatter(X,Z,label='Pr_eff from CYGNSS in dB')
#     plt.xlabel('Latitude in degree')
#     plt.ylabel('Surface Reflectivity and Soil Moisture')
#     plt.title(f'Correlation = {CR}')
#     plt.legend()
#     plt.show()    
    return CR

'''------------------------------------------Software for Correlating all the 2020 data------------------------------------------

Input  : m-n(Period of time for which we are correlating the data)
         a = method of interpolation

Output : Correlation between the Surface reflectivity of CYGNSS interpolated along SMAP co-ordinates and the Soil Moisture data of SMAP

'''

def CYGNSS_SR_SMAP_SM_Corr(m,n,a):
    for i in range(m,n+1):
        if i<=31:
            d = i
            if d <= 9:                
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_00{i}.csv'
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_01\SMAP_L3_SM_P_2020010{d}_R17000_001.h5'
                print(f'0{i} Jan 2020 Correlation =',Corr(path1,path2,a))
            else:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_01\SMAP_L3_SM_P_202001{d}_R17000_001.h5'
                print(f'{i} Jan 2020 Correlation =',Corr(path1,path2,a))
        elif i<=60:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
            d = i-31
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_02\SMAP_L3_SM_P_2020020{d}_R17000_001.h5'
                print(f'0{d} Feb 2020 Correlation =',Corr(path1,path2,a))
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_02\SMAP_L3_SM_P_202002{d}_R17000_001.h5'
                print(f'{d} Feb 2020 Correlation =',Corr(path1,path2,a))
        elif i<=91:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
            d = i-60
            if d <= 9:                
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_03\SMAP_L3_SM_P_2020030{d}_R17000_001.h5'
                print(f'0{d} March 2020 Correlation =',Corr(path1,path2,a))
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_03\SMAP_L3_SM_P_202003{d}_R17000_001.h5'
                print(f'{d} March 2020 Correlation =',Corr(path1,path2,a))
        elif i<=121:
            d = i-91
            if i <= 99:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_2020040{d}_R17000_001.h5'
                print(f'0{d} April 2020 Correlation =',Corr(path1,path2,a))
            else:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
                if d<=9:
                    path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_2020040{d}_R17000_001.h5'
                    print(f'0{d} April 2020 Correlation =',Corr(path1,path2,a))
                else:
                    path2 = f'D:\EG\Project Data\SMAP_DATA\Month_04\SMAP_L3_SM_P_202004{d}_R17000_001.h5'
                    print(f'{d} April 2020 Correlation =',Corr(path1,path2,a))
        elif i<=152:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-121
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_05\SMAP_L3_SM_P_2020050{d}_R17000_001.h5'
                print(f'0{d} May 2020 Correlation =',Corr(path1,path2,a))
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_05\SMAP_L3_SM_P_202005{d}_R17000_001.h5'
                print(f'{d} May 2020 Correlation =',Corr(path1,path2,a))
        elif i<=182:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-152
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_06\SMAP_L3_SM_P_2020060{d}_R17000_001.h5'
                print(f'0{d} June 2020 Correlation =',Corr(path1,path2,a))
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_06\SMAP_L3_SM_P_202006{d}_R17000_001.h5'
                print(f'{d} June 2020 Correlation =',Corr(path1,path2,a))
        elif i<=213:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-182
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_07\SMAP_L3_SM_P_2020070{d}_R17000_001.h5'
                print(f'0{d} July 2020 Correlation =',Corr(path1,path2,a))
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_07\SMAP_L3_SM_P_202007{d}_R17000_001.h5'
                print(f'{d} July 2020 Correlation =',Corr(path1,path2,a))
        elif i<=244:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-213
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_08\SMAP_L3_SM_P_2020080{d}_R17000_001.h5'
                print(f'0{d} Aug 2020 Correlation =',Corr(path1,path2,a))
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_08\SMAP_L3_SM_P_202008{d}_R17000_001.h5'
                print(f'{d} Aug 2020 Correlation =',Corr(path1,path2,a))
        elif i<=274:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-244
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_09\SMAP_L3_SM_P_2020090{d}_R17000_001.h5'
                print(f'0{d} Sept 2020 Correlation =',Corr(path1,path2,a))
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_09\SMAP_L3_SM_P_202009{d}_R17000_001.h5'
                print(f'{d} Sept 2020 Correlation =',Corr(path1,path2,a))
        elif i<=305:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-274
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_10\SMAP_L3_SM_P_2020100{d}_R17000_001.h5'
                print(f'0{d} Oct 2020 Correlation =',Corr(path1,path2,a))
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_10\SMAP_L3_SM_P_202010{d}_R17000_001.h5'
                print(f'{d} Oct 2020 Correlation =',Corr(path1,path2,a))
        elif i<=335:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-305
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_11\SMAP_L3_SM_P_2020110{d}_R17000_001.h5'
                print(f'0{d} Nov 2020 Correlation =',Corr(path1,path2,a))
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_11\SMAP_L3_SM_P_202011{d}_R17000_001.h5'
                print(f'{d} Nov 2020 Correlation =',Corr(path1,path2,a))
        elif i<=366:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-335
            if d <= 9:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_12\SMAP_L3_SM_P_2020120{d}_R17000_001.h5'
                print(f'0{d} Dec 2020 Correlation =',Corr(path1,path2,a))
            else:
                path2 = f'D:\EG\Project Data\SMAP_DATA\Month_12\SMAP_L3_SM_P_202012{d}_R17000_001.h5'
                print(f'{d} Dec 2020 Correlation =',Corr(path1,path2,a))                
                
'''-------------------------------------------------------------END----------------------------------------------------------------'''