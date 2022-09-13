'''
This script is for masking the CYGNSS raw daily data into 36 Km x 36 Km single pixel of SMAP '''

import numpy as np
import pandas as pd
from pandas import DataFrame

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

def SR_CYGNSS(path):
    Cygnss_csv1 = pd.read_csv(path)
    mask = ((Cygnss_csv1['ddm_snr']>2) & (Cygnss_csv1['sp_rx_gain']>0) & (Cygnss_csv1['sp_inc_angle']<40))
    Cygnss_csv = Cygnss_csv1[mask]
    Data = Cygnss_csv
    SR2 = Surface_Reflectivity2(Data)             # Calculating Effective Surface Reflectivity
    DF2 = pd.DataFrame(Cygnss_csv['sp_lat'])
    DF2.columns = ['sp_lat']
    DF2['sp_lon'] = Cygnss_csv['sp_lon']
    DF2['SP_I']   = Cygnss_csv['sp_inc_angle']
    DF2['SR_eff'] = SR2
    mask1 = ((DF2['SR_eff']<20) & (DF2['SR_eff']>0))
    DF3 = DF2[mask1]
    S_1 = DF3
    return S_1

def SR_1(path,lat,lon):
    S_1  = SR_CYGNSS(path)
    mask = ((S_1['sp_lat']>(lat-0.18)) & (S_1['sp_lat']<(lat+0.18)) & (S_1['sp_lon']>(lon-0.18)) & (S_1['sp_lon']<(lon+0.18)))
    S_1  = S_1[mask]
    SR   = pd.DataFrame(S_1['SR_eff'])
    Inc  = S_1['SP_I']
    SR.columns = ['SR_eff']
    SR['SP_I'] = Inc
    SR['sp_lat'] = S_1['sp_lat']
    SR['sp_lon'] = S_1['sp_lon']
    return SR

def CYGNSS_Data_availability(m,n,lat,lon):
    for i in range(m,n+1):
        if i<=9:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_00{i}.csv'
            SR1 = SR_1(path1,lat,lon)
            SR1.to_csv(f'CYGNSS_Data_Availability_Day_00{i}_{lat}_{lon}.csv',index = False)
        elif i<=99:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
            SR1 = SR_1(path1,lat,lon)
            SR1.to_csv(f'CYGNSS_Data_Availability_Day_0{i}_{lat}_{lon}.csv',index = False)
        elif i<=366:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            SR1 = SR_1(path1,lat,lon)
            SR1.to_csv(f'CYGNSS_Data_Availability_Day_{i}_{lat}_{lon}.csv',index = False)