''' -----------------Software for masking the CYGNSS Level1 version 3.0 science data for the Ganga Catchment region----------------'''

''' 
Latitude of the Ganga Catchment: 21-31N
Longitude of the Ganga Catchment: 73-89E

Input  : path of the Day

output : dataframe with variables;
                        v1 = 'sp_lat'
                        v2 = 'sp_lon'
                        v3 = 'sp_inc_angle'
                        v4 = 'sp_rx_gain'
                        v5 = 'gps_tx_power_db_w'
                        v6 = 'gps_ant_gain_db_i'
                        v7 = 'ddm_snr'
                        v8 = 'ddm_noise_floor'
                        v9 = 'rx_to_sp_range'
                        v10 = 'tx_to_sp_range'
                        v11 = 'quality_flags'
                        v12 = 'power_analog'   
                        
Libraries needed'''

import netCDF4 as nc
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from pandas import DataFrame

'''------------------------------------------------- Stacking Single Satellite Data----------------------------------------------- '''

''' 
Input : 
1. data of single satellite
2. variable name
Output :
Stacked single column of the data for single satellite
'''
def SD(data,variable):
    if variable=='power_analog':
        P1 =  data.variables['power_analog'][:]
        X1 = P1.reshape(len(P1)*4*17*11,)
        X2 = pd.DataFrame(X1)
        X3 = X2.replace(to_replace=-9999,value=np.nan)
        X4 = X3.to_numpy()
        P = X4.reshape(len(P1),4,17,11)
        
        df1 = []
        for i in range(len(P)):
            for j in range(4):
                a = np.max(P[i][j][:])
                df1.append(a)
                
        df2 = np.array(df1)
        DDM_peak = df2.reshape(len(P),4)
        Pk = pd.DataFrame(DDM_peak)
        Pk.columns = ['p1','p2','p3','p4']
        x1_s1 = Pk.p1
        x1_s2 = Pk.p2
        x1_s3 = Pk.p3
        x1_s4 = Pk.p4
        m1 = x1_s1.append(x1_s2)
        m2 = m1.append(x1_s3)
        m3 = m2.append(x1_s4)
        Stacked = pd.DataFrame(m3)
    else:
        v1  =  data.variables[variable][:]                   
        u1 = pd.DataFrame(v1)
        P1 = u1.replace(to_replace=-9999,value=np.nan)
        P1.columns = ['p1','p2','p3','p4']
        x1_s1 = P1.p1
        x1_s2 = P1.p2
        x1_s3 = P1.p3
        x1_s4 = P1.p4
        m1 = x1_s1.append(x1_s2)
        m2 = m1.append(x1_s3)
        m3 = m2.append(x1_s4)
        Stacked = pd.DataFrame(m3)
    return Stacked  


''' ----------------------------Stacking Eight Satellites Data for Single Day for a particular variable----------------------------'''

''' 

Input : 
1. data of single day of all the satellites 
2. variable name

Output :
Stacked single column of the data for single day

'''
def Data_Stacking_D(Data1,Data2,Data3,Data4,Data5,Data6,Data7,Data8,v):
    SR1,SR2,SR3,SR4 = SD(Data1,v),SD(Data2,v),SD(Data3,v),SD(Data4,v)
    SR5,SR6,SR7,SR8 = SD(Data5,v),SD(Data6,v),SD(Data7,v),SD(Data8,v)
    sr1 = SR1.append(SR2)
    sr2 = sr1.append(SR3)
    sr3 = sr2.append(SR4)
    sr4 = sr3.append(SR5)
    sr5 = sr4.append(SR6)
    sr6 = sr5.append(SR7)
    Stacked = sr6.append(SR8)
    return Stacked

'''------------------------------------------Creating File Name for the days--------------------------------------------------'''

''' 

Input : 
1. i = satellite order from 1 to 8
2. j = month order from 1 to 12
3. d = day order from 1 to 31

Output : Filename

'''
def file(i,j,d):
    
    if j == 1:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020010{str(d)}-000000-e2020010{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202001{str(d)}-000000-e202001{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
            
    elif j == 2:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020020{str(d)}-000000-e2020020{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202002{str(d)}-000000-e202002{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'

    elif j == 3:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020030{str(d)}-000000-e2020030{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202003{str(d)}-000000-e202003{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
            
    elif j == 4:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020040{str(d)}-000000-e2020040{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202004{str(d)}-000000-e202004{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'

    elif j == 5:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020050{str(d)}-000000-e2020050{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202005{str(d)}-000000-e202005{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'            

    elif j == 6:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020060{str(d)}-000000-e2020060{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202006{str(d)}-000000-e202006{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'

    elif j == 7:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020070{str(d)}-000000-e2020070{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202007{str(d)}-000000-e202007{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
            
    elif j == 8:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020080{str(d)}-000000-e2020080{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202008{str(d)}-000000-e202008{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'

    elif j == 9:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020090{str(d)}-000000-e2020090{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202009{str(d)}-000000-e202009{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc' 
            
    elif j == 10:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020100{str(d)}-000000-e2020100{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202010{str(d)}-000000-e202010{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'

    elif j == 11:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020110{str(d)}-000000-e2020110{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202011{str(d)}-000000-e202011{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
            
    elif j == 12:
        if (d>=1) and (d<=9):
            file1 = f'cyg0{i}.ddmi.s2020120{str(d)}-000000-e2020120{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
        else:
            file1 = f'cyg0{i}.ddmi.s202012{str(d)}-000000-e202012{str(d)}-235959.l1.power-brcs.a30.d31.nc.nc'
            
    return file1

''' ------------------------------------------Stacking all satellites data for a  single day------------------------------'''

def DF_Data(path,v,d,j):
    f1,f2,f3,f4,f5,f6,f7,f8 = file(1,j,d),file(2,j,d),file(3,j,d),file(4,j,d),file(5,j,d),file(6,j,d),file(7,j,d),file(8,j,d)
    Path1,Path2,Path3,Path4 = f'{path}\{f1}',f'{path}\{f2}',f'{path}\{f3}',f'{path}\{f4}'
    Path5,Path6,Path7,Path8 = f'{path}\{f5}',f'{path}\{f6}',f'{path}\{f7}',f'{path}\{f8}'
    Data1,Data2 = Dataset(Path1,'r'),Dataset(Path2,'r')
    Data3,Data4 = Dataset(Path3,'r'),Dataset(Path4,'r')
    Data5,Data6 = Dataset(Path5,'r'),Dataset(Path6,'r')
    Data7,Data8 = Dataset(Path7,'r'),Dataset(Path8,'r')
    df = Data_Stacking_D(Data1,Data2,Data3,Data4,Data5,Data6,Data7,Data8,v)
    return df 

''' ------------------------------------------Filtering Data within Ganga Catchment---------------------------------------------- '''

def df_Data(path,j,d):
    
    v1 = 'sp_lat'
    v2 = 'sp_lon'
    v3 = 'sp_inc_angle'
    v4 = 'sp_rx_gain'
    v5 = 'gps_tx_power_db_w'
    v6 = 'gps_ant_gain_db_i'
    v7 = 'ddm_snr'
    v8 = 'ddm_noise_floor'
    v9 = 'rx_to_sp_range'
    v10 = 'tx_to_sp_range'
    v11 = 'quality_flags'
    v12 = 'power_analog'
    
    df1,df2,df3,df4    = DF_Data(path,v1,d,j),DF_Data(path,v2,d,j),DF_Data(path,v3,d,j),DF_Data(path,v4,d,j)
    df5,df6,df7,df8    = DF_Data(path,v5,d,j),DF_Data(path,v6,d,j),DF_Data(path,v7,d,j),DF_Data(path,v8,d,j)
    df9,df10,df11,df12 = DF_Data(path,v9,d,j),DF_Data(path,v10,d,j),DF_Data(path,v11,d,j),DF_Data(path,v12,d,j)
    
    df1.columns = ['sp_lat']
    df1['sp_lon'] = df2
    df1['sp_inc_angle'] = df3
    df1['sp_rx_gain'] = df4
    df1['gps_tx_power_db_w'] = df5
    df1['gps_ant_gain_db_i'] = df6
    df1['ddm_snr'] = df7
    df1['ddm_noise_floor'] = df8
    df1['rx_to_sp_range'] = df9
    df1['tx_to_sp_range'] = df10
    df1['quality_flags'] = df11
    df1['peak of power_analog'] = df12 
    
    mask = ((df1['sp_lat']>21) & (df1['sp_lat']<31)) & ((df1['sp_lon']>73) & (df1['sp_lon']<89))
    Df = df1[mask]
    
    Df_D = Df.dropna()
    return Df_D


'''------------------------------------------------Converting DataFrame into CSV---------------------------------------------------'''

'''
Input : m = number of days initial
        n = upto the day to which data needed
        
output : csv_file of 2020 for given period m-n

'''
def Ganga_Catchment_CSV(m,n):
    for i in range(m,n+1):
        if i<=31:
            if i <= 9:
                path = f'D:\EG\Project Data\CYGNSS DATA\Day_00{i}'
                df = df_Data(path,1,i)
                df.to_csv(f'Day_00{i}.csv',index = False)
            else:
                path = f'D:\EG\Project Data\CYGNSS DATA\Day_0{i}'
                df = df_Data(path,1,i)
                df.to_csv(f'Day_0{i}.csv',index = False)
        elif i<=60:
            path = f'D:\EG\Project Data\CYGNSS DATA\Day_0{i}'
            df = df_Data(path,2,i-31)
            df.to_csv(f'Day_0{i}.csv',index = False)
        elif i<=91:
            path = f'D:\EG\Project Data\CYGNSS DATA\Day_0{i}'
            df = df_Data(path,3,i-60)
            df.to_csv(f'Day_0{i}.csv',index = False)
        elif i<=121:
            if i <= 99:
                path = f'D:\EG\Project Data\CYGNSS DATA\Day_0{i}'
                df = df_Data(path,4,i-91)
                df.to_csv(f'Day_0{i}.csv',index = False)
            else:
                path = f'D:\EG\Project Data\CYGNSS DATA\Day_{i}'
                df = df_Data(path,4,i-91)
                df.to_csv(f'Day_{i}.csv',index = False)
        elif i<=152:
            path = f'D:\EG\Project Data\CYGNSS DATA\Day_{i}'
            df = df_Data(path,5,i-121)
            df.to_csv(f'Day_{i}.csv',index = False)
        elif i<=182:
            path = f'D:\EG\Project Data\CYGNSS DATA\Day_{i}'
            df = df_Data(path,6,i-152)
            df.to_csv(f'Day_{i}.csv',index = False)
        elif i<=213:
            path = f'D:\EG\Project Data\CYGNSS DATA\Day_{i}'
            df = df_Data(path,7,i-182)
            df.to_csv(f'Day_{i}.csv',index = False)
        elif i<=244:
            path = f'D:\EG\Project Data\CYGNSS DATA\Day_{i}'
            df = df_Data(path,8,i-213)
            df.to_csv(f'Day_{i}.csv',index = False)
        elif i<=274:
            path = f'D:\EG\Project Data\CYGNSS DATA\Day_{i}'
            df = df_Data(path,9,i-244)
            df.to_csv(f'Day_{i}.csv',index = False)
        elif i<=305:
            path = f'D:\EG\Project Data\CYGNSS DATA\Day_{i}'
            df = df_Data(path,10,i-274)
            df.to_csv(f'Day_{i}.csv',index = False)
        elif i<=335:
            path = f'D:\EG\Project Data\CYGNSS DATA\Day_{i}'
            df = df_Data(path,11,i-305)
            df.to_csv(f'Day_{i}.csv',index = False)
        elif i<=366:
            path = f'D:\EG\Project Data\CYGNSS DATA\Day_{i}'
            df = df_Data(path,12,i-335)
            df.to_csv(f'Day_{i}.csv',index = False)
            
'''-------------------------------------------------------END----------------------------------------------------------------'''