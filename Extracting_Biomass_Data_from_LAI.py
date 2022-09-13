from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc

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

def SR_CYGNSS(path,lat,lon):
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
    S_1 = DF3
    mask = ((S_1['sp_lat']>(lat-0.18)) & (S_1['sp_lat']<(lat+0.18)) & (S_1['sp_lon']>(lon-0.18)) & (S_1['sp_lon']<(lon+0.18)))
    S_1 = S_1[mask]
    SR = S_1['SR_eff']
    if len(SR)>0:
        SR = (np.sum(SR))/len(SR)
    else:
        SR = np.nan
    return SR

def SPI_CYGNSS(path,lat,lon):
    Cygnss_csv1 = pd.read_csv(path)
    mask = ((Cygnss_csv1['ddm_snr']>2) & (Cygnss_csv1['sp_rx_gain']>0) & (Cygnss_csv1['sp_inc_angle']<40))
    Cygnss_csv = Cygnss_csv1[mask]
    Data = Cygnss_csv
    SR2 = Surface_Reflectivity2(Data)             # Calculating Effective Surface Reflectivity
    DF2 = pd.DataFrame(Cygnss_csv['sp_lat'])
    DF2.columns = ['sp_lat']
    DF2['sp_lon'] = Cygnss_csv['sp_lon']
    DF2['SR_eff'] = SR2
    DF2['sp_inc_angle'] = Cygnss_csv['sp_inc_angle']
    mask1 = ((DF2['SR_eff']<20) & (DF2['SR_eff']>0))
    DF3 = DF2[mask1]
    S_1 = DF3
    mask = ((S_1['sp_lat']>(lat-0.18)) & (S_1['sp_lat']<(lat+0.18)) & (S_1['sp_lon']>(lon-0.18)) & (S_1['sp_lon']<(lon+0.18)))
    S_1 = S_1[mask]
    SPI = S_1['sp_inc_angle']
    if len(SPI)>0:
        SPI = (np.sum(SPI))/len(SPI)
    else:
        SPI = np.nan
    return SPI

def GLDAS_Data(path,lat1,lon1):
    SM = Dataset(path,'r')                 # Importing the SM files  
    
    lon = np.array(SM['lon'])
    lat = np.array(SM['lat'])
    SoilMoist_S_tavg = pd.DataFrame(np.array(SM['SoilMoist_S_tavg'][0]))
    SoilMoist_S_tavg = np.array(SoilMoist_S_tavg.replace(to_replace=-9999,value=np.nan))

    Lon,Lat = np.meshgrid(lon,lat)  
    S_1 = pd.DataFrame(Lat.flatten())
    S_2 = Lon.flatten()
    SM  = SoilMoist_S_tavg.flatten() 
    
    # Creating dataframe of GLDAS_Lat, GLDAS_Lon and GLDAS_SM
    S_1.columns = ['Lat']
    S_1['lon'] = S_2
    S_1['SM'] = SM  
    mask = ((S_1['Lat']>(lat1-0.18)) & (S_1['Lat']<(lat1+0.18)) & (S_1['lon']>(lon1-0.18)) & (S_1['lon']<(lon1+0.18)) & (S_1['SM']>0))
    Df = S_1[mask]
    SM2 = np.array(Df['SM'])
    if len(SM2)>0:
        SM1 = float((np.sum(SM2))/len(SM2))
    else:
        SM1 = np.nan
    return SM1/10

def CYGNSS_GLDAS_Data_availability(m,n,lat,lon):
    CYGNSS_Backscatter = []
    CYGNSS_SPI = []
    GLDAS_SM = []
    
    Day_No = []
    for i in range(m,n+1):
        Day_No1 = i
        Day_No.append(Day_No1)
        if i<=31:
            d = i
            if d <= 9:                
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_00{i}.csv'
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020010{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202001{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
        elif i<=60:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
            d = i-31
            if d <= 9:                
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020020{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202002{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
        elif i<=91:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
            d = i-60
            if d <= 9:                
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020030{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202003{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
        elif i<=121:
            d = i-91
            if i <= 99:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020040{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
                if d<=9:
                    path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020040{d}.022.nc4.SUB.nc4'
                    
                    SR = SR_CYGNSS(path1,lat,lon)
                    SPI = SPI_CYGNSS(path1,lat,lon)
                    SM = GLDAS_Data(path2,lat,lon)
                    CYGNSS_Backscatter.append(SR)
                    CYGNSS_SPI.append(SPI)
                    GLDAS_SM.append(SM)
                    
                else:
                    path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202004{d}.022.nc4.SUB.nc4'
                    
                    SR = SR_CYGNSS(path1,lat,lon)
                    SPI = SPI_CYGNSS(path1,lat,lon)
                    SM = GLDAS_Data(path2,lat,lon)
                    CYGNSS_Backscatter.append(SR)
                    CYGNSS_SPI.append(SPI)
                    GLDAS_SM.append(SM)
                    
        elif i<=152:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-121
            if d <= 9:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020050{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202005{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
        elif i<=182:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-152
            if d <= 9:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020060{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202006{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
        elif i<=213:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-182
            if d <= 9:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020070{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202007{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
        elif i<=244:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-213
            if d <= 9:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020080{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202008{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
        elif i<=274:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-244
            if d <= 9:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020090{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202009{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
        elif i<=305:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-274
            if d <= 9:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020100{d}.022.nc4.SUB.nc4'
                
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202010{d}.022.nc4.SUB.nc4'               
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
        elif i<=335:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-305
            if d <= 9:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020110{d}.022.nc4.SUB.nc4'              
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202011{d}.022.nc4.SUB.nc4'              
                SSR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
        elif i<=366:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            d = i-335
            if d <= 9:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A2020120{d}.022.nc4.SUB.nc4'              
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
                
            else:
                path2 = f'D:\EG\Project Data\GLDAS_SM\GLDAS_SM_Ganga_Catchment\GLDAS_CLSM025_DA1_D.A202012{d}.022.nc4.SUB.nc4'              
                SR = SR_CYGNSS(path1,lat,lon)
                SPI = SPI_CYGNSS(path1,lat,lon)
                SM = GLDAS_Data(path2,lat,lon)
                CYGNSS_Backscatter.append(SR)
                CYGNSS_SPI.append(SPI)
                GLDAS_SM.append(SM)
    
    Day_Nodf = pd.DataFrame(Day_No)
    Day_Nodf.columns = ['Day_No']
    Day_Nodf['CYGNSS_Backscatter'] = pd.DataFrame(CYGNSS_Backscatter)
    Day_Nodf['CYGNSS_SPI'] = pd.DataFrame(CYGNSS_SPI)
    Day_Nodf['GLDAS_SM'] = pd.DataFrame(GLDAS_SM)
    Day_Nodf.to_csv(f'Annual_Variations_{lat}_{lon}.csv',index = False)
    

'''String definition for the month number'''
def Month(i):
    if i<=9:
        b = '0'+str(i)
        return b
    elif i>=10 and i <=12:
        b = str(i)
        return b  
    
'''String definition for the day number in month(29 or 30 or 31)'''
def Day1(i):
    if i<=9:
        b = '0'+str(i)
        return b
    elif i>=10 and i <=31:
        b = str(i)
        return b
    
''' Formatting Date '''    
def MMDD(i):
    if (i>=1) and (i<=31):
        b1 = Month(1)                            ## Month Number(Jan)
        b2 = Day1(i)                             ## Day number in month
        cw = b1+b2
        dw = 1
                
    elif (i>=32) and (i<=59):
        b1 = Month(2)                            ## Month Number(Feb)
        b2 = Day1(i-31)                          ## Day number in month
        cw = b1+b2
        dw = 2
                
    elif (i>=60) and (i<=90):
        b1 = Month(3)                           ## Month Number(Mar)
        b2 = Day1(i-59)                         ## Day number in month
        cw = b1+b2
        dw = 3
                
    elif (i>=91) and (i<=120):
        b1 = Month(4)                           ## Month Number(Apr)
        b2 = Day1(i-90)                         ## Day number in month
        cw = b1+b2
        dw = 4
        
    elif (i>=121) and (i<=151):
        b1 = Month(5)                           ## Month Number(May)
        b2 = Day1(i-120)                        ## Day number in month
        cw = b1+b2
        dw = 5
        
    elif (i>=152) and (i<=181):
        b1 = Month(6)                           ## Month Number(Jun)
        b2 = Day1(i-151)                        ## Day number in month
        cw = b1+b2
        dw = 6
        
    elif (i>=182) and (i<=212):
        b1 = Month(7)                           ## Month Number(Jul)
        b2 = Day1(i-181)                        ## Day number in month
        cw = b1+b2
        dw = 7
        
    elif (i>=213) and (i<=243):
        b1 = Month(8)                           ## Month Number(Aug)
        b2 = Day1(i-212)                        ## Day number in month
        cw = b1+b2
        dw = 8
        
    elif (i>=244) and (i<=273):
        b1 = Month(9)                           ## Month Number(Sept)
        b2 = Day1(i-243)                        ## Day number in month
        cw = b1+b2
        dw = 9
        
    elif (i>=274) and (i<=304):
        b1 = Month(10)                          ## Month Number(Oct)
        b2 = Day1(i-273)                        ## Day number in month 
        cw = b1+b2
        dw = 10
        
    elif (i>=305) and (i<=334):
        b1 = Month(11)                          ## Month Number(Nov)
        b2 = Day1(i-304)                        ## Day number in month
        cw = b1+b2
        dw = 11
        
    else:
        b1 = Month(12)                          ## Month Number(Dec)
        b2 = Day1(i-334)                        ## Day number in month
        cw = b1+b2
        dw = 12
    return cw,dw

### For 2018 LAI data

''' Formatting Date '''    
def MMDD_2018(i):
    if (i>=1) and (i<=31):
        b1 = Month(1)                            ## Month Number(Jan)
        b2 = Day1(i)                             ## Day number in month
        cw = b1+b2
        dw = 1
                
    elif (i>=32) and (i<=59):
        b1 = Month(2)                            ## Month Number(Feb)
        b2 = Day1(i-31)                          ## Day number in month
        cw = b1+b2
        dw = 2
                
    elif (i>=60) and (i<=90):
        b1 = Month(3)                           ## Month Number(Mar)
        b2 = Day1(i-59)                         ## Day number in month
        cw = b1+b2
        dw = 3
                
    elif (i>=91) and (i<=120):
        b1 = Month(4)                           ## Month Number(Apr)
        b2 = Day1(i-90)                         ## Day number in month
        cw = b1+b2
        dw = 4
        
    elif (i>=121) and (i<=151):
        b1 = Month(5)                           ## Month Number(May)
        b2 = Day1(i-120)                        ## Day number in month
        cw = b1+b2
        dw = 5
        
    elif (i>=152) and (i<=181):
        b1 = Month(6)                           ## Month Number(Jun)
        b2 = Day1(i-151)                        ## Day number in month
        cw = b1+b2
        dw = 6
        
    elif (i>=182) and (i<=212):
        b1 = Month(7)                           ## Month Number(Jul)
        b2 = Day1(i-181)                        ## Day number in month
        cw = b1+b2
        dw = 7
        
    elif (i>=213) and (i<=243):
        b1 = Month(8)                           ## Month Number(Aug)
        b2 = Day1(i-212)                        ## Day number in month
        cw = b1+b2
        dw = 8
        
    elif (i>=244) and (i<=273):
        b1 = Month(9)                           ## Month Number(Sept)
        b2 = Day1(i-243)                        ## Day number in month
        cw = b1+b2
        dw = 9
        
    elif (i>=274) and (i<=304):
        b1 = Month(10)                          ## Month Number(Oct)
        b2 = Day1(i-273)                        ## Day number in month 
        cw = b1+b2
        dw = 10
        
    elif (i>=305) and (i<=334):
        b1 = Month(11)                          ## Month Number(Nov)
        b2 = Day1(i-304)                        ## Day number in month
        cw = b1+b2
        dw = 11
        
    else:
        b1 = Month(12)                          ## Month Number(Dec)
        b2 = Day1(i-334)                        ## Day number in month
        cw = b1+b2
        dw = 12
    return cw,dw

def Daily_LAI_2020(m,n):
    for i in range(m,n+1):
        DD,Mon = MMDD(i)
        path   = f'D:\EG\Project Data\LAI Data\Month{Mon}'
        filename = f'\AVHRR-Land_v005-preliminary_AVH15C1_NOAA-19_2020{DD}.nc'
        Path = path + filename
        Data = nc.Dataset(Path,'r')

        LAI1 = np.array(Data['LAI'][0,1250:1360,5060:5185])               # For Ganga [0,1179:1380,5060:5380]
        LAI  = (pd.DataFrame(LAI1)).replace(to_replace=-100,value=np.nan)

        plt.figure(figsize = (10,4))
        plt.imshow(LAI,cmap = 'gist_ncar')
        cbar = plt.colorbar()
        cbar.set_label('LAI')
        plt.title(f'LAI in 50 m resolution {DD[0:2]}-{DD[2:]}-2020')
        
        
def masking(lat1,lon1,Df):
    mask = ((Df['Latitude']>(lat1-0.18)) & (Df['Latitude']<(lat1+0.18)) & (Df['Longitude']>(lon1-0.18)) & (Df['Longitude']<(lon1+0.18)))
    Df = Df[mask]
    return Df

def LAI_Daily_Data_2020(lat1,lon1):
    for i in range(1,366):
        DD,Mon = MMDD(i)
        path = f'D:\EG\Project Data\LAI Data\Month{Mon}'
        filename = f'\AVHRR-Land_v005-preliminary_AVH15C1_NOAA-19_2020{DD}.nc'
        Path = path+filename
        Data = nc.Dataset(Path,'r')
        
        LAI1 = np.array(Data['LAI'][0,1250:1360,5060:5185])
        LAI  = np.array((pd.DataFrame(LAI1)).replace(to_replace=-100,value=np.nan)).flatten()
        Lat = np.array(Data['latitude'])
        Lat = Lat[1250:1360]
        Lon = np.array(Data['longitude'])
        Lon = Lon[5060:5185]
        
        Lon,Lat = np.meshgrid(Lon,Lat)
        Lat = Lat.flatten()
        Lon = Lon.flatten()
        
        df  = pd.DataFrame(Lat)
        df.columns = ['Latitude']
        df['Longitude'] = Lon
        df['LAI'] = LAI
        
        DF = masking(lat1,lon1,df)
        DF.to_csv(f'LAI_Data_Within_36_Km_Resolution_Cell_Day_{i}_{lat1}_{lon1}.csv',index = False)
        
        
def LAI_Daily_2020_Ploting(lat1,lon1,j):
    Day_No = []
    LAI_Y = []
    for i in range(1,366):
        Day_No.append(i)
        path = f'D:\EG\Project Data\LAI Data\LAI_{lat1}_{lon1}'
        filename = f'\LAI_Data_Within_36_Km_Resolution_Cell_Day_{i}_{lat1}_{lon1}.csv'
        Path = path+filename
        Data = pd.read_csv(Path)
        Data = Data.replace(to_replace=np.nan,value=0)
        LAI3 = Data['LAI'][j]
        LAI_Y.append(LAI3)
    Lat = Data['Latitude'][j]
    Lon = Data['Longitude'][j]
    
    plt.figure(figsize=(40,10))
    plt.scatter(Day_No,LAI_Y)
    plt.bar(Day_No,LAI_Y)
    plt.ylabel('''LAI Observations''',fontsize=20)
    plt.xlabel('Days of Year 2020',fontsize=20)
    plt.title(f'LAI observations within 36 Km grid having centroid_Lat:{lat1}_Lon:{lon1}  Target Lat: {Lat} Lon: {Lon}',fontsize=20)

def Daily_LAI_2018(m,n):
    for i in range(m,n+1):
        DD,Mon = MMDD_2018(i)
        path = f'D:\EG\Project Data\LAI_2018\Month{Mon}'
        filename = f'\AVHRR-Land_v005-preliminary_AVH15C1_NOAA-19_2018{DD}.nc'
        Path = path+filename
        Data = nc.Dataset(Path,'r')

        LAI1 = np.array(Data['LAI'][0,1250:1360,5060:5185])
        LAI = (pd.DataFrame(LAI1)).replace(to_replace=-100,value=np.nan)

        plt.figure(figsize=(10,4))
        plt.imshow(LAI,cmap='gist_ncar')
        cbar = plt.colorbar()
        cbar.set_label('LAI')
        plt.title(f'LAI in 50 m resolution {DD[0:2]}-{DD[2:]}-2018')
    
    
def LAI_Daily_Data_2018(lat1,lon1):
    for i in range(1,366):
        DD,Mon = MMDD_2018(i)
        path = f'D:\EG\Project Data\LAI_2018\Month{Mon}'
        filename = f'\AVHRR-Land_v005-preliminary_AVH15C1_NOAA-19_2018{DD}.nc'
        Path = path+filename
        Data = nc.Dataset(Path,'r')
        
        LAI1 = np.array(Data['LAI'][0,1250:1360,5060:5185])
        LAI  = np.array((pd.DataFrame(LAI1)).replace(to_replace=-100,value=np.nan)).flatten()
        Lat = np.array(Data['latitude'])
        Lat = Lat[1250:1360]
        Lon = np.array(Data['longitude'])
        Lon = Lon[5060:5185]
        
        Lon,Lat = np.meshgrid(Lon,Lat)
        Lat = Lat.flatten()
        Lon = Lon.flatten()
        
        df  = pd.DataFrame(Lat)
        df.columns = ['Latitude']
        df['Longitude'] = Lon
        df['LAI'] = LAI
        
        DF = masking(lat1,lon1,df)
        DF.to_csv(f'LAI_Data_Within_36_Km_Resolution_Cell_Day_{i}_{lat1}_{lon1}.csv',index = False)
        
        
def LAI_Daily_Ploting_2018(lat1,lon1,j):
    Day_No = []
    LAI_Y = []
    for i in range(1,366):
        Day_No.append(i)
        path = f'D:\EG\Project Data\LAI_2018\LAI_{lat1}_{lon1}'
        filename = f'\LAI_Data_Within_36_Km_Resolution_Cell_Day_{i}_{lat1}_{lon1}.csv'
        Path = path+filename
        Data = pd.read_csv(Path)
        Data = Data.replace(to_replace=np.nan,value=0)
        LAI3 = Data['LAI'][j]
        LAI_Y.append(LAI3)
    Lat = Data['Latitude'][j]
    Lon = Data['Longitude'][j]
    
    plt.figure(figsize=(40,10))
    plt.scatter(Day_No,LAI_Y)
    plt.bar(Day_No,LAI_Y)
    plt.ylabel('''LAI Observations''',fontsize=20)
    plt.xlabel('Days of Year 2018',fontsize=20)
    plt.title(f'LAI observations within 36 Km grid having centroid_Lat:{lat1}_Lon:{lon1}  Target Lat: {Lat} Lon: {Lon}',fontsize=20)
    

def Averaging_LAI_2018(m,n):
    LAI = np.array([[0 for col in range(125)] for row in range(110)])

    for i in range(m,n):
        DD,Mon = MMDD_2018(i)
        path = f'D:\EG\Project Data\LAI_2018\Month{Mon}'
        filename = f'\AVHRR-Land_v005-preliminary_AVH15C1_NOAA-19_2018{DD}.nc'
        CYGNSS_LAI_Path = path+filename
        Data = nc.Dataset(CYGNSS_LAI_Path,'r')

        LAI1 = np.array(Data['LAI'][0,1250:1360,5060:5185])
        LAI2 = (pd.DataFrame(LAI1)).replace(to_replace=-100,value=0)
        LAI  = LAI+LAI2

    LAI = (pd.DataFrame(LAI)).replace(to_replace=0,value=np.nan)/365
    Lat = np.array(Data['latitude'])
    Lat = Lat[1250:1360]
    Lon = np.array(Data['longitude'])
    Lon = Lon[5060:5185]

    Lon,Lat = np.meshgrid(Lon,Lat)
    S_1 = pd.DataFrame(Lat.flatten())
    S_1.columns = ['Latitude']
    S_1['Longitude'] = Lon.flatten()
    S_1['LAI']  = np.array(LAI).flatten() 
    
    plt.figure(figsize=(10,4))
    plt.contourf(Lon,Lat,LAI,cmap='gist_ncar')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    cbar = plt.colorbar()
    cbar.set_label('LAI')
    plt.title('LAI in 50 m resolution in the Chambal Catchment')
    
    return S_1