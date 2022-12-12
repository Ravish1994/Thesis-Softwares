import numpy as np
import pandas as pd
from pandas import DataFrame

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
    SR_counts = len(SR)
    return SR_counts

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

def CYGNSS_Data_availability(m,n,lat,lon):
    Day_Number = []
    SR_counts  = []
    for i in range(m,n+1):
        Day_Number.append(i)
        if i<=9:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_00{i}.csv'
            SR1 = SR_1(path1,lat,lon)
            SR_counts.append(SR1)
        elif i<=99:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_0{i}.csv'
            SR1 = SR_1(path1,lat,lon)
            SR_counts.append(SR1)
        elif i<=366:
            path1 = f'D:\EG\Project Data\CYGNSS_CSV_Files\Day_{i}.csv'
            SR1 = SR_1(path1,lat,lon)
            SR_counts.append(SR1)
    Day_Number = pd.DataFrame(Day_Number)
    Day_Number.columns = ['Days']
    Day_Number['Number of CyGNSS Observations'] = SR_counts
    Day_Number.to_csv(rf'D:\EG\Project Data\CYGNSS_Data_in_0p36Dg\CYGNSS_ObsCounts_36KmGrid\Observations_Counts_of_CYGNSS_within_36km_Gridcell_{lat}_{lon}.csv',index = False)
    
def CYGNSS_Data_availability_plotting(lat,lon):
    import numpy as np
    import pandas as pd
    from pandas import DataFrame
    import matplotlib.pyplot as plt
    
    Day_Number = pd.read_csv(f'D:\EG\Project Data\CYGNSS_Data_in_0p36Dg\CYGNSS_ObsCounts_36KmGrid\Observations_Counts_of_CYGNSS_within_36km_Gridcell_{lat}_{lon}.csv')
    plt.figure(figsize=(30,10))
    plt.bar(Day_Number['Days'],Day_Number['Number of CyGNSS Observations'])
    N = len(Day_Number[Day_Number['Number of CyGNSS Observations']>0])
    plt.scatter([1],[0],label=f'''Number of Days Having CYGNSS Data Availability within 36 Km 
    grid having centroid as Lat {lat} Lon {lon}: {N}''')
    plt.legend()
    plt.ylabel('''Number of CyGNSS Observations''',fontsize=20)
    plt.xlabel('Days of Year 2020',fontsize=20)
    plt.title(f'Number Of CyGNSS observations within 36 Km grid having centroid_Lat:{lat}_Lon:{lon}',fontsize=20)
    
    
def CYGNSS_Data_availability_plotting2(lat,lon):
    import numpy as np
    import pandas as pd
    from pandas import DataFrame
    import matplotlib.pyplot as plt
         
    Counts_restricted = []
    for i in range(1,367):
        Counts_restricted.append(4)
        
    Day_Number = pd.read_csv(f'D:\EG\Project Data\CYGNSS_Data_in_0p36Dg\CYGNSS_ObsCounts_36KmGrid\Observations_Counts_of_CYGNSS_within_36km_Gridcell_{lat}_{lon}.csv')
    
    N = len(Day_Number[Day_Number['Number of CyGNSS Observations']>0])
    plt.figure(figsize=(30,10))
    plt.bar(Day_Number['Days'],Day_Number['Number of CyGNSS Observations'])
    plt.scatter([0],[0],color='white',s=5,label=f'''Number of Days Having CYGNSS Data Availability within 36 Km 
    grid having centroid as Lat {lat} Lon {lon}: {N}''')
    plt.legend(fontsize=30)
    plt.plot(Day_Number['Days'],Counts_restricted)
    plt.ylabel('''Number of CyGNSS Observations''',fontsize=20)
    plt.xlabel('Days of Year 2020',fontsize=20)
    plt.title(f'Number Of CyGNSS observations within 36 Km grid having centroid_Lat:{lat}_Lon:{lon}',fontsize=20)
    
def CYGNSS_Counts(path,lat,lon):
    Data = pd.read_csv(path)
    SR = Data['SR_eff']
    le = len(SR)
    if le>0:
        l = le
    else:
        l = 0
    return l

def CYGNSS_Data_availability_plotting3(m,n,lat,lon):
    path  = f'D:\EG\Project Data\CYGNSS_Obs_Chambal_{lat}_{lon}\Chambal_{lat}_{lon}'
    import numpy as np
    import pandas as pd
    from pandas import DataFrame
    import matplotlib.pyplot as plt
    Day_Number = []
    SR_counts  = []
    for i in range(m,n+1):
        Day_Number.append(i)
        if i<=9:
            path1 = path + f'\CYGNSS_Data_Availability_Day_00{i}_{lat}_{lon}.csv'
            SR1 = CYGNSS_Counts(path1,lat,lon)
            SR_counts.append(SR1)
        elif i<=99:
            path1 = path + f'\CYGNSS_Data_Availability_Day_0{i}_{lat}_{lon}.csv'
            SR1 = CYGNSS_Counts(path1,lat,lon)
            SR_counts.append(SR1)
        elif i<=366:
            path1 = path + f'\CYGNSS_Data_Availability_Day_{i}_{lat}_{lon}.csv'
            SR1 = CYGNSS_Counts(path1,lat,lon)
            SR_counts.append(SR1)
            
    Counts_restricted = []
    for i in range(1,367):
        Counts_restricted.append(4)
        
    Day_Number = pd.DataFrame(Day_Number)
    Day_Number.columns = ['Days']
    Day_Number['Number of CyGNSS Observations'] = SR_counts
    N = len(Day_Number[Day_Number['Number of CyGNSS Observations']>0])
    plt.figure(figsize=(30,10))
    plt.bar(Day_Number['Days'],Day_Number['Number of CyGNSS Observations'])
    plt.scatter([1],[0],label=f'''Number of Days Having CYGNSS Data Availability within 36 Km 
    grid having centroid as Lat {lat} Lon {lon}: {N}''')
    plt.legend()
    plt.plot(Day_Number['Days'],Counts_restricted)
    plt.ylabel('''Number of CyGNSS Observations''',fontsize=20)
    plt.xlabel('Days of Year 2020',fontsize=20)
    plt.title(f'Number Of CyGNSS observations within 36 Km grid having centroid_Lat:{lat}_Lon:{lon}',fontsize=20)

   