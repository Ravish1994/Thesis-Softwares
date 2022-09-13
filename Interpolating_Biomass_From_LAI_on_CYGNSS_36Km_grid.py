import numpy as np
import pandas as pd
from scipy.interpolate import griddata

'''
   Interpolating Biomass from daily LAI of 2020 
   On the CYGNSS satellite observation points 
   Within every single SMAP grid cell for 366 days  '''


# Interpolation of the Biomass data on the CYGNSS daily data points within 36 km grid
def fun1(i):
    if i<= 9:
        b = '00'+str(i)
        return b
    elif i>=10 and i<=99:
        b = '0'+str(i)
        return b
    elif i>=100 and i<=366:
        b = str(i)
        return b 

### Observations of CYGNSS with Biomass for a particular pixel 

def CYGNSS_LAI(lat,lon):
    for i in range(1,366):
        ## CYGNSS file Path
        day = fun1(i)
        p1 = f'D:\\EG\\Project Data\\CYGNSS_Obs_Chambal_{lat}_{lon}' 
        p2 = f'\Chambal_{lat}_{lon}\\CYGNSS_Data_Availability_Day_{day}_{str(lat)}_{str(lon)}.csv'
        cygpath = p1+p2
        CYGNSS = pd.read_csv(cygpath)
        
        ## LAI file path
        path = f'D:\EG\Project Data\LAI_2018\LAI_{lat}_{lon}'
        filename = f'\LAI_Data_Within_36_Km_Resolution_Cell_Day_{i}_{lat}_{lon}.csv'
        Path = path+filename
        LAI_DS = pd.read_csv(Path)
        
        Lat_LAI = LAI_DS['Latitude']
        Lon_LAI = LAI_DS['Longitude']
        LAI_DS1 = LAI_DS['LAI']
        
        Lat_CYGNSS = CYGNSS['sp_lat']
        Lon_CYGNSS = CYGNSS['sp_lon']
        LAI_CYGNSS = griddata((Lat_LAI,Lon_LAI),LAI_DS1,(Lat_CYGNSS,Lon_CYGNSS), method = 'nearest')
        CYGNSS['LAI'] = np.array(LAI_CYGNSS)
        CYGNSS = CYGNSS.replace(to_replace=np.nan,value=0)
        CYGNSS.to_csv(f'CYGNSS_Data_Availability_Day_{day}_{lat}_{lon}.csv',index=False)