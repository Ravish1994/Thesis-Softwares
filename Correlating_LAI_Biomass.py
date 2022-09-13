import numpy as np
import sympy
from sympy import *
import pandas as pd
import matplotlib.pyplot as plt

## Spatial averaging to get data on similar grid
## Averaging function to average the LAI and AGB in the SMAP grid cell of 36 x 36 Km
def Avg(lat1,lon1,Df,Var):
    mask = ((Df['Latitude']>(lat1-0.18)) & (Df['Latitude']<(lat1+0.18)) & (Df['Longitude']>(lon1-0.18)) & (Df['Longitude']<(lon1+0.18)))
    Df = Df[mask]
    
    T_Var = Df[f'{Var}']
    le = len(T_Var)
    if le>0:
        T_Var_avg = np.mean(T_Var)
    else:
        T_Var_avg = np.nan
    return T_Var_avg

def Correlate_LAI_Biomass(DS_LAI,DS_BM):
    ## Averaging operation is done in the Chambal catchment only
    ## Lat and Lon of the Centroid of 36 Km SMAP grid cells
    lat = [22.18,22.54,22.90,23.26,23.62,24.98,24.34,24.70,25.06,25.42,25.78,26.14,26.50,26.86,27.22]
    lon = [73.18,73.54,73.9,74.26,74.62,74.98,75.34,75.7,76.06,76.42,76.78,77.14,77.5,77.86,78.22,78.58,78.94,79.3,79.66]

    Lat_LAI = []
    Lon_LAI = []
    LAI_36  = []
    BM_36   = []
    for i in range(len(lat)):
        lat1 = lat[i]
        for j in range(len(lon)):
            lon1 = lon[j]
            LAI_avg_36 = Avg(lat1,lon1,DS_LAI,'LAI')
            BM_avg_36  = Avg(lat1,lon1,DS_BM,'AGB_MgpHa')
            LAI_36.append(LAI_avg_36)
            BM_36.append(BM_avg_36)
            Lat_LAI.append(lat1)
            Lon_LAI.append(lon1)

    DF = pd.DataFrame(Lat_LAI)
    DF.columns = ['Latitude']
    DF['Longitude'] = Lon_LAI
    DF['BM']  = BM_36
    DF['LAI'] = LAI_36
    DF = DF.dropna()

    DF_F = DF.drop(['Latitude','Longitude'],axis = 1)
    CR = (np.array(DF_F.corr()))*100
    CR1 = np.round(CR[0][1],2) 

    plt.figure(figsize=(20,5))
    plt.plot(DF_F['LAI'],label='Annual average LAI of 2018 in m2/m2')
    plt.plot(DF_F['BM'],label='Annual average LAI of 2018 Biomass in Mg/Hectare')
    plt.title(f'''Correlation: {np.abs(CR1)} %''') 
    plt.legend()
    
    return DF