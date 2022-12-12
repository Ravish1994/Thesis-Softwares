import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

import geopandas as gpd     
cmap = 'gist_ncar'
fp = r'indian_districts.shp'
map_df = gpd.read_file(fp) 
map_df_msk = map_df[((map_df['latitude']>22) & (map_df['latitude']<27)) & ((map_df['longitude']>72) & (map_df['longitude']<79))]
DF3 = map_df_msk

## Spatial averaging to get data on similar grid
## Averaging function to average the LAI and AGB in the SMAP grid cell of 36 x 36 Km
def Avg(lat1,lon1,Df,Var):
    mask = ((Df['Latitude']>(lat1-0.025)) & (Df['Latitude']<(lat1+0.025)) & (Df['Longitude']>(lon1-0.025)) & (Df['Longitude']<(lon1+0.025)))
    Df = Df[mask]
    
    T_Var = Df[f'{Var}']
    le = len(T_Var)
    if le>0:
        T_Var_avg = np.mean(T_Var)
    else:
        T_Var_avg = np.nan
    return T_Var_avg

def Resample_LAI_Biomass(DS_LAI,DS_BM):
    ## Averaging operation is done in the Chambal catchment only
    ## Creating Lat and Lon of the Centroid of 6 Km grid cells
    lat = np.arange(22.025,28,0.05)
    lon = np.arange(73.025,80,0.05)

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
    return DF,Lat_LAI,Lon_LAI,LAI_36,BM_36
    
def Correlate_LAI_BM(DF):
    DF2   = DF.copy()
    DF1   = DF2

    DF_F = DF1.drop(['Latitude','Longitude'],axis = 1)
    DF_F = DF_F[DF_F['BM']>55]
    CR = (np.array(DF_F.corr()))*100
    CR1 = np.round(CR[0][1],2) 
    
    DF = DF.replace(to_replace = np.nan,value = 0)
    lat = np.arange(22.025,27,0.05)
    lon = np.arange(73.025,79,0.05)
    Lat,Lon = np.meshgrid(lat,lon)
    
    Latitude  = DF['Latitude']
    Longitude = DF['Longitude']
    BM        = DF['BM']
    LAI       = DF['LAI']
    
    BM_mesh   = griddata((Latitude,Longitude),BM,(Lat,Lon), method = 'nearest')
    LAI_mesh  = griddata((Latitude,Longitude),LAI,(Lat,Lon), method = 'nearest')
    
    # Creating contour colour map of Biomass and LAI data on 5 Km grid
    fig , ax2 = plt.subplots(figsize=(20, 10))
    DF3.plot(color = 'white', edgecolor = 'black',axes=ax2)
    plt.pcolor(Lon,Lat,BM_mesh,cmap="YlGnBu")
    DF3.plot(color = "white", edgecolor = 'black',axes=ax2,alpha=0.1)
    cbar = plt.colorbar()
    cbar.set_label('Mg/ha')
    plt.xlabel('Longitude',fontsize=25)
    plt.ylabel('Lattitude',fontsize=25)
    plt.title(f'Biomass Data of 5 km spatial resolution',fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=15)
    
    fig , ax2 = plt.subplots(figsize=(20, 10))
    DF3.plot(color = 'white', edgecolor = 'black',axes=ax2)
    plt.pcolor(Lon,Lat,LAI_mesh,cmap="YlGnBu")
    DF3.plot(color = "white", edgecolor = 'black',axes=ax2,alpha=0.1)
    cbar = plt.colorbar()
    cbar.set_label('m2/m2')
    plt.xlabel('Longitude',fontsize=25)
    plt.ylabel('Lattitude',fontsize=25)
    plt.title(f'Correlation with Biomass: {CR1}% LAI Data of 5 km grid cell',fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=15)
    
def Spatial_Visualization_After_Normalization(DF):
    DF  = DF.replace(to_replace = np.nan,value = 0)
    DF['LAI'] = (DF['LAI']-np.min(np.array(DF['LAI'])))/(np.max(np.array(DF['LAI']))-np.min(np.array(DF['LAI'])))
    DF['BM']  = (DF['BM']-np.min(np.array(DF['BM'])))/(np.max(np.array(DF['BM']))-np.min(np.array(DF['BM'])))
    lat = np.arange(22.025,27,0.05)
    lon = np.arange(73.025,79,0.05)
    Lat,Lon = np.meshgrid(lat,lon)
    
    Latitude  = DF['Latitude']
    Longitude = DF['Longitude']
    BM        = DF['BM']
    LAI       = DF['LAI']
    
    BM_mesh   = griddata((Latitude,Longitude),BM,(Lat,Lon), method = 'nearest')
    LAI_mesh  = griddata((Latitude,Longitude),LAI,(Lat,Lon), method = 'nearest')
    
    # Creating contour colour map of Biomass and LAI data on 5 Km grid
    fig , ax2 = plt.subplots(figsize=(20, 10))
    DF3.plot(color = 'white', edgecolor = 'black',axes=ax2)
    plt.pcolor(Lon,Lat,BM_mesh,cmap="YlGnBu")
    DF3.plot(color = "white", edgecolor = 'black',axes=ax2,alpha=0.1)
    cbar = plt.colorbar()
    cbar.set_label('Mg/ha')
    plt.xlabel('Longitude',fontsize=25)
    plt.ylabel('Lattitude',fontsize=25)
    plt.title(f'Biomass Data of 5 km spatial resolution',fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=15)
    
    fig , ax2 = plt.subplots(figsize=(20, 10))
    DF3.plot(color = 'white', edgecolor = 'black',axes=ax2)
    plt.pcolor(Lon,Lat,LAI_mesh,cmap="YlGnBu")
    DF3.plot(color = "white", edgecolor = 'black',axes=ax2,alpha=0.1)
    cbar = plt.colorbar()
    cbar.set_label('m2/m2')
    plt.xlabel('Longitude',fontsize=25)
    plt.ylabel('Lattitude',fontsize=25)
    plt.title(f'LAI Data of 5 km grid cell',fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=15)
    
def Spatial_Visualization_LAI_BM(DF):
    DF  = DF.replace(to_replace = np.nan,value = 0)
    DF  = DF[DF['LAI']>0]
    
    lat = np.arange(22.025,27,0.05)
    lon = np.arange(73.025,79,0.05)
    Lat,Lon = np.meshgrid(lat,lon)
    
    Latitude  = DF['Latitude']
    Longitude = DF['Longitude']
    BM        = DF['BM']
    LAI       = DF['Calculated Biomass']
    LAI1      = DF['LAI']
    
    BM_mesh   = griddata((Latitude,Longitude),BM,(Lat,Lon), method = 'nearest')
    LAI_mesh  = griddata((Latitude,Longitude),LAI,(Lat,Lon), method = 'nearest')
    LAI_mesh1 = griddata((Latitude,Longitude),LAI1,(Lat,Lon), method = 'nearest')
    
    # Creating contour colour map of Biomass and LAI data on 5 Km grid
    # Creating contour colour map of Biomass and LAI data on 5 Km grid
    fig , ax2 = plt.subplots(figsize=(20, 10))
    DF3.plot(color = 'white', edgecolor = 'black',axes=ax2)
    plt.pcolor(Lon,Lat,BM_mesh,cmap="YlGnBu")
    DF3.plot(color = "white", edgecolor = 'black',axes=ax2,alpha=0.1)
    cbar = plt.colorbar()
    cbar.set_label('Kg/m2')
    plt.xlabel('Longitude',fontsize=25)
    plt.ylabel('Lattitude',fontsize=25)
    plt.title(f'Actual Biomass from ESA CCI',fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=15)
    
    fig , ax2 = plt.subplots(figsize=(20, 10))
    DF3.plot(color = 'white', edgecolor = 'black',axes=ax2)
    plt.pcolor(Lon,Lat,LAI_mesh,cmap="YlGnBu")
    DF3.plot(color = "white", edgecolor = 'black',axes=ax2,alpha=0.1)
    cbar = plt.colorbar()
    cbar.set_label('Kg/m2')
    plt.xlabel('Longitude',fontsize=25)
    plt.ylabel('Lattitude',fontsize=25)
    plt.title(f'Calculated Biomass from LAI',fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=15)
    
    fig , ax2 = plt.subplots(figsize=(20, 10))
    DF3.plot(color = 'white', edgecolor = 'black',axes=ax2)
    plt.pcolor(Lon,Lat,LAI_mesh1,cmap="YlGnBu")
    DF3.plot(color = "white", edgecolor = 'black',axes=ax2,alpha=0.1)
    cbar = plt.colorbar()
    cbar.set_label('m2/m2')
    plt.xlabel('Longitude',fontsize=25)
    plt.ylabel('Lattitude',fontsize=25)
    plt.title(f'LAI from NOAA',fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=15)