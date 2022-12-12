from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import griddata

import geopandas as gpd     
cmap = 'gist_ncar'
fp = r'indian_districts.shp'
map_df = gpd.read_file(fp) 
map_df_msk = map_df[((map_df['latitude']>22) & (map_df['latitude']<27)) & ((map_df['longitude']>72) & (map_df['longitude']<79))]
DF3 = map_df_msk

def Biomass_Data_In_Chambal(Big_Data_File_Path):
    ## Slicing the Big Biomass data inside the Ganga Catchment
    Biomass_BigData = Dataset(Big_Data_File_Path, 'r')
    AGB = Biomass_BigData.variables['agb']
    Lat = Biomass_BigData.variables['lat']
    Lon = Biomass_BigData.variables['lon']

    # Lat = np.array(Lat[55100:66375])
    # Lon = np.array(Lon[284630:302630])
    # sliced_agb = np.array(AGB[0,55100:66375,284630:302630])

    ### Extracting Data inside Chambal catchment small part in Ganga Catchment
    # The Chambal Basin lies between latitudes 22°27'N and 27°20'N and longitudes 73°20'E and 79°15'E
    Lat = np.array(Lat[59000:65250])
    Lon = np.array(Lon[284630:291730])
    sliced_agb = np.array(AGB[0,59000:65250,284630:291730])

    ## Creating meshgrid of latitude and longitude of spatial resolution of 100m
    Lon,Lat = np.meshgrid(Lon,Lat)

    # Creating contour colour map of Biomass data 
    fig , ax2 = plt.subplots(figsize=(20, 10))
    DF3.plot(color = 'white', edgecolor = 'black',axes=ax2)
    plt.pcolor(Lon,Lat,sliced_agb,cmap="YlGnBu")
    DF3.plot(color = None, edgecolor = 'black',axes=ax2,alpha=0.1)
    cbar = plt.colorbar()
    cbar.set_label(f'Mg/ha')
    plt.xlabel('Longitude',fontsize=25)
    plt.ylabel('Lattitude',fontsize=25)
    plt.title(f'Biomass Data of 100 m spatial resolution in the Chambal Catchment',fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=15)
    
    Df_Biomass = pd.DataFrame(Lat.flatten())
    Df_Biomass.columns = ['Latitude']
    Df_Biomass['Longitude'] = Lon.flatten()
    Df_Biomass['Biomass_Mg_per_ha'] = sliced_agb.flatten()
    lat = Df_Biomass['Latitude']
    lon = Df_Biomass['Longitude']
    AGB = Df_Biomass['Biomass_Mg_per_ha']
    
    ## Saving chunked file in the Chambal catchment
#     AGB_sorted_Df = pd.DataFrame(lat)
#     AGB_sorted_Df.columns = ['Latitude']
#     AGB_sorted_Df['Longitude'] = lon
#     AGB_sorted_Df['AGB_MgpHa'] = AGB
#     Mask = (AGB_sorted_Df['AGB_MgpHa']>0)
#     AGB_sorted_Df = AGB_sorted_Df[Mask]
#     AGB_sorted_Df.to_csv(f'Biomass_Data.csv',index = False)