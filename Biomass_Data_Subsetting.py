from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import griddata

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
    plt.figure(figsize=(10,4))
    plt.contourf(Lon,Lat,sliced_agb,cmap='gist_ncar')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    cbar = plt.colorbar()
    cbar.set_label('Mg/ha')
    plt.title('Biomass Data of 100 m spatial resolution in the Chambal Catchment')

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