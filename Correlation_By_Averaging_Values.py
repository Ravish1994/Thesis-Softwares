import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from termcolor import colored
import seaborn as sns

'''-----------------------------------------Sub Functions for Taking Average-----------------------------------------'''
'''
   Input: Lat : Latitude centroid of the SMAP grid cell
          Lon : Longitude centroid of the SMAP grid cell
          DF  : Daily DataFrame of CYGNSS
          Var : Variable to be masked 'Surface Reflectivity' 
          
   Output: Average value of Surface Reflectivity within SMAP 36 Km grid cell'''

def SR_gridded_Avg(Lat,Lon,DF,Var):
    mask = ((DF['sp_lat']>(Lat-0.18)) & (DF['sp_lat']<(Lat+0.18))) & ((DF['sp_lon']>(Lon-0.18)) & (DF['sp_lon']<(Lon+0.18)))
    Df1 = DF[mask]  
    SR  = np.array(Df1[f'{Var}'])
    SR1 = np.mean(SR)
    if SR1!=np.nan:
        SR2 = SR1
    else:
        SR2 = 0
    return SR2

'''
   Input: Lat : Latitude centroid of the SMAP grid cell
          Lon : Longitude centroid of the SMAP grid cell
          DF  : Daily DataFrame of SMAP
          Var : Variable to be masked 'SMAP_SM' 
          
   Output: Single value of Soil Moisture within SMAP 36 Km grid cell'''

def SM_gridded_Avg(Lat,Lon,DF,Var):
    mask = ((DF['SMAP_Latitude']>(Lat-0.18)) & (DF['SMAP_Latitude']<(Lat+0.18))) & ((DF['SMAP_Longitude']>(Lon-0.18)) & (DF['SMAP_Longitude']<(Lon+0.18)))
    Df1 = DF[mask]  
    SM  = np.array(Df1[f'{Var}'])
    
    ## Mean of single value will be that value only
    SM1 = np.mean(SM)
    if SM1!=np.nan:
        SM2 = SM1
    else:
        SM2 = 0
    return SM2

''' ----------------Averaging SMAP Soil Moisture and CYGNSS Surface reflectivity in 0.36 x 0.36 degree grid--------------------------'''

def Averaging_Gridd(DF1,DF2):
    Lat1 = np.arange(21,31,0.36)
    Lon1 = np.arange(73,89,0.36)
    
    SP_lat = []  # Latitude centroids of SMAP
    SP_lon = []  # Longitude centroids of SMAP
    SP_SR  = []  # Surface Reflectivity values within SMAP grid cell
    SP_SM  = []  # Soil Moisture values within SMAP grid cell
    for i in range(len(Lat1)):
        Lat = Lat1[i]
        for j in range(len(Lon1)):
            Lon = Lon1[j]
            SR  = SR_gridded_Avg(Lat,Lon,DF1,'SR_eff2')
            SM  = SM_gridded_Avg(Lat,Lon,DF2,'SMAP_SM')
            
            SP_lat.append(Lat)
            SP_lon.append(Lon)
            SP_SR.append(SR)
            SP_SM.append(SM)
            
    DF = pd.DataFrame(SP_lat)
    DF.columns = ['SP_lat']
    DF['SP_lon'] = SP_lon
    DF['SP_SR'] = SP_SR
    DF['SP_SM'] = SP_SM
    return DF,len(Lat1),len(Lon1)


'''-------------------------------------------------Correlating the Averaged values--------------------------------------------'''
def Correlating_CYGNSS_SR_SMAP_SM(path1,path2,d):
    
    DF1 = pd.read_csv(path1)
    DF2 = pd.read_csv(path2)
    
    # Since SMAP have 36Km x 36Km EASE gridded data
    DF,m,n  = Averaging_Gridd(DF1,DF2)
    
    # Finding correlation 
    DF1  = DF.drop(['SP_lat','SP_lon'],axis=1)
    Df_N = DF1.dropna()
    Corr = np.array(Df_N.corr())
    CR   = round(Corr[0][1]*100,2)
    
    X = np.array(DF['SP_lat'])
    X = X.reshape(m,n)
    Y = np.array(DF['SP_lon'])
    Y = Y.reshape(m,n) 
    Z = np.array(DF['SP_SR'])
    Z = Z.reshape(m,n)
    W = np.array(DF['SP_SM'])
    W = W.reshape(m,n)
    
    plt.figure(figsize = (20,5))
    plt.subplot(1,2,1)
    plt.contourf(X,Y,Z,cmap = 'gist_ncar')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    cbar = plt.colorbar()
    cbar.set_label('Surface Reflectivity in dB')
    plt.title(f'{d} : Correlation = {CR} %')
    
    plt.subplot(1,2,2)
    plt.contourf(X,Y,W,cmap = 'gist_ncar')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    cbar = plt.colorbar()
    cbar.set_label('Soil Moisture in cm3/cm3')
    plt.title(f'{d} : Correlation = {CR} %')
    
    Points = []
    for i in range(len(X.flatten())):
        Points.append(i)
        
    plt.figure(figsize=(40,15))
    plt.plot(Points[:1000],W.flatten()[:1000],label = 'SMAP Soil Moisture')
    plt.plot(Points[:1000],Z.flatten()[:1000]/20,label    = 'Normalized CYGNSS Derived Surface Reflectivity')
    plt.ylabel('''Surface reflectiviy and Soil Moisture''',fontsize=20)
    plt.xlabel('Points',fontsize=20)
    plt.legend(fontsize=20) 
    
    sns.lmplot(x='SP_SR',y='SP_SM',data=DF,aspect=4,height=6)
    plt.xlabel('Normalized CYGNSS Derived Surface Reflectivity')
    plt.ylabel('SMAP Soil Moisture in cm3/cm3')
    plt.title(f'''Fitting a Regression Line using lmplot''')
    

'''
Input : Day Range intended for correlation
Output: Correlation Coefficient Daywise'''

def Corr_SM_SMAP_CYGNSS(m,n):
    for i in range(m,n+1):
        if i<=9:
            path1 = f'D:\EG\Project Data\CYGNSS_features\CYGNSS_Features_Day_00{i}.csv'
            path2 = f'D:\EG\Project Data\CYGNSS_Featured_Data_Along_SMAP_Coordinate\SM_CYGNSS_Features_Day_00{i}.csv'
            d = f'Day_00{i}'
            Correlating_CYGNSS_SR_SMAP_SM(path1,path2,d)
            
        elif i<=99:
            path1 = f'D:\EG\Project Data\CYGNSS_features\CYGNSS_Features_Day_0{i}.csv'
            path2 = f'D:\EG\Project Data\CYGNSS_Featured_Data_Along_SMAP_Coordinate\SM_CYGNSS_Features_Day_0{i}.csv'
            d = f'Day_0{i}'
            Correlating_CYGNSS_SR_SMAP_SM(path1,path2,d)
            
        else:
            path1 = f'D:\EG\Project Data\CYGNSS_features\CYGNSS_Features_Day_{i}.csv'
            path2 = f'D:\EG\Project Data\CYGNSS_Featured_Data_Along_SMAP_Coordinate\SM_CYGNSS_Features_Day_{i}.csv'
            d = f'Day_{i}'
            Correlating_CYGNSS_SR_SMAP_SM(path1,path2,d)
                        
'''---------------------------------------------------------------END-----------------------------------------------------------------'''