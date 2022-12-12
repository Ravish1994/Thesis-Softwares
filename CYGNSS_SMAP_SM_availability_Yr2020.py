import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error

'''
   Input : path of all csv files saved previously
           n: numbers of files 
           
   Output : Visualizing Correlation with band plots of standard deviation

''' 

def Plotting_Variations(path,n):
    lat = [21,22,23,24,25,26,27,28,29,30,31]
    lon = [73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89]
    Lat = []
    Lon = []
    for i in range(len(lat)):
        lat1 = lat[i]
        for j in range(len(lon)):
            lon1 = lon[j]
            Lat.append(lat1)
            Lon.append(lon1)
    Lon = (np.array(Lon)).flatten()
    Lat = (np.array(Lat)).flatten()
    
    Path = []
    for i in range(len(Lat)):
        lat = str(Lat[i])
        lon = str(Lon[i])
        File_Name = lat+lon
        Path1 = path + '\Annual_Variations_' + lat + '_' + lon + '.csv'
        Path.append(Path1)
        
    Path = Path[0:n]
    for Path in Path:
        Path_1 = Path
        Lat_Lon = Path_1[-9:-4]
        lat = (Lat_Lon[0:2])
        lon = (Lat_Lon[3:])
        
        df          = pd.read_csv(Path_1)
        Day_No      = df['Day_No']
        
        CYGNSS_Data = df['Cygnss_SR']
        ## Applying change detection model on SR to retrieve Soil Moisture
        Min         = 0      ## SR corresponding to the driest condition
        Max         = 20     ## SR corresponding to the wettest condition
  
        SmN         = (CYGNSS_Data-Min)/(Max-Min)         ## Finding SM
       
        SMAP_Data1  = df['SMAP_SM']
        std1        = df['Cygnss_SD']/20
        
        df2              = pd.DataFrame(SmN)
        df2.columns      = ['CYGNSS_SM']
        df2['Day_No']    = pd.DataFrame(df['Day_No'])
        SMAP_Data1       = pd.DataFrame(df['SMAP_SM'])
        df2['SMAP_SM']   = pd.DataFrame(SMAP_Data1/100)
        df2['STD_SM']    = pd.DataFrame(df['Cygnss_SD']/20)
        df2['CYGNSS_SM'] = pd.DataFrame(SmN)
        
        df2            = df2.dropna()
        CR             = (np.array(df2.corr()))*100
        CR1            = np.round(CR[2][0],2)
        
        act   = pd.DataFrame(SMAP_Data1).replace(to_replace=np.nan,value=0)
        pred  = pd.DataFrame(CYGNSS_Data).replace(to_replace=np.nan,value=0)
        rmse  = mean_squared_error(act/100, pred/20, squared=False)
        
        a2 = np.array(act)/100
        a1 = np.array(pred)
        a1 = (a1-np.min(a1))/(np.max(a1)-np.min(a1))
        denominator = np.sum((a1 - np.mean(a1))**2)
        numerator   = np.sum((a2 - a1)**2)
        nse_val     = 1 - (numerator/denominator)
        
        Day_No = df2['Day_No']
        plt.figure(figsize=(30,10))
        plt.plot(Day_No,df2['CYGNSS_SM'],label  = 'CYGNSS Soil Moisture',color='blue')
        plt.fill_between(Day_No,df2['CYGNSS_SM']-df2['STD_SM'],df2['CYGNSS_SM']+df2['STD_SM'],color='grey')
        plt.scatter(Day_No,df2['SMAP_SM'],label = 'SMAP Soil Moisture',color='black')
        plt.ylabel('Soil Moisture',fontsize=20)
        plt.xlabel('Days of Year 2020',fontsize=20)
        plt.title(f'Latitude:{lat} Longitude:{lon} Correlation:{CR1}% RMSE:{round(rmse,4)} NSE:{np.round(nse_val,2)}',fontsize=20)
        plt.legend(fontsize=30) 

        plt.figure(figsize=(4,4))
        plt.scatter(df2['SMAP_SM'],df2['CYGNSS_SM'],s=20)
        plt.ylabel('CYGNSS Derived Soil Moisture')
        plt.xlabel('SMAP Soil Moisture')
        plt.plot([0,1],[0,1],c='gray')
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.xticks(np.arange(0, 1, 0.2))
        plt.yticks(np.arange(0, 1, 0.2))
        
        import seaborn as sns
        plt.figure(figsize=(2,2))
        sns.lmplot(x='SMAP_SM',y='CYGNSS_SM',data=df2,line_kws={'color': 'black'})
        plt.xlabel('SMAP Soil Moisture')
        plt.ylabel('CYGNSS Derived Soil Moisture')
        plt.plot([0,1],[0,1],c='gray')
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.xticks(np.arange(0, 1, 0.2))
        plt.yticks(np.arange(0, 1, 0.2))