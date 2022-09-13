import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt 
import seaborn as sns
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

'''
   Input : Path : Path of the file containing CYGNSS SR and GLDAS SM
           lat and lon: Centroid of latitude and longitude of the GLDAS single grid cell
           
   Output : Visualization of correlation

'''
        
def SR_GLDAS_SM_Corr(path,lat,lon):
    Dataset = pd.read_csv(path)
    Dataset = Dataset.dropna()
    
    y_pred = np.array((Dataset.drop(['Day_No','GLDAS_SM'],axis=1))/20)     # Independet variable
    y_act  = np.array((Dataset['GLDAS_SM'])/1000)                          # dependent variable
    
    # Calculation of RMSE
    act   = np.array(pd.DataFrame(y_act).replace(to_replace=np.nan,value=0))
    pred  = np.array(pd.DataFrame(y_pred).replace(to_replace=np.nan,value=0))
    
    if len(act)>0:
        rmse  = mean_squared_error(act, pred, squared=False)
    else:
        rmse  = np.nan
    
    # R_square
    if len(act)>0:
        R2_score = r2_score(act,pred)
    else:
        R2_score = np.nan
        
    Dataset['Cygnss_SR'] = Dataset['Cygnss_SR']/20
    Dataset['GLDAS_SM']  = Dataset['GLDAS_SM']/1000
    df1 = (Dataset.drop(['Day_No'],axis=1))
    df1  = df1.dropna()        
    CR = (np.array(df1.corr()))*100
    CR1 = np.round(CR[0][1],2)
    
    plt.figure(figsize=(30,8))
    plt.scatter(Dataset['Day_No'],Dataset['Cygnss_SR'],label = 'CYGNSS Backscatter Normalized between 0 and 1')
    plt.scatter(Dataset['Day_No'],Dataset['GLDAS_SM'],label = 'GLDAS Soil Moisture values between 0 and 1')
    plt.ylabel('''Soil Moisture and Normalized Backscatter''',fontsize=15)
    plt.xlabel('Days of Year 2020',fontsize=15)
    plt.title(f'Latitude:{lat} Longitude:{lon} Correlation:{CR1}%  RMSE : {round(rmse,2)}',fontsize=15)
    plt.legend(fontsize=15)
    
    sns.lmplot(x='Cygnss_SR',y='GLDAS_SM',data=Dataset,aspect=4,height=6)
    plt.xlabel('Cygnss Backscatter')
    plt.ylabel('GLDAS_SM')
    plt.title(f'''Fitting a Regression Line using lmplot''')