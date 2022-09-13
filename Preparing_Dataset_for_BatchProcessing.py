import numpy as np
import sympy
from sympy import *
import pandas as pd
import matplotlib.pyplot as plt

## Creating day number to read CYGNSS files
def Day(i):
    if i<= 9:
        b = '00'+str(i)
        return b
    elif i>=10 and i<=99:
        b = '0'+str(i)
        return b
    elif i>=100 and i<=366:
        b = str(i)
        return b 
    
def Data_Batch_2020(lat,lon):
    C1,C2,C3,C4,C6,C7 = [],[],[],[],[],[]
    Df  = pd.DataFrame(np.array(C1))
    Df.columns = ['SR_eff']
    Df['SP_I'] = C2
    Df['sp_lat'] = C3
    Df['sp_lon'] = C4
    Df['LAI'] = C6
    Df['Day_No'] = C7

    for i in range(1,366):
        p1 = f'D:\EG\Project Data\CYGNSS_Obs_Chambal_{lat}_{lon}\LAI_{lat}_{lon}'
        p2 = f'\CYGNSS_Data_Availability_Day_{Day(i)}_{lat}_{lon}.csv'
        path1  = p1+p2
        Data1 = pd.read_csv(path1)
        Data1['Day_No'] = np.array([i for row in range(len(Data1))])

        # Stack the DataFrames on top of each other
        Stacked_Data = pd.concat([Df,Data1], axis=0)
        Df = Stacked_Data
    mask  = Df['LAI']>0
    Df1   = Df[mask]
    return Df1