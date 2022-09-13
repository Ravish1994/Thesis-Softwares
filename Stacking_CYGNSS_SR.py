'''----------------------------------------------- Important Libraries-----------------------------------------------------------------'''
import numpy as np
import pandas as pd
from pandas import DataFrame
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import seaborn as sns

'''------------------------------------------------Creating String for the Month-------------------------------------------------'''
'''
   Input  : Month Number
   Output : String of Month number
   
   Example, 
   Input  = 1, 11, 10
   Output = '01', '11', '10'

'''
def Month_No(i):
    if i <=9:
        Month_No = '0'+str(i)
    else:
        Month_No = str(i)
    return Month_No 

'''------------------Function to spatially correlate the cygnss surface reflectivity and SMAP soil moisture ---------------------------'''

'''Masking SR in between 0 and 20 dB'''
def Masking_Data_SR(SR):
    m,n = SR.shape
    for i in range(m):
        for j in range(n):
            a = SR[i][j]
            if (a<0) or (a>20):
                SR[i][j] = 0
            else:
                SR[i][j] = SR[i][j]
    return SR

'''
   Input  : Path of the file containing daily Surface Reflectivity and Soil Moisture values
            i : Day Number in the Year 2020 ranges from 001 to 366
   Output : Spatial correlation between both Surface Reflectivity and Soil Moisture values and their respective image plots

'''               
def Correl_SR(path,i):
    Data      = pd.read_csv(path)
    
    DS        = Data.copy()
    Mask      = (DS['SR_SMAP']>0) & (DS['SR_SMAP']<20)
    New_SR_SM = (DS[Mask]).dropna()                                     ## Dropping NAN values before making correlation
    Corr      = np.array(New_SR_SM.corr())                      ## Finding Correlation matrix using DF.corr() function
    CR        = np.round((Corr[0][1])*100,2)                    ## Accessing the correlation value stored in correlation matrix array  
    
    ## This Code is only for the visualization purpose of Soil moisture and Surface Reflectivity variation in space
    SM        = np.array(Data['SM_SMAP'])
    SR        = np.array(Data['SR_SMAP'])
    
    SM1       = SM.reshape(34,47)                               ## This is the size of the SMAP grid
    SR1       = SR.reshape(34,47) 
    SR2       = Masking_Data_SR(SR1)
    
    ## Showing the image of the reshaped Surface reflectivity
    plt.figure(figsize=(20,20))
    plt.subplot(1,2,1)
    cmap = 'gist_ncar'                                           ## Colourmap of the colorbar
    img  = plt.imshow(SR2,interpolation = 'gaussian',cmap=cmap)  
    cbar = plt.colorbar(img,shrink=0.27)                         ## Taking interpolation as gaussian to fill the values in neighbourhood
    cbar.set_label('Surface Reflectivity in dB')
    plt.title(f'''Surface Reflectivity of CYGNSS 
    interpolated along SMAP gridpoints''')  
    
    ## Showing the image of the reshaped Soil Moisture
    plt.subplot(1,2,2)
    cmap = 'gist_ncar'
    img = plt.imshow(SM1,interpolation = 'gaussian',cmap=cmap)   
    cbar = plt.colorbar(img,shrink=0.27)
    cbar.set_label('Soil Moisture in cm3/cm3')
    plt.title(f'''Soil Moisture of SMAP gridpoints
    Time: Day{i} 2020 Correlation: {CR} %''') 
    
    Points = []
    for i in range(len(SR2.flatten())):
        Points.append(i)
        
    plt.figure(figsize=(40,15))
    plt.plot(Points[:1000],SM1.flatten()[:1000],label = 'SMAP Soil Moisture')
    plt.plot(Points[:1000],SR2.flatten()[:1000]/20,label    = 'Normalized CYGNSS Derived Surface Reflectivity')
    plt.ylabel('''Surface reflectiviy and Soil Moisture''',fontsize=20)
    plt.xlabel('Points',fontsize=20)
    plt.legend(fontsize=20) 
    
    sns.lmplot(x='SR_SMAP',y='SM_SMAP',data=DS,aspect=4,height=6)
    plt.xlabel('Normalized CYGNSS Derived Surface Reflectivity')
    plt.ylabel('SMAP Soil Moisture in cm3/cm3')
    plt.title(f'''Fitting a Regression Line using lmplot''')
    
'''---------------------------------Function to Generate path and to call previous function---------------------------------------------'''
'''
   Input  : m: Initial number of day
            n: Last number of day 
   Output : Spatial correlation between both Surface Reflectivity and Soil Moisture values and their respective image plots

''' 
def Correlating_SM_SR(m,n):
    str1 = 'D:\EG\Project Data\Plotting_Data_SRSM'           ## User can change this path according to their files location
    f1   = '\Day_'
    f3   = '.csv'
    
    for i in range(m,n+1):
        if i<=9:
            f2   = '00' + str(i)                             ## Here f2 is generating the Day_number from 001 - 366
            path = str1 + f1 + f2 + f3
            Correl_SR(path,i)
        elif i<=99:
            f2   = '0' + str(i)
            path = str1 + f1 + f2 + f3
            Correl_SR(path,i)
        else:
            f2   = str(i)
            path = str1 + f1 + f2 + f3
            Correl_SR(path,i)
            
            
'''-------------------------------------------------End of Daily Spatial Correlation--------------------------------------------------'''

'''
   1. Since above big data gaps were found in daily SMAP and SR data 
   2. So in order to visualize the variations of data at every points inside the catchment 
   3. We need to stack the data of every single week '''

'''----------------------------------------------Stacking and Averaging of weekly data---------------------------------------------'''

'''
   Generating String of Day number using str2(i) function:
   Input : Day Number
   Ouput : String format of Day Number

'''
            
def str2(i):
    if i<=9:
        f2 = '00'+str(i)
    elif i<=99:
        f2 = '0'+str(i)
    else:
        f2 = str(i)
    return f2


'''
   Generating Weekly paths of the data:
   Input : Day Number
   Ouput : Weekly paths of the data according to the day number as starting day of the week

'''           
def CYGNSS_weekly(i):
    str1 = f'D:\EG\Project Data\Plotting_Data_SRSM\Day'+'_'
    str3 = '.csv'
    path1 = str1+str2(i)+str3
    path2 = str1+str2(i+1)+str3
    path3 = str1+str2(i+2)+str3
    path4 = str1+str2(i+3)+str3
    path5 = str1+str2(i+4)+str3
    path6 = str1+str2(i+5)+str3
    path7 = str1+str2(i+6)+str3
    return path1,path2,path3,path4,path5,path6,path7


'''
   Input :  Week Number in a Particular Month, Starting Day number of that particular Week in that particular month
   Output:  Generating weekly paths according to the week number
   
'''
def Weekly_Path(Week,m):
    if Week==1:
        i = m
        path1,path2,path3,path4,path5,path6,path7 = CYGNSS_weekly(i)        ## Paths of the 1st week of the month
    elif Week ==2:
        i = m+7
        path1,path2,path3,path4,path5,path6,path7 = CYGNSS_weekly(i)        ## Paths of the 2nd week of the month
    elif Week ==3:
        i = m+14
        path1,path2,path3,path4,path5,path6,path7 = CYGNSS_weekly(i)        ## Paths of the 3rd week of the month
    else:
        i = m+21
        path1,path2,path3,path4,path5,path6,path7 = CYGNSS_weekly(i)        ## Paths of the 4rth week of the month
    return path1,path2,path3,path4,path5,path6,path7


'''
   Input :  Month Number of the year 2020, Week Number in a Particular Month
   Output:  Generating weekly paths according to the week number

'''
    
def Weekly_path_Generating_CYGNSS(Month_No,Week_No): 
    if Month_No ==1:
        m = 1
    elif Month_No ==2:
        m = 32
    elif Month_No ==3:
        m = 61       
    elif Month_No ==4:
        m = 92
    elif Month_No ==5:
        m = 122
    elif Month_No ==6:
        m = 153      
    elif Month_No ==7:
        m = 183       
    elif Month_No ==8:
        m = 214       
    elif Month_No ==9:
        m = 245
    elif Month_No ==10:
        m = 275
    elif Month_No ==11:
        m = 306     
    else:
        m = 336
        
    ## Calling to previous function to generate the paths according to starting day number of the week
    path1,path2,path3,path4,path5,path6,path7 = Weekly_Path(Week_No,m)   
    return path1,path2,path3,path4,path5,path6,path7 


'''
    Input : Path of a single day file
    Ouput : Gridded values of surface reflectivity on SMAP grid points

'''
def SMSR_Daily(path):
    Data = pd.read_csv(path)
    SM   = np.array(Data['SM_SMAP'])
    SR   = np.array(Data['SR_SMAP'])
    SM1  = pd.DataFrame(SM.reshape(34,47))          ## Reshaping it into its original shape of the gridded form
    SR1  = pd.DataFrame(SR.reshape(34,47))          ## Reshaping surface reflectivity interpolated along SMAP grid points 
    SR2  = pd.DataFrame(Masking_Data_SR(np.array(SR1)))
    
    SR2  = np.array(SR2.replace(to_replace=np.nan,value=0))  ## Replacing the soil moisture and surface reflectivity values of NAN by 0
    SM1  = np.array(SM1.replace(to_replace=np.nan,value=0))  ## So that stacking or averaging operation will be easy
    return SM1,SR2

'''
    Input : Weekly Paths 
    Ouput : Weekly Gridded values of surface reflectivity on SMAP grid points

'''
def Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7):
    SM1,SR1 = SMSR_Daily(path1)
    SM2,SR2 = SMSR_Daily(path2)
    SM3,SR3 = SMSR_Daily(path3)
    SM4,SR4 = SMSR_Daily(path4)
    SM5,SR5 = SMSR_Daily(path5)
    SM6,SR6 = SMSR_Daily(path6)
    SM7,SR7 = SMSR_Daily(path7)
    return SM1,SM2,SM3,SM4,SM5,SM6,SM7,SR1,SR2,SR3,SR4,SR5,SR6,SR7


'''
    Input : Grids of the two days surface reflectivity
    Ouput : Stacked values of surface reflectivity on SMAP grid points

'''
def Fill_Gap(SR1,SR2):
    m,n = SR1.shape
    for i in range(m):
        for j in range(n):
            if SR1[i][j]==0:            ## If surface reflectivity value available in previous grid is zero at any grid point
                SR1[i][j] = SR2[i][j]   ## Replacing previous zero surface reflectivity values by next day value at that point
    return SR1



'''
   Input : Month number and week number
   Ouput : Spatial correlation between weekly stacked Surface Reflectivity and Soil Moisture values and their respective image plots
   
'''
def Weekly_Stacking_SM(Month_No,Week_No):
    
    ## Creating weekly paths of the files
    path1,path2,path3,path4,path5,path6,path7               = Weekly_path_Generating_CYGNSS(Month_No,Week_No)
    
    ## Weekly Gridded values of surface reflectivity and Soil Moisture 
    SM1,SM2,SM3,SM4,SM5,SM6,SM7,SR1,SR2,SR3,SR4,SR5,SR6,SR7 = Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7)
    
    ## Weekly Stacking Soil Moisture 
    SM_u1 = Fill_Gap(SM1,SM2)                  # Stacking Soil Moisture Values of 1st and 2nd day of soil moisture of 1st week
    SM_u2 = Fill_Gap(SM_u1,SM3)                # Stacking Soil Moisture Values of 2nd and 3rd day of soil moisture of 1st week
    SM_u3 = Fill_Gap(SM_u2,SM4)                # Stacking Soil Moisture Values of 3rd and 4th day of soil moisture of 1st week
    SM_u4 = Fill_Gap(SM_u3,SM5)                # Stacking Soil Moisture Values of 4th and 5th day of soil moisture of 1st week
    SM_u5 = Fill_Gap(SM_u4,SM6)                # Stacking Soil Moisture Values of 5th and 6th day of soil moisture of 1st week
    SM_u6 = Fill_Gap(SM_u5,SM7)                # Stacking Soil Moisture Values of 6th and 7th day of soil moisture of 1st week
    
    SM_u6 = pd.DataFrame(SM_u6)
    SM_u6 = np.array(SM_u6)
    
    
    ## Weekly Stacking surface reflectivity 
    SR_u1 =   Fill_Gap(SR1,SR2)            # Stacking surface reflectivity Values of 1st and 2nd day of surface reflectivity of 1st week
    SR_u2 = Fill_Gap(SR_u1,SR3)            # Stacking surface reflectivity Values of 2nd and 3rd day of surface reflectivity of 1st week
    SR_u3 = Fill_Gap(SR_u2,SR4)            # Stacking surface reflectivity Values of 3rd and 4th day of surface reflectivity of 1st week
    SR_u4 = Fill_Gap(SR_u3,SR5)            # Stacking surface reflectivity Values of 4th and 5th day of surface reflectivity of 1st week
    SR_u5 = Fill_Gap(SR_u4,SR6)            # Stacking surface reflectivity Values of 5th and 6th day of surface reflectivity of 1st week
    SR_u6 = Fill_Gap(SR_u5,SR7)            # Stacking surface reflectivity Values of 6th and 7th day of surface reflectivity of 1st week
    
    SR_u6 = pd.DataFrame(SR_u6)
    SR_u6 = np.array(SR_u6)
    
    ## Making correlation of weekly stacked data
    SR_SM            = pd.DataFrame(SR_u6.flatten())               # Flattening the weekly stacked Surface Reflectivity data points
    SR_SM.columns    = ['SR_SMAP']
    SR_SM['SM_SMAP'] = np.array((np.array(SM_u6)).flatten())       # Flattening the weekly stacked SOIL Moisture data points
    
    DS               = SR_SM.copy()
    Mask             = (DS['SR_SMAP']>0) & (DS['SR_SMAP']<20)
    New_SR_SM        = (DS[Mask]).dropna()
    
    Corr             = np.array(New_SR_SM.corr())
    CR               = np.round((Corr[0][1])*100,2)
    
    ## Visualizing the weekly stacked data of soil moisture and surface reflectivity
    plt.figure(figsize=(20,20))
    plt.subplot(1,2,1)
    cmap = 'gist_ncar'
    img = plt.imshow(SR_u6,interpolation = 'gaussian',cmap=cmap)
    cbar = plt.colorbar(img,shrink=0.27)
    cbar.set_label('Surface Reflectivity in dB')
    plt.title(f'''Month No. {Month_No} Week No. {Week_No}
    Surface Reflectivity of CYGNSS interpolated along SMAP gridpoints''') 
    
    plt.subplot(1,2,2)
    cmap = 'gist_ncar'
    img = plt.imshow(SM_u6,interpolation = 'gaussian',cmap=cmap)
    cbar = plt.colorbar(img,shrink=0.27)
    cbar.set_label('Soil Moisture in cm3/cm3')
    plt.title(f'''Soil Moisture of SMAP gridpoints: Month No. {Month_No} Week No. {Week_No} 
    Correlation: {CR} %''') 
    
    
    Points = []
    for i in range(len(SM_u6.flatten())):
        Points.append(i)
        
    plt.figure(figsize=(40,15))
    plt.plot(Points[:1000],SM_u6.flatten()[:1000],label = 'SMAP Soil Moisture')
    plt.plot(Points[:1000],SR_u6.flatten()[:1000]/20,label    = 'Normalized CYGNSS Derived Surface Reflectivity')
    plt.ylabel('''Surface reflectiviy and Soil Moisture''',fontsize=20)
    plt.xlabel('Points',fontsize=20)
    plt.legend(fontsize=20) 
    
    sns.lmplot(x='SR_SMAP',y='SM_SMAP',data=DS,aspect=4,height=6)
    plt.xlabel('Normalized CYGNSS Derived Surface Reflectivity')
    plt.ylabel('SMAP Soil Moisture in cm3/cm3')
    plt.title(f'''Fitting a Regression Line using lmplot''')

'''
   Input : Month number and week number
   Ouput : Spatial correlation between weekly averaged Surface Reflectivity and Soil Moisture values and their respective image plots
   
'''      
def WeeklyAvg_SM(Month_No,Week_No):
    ## Creating weekly paths of the files
    path1,path2,path3,path4,path5,path6,path7 = Weekly_path_Generating_CYGNSS(Month_No,Week_No)
    
    ## Weekly Gridded values of surface reflectivity and Soil Moisture
    SM1,SM2,SM3,SM4,SM5,SM6,SM7,SR1,SR2,SR3,SR4,SR5,SR6,SR7 = Read_Weekly_SM(path1,path2,path3,path4,path5,path6,path7)
    
    ## Taking weekly average of the surface reflectivity and soil moisture
    SM_avg = np.array(pd.DataFrame((SM1+SM2+SM3+SM4+SM5+SM6+SM7)/7))    
    SR_avg = np.array(pd.DataFrame((SR1+SR2+SR3+SR4+SR5+SR6+SR7)/7))
    
    ## Making correlation of weekly averaged data
    SR_SM            = pd.DataFrame(SR_avg.flatten())
    SR_SM.columns    = ['SR_SMAP']
    SR_SM['SM_SMAP'] = np.array((np.array(SM_avg)).flatten())
    
    DS               = SR_SM.copy()
    Mask             = (DS['SR_SMAP']>0) & (DS['SR_SMAP']<20)
    New_SR_SM        = (DS[Mask]).dropna()
    Corr             = np.array(New_SR_SM.corr())
    CR               = np.round((Corr[0][1])*100,2)
    
    ## Visualizing the weekly averaged data of soil moisture and surface reflectivity
    plt.figure(figsize=(20,20))
    plt.subplot(1,2,1)
    cmap = 'gist_ncar'
    img = plt.imshow(SR_avg,interpolation = 'gaussian',cmap=cmap)
    cbar = plt.colorbar(img,shrink=0.27)
    cbar.set_label('Avg Surface Reflectivity in dB')
    plt.title(f'''Month No. {Month_No} Week No. {Week_No}
    Avg Surface Reflectivity of CYGNSS interpolated along SMAP gridpoints''') 
    
    plt.subplot(1,2,2)
    cmap = 'gist_ncar'
    img = plt.imshow(SM_avg,interpolation = 'gaussian',cmap=cmap)
    cbar = plt.colorbar(img,shrink=0.27)
    cbar.set_label('Avg Soil Moisture in cm3/cm3')
    plt.title(f'''Avg Soil Moisture of SMAP gridpoints: Month No. {Month_No} Week No. {Week_No} 
    Correlation: {CR} %''') 
    
    Points = []
    for i in range(len(SM_avg.flatten())):
        Points.append(i)
        
    plt.figure(figsize=(40,15))
    plt.plot(Points[:1000],SM_avg.flatten()[:1000],label = 'SMAP Soil Moisture')
    plt.plot(Points[:1000],SR_avg.flatten()[:1000]/20,label    = 'Normalized CYGNSS Derived Surface Reflectivity')
    plt.ylabel('''Surface reflectiviy and Soil Moisture''',fontsize=20)
    plt.xlabel('Points',fontsize=20)
    plt.legend(fontsize=20) 
    
    sns.lmplot(x='SR_SMAP',y='SM_SMAP',data=DS,aspect=4,height=6)
    plt.xlabel('Normalized CYGNSS Derived Surface Reflectivity')
    plt.ylabel('SMAP Soil Moisture in cm3/cm3')
    plt.title(f'''Fitting a Regression Line using lmplot''')