import numpy as np
import pandas as pd
from pandas import DataFrame
import seaborn as sns
import sympy
from sympy import *

# Calculations of Biomass using LAI
def LAI_L(LAI):
    LAI_F     = 2.1
    SLA       = 30
    wd_fract  = 0.25
    LAI_MX    = 2.5
    foliar_W  = 1.30
    wood_W    = 0.54
    foliar_C  = LAI*LAI_F/SLA
    wood_C    = wd_fract*LAI_MX*1.25
    veg_water = foliar_C*2.22*foliar_W + wood_C*2.22*wood_W
    L         = 2.22*foliar_C + 2.22*wood_C + veg_water
    return L

# Water Cloud Model Equation
a,b,c,d,e1,e2,Mvt,Sgma,theeta,l = symbols('a b c d e1 e2 Mvt Sgma theeta l')
f   = Sgma-(a*((l)**e1)*cos(theeta)*(1-exp(-2*b*((l)**e2)/cos(theeta)))+c*Mvt*exp(-2*b*((l)**e2)/cos(theeta))+d*exp(-2*b*((l)**e2)/cos(theeta)))

# Model Parameters (n, u, r, P)
def Adjustment_Parameters1(Sgma0,DD):
    P = np.identity(len(Sgma0))         # Weight matrix (Precision matrix)
    n = len(Sgma0)                      # As total n erroneous observations are there  
    
    # We have to calculate the model parameters and soil moisture 
    u = 6+len(DD)             # A, B, C, D, E1, E2, Mvt6, Mvt8, ....., Mvt365 Assume those day soil moisture only having Biomass
    r = n-u                   # Redundancy
    return P,n,u,r

# Generating Design Matrix
def Der(f,var):
    A1 = (Derivative(f,var)).doit()
    return A1

def Design_Matrix(Df1,a,b,c,d,e1,e2,Mvt,DD):
    C1 = Der(f,a)
    C2 = Der(f,b)
    C3 = Der(f,c)
    C4 = Der(f,d)
    C5 = Der(f,e1)
    C6 = Der(f,e2)
    C7 = Der(f,Mvt)
    A1 = []
    for j in range(len(Df1)):
        Single_Row_A = [C1,C2,C3,C4,C5,C6]
        for i in range(len(DD)):
            d1 = np.array(DD)[i]
            d2 = int(np.array(Df1['Day_No'])[j])
            if d1 == d2:
                Single_Row_A += [C7]
            else:
                Single_Row_A += [0]
        A1 += Single_Row_A
    A2 = np.array(A1)
    A1 = A2.reshape(len(Df1),(len(DD)+6))
    return A1

# Substituting Values in Design Matrix
def Resubstituting_Values(A1,Sgma0,L,Theta,SM_GLDAS,a0,b0,c0,d0,e10,e20):
    n1,n2 = A1.shape
    A_1   = []
    for i in range(n1):
        Sgma1  = np.array(Sgma0)[i]
        L1     = np.array(L)[i]
        Theta1 = np.array(Theta)[i]
        Mv0    = np.array(SM_GLDAS)[i]
        for j in range(n2):
            a1 = A1[i][j]
            if a1!=0:
                aij = a1.subs({Sgma:Sgma1,l:L1,theeta:Theta1,Mvt:Mv0,a:a0,b:b0,c:c0,d:d0,e1:e10,e2:e20})
                A_1.append(aij)
            else:
                A_1.append(0)
    A_1 = np.array(A_1)    
    A1  = A_1.reshape(n1,n2)
    A   = np.array(A1)
    return A

## Residual vector
def Misclosure(Sgma0,L,Theta,SM_GLDAS,a0,b0,c0,d0,e10,e20,f):
    W_0 = []
    for i in range(len(Sgma0)):
        Sgma1  = np.array(Sgma0)[i]
        L1     = np.array(L)[i]
        Theta1 = np.array(Theta)[i]
        Mv0    = np.array(SM_GLDAS)[i]
        f5 = f.subs({Sgma:Sgma1,l:L1,theeta:Theta1,Mvt:Mv0,a:a0,b:b0,c:c0,d:d0,e1:e10,e2:e20})
        W_0.append(f5)
    W_0 = np.array(W_0)
    W1 = W_0.reshape(len(Sgma0),1)
    return W1

# Converting the elements of the matrix from the object to the float
def Object_toFloat(Matrix):
    m,n = Matrix.shape
    M1 = []
    for i in range(m):
        for j in range(n):
            M1.append(float(Matrix[i,j]))
    M2 = np.array(M1)
    M = M2.reshape(m,n)
    return M

# Updating Soil Moisture
def New_SM(DF_batch,X_SM,idx):
    Mv0 = X_SM
    DF_batch = DF_batch.drop(['GLDAS_SM'], axis=1)
    SM_initialized = []
    for i in range(len(DF_batch)):
        D  = int(np.array(idx)[i])
        SM = Mv0[D]
        SM_initialized.append(SM)
    SM_initialized = np.array(SM_initialized).reshape(len(SM_initialized),) 
    return SM_initialized

# Performing more Iteration
def Adjustment_Parameters(DF_batch,idx,A1,Sgma0,L,Theta,X,f,P):
    X_SM     = X[6:]
    GLDAS_SM = New_SM(DF_batch,X_SM,idx)
    A        = Resubstituting_Values(A1,Sgma0,L,Theta,GLDAS_SM,X[0][0],X[1][0],X[2][0],X[3][0],X[4][0],X[5][0])
    W        = Misclosure(Sgma0,L,Theta,GLDAS_SM,X[0][0],X[1][0],X[2][0],X[3][0],X[4][0],X[5][0],f)
    N        = np.array(np.transpose(A)@P@A)
    U        = (np.transpose(A))@P@W
    N        = Object_toFloat(N)
    Part1    = np.linalg.pinv(N)       ## Inverse of N
    Part2    = Object_toFloat(Part1)   ## Converting elements to float
    dX       = -Part2@U                ## Change in X
    return A,W,N,U,dX

def Plotting_Variations(Df,Var1,Var2,label1,label2,lat,lon,RMSE,CR):
    import numpy as np
    import matplotlib.pyplot as plt
    plt.figure(figsize=(30,8))
    plt.scatter(Df['Day_No'],Df[f'{Var1}'],label = label1)
    plt.scatter(Df['Day_No'],Df[f'{Var2}'],label = label2)
    
    plt.title(f'''Latitude: {lat} Longitude: {lon} RMSE: {round(RMSE,4)} Correlation: {CR} %''',size=20) 
    plt.xlabel('Day number of the year 2020',size=20)
    plt.ylabel('Volumetric Soil Moisture',size=20)
    plt.ylim(0,1)
    plt.xticks(np.arange(1, 370, 10),size=15)
    plt.yticks(np.arange(0, 1, 0.1),size=15)
    plt.legend(fontsize=30)

    plt.figure(figsize=(3,3))
    plt.scatter(Df[f'{Var1}'],Df[f'{Var2}'],s=40)
    plt.xlabel(f'{Var1}')
    plt.ylabel(f'{Var2}')
    plt.plot([0,1],[0,1],c='gray')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xticks(np.arange(0, 1, 0.2))
    plt.yticks(np.arange(0, 1, 0.2))
    
    import seaborn as sns
    plt.figure(figsize=(2,2))
    sns.lmplot(x=f'{Var1}',y=f'{Var2}',data=Df,line_kws={'color': 'black'})
    plt.xlabel(f'{Var1}')
    plt.ylabel(f'{Var2}')
    plt.plot([0,1],[0,1],c='gray')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xticks(np.arange(0, 1, 0.2))
    plt.yticks(np.arange(0, 1, 0.2))

# Correlation of WCM with GLDAS soil moisture
def New_SM1(DF_batch,X_SM,idx):
    Mv0 = X_SM
    SM_initialized = []
    for i in range(len(DF_batch)):
        D  = int(np.array(idx)[i])
        SM = float(Mv0[D])
        SM_initialized.append(SM)
    SM_initialized = np.array(SM_initialized).reshape(len(SM_initialized),) 
    DF_batch['WCM_SM'] = SM_initialized
    return DF_batch

def Bach_Adjustment_For_366_Days_WCM(lat,lon):
    # Collecting whole year data for a single Pixel Point
    import Preparing_Dataset_for_BatchProcessing
    import imp
    imp.reload(Preparing_Dataset_for_BatchProcessing)
    import Preparing_Dataset_for_BatchProcessing as r5

    DF_batch = r5.Data_Batch_2020(lat,lon)

    # GLDAS Soil Moisture Data for initialization purpose of 366 days of 2020 within a pixel of 36 x 36 Km of SMAP
    GLDAS_Path = f'D:\EG\Project Data\CYGNSS_Obs_Chambal_{lat}_{lon}'+f'\GLDAS_SM_{lat}_{lon}_Day_1-366.csv'
    GLDAS_SM1  = pd.read_csv(GLDAS_Path)
    GLDAS_SM1  = GLDAS_SM1.drop(['CYGNSS_Backscatter','CYGNSS_SPI','AGB'],axis=1)
    GLDAS_SM   = np.array(GLDAS_SM1['GLDAS_SM'])
    Mv0        = GLDAS_SM

    SM_initialized = []
    for i in range(len(DF_batch)):
        D = np.array(DF_batch['Day_No'])[i]-1
        SM_initialized.append(Mv0[int(D)])
    DF_batch['GLDAS_SM'] = SM_initialized

    # Unique days
    Df1 = DF_batch
    DD  = (Df1['Day_No'])
    DD  = DD.drop_duplicates()

    # Collecting observations
    Sgma0    = Df1['SR_eff']
    L        = LAI_L(Df1['LAI'])
    Theta    = Df1['SP_I']
    GLDAS_SM = Df1['GLDAS_SM']
    Day_No   = Df1['Day_No']

    P,n,u,r = Adjustment_Parameters1(Sgma0,DD)
    A1 = Design_Matrix(Df1,a,b,c,d,e1,e2,Mvt,DD)
    a0    = 0.0094645
    b0    = 0.09645
    c0    = 13.92
    d0    = 0
    e10   = -5
    e20   = -190       ## Should be less than or -190

    SM_initial = np.array(DF_batch['GLDAS_SM'].drop_duplicates())
    X0         = [a0,b0,c0,d0,e10,e20]  # X0 is the Initial Parameters
    for i in range(len(SM_initial)):
        X0.append(SM_initial[i])

    A = Resubstituting_Values(A1,Sgma0,L,Theta,GLDAS_SM,a0,b0,c0,d0,e10,e20)

    W = Misclosure(Sgma0,L,Theta,GLDAS_SM,a0,b0,c0,d0,e10,e20,f)

    # Finding Normal Matrix
    N = np.array(np.transpose(A)@P@A)

    N = Object_toFloat(N)
    U = (np.transpose(A))@P@W

    # Change in Parameters
    Part1  = np.linalg.pinv(N)       ## Inverse of N
    Part2  = Object_toFloat(Part1)   ## Converting elements to float
    dX     = -Part2@U                ## Change in X
    X0 = np.array(X0)
    X0 = X0.reshape(len(X0),1)
    X = X0+dX

    # Creating Index for updating soil moisture
    Df = pd.DataFrame(DF_batch['Day_No'].value_counts(sort=False))['Day_No']
    DD = Df.values
    idx = []
    m = 0
    for i in range(len(DD)):
        j = DD[i]
        for k in range(j):
            idx.append(m)
        m = m+1

    # 1st Iteration
    X_SM = X[6:]
    GLDAS_SM = New_SM(DF_batch,X_SM,idx)
    A        = Resubstituting_Values(A1,Sgma0,L,Theta,GLDAS_SM,X[0][0],X[1][0],X[2][0],X[3][0],X[4][0],X[5][0])
    W        = Misclosure(Sgma0,L,Theta,GLDAS_SM,X[0][0],X[1][0],X[2][0],X[3][0],X[4][0],X[5][0],f)
    N        = np.array(np.transpose(A)@P@A)
    U        = (np.transpose(A))@P@W
    N        = Object_toFloat(N)
    Part1    = np.linalg.pinv(N)       ## Inverse of N
    Part2    = Object_toFloat(Part1)   ## Converting elements to float
    dX1      = -Part2@U                ## Change in X

    X = X+dX1
    dX = dX1
    for i in range(1,200):
        A,W,N,U,dX1 = Adjustment_Parameters(DF_batch,idx,A1,Sgma0,L,Theta,X,f,P)
        Thres       = np.max(abs(dX-dX1))
        if Thres<10**(-5):
            A_adjusted  = A
            W_adjusted  = W
            N_adjusted  = N
            U_adjusted  = U
            X_adjusted  = X
            dX_adjusted = dX
            Number_Iteration = i
            break
        else:
            dX = dX1
            X  = X+dX
    A = pd.DataFrame(A_adjusted)
    W = pd.DataFrame(W_adjusted) 
    N = pd.DataFrame(N_adjusted)
    X = pd.DataFrame(X_adjusted)
    dX = pd.DataFrame(dX_adjusted)
    DD1 = DF_batch['Day_No'].drop_duplicates()
    DD  = DD1.values
    SM_adjusted = X[6:].values

    # Correlation of WCM with SMAP soil moisture
    SMAP_SM = pd.read_csv(f'D:\EG\Project Data\CYGNSS_Data_in_0p36Dg\SMAP_RF_SM\SMAP_SM_Variations_{lat}_{lon}.csv')  
    SMAP_SM['SMAP_SM'] = SMAP_SM['SMAP_SM']/100
    x2  = np.array(SMAP_SM['SMAP_SM'])
    SMAP_SM1 = SMAP_SM
    for j in range(len(DD)):
        D2 = int(DD[j])
        SM = float(SM_adjusted[j][0])
        SMAP_SM1['SMAP_SM'] = SMAP_SM1['SMAP_SM'].replace(np.array(SMAP_SM1['SMAP_SM'])[D2-1],SM)  
    x1  = np.array(SMAP_SM1['SMAP_SM'])
    DF1 = pd.DataFrame(np.array(SMAP_SM1['Day_No']))
    DF1.columns = ['Day_No']
    DF1['SMAP_SM'] = x2
    DF1['Improved_SM'] = x1
    SMAP_SM2 = DF1[DF1['Improved_SM']<0.9]

    CR = np.array(SMAP_SM2.corr())
    CR = round((CR[1][2])*100,3)
    RMSE = round(np.sum((SMAP_SM2['SMAP_SM'] - SMAP_SM2['Improved_SM'])**2)/len(SMAP_SM2['SMAP_SM']),3)
    Plotting_Variations(SMAP_SM2
                        ,'SMAP_SM'
                        ,'Improved_SM'
                        ,'SMAP Soil Moisture on Barren Land'
                        ,'After removing vegetation attenuation WCM_SM',lat,lon,RMSE,CR)

    X_SM     = X[6:].values
    GLDAS_WCM_SM = New_SM1(DF_batch,X_SM,idx)
    GLDAS_WCM_SM1 = GLDAS_WCM_SM.drop(['SR_eff','SP_I','sp_lat','sp_lon','LAI'], axis=1)
    GLDAS_WCM_SM1 = GLDAS_WCM_SM1.drop_duplicates()
    GLDAS_WCM_SM1 = GLDAS_WCM_SM1[GLDAS_WCM_SM1['WCM_SM']<1]
    CR = np.array(GLDAS_WCM_SM1.corr())
    CR = round((CR[1][2])*100,3)
    RMSE = round(np.sum((GLDAS_WCM_SM1['GLDAS_SM'] - GLDAS_WCM_SM1['WCM_SM'])**2)/len(GLDAS_WCM_SM1['GLDAS_SM']),3)
    Plotting_Variations(GLDAS_WCM_SM1
                        ,'GLDAS_SM'
                        ,'WCM_SM'
                        ,'GLDAS Soil Moisture on Barren Land'
                        ,'After removing vegetation attenuation WCM soil moisture',lat,lon,RMSE,CR)

    # Correlation of SMAP with GLDAS soil moisture
    GLDAS_SM1['SMAP_SM'] = SMAP_SM2['SMAP_SM']
    CR = np.array(GLDAS_SM1.corr())
    CR = round((CR[1][2])*100,3)
    RMSE = round(np.sum((GLDAS_SM1['GLDAS_SM'] - GLDAS_SM1['SMAP_SM'])**2)/len(GLDAS_SM1['GLDAS_SM']),3)

    Plotting_Variations(GLDAS_SM1
                        ,'GLDAS_SM'
                        ,'SMAP_SM'
                        ,'GLDAS Soil Moisture on Barren Land'
                        ,'SMAP Soil Moisture on Barren Land',lat,lon,RMSE,CR)