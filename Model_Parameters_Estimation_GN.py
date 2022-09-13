import numpy as np
import sympy
from sympy import *
import pandas as pd
import matplotlib.pyplot as plt

# Creating Variables
a,b,c,d,e1,e2,Mvt,Sgma,theeta,l = symbols('a b c d e1 e2 Mvt Sgma theeta l')

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
 
''' Function to calculate : 
    P : weight matrix, 
    n : number of equations, 
    r : redundancy, 
    u : number of parameters'''

def Adjustment_Parameters(Sgma0):
    P = np.identity(len(Sgma0))         # Weight matrix (Precision matrix)
    n = len(Sgma0)                      # As total n erroneous observations are there  
    
    # We have to calculate the model parameters and soil moisture 
    u = 6+1             # A, B, C, D, E1, E2, Mvt
    r = n-u             # Redundancy
    return P,n,u,r

'''Function for Jacobian of non-linear expression with respect to a single variable'''
def Der(f,var):
    A1 = (Derivative(f,var)).doit()
    return A1

'''Function for Design Matrix'''
def Jacobian(f,X,n):
    A1 = []
    m   = len(X)
    for i in range(n):
        for j in range(m):
            X1 = X[j]
            a1 = Der(f,X1)
            A1.append(a1)
    A1 = np.array(A1)
    A = A1.reshape(n,7)
    return A

def Resubstituting_Values(A1,Sgma0,L,Theta,a0,b0,c0,d0,e10,e20,Mv0):
    n1,n2 = A1.shape
    A_1 = []
    for i in range(n1):
        Sgma1 = Sgma0[i]
        L1 = L[i]
        Theta1 = Theta[i]
        for j in range(n2):
            a1 = A1[i][j]
            a1 = a1.subs({Sgma:Sgma1,l:L1,theeta:Theta1,Mvt:Mv0,a:a0,b:b0,c:c0,d:d0,e1:e10,e2:e20})
            A_1.append(a1)      
    A_1 = np.array(A_1)    
    A1  = A_1.reshape(n1,n2)
    return A1

## Calculations of Biomass using LAI
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

## Converting the elements of the matrix from the object to the float
def Object_toFloat(Matrix):
    m,n = Matrix.shape
    M1 = []
    for i in range(m):
        for j in range(n):
            M1.append(float(Matrix[i,j]))
    M2 = np.array(M1)
    M = M2.reshape(m,n)
    return M  

## Updates in Parameters, N and U matrix 
def delX(A,P,Sgma0):
    N = (np.transpose(A))@P@A
    N = Object_toFloat(N)
    U = (np.transpose(A))@P@Sgma0
    Part1 = np.linalg.pinv(N)
    dX1   = Part1@U
    return dX1,N,U

'''Main Iteration Function'''

def WCM_Daily_SM(lat,lon):
    Day_No     = []
    WCM_SM     = []
    GLDAS_InSM = []
    a_var      = []
    b_var      = []
    c_var      = []
    d_var      = []
    e1_var     = []
    e2_var     = []
    for i in range(1,365):
        D = i
        Day_No.append(D)
        # Path of Dth Day Observations of CYGNSS and LAI 
        P1  = f'D:\EG\Project Data\CYGNSS_Obs_Chambal_{lat}_{lon}\LAI_{lat}_{lon}'
        P2  = f'\CYGNSS_Data_Availability_Day_{Day(D)}_{lat}_{lon}.csv'
        CYGNSS_LAI_Path = P1+P2
        
        # Path for the GLDAS SM data for initialization purpose of Soil Moisture Value
        GLDAS_Path = f'D:\EG\Project Data\CYGNSS_Obs_Chambal_{lat}_{lon}'+f'\GLDAS_SM_{lat}_{lon}_Day_1-366.csv'
        
        # Reading the CYGNSS and LAI data
        CYGNSS_LAI_Data = pd.read_csv(CYGNSS_LAI_Path)
        GLDAS_SM = pd.read_csv(GLDAS_Path)
        
        # Backscatter in dB
        Sgma0 = CYGNSS_LAI_Data['SR_eff']
        
        # LAI to Biomass
        LAI = CYGNSS_LAI_Data['LAI']
        L = LAI_L(LAI)
        
        # Incidence angle
        Theta = (CYGNSS_LAI_Data['SP_I']/180)*np.pi
        
        # GLDAS SM
        GLDAS_InSM.append(np.round(float(GLDAS_SM['GLDAS_SM'][D]),3))
    
        LAI1 = pd.DataFrame(LAI)
        LAI1 = LAI1.dropna()
        if len(LAI1)>6:
            LAI = LAI1
            # Modified Water Cloud Model Equation
            f   = Sgma-a*((l/100)**e1)*cos(theeta)*(1-exp(-2*b*((l/100)**e2)/cos(theeta)))+((c)*(Mvt)+d)*exp(-2*b*((l/100)**e2)/cos(theeta))
            
            # Initialising Parameters
            a0    = 0.1836
            b0    = 0.0862
            c0    = 13.92
            d0    = 0.195
            e10   = -5
            e20   = 0
            Mv0   = GLDAS_SM['GLDAS_SM'][D]
            X     = [a0,b0,c0,d0,e10,e20,Mv0]
            
            # Weight matrix and number of equations and number of parameters
            P,n,u,r = Adjustment_Parameters(Sgma0)
            
            # Design Matrix
            X1  = [a,b,c,d,e1,e2,Mvt]
            A1  = Jacobian(f,X1,n)
            A   = Resubstituting_Values(A1,Sgma0,L,Theta,a0,b0,c0,d0,e10,e20,Mv0)
            
            # Updates in Parameters
            dX,N,U = delX(A,P,Sgma0)
            
            # Updated initial parameters
            X  = X+dX
            
            # Residuals in Backscatter
            V  = A@dX-Sgma0   
            
            # Updated Backscatter
            Sgma0 = Sgma0+V
            
            for i in range(20):
                a0,b0,c0,d0,e10,e20,Mv0 = X[0],X[1],X[2],X[3],X[4],X[5],X[6]
                A                       = Resubstituting_Values(A1,Sgma0,L,Theta,a0,b0,c0,d0,e10,e20,Mv0)
                dX1,N,U                 = delX(A,P,Sgma0) 
                X1                      = X+dX1
                Thres1                  = (dX1-dX)**2
                Thres                   = np.max(Thres1**(0.5))
                if (Thres <(10**(-2))):
                    Xf = X
                    WCM_SM.append(np.round(float(Xf[6]),3))
                    a_var.append(np.round(float(Xf[0]),3))
                    b_var.append(np.round(float(Xf[1]),3))
                    c_var.append(np.round(float(Xf[2]),3))
                    d_var.append(np.round(float(Xf[3]),3))
                    e1_var.append(np.round(float(Xf[4]),3))
                    e2_var.append(np.round(float(Xf[5]),3))
                    break
                else:
                    X  = X1
                    dX = dX1
                    V  = A@dX-Sgma0   
                    Sgma0 = Sgma0+V                 
        else:
            WCM_SM.append(np.nan)
            a_var.append(np.nan)
            b_var.append(np.nan)
            c_var.append(np.nan)
            d_var.append(np.nan)
            e1_var.append(np.nan)
            e2_var.append(np.nan)
    DF = pd.DataFrame(Day_No)
    DF.columns = ['Day_No']
    DF['A']  = a_var
    DF['B']  = b_var
    DF['C']  = c_var
    DF['D']  = d_var
    DF['E1'] = e1_var
    DF['E2'] = e2_var
    DF['WCM_SM']  = WCM_SM
    return DF            

'''Visualization and correlation of WCM Soil Moisture with SMAP soil moisture'''
def WCM_SMAP_SM(WCM_SM,lat,lon):
    SMAP_SM = pd.read_csv(f'D:\EG\Project Data\CYGNSS_Data_in_0p36Dg\SMAP_RF_SM\SMAP_SM_Variations_{lat}_{lon}.csv')
    SMAP_SM['WCM_SM'] = WCM_SM['WCM_SM']
    SMAP_SMf = SMAP_SM.dropna()
    SMAP_SM = SMAP_SM.drop(['Day_No'],axis = 1)
    
    CR = (np.array(SMAP_SM.corr()))*100
    CR1 = np.round(CR[0][1],2) 
    RMSE = np.sum((SMAP_SMf['SMAP_SM']/100 - SMAP_SMf['WCM_SM'])**2)/len(SMAP_SMf['SMAP_SM'])
    plt.figure(figsize=(30,5))
    plt.scatter(SMAP_SMf['Day_No'],SMAP_SMf['SMAP_SM']/100,label='SMAP_SM')
    plt.scatter(SMAP_SMf['Day_No'],SMAP_SMf['WCM_SM'],label='WCM_SM')
    plt.title(f'''RMSE: {np.round(RMSE,4)} Correlation: {np.abs(CR1)} %''') 
    plt.legend() 
    
    plt.figure(figsize=(30,30))
    plt.subplot(6,1,1)
    plt.scatter(WCM_SM['Day_No'],WCM_SM['A'])
    plt.title(f'''Vegetation Parameter A at Latitude: {lat} Longitude: {lon} Mean: {np.mean(WCM_SM['A'])}''') 
    
    plt.subplot(6,1,2)
    plt.scatter(WCM_SM['Day_No'],WCM_SM['B'])
    plt.title(f'''Vegetation Parameter B at Latitude: {lat} Longitude: {lon} Mean: {np.mean(WCM_SM['B'])}''') 
    
    plt.subplot(6,1,3)
    plt.scatter(WCM_SM['Day_No'],WCM_SM['C'])
    plt.title(f'''Soil Parameter C at Latitude: {lat} Longitude: {lon} Mean: {np.mean(WCM_SM['C'])}''') 
    
    plt.subplot(6,1,4)
    plt.scatter(WCM_SM['Day_No'],WCM_SM['D'])
    plt.title(f'''Soil Parameter D at Latitude: {lat} Longitude: {lon} Mean: {np.mean(WCM_SM['D'])}''') 
    
    plt.subplot(6,1,5)
    plt.scatter(WCM_SM['Day_No'],WCM_SM['E1'])
    plt.title(f'''Canopy Geometry Parameter E1 at Latitude: {lat} Longitude: {lon} Mean: {np.mean(WCM_SM['E1'])}''') 
    
    plt.subplot(6,1,6)
    plt.scatter(WCM_SM['Day_No'],WCM_SM['E2'])
    plt.title(f'''Canopy Geometry Parameter E2 at Latitude: {lat} Longitude: {lon} Mean: {np.mean(WCM_SM['E2'])}''') 