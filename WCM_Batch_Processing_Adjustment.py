import numpy as np
import pandas as pd
from pandas import DataFrame
import seaborn as sns
import sympy
from sympy import *

# Creating Variables
a,b,c,d,e1,e2,Mvt,Sgma,theeta,l = symbols('a b c d e1 e2 Mvt Sgma theeta l')

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

# Jacobian of function wrt a single variable
def Der(f,var):
    A1 = (Derivative(f,var)).doit()
    return A1

def Adjustment_Parameters(Sgma0,DD):
    P = np.identity(len(Sgma0))         # Weight matrix (Precision matrix)
    n = len(Sgma0)                      # As total n erroneous observations are there  
    
    # We have to calculate the model parameters and soil moisture 
    u = 6+len(DD)             # A, B, C, D, E1, E2, Mvt6, Mvt8, ....., Mvt365 Assume those day soil moisture only having Biomass
    r = n-u                   # Redundancy
    return P,n,u,r

def N_Day(a1):
    a2 = a1[6:]
    for i in range(len(a2)):
        a3 = a2[i]
        if a3!=0:
            a4 = i
    return a4

def Resubstituting_Values(A1,Sgma0,L,Theta,a0,b0,c0,d0,e10,e20,Mv0,DD):
    n1,n2 = A1.shape
    A_1 = []
    for i in range(n1):
        Sgma1 = np.array(Sgma0)[i]
        L1 = np.array(L)[i]
        Theta1 = np.array(Theta)[i]
        
        a2 = A1[i]
        k  = N_Day(a2)
        d1 = int(np.array(DD)[k])
        SM = Mv0[d1-1]
        for j in range(n2):
            a1 = A1[i][j]
            if a1!=0:
                aij = a1.subs({Sgma:Sgma1,l:L1,theeta:Theta1,Mvt:SM,a:a0,b:b0,c:c0,d:d0,e1:e10,e2:e20})
                A_1.append(aij)
            else:
                A_1.append(0)
    A_1 = np.array(A_1)    
    A1  = A_1.reshape(n1,n2)
    A   = np.array(A1)
    return A

## From 2nd iteration since size of Mv0 reduces to number of days available
def Resubstituting_Values2(A1,Sgma0,L,Theta,a0,b0,c0,d0,e10,e20,X0):
    n1,n2 = A1.shape
    A_1 = []
    for i in range(n1):
        Sgma1 = np.array(Sgma0)[i]
        L1 = np.array(L)[i]
        Theta1 = np.array(Theta)[i]
        
        a2 = A1[i]
        k  = N_Day(a2)
        Xk = X0[6:]
        SM = Xk[k]
        for j in range(n2):
            a1 = A1[i][j]
            if a1!=0:
                aij = a1.subs({Sgma:Sgma1,l:L1,theeta:Theta1,Mvt:SM,a:a0,b:b0,c:c0,d:d0,e1:e10,e2:e20})
                A_1.append(aij)
            else:
                A_1.append(0)
    A_1 = np.array(A_1)    
    A1  = A_1.reshape(n1,n2)
    A   = np.array(A1)
    return A

# Updated Parameters
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

## Batch Processing with 1 year data
def Bach_processing(Df1,lat,lon):
    DD = (Df1['Day_No'])
    DD = DD.drop_duplicates()
    
    Sgma0  = Df1['SR_eff']
    L      = LAI_L(Df1['LAI'])
    Theta  = Df1['SP_I']
    Day_No = Df1['Day_No']

    # GLDAS Soil Moisture Data for initialization purpose of 366 days of 2020 within a pixel of 36 x 36 Km of SMAP
    GLDAS_Path = f'D:\EG\Project Data\CYGNSS_Obs_Chambal_{lat}_{lon}'+f'\GLDAS_SM_{lat}_{lon}_Day_1-366.csv'
    GLDAS_SM1  = pd.read_csv(GLDAS_Path)
    GLDAS_SM1  = GLDAS_SM1.drop(['CYGNSS_Backscatter','CYGNSS_SPI','AGB'],axis=1)
    GLDAS_SM   = np.array(GLDAS_SM1['GLDAS_SM'])
    Mv0        = GLDAS_SM

    # Modified Water Cloud Model Equation
    f   = Sgma-a*((l)**e1)*cos(theeta)*(1-exp(-2*b*((l)**e2)/cos(theeta)))+c*Mvt*exp(-2*b*((l)**e2)/cos(theeta))+d*exp(-2*b*((l)**e2)/cos(theeta))

    # Model Parameters (n, u, r, P)
    P,n,u,r = Adjustment_Parameters(Sgma0,DD)

    # Generating Design Matrix
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

    ## Substituting Values in Design Matrix
    # Initializing Parameters
    a0    = 0.1836
    b0    = 0.0862
    c0    = 13.92
    d0    = 0.195
    e10   = -5
    e20   = 0
    X0 = [a0,b0,c0,d0,e10,e20]
    for i in range(len(DD)):
        dd = int(np.array(DD)[i])
        SM = np.round(Mv0[dd-1],4)
        X0.append(SM)

    # Substituting Inialized Parameters and the Observations of CYGNSS and all
    A = Resubstituting_Values(A1,Sgma0,L,Theta,a0,b0,c0,d0,e10,e20,Mv0,DD)
    
    # Updates in Parameters
    dX,N,U = delX(A,P,Sgma0)
    X = X0+dX

    # Residuals in Backscatter
    V  = A@dX-Sgma0   

    # Updated Backscatter
    Sgma0 = Sgma0+V

    for i in range(20):
        a0,b0,c0,d0,e10,e20,Mv0 = X[0],X[1],X[2],X[3],X[4],X[5],X[6:]
        A                       = Resubstituting_Values2(A1,Sgma0,L,Theta,a0,b0,c0,d0,e10,e20,X)
        dX1,N,U                 = delX(A,P,Sgma0) 
        X1                      = X+dX1
        Thres1                  = (dX1-dX)**2
        Thres                   = np.max(Thres1**(0.5))
        if (Thres <(10**(-5))):
            Xf = X
            break
        else:
            X  = X1
            dX = dX1
            V  = A@dX-Sgma0   
            Sgma0 = Sgma0+V  
    return Sgma0,Sgma0_updated,V,A,n,u,r,Xf,dX,DD

## Batch Processing with 1 year data
def Bach_processing1(Df1,lat,lon):
    DD = (Df1['Day_No'])
    DD = DD.drop_duplicates()
    
    Sgma0  = Df1['SR_eff']
    L      = LAI_L(Df1['LAI'])
    Theta  = Df1['SP_I']
    Day_No = Df1['Day_No']

    # GLDAS Soil Moisture Data for initialization purpose of 366 days of 2020 within a pixel of 36 x 36 Km of SMAP
    GLDAS_Path = f'D:\EG\Project Data\CYGNSS_Obs_Chambal_{lat}_{lon}'+f'\GLDAS_SM_{lat}_{lon}_Day_1-366.csv'
    GLDAS_SM1  = pd.read_csv(GLDAS_Path)
    GLDAS_SM1  = GLDAS_SM1.drop(['CYGNSS_Backscatter','CYGNSS_SPI','AGB'],axis=1)
    GLDAS_SM   = np.array(GLDAS_SM1['GLDAS_SM'])
    Mv0        = GLDAS_SM

    # Modified Water Cloud Model Equation
    f   = Sgma-a*((l)**e1)*cos(theeta)*(1-exp(-2*b*((l)**e2)/cos(theeta)))+c*Mvt*exp(-2*b*((l)**e2)/cos(theeta))+d*exp(-2*b*((l)**e2)/cos(theeta))

    # Model Parameters (n, u, r, P)
    P,n,u,r = Adjustment_Parameters(Sgma0,DD)

    # Generating Design Matrix
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

    ## Substituting Values in Design Matrix
    # Initializing Parameters
    a0    = 0.1836
    b0    = 0.0862
    c0    = 13.92
    d0    = 0.195
    e10   = -5
    e20   = 0
    X0 = [a0,b0,c0,d0,e10,e20]
    for i in range(len(DD)):
        dd = int(np.array(DD)[i])
        SM = np.round(Mv0[dd-1],4)
        X0.append(SM)

    # Substituting Inialized Parameters and the Observations of CYGNSS and all
    A = Resubstituting_Values(A1,Sgma0,L,Theta,a0,b0,c0,d0,e10,e20,Mv0,DD)
    
    # Updates in Parameters
    dX,N,U = delX(A,P,Sgma0)
    X = X0+dX

    # Residuals in Backscatter
    V  = A@dX-Sgma0   

    # Updated Backscatter
    Sgma0_updated = Sgma0+V
    
    return Sgma0,Sgma0_updated,V,A,n,u,r,X,dX,DD