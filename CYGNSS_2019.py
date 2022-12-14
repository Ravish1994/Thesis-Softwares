import pandas as pd
import numpy as np
from netCDF4 import Dataset
import xarray
import requests
import urllib.request
import xarray as xr
import io

''' 
CYGNSS Level 1 Science Data Record Version3.1
Publisher: PO.DAAC
Release Date: 2021-Dec-23
https://podaac-opendap.jpl.nasa.gov/opendap/hyrax/allData/cygnss/cygnss_v3.1/L1/v3.1/2019/
Basic inputs:
m = starting Day number 
n = ending Day number 
                      ***********Updating functions:
Update(Add) your required variables as per your need
---------------------------------------------------SUB FUNCTIONS--------String definition for the day number in year(365 or 366)'''
def fun1(i):
    if i<= 9:
        b = '00'+str(i)
        return b
    elif i>=10 and i<=99:
        b = '0'+str(i)
        return b
    elif i>=100 and i<=365:
        b = str(i)
        return b   
    
'''String definition for the satellite number'''
def fun2(i):
    if i>=0 and i<=8:
        b = '0'+str(i)
        return b
    
'''String definition for the month number'''
def fun3(i):
    if i<=9:
        b = '0'+str(i)
        return b
    elif i>=10 and i <=12:
        b = str(i)
        return b    
    
'''String definition for the day number in month(28 or 30 or 31)'''
def fun4(i):
    if i<=9:
        b = '0'+str(i)
        return b
    elif i>=10 and i <=31:
        b = str(i)
        return b

''' Formatting Date '''    
def fun5(i):
    if (i>=1) and (i<=31):
        b1 = fun3(1)                           ## Month Number(Jan)
        b2 = fun4(i)                           ## Day number in month
        c1 = b1+b2
                
    elif (i>=32) and (i<=59):
        b1 = fun3(2)                           ## Month Number(Feb)
        b2 = fun4(i-31)                           ## Day number in month
        c1 = b1+b2
                
    elif (i>=60) and (i<=90):
        b1 = fun3(3)                           ## Month Number(Mar)
        b2 = fun4(i-59)                        ## Day number in month
        c1 = b1+b2
                
    elif (i>=91) and (i<=120):
        b1 = fun3(4)                           ## Month Number(Apr)
        b2 = fun4(i-90)                        ## Day number in month
        c1 = b1+b2
        
    elif (i>=121) and (i<=151):
        b1 = fun3(5)                           ## Month Number(May)
        b2 = fun4(i-120)                       ## Day number in month
        c1 = b1+b2
        
    elif (i>=152) and (i<=181):
        b1 = fun3(6)                           ## Month Number(Jun)
        b2 = fun4(i-151)                       ## Day number in month
        c1 = b1+b2
        
    elif (i>=182) and (i<=212):
        b1 = fun3(7)                           ## Month Number(Jul)
        b2 = fun4(i-181)                       ## Day number in month
        c1 = b1+b2
        
    elif (i>=213) and (i<=243):
        b1 = fun3(8)                           ## Month Number(Aug)
        b2 = fun4(i-212)                       ## Day number in month
        c1 = b1+b2
        
    elif (i>=244) and (i<=273):
        b1 = fun3(9)                           ## Month Number(Sept)
        b2 = fun4(i-243)                       ## Day number in month
        c1 = b1+b2
        
    elif (i>=274) and (i<=304):
        b1 = fun3(10)                          ## Month Number(Oct)
        b2 = fun4(i-273)                       ## Day number in month 
        c1 = b1+b2
        
    elif (i>=305) and (i<=334):
        b1 = fun3(11)                          ## Month Number(Nov)
        b2 = fun4(i-304)                       ## Day number in month
        c1 = b1+b2
        
    elif (i>=335) and (i<=365):
        b1 = fun3(12)                          ## Month Number(Dec)
        b2 = fun4(i-334)                       ## Day number in month
        c1 = b1+b2
    return c1

''' Formatting Satellite IDs'''
def fun6(i,j):   
    a4 = '-000000-e2019'       ## +b1+b2
    return a4

def fun7(i,j):
    if (i==12) and (j==6):
        a5 = '-205725.l1.power-brcs.a31.d32.nc.nc?'
    else:
        a5 = '-235959.l1.power-brcs.a31.d32.nc.nc?'
    return a5


'''--------------------------------Collecting Samples for 2020 for every satellites every day--------------------------------------'''

       # Cyg01   Cyg02   Cyg03   Cyg04   Cyg05   Cyg06   Cyg07   Cyg08
TS1 = [[85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_001          
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_002           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_003             
       [85000, 85000, 85000, 85000, 85000, 84000, 85000, 85000],     # Day_004          
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_005
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_006           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_007
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_008           
       [85000, 85000, 85000, 85000, 85000, 85000, 82316, 85000],     # Day_009           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_010
       [80000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_011
       [85000, 85000, 85000, 85000, 85000, 76000, 85000, 85000],     # Day_012            
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_013
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_014           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_015           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_016        January(01/2020)      
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_017
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_018           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_019           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_020          
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_021
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_022           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_023
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_024           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_025           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_026
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_027           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_028           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_029           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_030            
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_031
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_032
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_033       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_034       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_035
       [85000, 85000, 84000, 85000, 85000, 85000, 85000, 85000],     # Day_036      
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_037       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_038       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_039
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_040       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_041       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_042   
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_043
       [80000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_044
       [84000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_045
       [84000, 85000, 85000, 85000, 85000, 78000, 85000, 85000],     # Day_046        February(02/2020)
       [85000, 85000, 85000, 85000, 78000, 85000, 85000, 85000],     # Day_047
       [85000, 85000, 85000, 85000, 85000, 85000, 81000, 85000],     # Day_048       
       [85000, 84000, 85000, 85000, 85000, 85000, 85000, 82000],     # Day_049       
       [85000, 80000, 85000, 85000, 85000, 85000, 85000, 84000],     # Day_050      
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_051       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_052       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_053       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_054       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_055       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_056       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_057       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_058
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_059       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_060 
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_061
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_062
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_063
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_064       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_065
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_066
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_067
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_068       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_069       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_070       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_071
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_072
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_073
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_074
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_075        March(03/2020)
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_076
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_077
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_078
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_079       
       [85000, 85000, 85000, 84000, 85000, 85000, 85000, 85000],     # Day_080
       [85000, 85000, 85000, 80000, 85000, 85000, 85000, 85000],     # Day_081
       [85000, 83000, 85000, 85000, 85000, 83000, 85000, 85000],     # Day_082
       [80000, 85000, 80000, 85000, 85000, 85000, 85000, 85000],     # Day_083
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 80000],     # Day_084
       [85000, 85000, 85000, 80000, 85000, 85000, 85000, 85000],     # Day_085       
       [85000, 85000, 85000, 85000, 80000, 85000, 80000, 85000],     # Day_086
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_087
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_088
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_089       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_090       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_091
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_092
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_093
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_094       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_095
       [85000, 85000, 85000, 85000, 80000, 85000, 85000, 85000],     # Day_096       
       [85000, 85000, 85000, 80000, 85000, 85000, 85000, 85000],     # Day_097       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_098
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_099       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_100
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_101
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_102
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_103       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_104
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_105       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_106         April(04/2020)       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_107        
       [85000, 85000, 85000, 85000, 85000, 85000, 80000, 85000],     # Day_108
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_109       
       [85000, 85000, 85000, 85000, 85000, 85000, 80000, 85000],     # Day_110       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_111
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_112
       [80000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_113
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 80000],     # Day_114       
       [85000, 80000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_115
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_116       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_117       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_118       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_119       
       [80000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_120       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_121 
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_122
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_123
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_124       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_125
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_126       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_127       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_128
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_129       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_130        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_131
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_132
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_133       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_134       
       [85000, 85000, 85000, 85000, 85000, 85000, 80000, 85000],     # Day_135       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_136          May(05/2020)     
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_137       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_138
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_139       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_140       
       [85000, 85000, 85000, 85000, 85000, 80000, 85000, 85000],     # Day_141         
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_142
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_143
       [85000, 85000, 85000, 85000, 85000, 80000, 85000, 85000],     # Day_144        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_145
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_146        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_147       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_148       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_149       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_150        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_151
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_152
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_153
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_154       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_155
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_156       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_157       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_158
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_159        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_160        
       [80000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_161       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_162
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 80000],     # Day_163
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_164        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_165          June(06/2020)
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_166       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_167       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_168
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_169       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_170       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_171      
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 80000],     # Day_172
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_173       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_174       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_175
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_176        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_177       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_178
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_179       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_180       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_181       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_182
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_183       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_184        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_185
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_186       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_187       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_188
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_189       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_190        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_191       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_192
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_193       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_194        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_195
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_196          July (07/2020)      
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_197       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_198       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_199        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_200           
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_201
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_202
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_203       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_204       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_205       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_206       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_207        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_208       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_209       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_210        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_211       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_212       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_213 
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_214         
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_215       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_216       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_217       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_218       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_219        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_220       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_221       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_222       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_223       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_224         
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_225       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_226       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_227       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_228       August (08/2020)     
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_229       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_230       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_231       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_232       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_233       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_234         
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_235       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_236         
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_237        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_238       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_239       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_240       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_241        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_242       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_243       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_244
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_245       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_246       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_247        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_248
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_249        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_250        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_251       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_252       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_253
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_254        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_255       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_256       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_257       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_258       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_259       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_260        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_261        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_262      September(09/2020)
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_263
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_264       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_265      
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_266      
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_267       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_268
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_269       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_270       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_271
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_272
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_273
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_274
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_275
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_276       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_277        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_278       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_279        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_280       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_281      
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_282       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_283       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_284        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_285
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_286       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_287       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_288       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_289       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_290        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_291
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_292       October(10/2020)   
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_293       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_294       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_295       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_296       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_297       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_298       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_299       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_300        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_301       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_302       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_303       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_304       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_305
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_306       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_307       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_308       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_309        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_310        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_311       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_312       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_313       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_314       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_315       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_316       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_317       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_318       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_319       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_320       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_321        November(11/2020)       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_322       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_323       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_324       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_325       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_326       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_327       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_328       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_329        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_330       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_331       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_332       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_333       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_334       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_335
       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_336        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_337        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_338
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_339       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_340       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_341        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_342       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_343       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_344       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_345
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_346       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_347       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_348       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_349       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_350        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_351       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_352      
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_353 --------> December (12/2020)      
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_354       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_355       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_356       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_357       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_358       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_359        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_360       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_361       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_362       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_363       
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_364        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000],     # Day_365        
       [85000, 85000, 85000, 85000, 85000, 85000, 85000, 85000]]     # Day_366
       

    
'''-------------------------------------------------Generating Link--------------------------------------------------------'''

def Newly_Added_Variables(List_Variables,TS):
    n = len(List_Variables)
    s3 = '%5B0:1:16%5D%5B0:1:10%5D'
    S = f'%5B0:1:{TS}%5D%5B0:1:3%5D'                         ## Sample range
    L = ''
    for V in List_Variables:
        if V=='power_analog':
            L = L+f',{V}'+S+s3
        else:
            L = L+f',{V}'+S
    return L

def Creating_Array_Links(Variables):
    a1 = 'https://podaac-opendap.jpl.nasa.gov/opendap/hyrax/allData/cygnss/cygnss_v3.1/L1/v3.1/2019/'  ## +a, Day Number in Year
    a2 = '/cyg'                                                          ## +b, Satellite number
    a3 = '.ddmi.s2019'                                                   ## +b1, Month number +b2, Day Number in month
    v3  = 'sp_lat' ## +S
    v4  = ',sp_lon' ## +S
    url1 = []
    for i in range(1,366):
        a = fun1(i)                            ## Day number in year
        B = fun5(i)
        for j in range(1,9):                   ## Satellite number
            a4 = fun6(i,j)
            a5 = fun7(i,j)
            TS = TS1[i-1][j-1]
            s3 = '%5B0:1:16%5D%5B0:1:10%5D'
            S = f'%5B0:1:{TS}%5D%5B0:1:3%5D'                         ## Sample range
            b = fun2(j)                                              ## Satellite number
            url_1 = a1+a+a2+b+a3+B+a4+B+a5+v3+S+v4+S
            url_2 = Newly_Added_Variables(Variables,TS)  ## User Variables
            url_0 = url_1+url_2
            url1.append(url_0)
    url2 = np.array(url1)
    url = url2.reshape(365,8)
    
    '''Handling missing satellite's data'''
    L1 = [12,13,14,20,21,57,58,58,59,70,71,96,97,192,193,208,208,209,209,210,211,212,227,228,230,231,239,240,249,249,250,251,253,254,260,
          261,262,263,287,288,289,290,317,347,348,349]
    L2 = [5,5,5,6,6,6,3,6,3,7,7,3,1,7,7,2,6,2,6,4,4,4,2,2,3,3,1,1,2,4,4,4,2,2,0,0,0,0,7,7,7,7,6,6,6,6]
    for i,j in zip(L1,L2):
        url[i][j] = url[57][5]
    return url

'''Code for subsetting'''
def masking_data_in_Region(ds,lat_lowerLim,lat_upperLim,lon_lowerLim,lon_upperLim):
    Keys = ds.variables.keys()
    Keys = [K for K in Keys]
    K1 = ['sample', 'ddm', 'sp_lat', 'sp_lon']
    for K1 in K1:
        Keys.remove(f'{K1}')
    Lat1 = (np.array(ds['sp_lat'])).flatten()
    Lat1 = (pd.DataFrame(Lat1)).replace(to_replace=-9999,value=np.nan)
    Lat1.columns = ['sp_lat']
    Lat1['sp_lon'] = (pd.DataFrame((np.array(ds['sp_lon'])).flatten())).replace(to_replace=-9999,value=np.nan)
    for K in Keys:
        if K=='power_analog':
            Pi = np.array(ds[f'{K}'])
            l,m,n,o = Pi.shape
            P1 =  (pd.DataFrame((np.array(ds[f'{K}'])).flatten())).replace(to_replace=-9999,value=np.nan)
            P1 = P1.to_numpy()
            P1 = P1.reshape(l,4,17,11)
            df1 = []
            for i in range(l):
                for j in range(4):
                    a = np.max(P1[i][j][:])
                    df1.append(a)      
            df2 = np.array(df1)
            DDM_peak = df2.reshape(l,4)
            Lat1['ddm_peak'] = DDM_peak.flatten()
        else:
            Lat1[f'{K}'] = (pd.DataFrame((np.array(ds[f'{K}'])).flatten())).replace(to_replace=-9999,value=np.nan)
    mask = ((Lat1['sp_lat']>lat_lowerLim) & (Lat1['sp_lat']<lat_upperLim)) & ((Lat1['sp_lon']>lon_lowerLim) & (Lat1['sp_lon']<lon_upperLim))
    Lat1 = Lat1[mask]
    return Lat1
 
                
def Subsetting_CYGNSS_Data(m,n,lat1,lon1,lat2,lon2,Variables):
    url = Creating_Array_Links(Variables)
    for i in range(m-1,n):
        B = fun5(i+1)
        for j in range(7,8):
            urls = url[i][j]
            r = requests.get(urls)
            Path = rf"D:\EG\Project Data\CYGNSS_Raw_Ganga_Catchment_2019\Day_{i+1}\cyg0{j+1}.ddmi.s2019{B}-000000-e2019{B}-235959.l1.power-brcs.a30.d31.nc.nc"
            with open(Path,'wb') as f:
                f.write(r.content) 
            dataset = Dataset(Path,'r')
            Subset_Data = masking_data_in_Region(dataset,lat_lowerLim=lat1,lat_upperLim=lat2,lon_lowerLim=lon1,lon_upperLim=lon2)
            xr = xarray.Dataset.from_dataframe(Subset_Data)
            xr.to_netcdf(Path)
                               
 # ---------------------------------------------------END OF CODE-----------------------------------------------------------------