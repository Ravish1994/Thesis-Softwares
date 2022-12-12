''' -------------------------------------CYGNSS Data Sub-setting software for 2019---------------------------------'''
import CYGNSS_2019
import imp
imp.reload(CYGNSS_2019)
import CYGNSS_2019 as SC1

''' -------------------------------------CYGNSS Data Sub-setting software for 2020---------------------------------'''
import CYGNSS_2020
import imp
imp.reload(CYGNSS_2020)
import CYGNSS_2020 as SC2

''' -------------------------------------CYGNSS Data Sub-setting software for 2021---------------------------------'''
import CYGNSS_2021
import imp
imp.reload(CYGNSS_2021)
import CYGNSS_2021 as SC3


''' -------------------------------------CYGNSS Data Sub-setting software for 2019, 2020, 2021---------------------------------'''

def Subsetting_CYGNSS_Data(Year,Initial_Day_Number,Ending_Day_Number,lat1,lon1,lat2,lon2,Variables):
    if Year==2019:
        SC1.Subsetting_CYGNSS_Data(Initial_Day_Number,Ending_Day_Number,lat1,lon1,lat2,lon2,Variables)
    elif Year==2020:
        SC2.Subsetting_CYGNSS_Data(Initial_Day_Number,Ending_Day_Number,lat1,lon1,lat2,lon2,Variables)
    elif Year==2021:
        SC3.Subsetting_CYGNSS_Data(Initial_Day_Number,Ending_Day_Number,lat1,lon1,lat2,lon2,Variables)
    else:
        print('Enter the Year 2019 or 2020 or 2021')
        
'''-------------------------------------------------------End of The Code--------------------------------------------------------'''