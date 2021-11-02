'''
Created on Feb 19, 2021

Creating a function to read depth and output temperature,light, and oxygen for specific depth. 

@author: garmin
'''
from numpy import *
from pylab import *
from matplotlib import *

#Read csv file into code
data=genfromtxt('..\\Data\\Depth_T_O2_L.csv',delimiter=',')

#want all rows except for the first for all four parameters from the data
depth_d=data[1:,0]      #depth
temp_d=data[1:,1]       #temperature
o2_d=data[1:,2]         #oxygen
lit_d=data[1:,3]        #light
    
#create a function to output temperature and oxygen as a function of depth
def TO2L(depth):
    if depth >= 1803:
        return(o2_d[1799],temp_d[1799],lit_d[1799])
        
    if round(depth)<=depth:
        i=round(depth,0)-4  #want this i to give us the element that corresponds with depth if the depth rounds down
    elif round(depth)>depth:
        i=round(depth,0)-5  #want this i to give us the element that corresponds with depth if the depth will round up
    
    i=int(i)            #need i to be an integer, not float, so that it can be used to identify array elements
    o2=o2_d[i]+(depth-depth_d[i])*((o2_d[i+1]-o2_d[i])/(depth_d[i+1]-depth_d[i]))           #linear interpolation equation for o2
    temp=temp_d[i]+(depth-depth_d[i])*((temp_d[i+1]-temp_d[i])/(depth_d[i+1]-depth_d[i]))    #linear interpolation equation for temp
    lit=lit_d[i]+(depth-depth_d[i])*((lit_d[i+1]-lit_d[i])/(depth_d[i+1]-depth_d[i]))        #linear interpolation equation for light
               
    return(o2,temp,lit)

def Arr(Temp):
    A = 10000
    TrefK = 298.15
    TempK = Temp + 273.15
    return exp(-A*((1/TempK)-(1/TrefK))) #arrehenius equation (Geider, 1997) function of temperature

#########
O2,temp,lit=TO2L(1802)     
#print(O2,temp,lit)

        
