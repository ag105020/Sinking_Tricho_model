'''
Created on Oct 11, 2018
This function gives viscusity mu at a certain temperature
It is based on Kestin 1978 and "15 Computation of diffusion coefficient for different temperature.xlsx"
learned from 633 00 07
@author: Keisuke
'''

from pylab import *

def TemperatureViscosity(Temperature):
    mu10 = 1306.9
    mu15 = 1138.2
    mu20 = 1002
    mu25 = 890.3
    mu30 = 797.5
    mu35 = 719.5
    Kelvin = 273.15
    
    TC = Temperature - Kelvin
    
    if TC >= 10 and TC <15:
        Tlow = 10
        Thigh = 15
        mu = mu10 + (mu15-mu10)/(Thigh-Tlow)*(Temperature-Tlow-Kelvin)
    elif TC >= 15 and TC <20:
        Tlow = 15
        Thigh = 20
        mu = mu15 + (mu20-mu15)/(Thigh-Tlow)*(Temperature-Tlow-Kelvin)
    elif TC >= 20 and TC <25:
        Tlow = 20
        Thigh = 25
        mu = mu20 + (mu25-mu20)/(Thigh-Tlow)*(Temperature-Tlow-Kelvin)
    elif TC >= 25 and TC <30:
        Tlow = 25
        Thigh = 30
        mu = mu25 + (mu30-mu25)/(Thigh-Tlow)*(Temperature-Tlow-Kelvin)
    elif TC >= 30 and TC <35:
        Tlow = 30
        Thigh = 35
        mu = mu30 + (mu35-mu30)/(Thigh-Tlow)*(Temperature-Tlow-Kelvin)
    else: 
        #mu = -475.9*log(Temperature-Kelvin) + 2415.9
        mu = 51.049e-10*TC**6 - 18041e-10*TC**5 + 27.597e-5*TC**4 - 2440.4e-5*TC**3 + 1.42426*TC**2 - 59.753*TC + 1783.8 #from "15 Computation of diffusion coefficient for different temperature.xlsx"
        #mu = 0.2236*TC - 32.024*TC + 1572.6
    return mu