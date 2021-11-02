'''
Created on Apr 5, 2021
needs to be non-linear before 500 meters. At 500 meters, velocity is unchanging.
@author: garmin
'''

from pylab import *

#before 500 meters: changing velocity with depth
def SinkV(time,a,b):
    time=time/86400             #convert time to days

    c=log(40*a+b)
    depth=(1/a)*(e**((a*time)+c)-b)

    return(depth)


