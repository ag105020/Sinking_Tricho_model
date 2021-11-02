'''
Created on May 18, 2014

@author: Keisuke
'''
from pylab import * 

def sf(fignumber):
    First_part="C:\\Users\\Keiin\\OneDrive\\Desktop\\figures\\02\\11 Mar\\"
    Figure_number=str(fignumber)
    Last_part=".png"
    savefig(First_part+Figure_number+Last_part,dpi=300)
