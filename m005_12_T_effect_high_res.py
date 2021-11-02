'''
Created on Mar 21, 2018
From t000 41 16
@author: Keisuke
''
Here I add data for O2
'''

from pylab import *
from McCarty02 import *
from McCarty01 import *
from Depth_O2_Temp_FINAL import *
from sf import *
from TemperatureViscosity import *
from Sinking_vel8 import *

def tricho(a1,b1,CsFac):
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Time setting 
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    mu25 = TemperatureViscosity(25+273.15)
    days = 7
    Tmax = 86400*days #(s) maximum time step (177-32)
    dT = 50
    T = arange(0.,Tmax+10**(-10),dT)  #(s) time array (177-32)
    U = (T/dT).astype(int)                   #(s) U array for for loop (from a704_16)
#     Depth = 40 + T*w/86400  #(m) Depth: start sinking at 40(m)
    Depth = SinkV(T,a1,b1) #(m)
    LD = ones(size(T))  #1 = light, 0 = dark
    dl=46740     #(hours) Day length. From "13 Time and location.xlsx"
#    LD[0] = 0
    #LD[T>=dl]=0
    for i in arange(days):
        LD[(T>=dl+i*86400) & (T<86400+i*86400)]=0  #applying it to each day with for loop
    
    #LD[0]=0
    #-----------------------------------
    #photosynthesis related parameters
    #-----------------------------------
    Pmax0=0.002681321
    O=0.01  #(m2 / umol photon) absorption cross section (from a704_16)
    Tau=1    #(s) handling time (from a704_16) (Irina: half life of RNA is less than 1 minute (30 (s)) or so)
    Rd=0.2*ones(size(T)) #* sin(pi*T/dl)     #(dimensionless) The ratio of diazocyte
#     Rd=0.325*ones(size(T)) #Average through the day 
    for i in arange(days):
        Rd[(T>10800+i*86400) & (T<32400+i*86400)]=0.45  #applying it to each day with for loop
    Rp=1-Rd     #(dimensionless) The ratio of photosynthetic cells
    Chlmaxg=0.048     #(gchl gC-1) Chlorophyll max of the cell(186-40)
    Chlmaxmol=Chlmaxg*12*55/868      #(molCchl molC-1) Chlorophyll max of the cell in mol (186-40)
    ro=18333    #(molC / m3) carbon density in the cell
    Chl=Chlmaxmol*0.8*ro        #(molCchl m-3) Chlorophyll concentration
    Chl=Chl*Rp
    Ync=1/6.3   #(molN / molC) ratio of N:C (LaRoche 2005)  
    Mumax=0.25/86400  #(s-1) maximum growth rate
    lmax0=ro*Mumax  #(molC m-3 s-1)
    Kc=ro/5     #(molC / m3) half saturation constant of carbon storage use
    Kn=Kc*Ync   #(molN / m3) half saturation constant of nitrogen storage use

    Fnitroge = 0.4
    
    NstoMax=ro*Ync      #(molN m-3) Maximum nitrogen storage
    CstoMax=ro*2
    NtoC=McCarty02(0.6)
    E = McCarty01(0.6)
    
    #---------------------------------------
    #Array preparation for each prarameter
    #---------------------------------------
    o=copy(T)*0             # this creates zero array for the right size for the time steps
    dCsto=copy(o)
    dNsto=copy(o)
    PI=copy(o)
    Res=copy(o)
    ls=copy(o)
    Csto=copy(o)
    Nfix=copy(o)
    Res2=copy(o)
    Nsto=copy(o)
    Res1=copy(o)
    Cn2fix=copy(o)
    Ox1=copy(o)
    Ox2=copy(o)
    Ox3=copy(o)
    Csto2=copy(o)
    CcN2fix=copy(o)
    Res2c=copy(o)
    which=copy(o)
    Nf1=copy(o)
    Nf2=copy(o)
    Nf3=copy(o)
    ResN2fix=copy(o)
    ResPro=copy(o)
    Iarray=copy(o)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Initial condition
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #=======================
    # C store factor array
    #=======================
    #CsFac = arange(0,2+1e-10,0.2)  #(s) time array (177-32)
    #CsFac = array([0.5,2])
    U2 = arange(size(CsFac)).astype(int)         #(s) U array for for loop (from a704_16)
    Uij =copy(U2)*0
    
    CstoBox= zeros((size(Uij),size(U)))

    #=================
    
    Ox1[0]=0.213
    Ox2[0]=0.213
    Ox3[0]=0.213
    
    Nsto[0]=1000
        
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For statement
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    for j in U2:
        Csto[0]=ro*CsFac[j]
        Clim2 = 10
        Ui = 0
        for i in U:
            Ox4,Temp,I = TO2L(Depth[i])  
            if LD[i] == 0:
               I = 0
            #---T effect----
            lmax = lmax0*Arr(Temp)
            ResMax=lmax*200  #(molC m-3 s-1) Maximum respiration rate  
            NfixFull=lmax*Ync*25       #(molN m-3 s-1) Maximmum nitrogen fixation rate
            Pmax = Pmax0*Arr(Temp)
            NfixMax=NfixFull*Rd*Fnitroge
            K = 273.15
            mu = TemperatureViscosity(Temp+K)
            fTd = (Temp+K)/(25+K)/mu*mu25 #factor for tempearture dependence on diffusion
            
            #---------------
            
            CaFactor = 1265
            CaFactor = 1
            Ca= 1.8*fTd    #Constant for A
            Ca0 = 1.8*fTd
            #Ca=100
            #E=0.7
            
            A13=2     *Ca*CaFactor
            A12=0.0345    *Ca*CaFactor
            A32=2    *Ca*CaFactor
            A34=2      *Ca0
            
            A13=A13*Rp[i]
            A32=A32*Rd[i]
            
            
            OxCri=0.1
            
            #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            #Preparing parts
            #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            PI[i] = Pmax * (1 - exp(-O * I * Tau)) * (CstoMax-Csto[i])/CstoMax     #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)
            
            Clim=Csto[i] / (Csto[i] + Kc)
            Nlim=Nsto[i] / (Nsto[i] + Kn)
            ls[i] = lmax * min(Clim, Nlim)    #(C mol m-3 s-1) biomass synthesis rate
            
            if Clim>Nlim:
                which[i]=1
            else:
                which[i]=0
            
            Res1[i] = ls[i] * E
            Csto2[i] = Csto[i] / Rd[i]
            Res2c[i] = ResMax*Csto2[i] / (Csto2[i] + Kc) * Rd[i]
            
            #=========================================================================================================
            # O2 flux part Kei210-33 + 02EquationSolving02simplified.nb and 03EquationSolving03different.nb
            #=========================================================================================================
            
            #-------------------------------------------------------
            # When Ox2 = 0 "03EquationSolving03different.nb"
            #-------------------------------------------------------
            a = -Res1[i]  + PI[i] * Chl[i]
            Ox1[i] = -((-a*A13-a*A32-a*A34-A13*A34*Ox4)/(A12*A13+A12*A32+A13*A32+A12*A34+A13*A34))
            Ox3[i] = -((-a*A13-A12*A34*Ox4-A13*A34*Ox4)/(A12*A13+A12*A32+A13*A32+A12*A34+A13*A34))
            b   = -(a*A12*A13+a*A12*A32+a*A13*A32+a*A12*A34+A12*A13*A34*Ox4+A12*A32*A34*Ox4+A13*A32*A34*Ox4)\
                   /(A12*A13+A12*A32+A13*A32+A12*A34+A13*A34)
            Ox2[i] = 0
            
            Res2[i] = -b
            
            #-------------------------------------------------------
            # Testing b
            #-------------------------------------------------------
            
            if Res2[i] > Res2c[i]:
            
                #-------------------------------------------------------
                # When Ox2 =/ 0 "02EquationSolving02simplified.nb"
                #-------------------------------------------------------
                b= -Res2c[i]
                Ox1[i] = -((-a*A12*A13-a*A12*A32-a*A13*A32-a*A12*A34-a*A32*A34-A12*A13*b-A12*A32*b-A13*A32*b\
                         -A12*A34*b-A12*A13*A34*Ox4-A12*A32*A34*Ox4-A13*A32*A34*Ox4)/((A12*A13+A12*A32+A13*A32)*A34))
                Ox2[i] = -((-a*A12*A13-a*A12*A32-a*A13*A32-a*A12*A34-A12*A13*b-A12*A32*b-A13*A32*b-A12*A34*b\
                         -A13*A34*b-A12*A13*A34*Ox4-A12*A32*A34*Ox4-A13*A32*A34*Ox4)/((A12*A13+A12*A32+A13*A32)*A34))
                Ox3[i] = (a+b+A34*Ox4)/A34
            
                Res2[i] = Res2c[i]
               
            #=========================================================================================================
            if LD[i] == 0:
                Res2[i] = 0
            
            Nf1[i] = max((OxCri - Ox2[i]) / OxCri, 0)
            Nf2[i] = max((NstoMax - Nsto[i]) / NstoMax, 0)
            Nf3[i] = Csto2[i] / (Csto2[i] + Kc)
            
            Nfix[i] = NfixMax[i] * Nf1[i] * Nf2[i] * Nf3[i]   #(molN m-3 s-1) Nitrogen fixation rate per biomass carbon
            Nfix[LD==0] = 0
            CcN2fix[i] = Nfix[i]         #(molC m-3 s-1) carbon consumption for nitrogen fixation 
            ResN2fix[i] = Nfix[i] * NtoC            #(molC m-3 s-1) Respiration for nitrogen fixation (in carbon)
            ResPro[i] = Res2[i] - ResN2fix[i] 
            Res[i] = Res1[i] + Res2[i] 
             
            #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            #time step change
            #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            dCsto[i] = PI[i] * Chl[i] - Res[i] - ls[i] - Cn2fix[i] - CcN2fix[i]
            dNsto[i] = Nfix[i] - ls[i] * Ync

            #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            #obtaining next time step values
            #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            if i<U[-1]:
                Csto[i+1] = Csto[i] + dCsto[i] * dT
                Nsto[i+1] = Nsto[i] + dNsto[i] * dT
            
            Ui = i

            if i != 0:  #To avoid getting out at the initial point when the cells are still fixing caron
                if Csto[i]<Clim2:
                    break
        Uij[j] = Ui 
        
        CstoBox[j,:] = Csto
    
    Tdep = Uij*dT #(s) C storage depleted time
    md = 40 #(m) Depth of the bottom of the mixed layer
    DD = SinkV(Tdep,a1,b1) #(m) Depletion depth
    
    return DD

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Plot
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

rcParams.update({'axes.linewidth':2})
rcParams.update({'xtick.major.width':1.5})
rcParams.update({'ytick.major.width':1.5})

rcParams.update({'lines.markersize': 10})
rcParams.update({'lines.markeredgewidth': 0.5})
rcParams.update({'font.size': 22})
rcParams.update({'lines.linewidth': 3.5})
rcParams.update({'figure.autolayout': True})
rcParams['figure.figsize']=8,6.5
rcParams.update({'figure.facecolor':'W'})
     
rcParams.update({'mathtext.default': 'regular' })
rcParams.update({'patch.edgecolor':'none'})
rcParams.update({'xtick.major.pad': 10})
rcParams.update({'ytick.major.pad': 10})

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

CsFac = arange(0,2+1e-10,0.01)  #(s) time array (177-32)
a = ('',0.0604,0.0604,0.0604) #Coefficient
b = ('',86.8,332.86,578.92) #Coefficient
DD1 = tricho(a[1],b[1],CsFac)
DD2 = tricho(a[2],b[2],CsFac)
DD3 = tricho(a[3],b[3],CsFac)

DD10 = copy(DD1)
DD20 = copy(DD2)
DD30 = copy(DD3)

DD1[(DD30<DD10) & (DD30<DD20)] = DD30[(DD30<DD10) & (DD30<DD20)]
DD1[(DD20<DD10) & (DD20<DD30)] = DD20[(DD20<DD10) & (DD20<DD30)]
DD3[(DD30<DD10) & (DD20<DD10)] = DD10[(DD30<DD10) & (DD20<DD10)]  
DD3[(DD30<DD20) & (DD10<DD20)] = DD20[(DD30<DD20) & (DD10<DD20)]

figure(1) 
plot([0,5],[1000,1000],'--',color='k',linewidth = 1.5)
plot(CsFac,DD2,'-',zorder=5)
fill_between(CsFac, DD1, DD3, where=DD3 >= DD1, facecolor='#9BC2E6',edgecolor = "none",alpha = 0.7,zorder=4)
xlabel('Initial $\mathit{R_{Sto}}$')
ylabel('Depletion depth (m)')
gca().invert_yaxis()
xlim(0,2)
ylim(top=0,bottom=1250)
sf('DD_high_res')
    
show()
    
    
    
    
    