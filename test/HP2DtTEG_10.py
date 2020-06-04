import time

start_time = time.time()

import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import style
import matplotlib.animation as animation
import matplotlib.ticker as ticker
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from tkinter import * 

from pylab import *

import numpy as np
import pandas as pd
import scipy.linalg
import scipy.sparse.linalg as spla

import sys
import pickle
import inspect
import datetime
import os, shutil

timeElapsed1_old = time.time() - start_time

def init_n_print():
    print('\n','\n')
    print('=====================================================================================')
    print('Rui Vieira, 2019')
    print('HP_HT2Dt_TE1Dt')
    print('version 8')
    print('','2D transient simulation of heat transfer with thermoelectric effects ','\n')

    print('Date: ', datetime.datetime.now())
    print('=====================================================================================', '\n')
    
    # deletes previous results
    if not os.listdir('results'):
        return

    delete_ = input('Delete previous results? ')
    if str(delete_) == str(1):
        print('\nDeleting old results!')
        
        folder = 'results'
        for the_file in os.listdir(folder):
            file_path = os.path.join(folder, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                #elif os.path.isdir(file_path): shutil.rmtree(file_path)
            except Exception as e:
                print(e)

def ETA(start_time,time_par,final_time_par,error_par,it):
    
    global timeElapsed1_old
    
    elapsed_time = time.time() - start_time
    
    print('')
    print('==========================================================================')
    print('Writing step: ', str(time_par), '/',str(final_time_par), ' [s]', '\n')
    print('Error = ',str(error_par))

    if write_increment == time_par+1:
        print('in ', str(it),' iterations; ',
              str( -((elapsed_time)*( 1/write_increment)*(time_par-final_time_par))/(60*60) ),
              'hrs remaining (estimated)','\n')
    else:
        print('in ', str(it),' iterations for ',
              str(elapsed_time-timeElapsed1_old),'s; ETA: ',
              str( -((elapsed_time-timeElapsed1_old)*( 1/write_increment)*(time_par-final_time_par))/(60*60) )
              ,'hrs','\n')
    
    
    print('Clock Time = ', str(elapsed_time),' s')
    print('==========================================================================')
    print('')
    
    timeElapsed1_old = elapsed_time

''' PROPERTIES AND HT COEFFICIENTS '''

def Water_Properties(Temp):
    
    ## water properties
    T = Temp
 
    # [Pa.s] LIQUID water viscosity (ZGGRAFOS,1987)
    visc_liq = 3.8208*(10**-2)*(T-252.33)**-1 

    # [kg/m3] LIQUID water volumic mass 1atm (ZGGRAFOS,1987)
    rho_liq = -3.0115*(10**-6)*(T**3)+9.6272*(10**-4)*(T**2)-0.11052*T+1022.4 

    # [W/mK] water thermal conductivity
    k_water = 4.2365*(10**-9)*(T**3)-1.144*(10**-5)*(T**2)+7.1959*(10**-3)*T-0.63262 

    # [J/kgK] water cp (Zografos, 1986)
    Cp_water = 1000*(1.785*(10**-7)*(T**3)-1.9149*(10**-4)*(T**2)+6.7953*(10**-2)*T-3.7559)
    
    properties = [visc_liq, rho_liq, k_water,Cp_water]

    return properties

def R_cooler_water (rho_water, k_water,k_Al, Cp_water,visc_water, dx):
    
    #WATER COOLER MODULE
    #the cooling plates have small channels where the water flows 
    #the correlation for laminar flow can be found in tese: Rui Vieira, 2017


    ## cooler geometry #################################################################
    global m_cooler
    m_cooler = 0.0207916666666667 ## [kg/s] cooling water mass flow rate per plate
    Nc = 70  # number of channels
    channel_widt = 0.001  # channel width [m]
    fin_widt = 0.001  # fin widt [m]
    fin_heigt = 0.01  # fin heigt [m]
    channel_length = 0.13  # channel length [m]
    lateral_length = 0.14  # lateral length [m] comprimento transversal
    Metal_tick = 0.01  # [m] metal thickness joined by the fins

    fin_area_coller = 2 * fin_heigt * channel_length * Nc  # [m2] Af finned area
    fin_perimeter_coller = 2 * (fin_widt + channel_length)  # [m] fin tip perimeter
    non_fin_area = Nc * fin_widt * channel_length  # [m2] Ab area of the base with no fin
    total_HT_area = fin_area_coller + non_fin_area  # [m2] At=Ab+Af total heat tranfer area
    flow_area = fin_widt * fin_heigt  # [m2] single channel flow area
    Dh_cooler = (4 * flow_area) / (2 * (fin_widt + fin_heigt)
                                )  # [m] hidraulic diameter 4A/P
    D_charact_cooler = (flow_area)**0.5  # [m] characteristic lenght A**1/2
    ###################################################################################


    #geometry definition
    fin_area_coller=2*fin_heigt*channel_length*Nc # [m2] Af finned area
    fin_perimeter_coller=2*(fin_widt+channel_length) # [m] fin tip perimeter
    non_fin_area=Nc*fin_widt*channel_length # [m2] Ab area of the base with no fin
    total_HT_area=fin_area_coller+non_fin_area # [m2] At=Ab+Af total heat tranfer area
    flow_area=fin_widt*fin_heigt # [m2] single channel flow area
    D_charact_cooler=(flow_area)**0.5 # [m] characteristic lenght A**1/2

    #Cooler convection model constants, check tese Rui Vieira 2017

    Pr_water =(Cp_water *visc_water )/k_water 

    C1=3.24 # Uniform heat power
    C2=1.5
    C3=0.409 # Uniform heat power
    C4=2
    Par_forma=0.1 # Par de forma (>=90)
    ee=fin_widt/fin_heigt # channel aspect ratio e=a/b of a rectangule

    fPr =0.564/((1+(1.664*Pr_water **(1/6))**(9/2))**(2/9)) # Uniform heat power
    m_coef_cool =2.27+1.65*Pr_water **(1/3) # m coefficient 


    # Flow conditions
    #            m_cooler_channel=m_cooler/Nc # [kg/s] cooling water mass flow per channel

    # [m3/s] cooling water volumic flow per channel
    Q_water_channel =(m_cooler/rho_water )/Nc 

    # [m/s] water velocity in the channels
    u_water =Q_water_channel /flow_area 


    #water Reynolds number
    Re_water =(rho_water *u_water *D_charact_cooler)/visc_water  

    L_coef =0.058*D_charact_cooler*Re_water  # [m] L' entry region
    Z_plus =((L_coef /channel_length)**-1)/(Re_water ) # Z+ coefficient
    Z_star =((channel_length/D_charact_cooler)/(Re_water *Pr_water )) # Z* coefficient

    # fRe(A**0.5) coefficient
    fRe =((12/((ee**0.5)*(1+ee)*(1-((192*ee)/(np.pi**5))*np.tanh(np.pi/(2*ee)))))**2+(3.44/Z_plus **0.5)**2)**0.5 

    #Nusselt number
    Nu_water =((((C4*fPr )/(Z_star **0.5))**m_coef_cool )+((((C2*C3*(fRe /Z_star )**(1/3))**5)+(C1*(fRe /(8*(np.pi**0.5)*(ee**Par_forma))))**5)**(m_coef_cool /5)))**(1/m_coef_cool ) 
    # [Wm2/K] heat transfer coefficient
    h_water =(Nu_water *k_water )/D_charact_cooler 

    # m fin efficiency coeficient
    m_coef_water = ((h_water *fin_perimeter_coller)/(k_Al*non_fin_area))**0.5 
    # fin efficiency
    fin_eff_water =(np.tanh(m_coef_water *fin_heigt))/(m_coef_water *fin_heigt) 
    # group fin efficiency
    fin_eff_water_group =1-(fin_area_coller/total_HT_area)*(1-fin_eff_water ) 


    #Cooler thermal resistance [K/W]
    R_cooler = (1/(h_water *(total_HT_area/(channel_length/dx))*fin_eff_water_group ))/2
    
    return R_cooler

def Prop_air(Temp):
    
    T=Temp

    #AIR PROPERTIES IN S.I. (Zografos, 1986) (equações na tese Rui Vieira, 2017 no ultimo anexo)
    #Calculates the air resistance for Th 

    #viscosity as a function of hot inlet temperature [Pa.s]
    air_viscosity=((2.5914*10**(-15))*T**3)-((1.4346*10**(-11))*T**2)+((5.0523*10**(-8))*T)+4.113*10**(-6) #[Pa.s]
            
    #density as a function of hot inlet temperature [kg/m3]
    air_density=101325/((8314.4598/28.97)*T) #[kg/m3]
            
    #air conductivity [W/mK]
    air_conductivity=((0.000000000015207)*T**3) - ((0.000000048574)*T**2 )+( (0.00010184)*T)- 0.00039333
            
    #Specific Heat [J/kg.K]
    air_specific_heat=(((1.3864*10**-13)*T**4)-((6.4747*10**-10)*T**3)+((1.0234*10**-6)*T**2)-(4.3282*10**-4)*T+1.0613)*1000
    
    properties = [air_viscosity,air_density,air_conductivity,air_specific_heat]
    
    return properties

def R_Exhaust (air_viscosity,air_conductivity,air_specific_heat,m_air, N_HE):
    

    #geometry definition

    s=0.003   # [m] fin spacing
    h=0.0189 # [m] fin heigth
    tick=0.0005 # [m] fin width
    l=0.0531 # [m] offset length

    number_of_channels=100 # [-] number of channels 
    number_of_offsets=10 # [-] number of offsets
    
    global L_lat
    L_lat = number_of_channels*(s+tick)

    N_slices = N_HE #1+HE_length/dx # number of slices 

    k_ss=20 # [W/mK] thermal conductivity for the offset fins
    #20 W/mK average value for the thermal conductivity of the stainless steel

    #estas areas são para o PC todo, para os calculos das alhetas,ao
    # dividir o PC a meio, é necessário dividir por 2!
    total_area=2*(number_of_channels*number_of_offsets)*(s*l+h*l) #[m2] total heat transfer area 

    #hidraulic diameter for offset fins

    alfa_geo=s/h
    sigma_geo=tick/l
    gama_geo=tick/s

    D_hidraulic=(4*s*h*l)/(2*(s*l+h*l+tick*h)+tick*s) # [m] hydrauslic dyameter

    #Dynamic and themal parameters 
    #Reynolds Number
    Re=(4/np.pi)*(((m_air/number_of_channels))/(air_viscosity*D_hidraulic))

    #Colburn factor (j)
    j=(0.6522*(Re**-0.5403)*(alfa_geo**-0.1541)*(sigma_geo**0.1499)*(gama_geo**-0.0678))*(1+(5.269*10**-5)*(Re**1.34)*(alfa_geo**0.504)*(sigma_geo**0.456)*(gama_geo**-1.055))**0.1

    #Frication factor (f)
    #f=(9.6243*(Re**(-0.7422))*(alfa_geo**(-0.1856))*(sigma_geo**(0.3053))*(gama_geo**(-0.2659)))*(1+(7.669*(10**-8))*(Re**(4.429))*(alfa_geo**(0.92))*(sigma_geo**(3.767))*(gama_geo**(0.236)))**(0.1)
    #Performance factor j/f
    #performance_factor=j/f

    #Prandtl number
    Pr=(air_specific_heat*air_viscosity)/air_conductivity

    #Nusselt number
    Nu=j*Re*(Pr)**(1/3)

    #Heat tranfer coefficient [W/m2K]
    h_air=(air_conductivity*Nu)/D_hidraulic

    #air velocity [m/s]
    #u_air(x)=((m(1,1))/number_of_channels)/(air_density*h*s)
    #Pressure drop (dP) [Pa]
    #dP(x)=(0.5*air_density*(u_air(x)*u_air(x))*dx*f(x))/D_hidraulic

    #Fin efficiency coefficient

    m_fin=((h_air*2*(tick+l))/(k_ss*tick*l))**0.5

    #Fin efficiency 
    fin_eff=(np.tanh(m_fin*h/2))/(m_fin*h/2)

    #Fin group efficiency  
    fin_eff_group=1-((2*h/2)/(2*h/2+s))*(1-fin_eff)

    #Exhaust heat exchanger thermal resistance [K/W]
    R_air=(1/(((total_area/2)/N_slices)*h_air*fin_eff_group))
    
    return R_air

def Prop_TEG_Cp(self):
    
    T = self
    return (36.858*(10**-3)*T + 117.02 - 161744*T**-2)/(800.761*10**-3)

def Prop_TEG_N(self):

    select = 1


    T = self # all temperatures must be in [K]



    if select == 1:
        
        ## Hi-z data
        
        # Seebeck coeficient [V/K]
        alfa_Bi2Te3N = 0.00007423215 - 0.0000015018 * T + 0.0000000029361 * T**2 - 0.000000000002499 * T**3 + 0.000000000000001361 * T**4 # [V/K]
        
        # electrical resistivity [ohm.m]
        rho_Bi2Te3N = -0.00195922 + 0.00001791526 * T - 0.00000003818 * T**2 + 0.000000000049186 * T**3 - 0.0000000000000298 * T**4# ohm.cm
        rho_Bi2Te3N = rho_Bi2Te3N/100 # ohm.m
        
        # thermal conductivity [W/mK]
        k_Bi2Te3N = 1.425785 + 0.006514882 * T - 0.00005162 * T**2 + 0.00000011246 * T**3 - 0.000000000076 * T**4 # W/mK
        
    else:
        
        # Cyprus data
        
        alfa_Bi2Te3N = -(3.8273*(10**-11)*(T**5) - 9.3841*(10**-8)*(T**4) + 9.1048*(10**-5)*(T**3) - 4.3783*(10**-2)*(T**2) + 1.0608*(10**+1)*T - 8.7621*(10**+2))*(10**-6) #R2 = 9.981710**-01
        
        rho_Bi2Te3N = 1.3292*(10**-15)*(T**4) - 2.8350*(10**-12)*(T**3) + 2.3015*(10**-9)*(T**2) - 8.3722*(10**-7)*T + 1.7220*(10**-4) #R2 = 9.9562E-01
        
        k_Bi2Te3N = -7.3633*(10**-16)*(T**6) + 2.4093*(10**-12)*(T**5) - 3.1602*(10**-9)*(T**4) + 2.1370*(10**-6)*(T**3) - 7.8786*(10**-4)*(T**2) + 1.5046*(10**-1)*T - 1.1021*(10**+1) #R2 = 9.9508E-01
        
    
    Prop_N = [alfa_Bi2Te3N, rho_Bi2Te3N, k_Bi2Te3N]
    
    return Prop_N

def Prop_TEG_P(self):
    
    
    T = self # all temperatures must be in [K]

    select = 1

    if select == 1:
        ## Hi-z data
        
        # Seebeck coeficient [V/K]
        alfa_Bi2Te3P = -0.0002559037 + 0.0000023184 * T - 0.000000003181 * T**2 + 0.0000000000009173 * T**3 - 0.000000000000000488 * T**4 # [V/K]
        
        # electrical resistivity [ohm.m]
        rho_Bi2Te3P = -0.002849603 + 0.00001967684 * T - 0.00000003317 * T**2 + 0.000000000034733 * T**3 - 0.000000000000019 * T**4 # ohm.cm
        rho_Bi2Te3P = rho_Bi2Te3P/100 # ohm.m
        
        # thermal conductivity [W/mK]
        k_Bi2Te3P = 6.9245746 - 0.05118914 * T + 0.000199588 * T**2 - 0.0000003891 * T**3 + 0.00000000030382 * T**4 # W/mK
        
    else:
        
        # Cyprus data
        
        alfa_Bi2Te3P = (3.8273*(10**-11)*(T**5) - 9.3841*(10**-8)*(T**4) + 9.1048*(10**-5)*(T**3) - 4.3783*(10**-2)*(T**2) + 1.0608*(10**+1)*T - 8.7621*(10**+2))*(10**-6) #R² = 9.981710**-01
        
        rho_Bi2Te3P = 1.3292*(10**-15)*(T**4) - 2.8350*(10**-12)*(T**3) + 2.3015*(10**-9)*(T**2) - 8.3722*(10**-7)*T + 1.7220*(10**-4) #R² = 9.9562E-01
        
        
        k_Bi2Te3P = -7.3633*(10**-16)*(T**6) + 2.4093*(10**-12)*(T**5) - 3.1602*(10**-9)*(T**4) + 2.1370*(10**-6)*(T**3) - 7.8786*(10**-4)*(T**2) + 1.5046*(10**-1)*T - 1.1021*(10**+1) #R² = 9.9508E-01


    Prop_P = [alfa_Bi2Te3P, rho_Bi2Te3P, k_Bi2Te3P] 
    return Prop_P

''' END PROPERTIES '''

def fun_TEG(Q1TEG,Q2TEG,u_TEG_n_P,u_TEG_n_N,u_TEG_N,u_TEG_P,dt,I_current):

    ## TEG Geometry
    global A_crosssec_N
    global A_crosssec_P
    global leg_side_N
    global leg_side_P

    # leg geometry
    n_leg = 98 # number of legs per TEG module
    leg_side_N = 3.8/1000 # [m] A_leg dimensions
    leg_side_P = 3.8/1000 # [m] A_leg dimensions
    A_crosssec_N = leg_side_N**2  # [m2]
    A_crosssec_P = leg_side_P**2  # [m2]
    leg_length = 0.16/100 # [m] TEG length in y direction

    # Al connector geometry
    L_Al_connector = 0 #0.15/100 # [m] connector thickkness
    #width_connector = 0.5/100# [m] connector width

    # total TEG module length
    L_total = leg_length + 2*L_Al_connector
    # module_widt_1 = 5.81/100 # TEG module widt 1
    # module_widt_2 = 5.81/100 # TEG module widt 2

    # resistance properties
    # R_connector = 2*10**-5 # [ohm] connector resistivity (simplificação...)
    R_contact = 0.0012 # [ohm] contact resistance

    # discretization
    global delta_x
    delta_x = (0.16/100)/3 #2*0.2/1000 # [m]

    global N_legs
    global N_connector

    N_legs = leg_length/delta_x # number of slices in the legs
    N_connector = round(L_Al_connector/delta_x) # number of slices in each connector
    global Nx_TEG
    Nx_TEG = int(N_legs+2*N_connector) # total number of slices

    rhoDense_TEG = 7700 # kg/m3

    Tamb = 20 + 273

    global Area_TEG_N
    global Area_TEG_P


    Area_TEG_N = leg_side_N*leg_side_N
    Area_TEG_P = leg_side_N*leg_side_N


    A_TEG_N = np.zeros([Nx_TEG+1,Nx_TEG+1]) # numerical matrix
    A_TEG_P = np.zeros([Nx_TEG+1,Nx_TEG+1]) # numerical matrix

    b_TEG_N = np.zeros([Nx_TEG+1]) # numerical matrix
    b_TEG_P = np.zeros([Nx_TEG+1]) # numerical matrix


    Tau_P = 0.00027 # [V/K] thomson coefficient
    Tau_N = -0.000156 # [V/K] thomson coefficient

    # Aluminum properties
    k_Al = 200
    cp_metal_HP = 910
    density_Al = 2712

    Temp_N = Tamb # ref temp is the Tamb
    Temp_P = Tamb # ref temp is the Tamb
    cp_Bi2Te3 = Prop_TEG_Cp (Temp_P)
    [alfa_Bi2Te3N, rho_Bi2Te3N, k_Bi2Te3N] = Prop_TEG_N(Temp_N)
    [alfa_Bi2Te3P, rho_Bi2Te3P, k_Bi2Te3P] = Prop_TEG_P(Temp_P)
    alfa_N_ref = alfa_Bi2Te3N
    alfa_P_ref = alfa_Bi2Te3P

    alfa_N_now = np.zeros([Nx_TEG+1])
    alfa_P_now = np.zeros([Nx_TEG+1])

    k_N_now = np.zeros([Nx_TEG+1])
    k_P_now = np.zeros([Nx_TEG+1])

    Cp_N_now = np.zeros([Nx_TEG+1])
    Cp_P_now = np.zeros([Nx_TEG+1])

    rho_N_before = np.zeros([Nx_TEG+1])
    rho_P_before = np.zeros([Nx_TEG+1])
    R_y = np.zeros([Nx_TEG+1])
    Vemf = np.zeros([Nx_TEG+1])

    density_before = np.zeros([Nx_TEG+1])

    if Q1TEG == 0: # delta_x 1st output condition 
        return delta_x

    # BC condition constant heat flux
    QN2in = k_Al*Area_TEG_N*(-u_TEG_N[Nx_TEG] + u_TEG_N[Nx_TEG-1])/delta_x
    QP2in = k_Al*Area_TEG_P*(-u_TEG_P[Nx_TEG] + u_TEG_P[Nx_TEG-1])/delta_x
    
    ff = 1
    

    VTEG = delta_x*dx*L_lat
    VTEG1 = VTEG #dy*dx*L_lat

    TQ1N = ff*(Q1TEG*dt)/(density_Al*VTEG1*cp_metal_HP)
    TQ2N = ff*( (QN2in - Q2TEG) *dt)/(density_Al*VTEG*cp_metal_HP)
    TQ1P = ff*(Q1TEG*dt)/(density_Al*VTEG1*cp_metal_HP)
    TQ2P = ff*( (QP2in - Q2TEG) *dt)/(density_Al*VTEG*cp_metal_HP)

    for y in range(0, Nx_TEG+1): # TEG delta_x iterator
        
        ## t time
        if y >= N_connector and y <= N_connector + N_legs: # material properties
            
            Temp_N = u_TEG_n_N[y] #u_TEG_n_N[Nx_HE-1][y]
            Temp_P = u_TEG_n_P[y] #u_TEG_n_P[Nx_HE-1][y]
            Temp = (Temp_N + Temp_P)/2
            
            [alfa_Bi2Te3N, rho_Bi2Te3N, k_Bi2Te3N] = Prop_TEG_N(Temp_N)
            [alfa_Bi2Te3P, rho_Bi2Te3P, k_Bi2Te3P] = Prop_TEG_P(Temp_P)
            
            rho_N_before[y] = rho_Bi2Te3N
            rho_P_before[y] = rho_Bi2Te3P
            
            density_before[y] = 7700
            
        else:
            
            rho_N_before[y] = 2.82e-08
            rho_P_before[y] = 2.82e-08
            
            density_before[y] = 2712
            
        # end material properties
        
        ## t+1 time
        if y >= N_connector and y <= N_connector + N_legs: # material properties
            
            Temp_N = u_TEG_N[y]
            Temp_P = u_TEG_P[y]
            Temp = (Temp_N + Temp_P)/2
            
            cp_Bi2Te3 = Prop_TEG_Cp (Temp)
            [alfa_Bi2Te3N, rho_Bi2Te3N, k_Bi2Te3N] = Prop_TEG_N(Temp_N)
            [alfa_Bi2Te3P, rho_Bi2Te3P, k_Bi2Te3P] = Prop_TEG_P(Temp_P)
            
            alfa_N_now[y] = alfa_Bi2Te3N
            alfa_P_now[y] = alfa_Bi2Te3P
            
            k_N_now[y] = k_Bi2Te3N
            k_P_now[y] = k_Bi2Te3P
            
            Cp_N_now[y] = cp_Bi2Te3
            Cp_P_now[y] = cp_Bi2Te3
            
        else:
            alfa_N_now[y] = 3.5*10**-6
            alfa_P_now[y] = 3.5*10**-6
            
            k_N_now[y] = k_Al
            k_P_now[y] = k_Al
            
            Cp_N_now[y] = cp_metal_HP
            Cp_P_now[y] = cp_metal_HP
            
        # end material properties
    
    for y in range(1,Nx_TEG+1):
        Vemf[y] = abs( alfa_N_now[y]*(u_TEG_N[y-1] - u_TEG_N[y]) ) + abs( alfa_P_now[y]*(u_TEG_P[y-1] - u_TEG_P[y]) )

    for y in range(1, Nx_TEG):    

        R_y[y] = (rho_P_before[y]*delta_x)/A_crosssec_P + (rho_N_before[y]*delta_x)/A_crosssec_N

        
        #alpha_PN_0 = (alfa_N_ref - alfa_P_ref) + (Tau_P - Tau_N)*np.log( u_TEG_N[0]/Tamb )
        #alpha_PN_L = (alfa_N_ref - alfa_P_ref) + (Tau_P - Tau_N)*np.log( ((u_TEG_N[Nx_TEG] + u_TEG_P[Nx_TEG])/2)/Tamb )
                
        
        E_T_N = (Tau_N * I_current * dt) / (Cp_N_now[y] * density_before[y] * Area_TEG_N * delta_x)
        E_T_P = (Tau_P * I_current * dt) / (Cp_P_now[y] * density_before[y] * Area_TEG_P * delta_x)
        
        E_N = (rho_N_before[y]*(I_current**2)*dt)/(Cp_N_now[y]*density_before[y]*(Area_TEG_N**2))
        E_P = (rho_P_before[y]*(I_current**2)*dt)/(Cp_P_now[y]*density_before[y]*(Area_TEG_P**2))
        
        F_N = ( (k_N_now[y]/(density_before[y]*Cp_N_now[y]))*dt )/(delta_x**2)
        F_P = ( (k_P_now[y]/(density_before[y]*Cp_P_now[y]))*dt )/(delta_x**2)
        
        A_TEG_N[y,y-1] = - E_T_N - F_N
        A_TEG_N[y,y+1] = E_T_N - F_N
        A_TEG_N[y,y] = 1 + 2*F_N
        
        b_TEG_N[y] = u_TEG_n_N[y] + E_N #u_TEG_n_N[Nx_HE-1][y] + E_N
        
        A_TEG_P[y,y-1] = - E_T_P - F_P
        A_TEG_P[y,y+1] = E_T_P - F_P
        A_TEG_P[y,y] = 1 + 2*F_P
        
        b_TEG_P[y] = u_TEG_n_P[y] + E_P #u_TEG_n_P[Nx_HE-1][y] + E_P
    
    R = sum(R_y)
    R0 = R

    
    A_TEG_N[0,0] = A_TEG_N[Nx_TEG,Nx_TEG] = 1
    b_TEG_N[0] = u_TEG_n_N[0] + E_N + TQ1N #u_TEG_n_N[Nx_HE-1][0] + E_N + TQ1N
    b_TEG_N[Nx_TEG] = u_TEG_n_N[Nx_TEG] + E_N + TQ2N # u_TEG_n_N[Nx_HE-1][Nx_TEG + 1] + E_N + TQ2N
    
    A_TEG_P[0,0] = A_TEG_P[Nx_TEG,Nx_TEG] = 1
    b_TEG_P[0] = u_TEG_n_P[0] + E_P + TQ1P #u_TEG_n_P[Nx_HE-1][0] + E_P + TQ1P
    b_TEG_P[Nx_TEG] = u_TEG_n_P[Nx_TEG] + E_P + TQ2P # u_TEG_n_P[Nx_HE-1][Nx_TEG + 1] + E_P + TQ2P
    
            
    c_N = spla.bicgstab(A_TEG_N, b_TEG_N, tol=1e-6)
    c_P = spla.bicgstab(A_TEG_P, b_TEG_P, tol=1e-6)

    # Fill u with vector c
    for y in range(0, Nx_TEG+1):
        u_TEG_N[y] = c_N[0][y].copy()
        u_TEG_P[y] = c_P[0][y].copy()
        
    # output 
    V = sum(Vemf) #- I_current*R
    return [u_TEG_N, u_TEG_P, R, V,delta_x,Area_TEG_N]

def Solver():

    global write_increment
    write_increment = 0.5
    write = write_increment
    dt = 0.01

    debug = 0
    select = str(2) #str(1) #input('AE [1] or WLTP [2]? ')
    
    if select == str(1):
        data = pd.read_csv(
            "AE.csv",
            names=['Exhaust Power [kW]', 'Exhaust Mass [kg/s]', 'T_EVO [K]'])
        N_time_steps = 1147
    elif select == str(2):
        data = pd.read_csv(
            "WLTP.csv",
            names=['Exhaust Power [kW]', 'Exhaust Mass [kg/s]', 'T_EVO [K]'])
        N_time_steps = 1800

    if debug == 0:
        data['Exhaust Mass [kg/s]'] *= 1e-3
        M_air = data['Exhaust Mass [kg/s]']
        Th = np.array(data['T_EVO [K]'])

    if debug == 1:
        
        N_time_steps = 2500
        M_air = np.zeros([N_time_steps + 1])
        Th = np.zeros([N_time_steps + 1])
        M_air[:] = 803/3600
        Th[:int(N_time_steps/3)]=1000
        Th[int(N_time_steps/3):]=273.15 + 20     

        # Aluminum properties
    rho_Al = 2700/2; Cp_Al = 910/2; k_Al = 204 #Cp_Al = 910
    cp_metal_HP = Cp_Al

    Tamb = 273.15 + 20

    Lx = 0.6; Ly = 0.15
    Nx = 60
    Ny = 15
    Nt = N_time_steps/dt
    global dx
    dx = Lx/Nx
    dy = Ly/Ny
    L_HP = 0.10
    n_HP = int(L_HP/dy) #- 1

    theta = 1 # fully implicit

    if theta == 1:
        print('\n','Fully implicit formulation','\n')
    elif theta == 0.5:
        print('\n','Crank-Nicholson formulation','\n')
    elif theta == 0:
        print('\n','Implicit formulation','\n')

    print( 'Delta t = ', dt,' [s]', '\n')

    print('\n','initializing... ')

    ## intial parameters
    number_of_coolers = 5  # number of coolers
    T_TEG_MAX = 230 + 273.15  # [K] maximum allowed temperature on the hot TEG face
    Twater_in = Tamb  # [k] water inlet temperature

    ## AIR HEAT EXCHANGER MODULE (OFFSET FINS) -- inlet conditions ####################
    N_slices_cooler = Nx / number_of_coolers  # number of dx slices per cooler length

    k_ss = 20  # [W/mK] thermal conductivity for the offset fins
    #20 W/mK average value for the thermal conductivity of the stainless steel

    delta_x = fun_TEG(0,0,0,0,0,0,0,0) # to ouput the delta_x

    # DECLARATIONS
    alpha = np.zeros([Ny+1])
    Fy = np.zeros([Ny+1])
    Fx = np.zeros(Ny+1)
    C = np.zeros([Ny+1])
    u   = np.zeros([Nx+1, Ny+1])      # unknown u at new time level
    u_n = np.zeros([Nx+1, Ny+1])      # u at the previous time level
    Q = np.zeros([Nx+1, Ny+1])
    N = (Nx+1)*(Ny+1)  # no of unknowns
    uoo = np.zeros([Nx+1])
    ucc = np.zeros([Nx+1])
    R_air = np.zeros([Nx+1])
    R_water = np.zeros([Nx+1])
    Qoo = np.zeros([Nx+1])
    Qcc = np.zeros([Nx+1])
    uLoo = np.zeros([Nx+1])
    uLcc = np.zeros([Nx+1])
    U = np.zeros([Nx+1])
    Rtot = np.zeros([Nx+1])

    NNTEG = Nx_TEG + 1

    u_TEG_N_main = np.zeros([Nx+1,NNTEG])
    u_TEG_P_main = np.zeros([Nx+1,NNTEG])
    u_TEG_N_main[:,:] = Tamb
    u_TEG_P_main[:,:] = Tamb
    u_TEG_N_old = u_TEG_N_main.copy()
    u_TEG_P_old = u_TEG_P_main.copy()
    u_TEG_N_n = u_TEG_N_main.copy()
    u_TEG_P_n = u_TEG_N_main.copy()

    Q1 = np.zeros([Nx+1])
    Q2 = np.zeros([Nx+1])
    Qin = np.zeros([Nx+1])
    Q_excess = np.zeros([Nx+1])

    u_TEG_n_N = np.zeros([NNTEG])
    u_TEG_n_P = np.zeros([NNTEG])

    u_TEG_n_N[:] = Tamb
    u_TEG_n_P[:] = Tamb

    A_HP = np.zeros([(Nx+1)*(n_HP+1), (Nx+1)*(n_HP+1)])
    b_HP = np.zeros([(Nx+1)*(n_HP+1)])
    #A_cool = np.zeros([(Nx+1)*(Ny - n_HP+1),(Nx+1)*(Ny - n_HP+1)])
    #b_cool = np.zeros([(Nx+1)*(Ny - n_HP+1)]) # np.zeros([int(Ny + 1 - n_HP)])

    A_cool = np.zeros([(Nx+1)*(Ny - n_HP),(Nx+1)*(Ny - n_HP)])
    b_cool = np.zeros([(Nx+1)*(Ny - n_HP)]) # np.zeros([int(Ny + 1 - n_HP)])

    u_n[:,:] = Tamb
    u[:,:] = Tamb

    if debug == 1:
        u[:,0:n_HP] = 273 +230
        u_n[:,0:n_HP] = 273 +230

    u_old = u_n.copy()

    Ix = range(0, Nx+1)
    Iy = range(0, Ny+1)
    IyHP = range(0, int(n_HP + 1))
    Icool = range(0,Ny - n_HP) #Icool = range(0,Ny - n_HP + 1)

    empty = 1 # void volume on the HP

    for j in Iy:
        if j > 0 and j < n_HP - 1:
            alpha[j] = empty * k_Al/(rho_Al*Cp_Al)
            C[j] = empty * rho_Al*Cp_Al
        else:
            alpha[j] = k_Al/(rho_Al*Cp_Al)
            C[j] = rho_Al*Cp_Al

    for j in Iy:
        Fy[j] = (alpha[j]*dt)/dy**2
        Fx[j] = (alpha[j]*dt)/dx**2

    print('Stabillity coef, Fy.max() + Fx.max() < 0.5 = ', Fy.max() + Fx.max())
    if Fy.max() + Fx.max() > 0.5:
        print('WARNING not stable, proced?')
        time.sleep(3)


    m = lambda i, j: j*(Nx+1) + i

        
    ''' cool part '''
    j = 0
    for i in Ix:
        p = m(i,j);  A_cool[p, p] = 1
    # Loop over all internal mesh points in y diretion
    # and all mesh points in x direction
    for j in Icool[1:-1]:
        i = 0;  p = m(i,j);  A_cool[p, p] = 1 # boundary
        for i in Ix[1:-1]:                 # interior points
            p = m(i,j)
            A_cool[p, m(i,j-1)] = - theta*Fy[j]
            A_cool[p, m(i-1,j)] = - theta*Fx[j]
            A_cool[p, p]        = 1 + 2*theta*(Fx[j]+Fy[j])
            A_cool[p, m(i+1,j)] = - theta*Fx[j]
            A_cool[p, m(i,j+1)] = - theta*Fy[j]
        i = Nx;  p = m(i,j);  A_cool[p, p] = 1  # boundary
    # Equations corresponding to j=Ny, i=0,1,... (u known)
    j = Ny - (n_HP + 1)
    for i in Ix:
        p = m(i,j);  A_cool[p, p] = 1

    N_TEG_legs = 1 #8 * 20
    fff = 1
    I_current = 0
    I_current_old = 0

    tol = 10**-4; max_iter = 200
    Erro_tot = 1
    #Thot = 800; Tcold = Tamb
    converged = False

    n = 1
    print('=====================================================================================')
    print('=====================================================================================')
    print('\n','solving')

    for t in range(int(Nt)):

        it = 0
        converged = False
        
        ## cycle information
        if t*dt > N_time_steps-1:
            uoo[0] = Th[int(np.ceil(t*dt))]
            m_air = M_air[int(np.ceil(t*dt))]/2
        else:
            uoo[0] = (int(np.ceil(t*dt)) - t*dt)*Th[int(np.ceil(t*dt))] + \
            (1-(int(np.ceil(t*dt)) - t*dt))*Th[int(np.ceil(t*dt))+1]
            m_air = ( (int(np.ceil(t*dt)) - t*dt)*M_air[int(np.ceil(t*dt))] + \
                    (1-(int(np.ceil(t*dt)) - t*dt))*M_air[int(np.ceil(t*dt))+1] )/2
        
        ''' heat spreading '''
        for i in range(0, Nx+1):
            if Q_excess[i] == 0:
                Q[i,n_HP] = sum(Q_excess)
            elif sum(Q_excess) > 0 and i == Nx+1:
                Q_reject = sum(Q_excess)
        ''' heat spreading '''

        while not converged:
                    
                    
            for i in range(0, Nx+1):
                if i == Nx or i == N_slices_cooler or i == (2*N_slices_cooler) or i == (3*N_slices_cooler) or i == (4*N_slices_cooler):
                    ucc[i] = Twater_in # cooler inlet
                
            ## BC convection
            for i in range(0, Nx): #Ix[:-1]:
                ## HE conv
                [air_viscosity,air_density,air_conductivity,air_specific_heat] = Prop_air(uoo[i])
                R_air[i] = R_Exhaust(air_viscosity,air_conductivity,air_specific_heat,m_air, Nx+1)
                uoo[i+1] = ((uoo[i]-u[i][0])/R_air[i] - m_air*air_specific_heat*uoo[i])/(- m_air*air_specific_heat)
                Qoo[i] = (uoo[i]-uoo[i+1])*m_air*air_specific_heat #(uoo[i] - u[i][0])/R_air[i] 
                uLoo[i] = u_n[i][0] + (Qoo[i]*dt)/(rho_Al*dx*dy*L_lat*Cp_Al) - ( ( (L_lat*dx)*k_Al*(u[i,0] - u[i,1])/dy ) *dt ) / (rho_Al*dx*dy*L_lat*Cp_Al)
                
                ## cooler conv
                [visc_liq, rho_liq, k_water,Cp_water] = Water_Properties(ucc[Nx-i])
                R_water[Nx-i] = R_cooler_water(rho_liq, k_water,k_Al, Cp_water,visc_liq, dx)
                
                if i != N_slices_cooler-1 and i != (2*N_slices_cooler-1) and i != (3*N_slices_cooler-1) and i != (4*N_slices_cooler-1): # and i != (5*N_slices_cooler-1): # i != Nx-1 and 
                    ucc[Nx-1-i] = ((u[Nx-i][Ny] - ucc[Nx-i])/R_water[Nx-i] + 2*m_cooler*Cp_water*ucc[Nx-i])/(2*m_cooler*Cp_water)
                    Qcc[Nx-i] = (ucc[Nx-1-i]-ucc[Nx-i])*2*m_cooler*Cp_water
                    uLcc[Nx-i] = u_n[Nx-i][Ny] - (Qcc[Nx-i]*dt)/(rho_Al*dx*dy*L_lat*Cp_Al) + ( ( (L_lat*dx)*k_Al*(u[Nx-i,Ny-1] - u[Nx-i,Ny])/dy ) *dt ) / (rho_Al*dx*dy*L_lat*Cp_Al)
            ## outflows  
            i = Nx
            [air_viscosity,air_density,air_conductivity,air_specific_heat] = Prop_air(uoo[i-1])
            R_air[i] = R_Exhaust(air_viscosity,air_conductivity,air_specific_heat,m_air, Nx+1)
            uoo_out = ((uoo[i]-u[i][0])/R_air[i] - m_air*air_specific_heat*uoo[i])/(- m_air*air_specific_heat)
            Qoo[i] = (uoo[i]-uoo_out)*m_air*air_specific_heat
            uLoo[i] = u_n[i][0] + (Qoo[i]*dt)/(rho_Al*dx*dy*L_lat*Cp_Al) - ( ( (L_lat*dx)*k_Al*(u[i,0] - u[i,1])/dy ) *dt ) / (rho_Al*dx*dy*L_lat*Cp_Al)
            
            [visc_liq, rho_liq, k_water,Cp_water] = Water_Properties(ucc[Nx-i])
            R_water[Nx-i] = R_cooler_water(rho_liq, k_water,k_Al, Cp_water,visc_liq, dx)
            ucc_out = ((u[Nx-i][Ny] - ucc[Nx-i])/R_water[Nx-i] + 2*m_cooler*Cp_water*ucc[Nx-i])/(2*m_cooler*Cp_water)
            Qcc[Nx-i] = (ucc_out-ucc[Nx-i])*2*m_cooler*Cp_water
            uLcc[Nx-i] = u_n[Nx-i][Ny] - (Qcc[Nx-i]*dt)/(rho_Al*dx*dy*L_lat*Cp_Al) + ( ( (L_lat*dx)*k_Al*(u[Nx-i,Ny-1] - u[Nx-i,Ny])/dy ) *dt) / (rho_Al*dx*dy*L_lat*Cp_Al)

            for i in [int(N_slices_cooler-1), int(2*N_slices_cooler-1), int(3*N_slices_cooler-1), int(4*N_slices_cooler-1), int(5*N_slices_cooler-1)]:
                [visc_liq, rho_liq, k_water,Cp_water] = Water_Properties(ucc[(Nx-i)])
                R_water[(Nx-i)] = R_cooler_water(rho_liq, k_water,k_Al, Cp_water,visc_liq, dx)
                ucc_out = ((u[Nx-i][Ny] - ucc[Nx-i])/R_water[Nx-i] + 2*m_cooler*Cp_water*ucc[Nx-i])/(2*m_cooler*Cp_water)
                Qcc[Nx-i] = (ucc_out-ucc[Nx-i])*2*m_cooler*Cp_water
                uLcc[Nx-i] = u_n[Nx-i][Ny] - (Qcc[Nx-i]*dt)/(rho_Al*dx*dy*L_lat*Cp_Al) + ( ( (L_lat*dx)*k_Al*(u[Nx-i,Ny-1] - u[Nx-i,Ny])/dy ) *dt) / (rho_Al*dx*dy*L_lat*Cp_Al)
                
                    
            ## End BC convection

                    
            if it > max_iter or np.isnan(Erro_tot) == 1:# or Erro_tot > 10:
                
                print('\n','=====================================================================================')
                print('================================= DOES NOT CONVERGE =================================')
                print('=====================================================================================')
                print('\n','time: ',t*dt,'time step: ', t,'[s], Error: ',Erro_tot,' with ',it, ' iterations')
                print('Error: ',Erro_tot)
                if Erro_tot < tol*10:
                    print('break the cycle, error = ', Erro_tot)
                    break
                
                sys.exit('does not converge')
                
            '''================================= HP section =================================''' 
            
            ''' HP part '''
            
            for j in Iy:
                if j > 0 and j < n_HP-1:
                    alpha[j] = 0.25 * k_Al/(rho_Al*Cp_Al)
                    C[j] = rho_Al/4*Cp_Al/4
                else:
                    pass
                
            for j in Iy:
                Fy[j] = (alpha[j]*dt)/dy**2
                Fx[j] = (alpha[j]*dt)/dx**2

            
            ## if A is time independent!!
            # Equation corresponding to mesh point (i,j) has number
            # j*(Nx+1)+i and will contribute to row j*(Nx+1)+i
            # in the matrix.

            # Equations corresponding to j=0, i=0,1,... (u known)
            j = 0
            for i in Ix:
                p = m(i,j);  A_HP[p, p] = 1
            # Loop over all internal mesh points in y diretion
            # and all mesh points in x direction
            for j in IyHP[1:-1]:
                i = 0;  p = m(i,j);  A_HP[p, p] = 1 # boundary
                for i in Ix[1:-1]:                 # interior points
                    
                    if u[i,n_HP] >= T_TEG_MAX: # active 
                        ''' HP F modifications '''
                        alpha[j] = empty * k_Al/(rho_Al*Cp_Al)
                        C[j] = empty * rho_Al*Cp_Al
                        Fy[j] = (alpha[j]*dt)/dy**2
                        Fx[j] = (alpha[j]*dt)/dx**2
                    else: # inactive 
                        alpha[j] = empty * k_Al/(rho_Al*Cp_Al)
                        C[j] = empty * rho_Al*Cp_Al
                        Fy[j] = (alpha[j]*dt)/dy**2
                        Fx[j] = (alpha[j]*dt)/dx**2

                    
                    p = m(i,j)
                    A_HP[p, m(i,j-1)] = - theta*Fy[j]
                    A_HP[p, m(i-1,j)] = - theta*Fx[j]
                    A_HP[p, p]        = 1 + 2*theta*(Fx[j]+Fy[j])
                    A_HP[p, m(i+1,j)] = - theta*Fx[j]
                    A_HP[p, m(i,j+1)] = - theta*Fy[j]
                i = Nx;  p = m(i,j);  A_HP[p, p] = 1  # boundary
            # Equations corresponding to j=Ny, i=0,1,... (u known)
            j = n_HP
            for i in Ix:
                p = m(i,j);  A_HP[p, p] = 1
            
            # Compute b_HP
            j = 0
            for i in Ix:
                b_HP[m(i,j)] = uLoo[i] # convective boundary
            for j in IyHP[1:-1]:
                i = 0
                b_HP[m(i,j)] = u[i+1,j] + (1-theta)*(Fx[j]*(u_n[i+1,j] - u_n[i,j] ) +\
                Fy[j]*(u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1])) +\
                theta*dt*Q[i][j]/C[j] + \
                (1-theta)*dt*Q[i][j]/C[j] # # adiabatic boundary

                for i in Ix[1:-1]:
                    p = m(i,j)                                # interior
                    b_HP[p] = u_n[i,j] + \
                    (1-theta)*(Fx[j]*(u_n[i+1,j] - 2*u_n[i,j] + u_n[i-1,j]) +\
                    Fy[j]*(u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1])) +\
                    theta*dt*Q[i][j]/C[j] + \
                    (1-theta)*dt*Q[i][j]/C[j]
                    
                i = Nx
                b_HP[m(i,j)] = u[i-1,j] + (1-theta)*(dx*Fx[j]*(u_n[i,j] + u_n[i-1,j]) +\
                dy*Fy[j]*(u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1])) +\
                theta*dt*Q[i][j]/C[j] + \
                (1-theta)*dt*Q[i][j]/C[j] # # adiabatic boundary
            # boundary
            j = n_HP
            for i in Ix:
                Q_qout = Q1[i] #(L_lat*dx)*k_Al*(u_n[i,n_HP-1]-u_n[i,n_HP])/dy
                Q_in = (L_lat*dx)*k_Al*(u_n[i,n_HP-1]-u_n[i,n_HP])/dy
                b_HP[m(i,j)] = u_n[i,j] + ((Q_in-Q_qout)*dt)/(rho_Al*dx*dy*L_lat*Cp_Al) #- k_Al*(L_lat*dx) * (u_n[i,j] - (u_TEG_N_main[i,1] + u_TEG_P_main[i,1] ) /2 )/delta_x # Q1 boundary + (k_Al * dx * lateral_length * (u_n[i,j] - (u_TEG_N[i,1] + u_TEG_P_main[i,1])/2))/dx 
            #SOLVER
            
            for i in Ix: # temperature control
                if u[i,n_HP] > T_TEG_MAX:
                    b_HP[m(i,j)] = T_TEG_MAX            
                    u[i,j] = T_TEG_MAX
                    
            c = spla.bicgstab(A_HP,b_HP,tol=1e-6)

            # Fill u with vector c
            for i in Ix:
                for j in IyHP:
                    u[i,j] = c[0][m(i,j)].copy()
                    
            '''================================= HP vapor =================================''' 
            for i in Ix:
                if u[i,n_HP] > T_TEG_MAX:
                    Q_excess[i] = ( (u[i,n_HP] - T_TEG_MAX)*rho_Al*(dx*dy*L_lat)*cp_metal_HP )/dt
                    b_HP[m(i,j)] = T_TEG_MAX            
                    u[i,j] = T_TEG_MAX
                    
                    

            '''================================= Cooler section =================================''' 
            
            # Compute b
            j = n_HP + 1
            for i in Ix:
                Q_qout = (L_lat*dx)*k_Al*(u_n[i,n_HP+1]-u_n[i,n_HP+2])/dy
                Q_in = Q2[i] #k_Al * (dx*L_lat) * ( (u_TEG_N_n[i,NNTEG-1] + u_TEG_P_n[i,NNTEG-1])/2 - u_n[i,n_HP+2])/dy
                b_cool[m(i,j - (n_HP+1))] = u_n[i,j] + ((Q_in - Q_qout)*dt)/(rho_Al*dx*dy*L_lat*Cp_Al) #- k_Al*(L_lat*dx) * (u_n[i,j] - (u_TEG_N_main[i,1] + u_TEG_P_main[i,1] ) /2 )/delta_x # Q2 boundary
                #b_cool[m(i,j - (n_HP+1))] = (u_TEG_N_n[i,NNTEG-1] + u_TEG_P_n[i,NNTEG-1])/2
            for j in range(n_HP+2,Ny):
                i = 0
                b_cool[m(i,j - (n_HP+1))] = u[i+1,j] + (1-theta)*(Fx[j]*(u_n[i+1,j] - u_n[i,j] ) +\
                Fy[j]*(u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1])) +\
                theta*dt*Q[i][j]/C[j] + \
                (1-theta)*dt*Q[i][j]/C[j] # # adiabatic boundary

                for i in Ix[1:-1]:
                    p = m(i,j - (n_HP+1))                                # interior
                    b_cool[p] = u_n[i,j] + \
                    (1-theta)*(Fx[j]*(u_n[i+1,j] - 2*u_n[i,j] + u_n[i-1,j]) +\
                    Fy[j]*(u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1])) +\
                    theta*dt*Q[i][j]/C[j] + \
                    (1-theta)*dt*Q[i][j]/C[j]
                    
                i = Nx
                b_cool[m(i,j - (n_HP+1))] = u[i-1,j] + (1-theta)*(dx*Fx[j]*(u_n[i,j] + u_n[i-1,j]) +\
                dy*Fy[j]*(u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1])) +\
                theta*dt*Q[i][j]/C[j] + \
                (1-theta)*dt*Q[i][j]/C[j] # # adiabatic boundary
            # boundary
            j = Ny
            for i in Ix:
                b_cool[m(i,j - (n_HP+1))] = uLcc[i] # convective boundary
            
            #SOLVER 
            c = spla.bicgstab(A_cool,b_cool, tol=1e-6)

            # Fill u with vector c
            for i in Ix:
                for j in range(n_HP + 1, Ny+1):
                    u[i,j] = c[0][m(i,j-(n_HP+1))].copy()


            '''================================= TEG section =================================''' 
            for i in Ix: # TEG calculations
                
                R_cont = 1/(10000*L_lat*dx)

                Q1[i] = ( u[i,n_HP] - (u_TEG_P_main[i,0] + u_TEG_N_main[i,0])/2)/R_cont # (L_lat*dx)*k_Al*(u_n[i,n_HP-1]-u_n[i,n_HP])/dy # k_Al * (dx*L_lat) * ( u_n[i,n_HP-1]-u_n[i,n_HP])/dy # [W/m2] heat flux out of TEG
                Q2[i] = ((u_TEG_N_main[i,NNTEG-1] + u_TEG_P_main[i,NNTEG-1])/2 - u[i,n_HP+1])/R_cont # k_Al * (dx*L_lat) * ( (u_TEG_N_n[i,NNTEG-1] + u_TEG_P_n[i,NNTEG-1])/2 - u_n[i,n_HP+1])/dy # [W/m2] heat flux out of TEG

                if t > 1:
                    [u_TEG_n_N, u_TEG_n_P, R, vemf,delta_x,Area_TEG_N] = fun_TEG(Q1[i],Q2[i],u_TEG_P_n[i,:],u_TEG_N_n[i,:],u_TEG_N_main[i,:],u_TEG_P_main[i,:],dt, I_current)
                    
                    Avoid = 0.2*(A_crosssec_N + A_crosssec_P)
                    Rtot[i] = ( R * L_lat * dx ) / ( Avoid + A_crosssec_N + A_crosssec_P )
                    U[i] = ( ( vemf * L_lat * dx ) / ( Avoid + A_crosssec_N + A_crosssec_P ) )/2 # vemf -> OCV for one pair
                    u_TEG_N_main[i,:] = u_TEG_n_N.copy()
                    u_TEG_P_main[i,:] = u_TEG_n_P.copy()
            

            if sum(Rtot) != 0:
                I_current = (sum(U))/(2*sum(Rtot)) # 2*Rtot -> Rload = Rinterno

                U_matchedLoad = U-I_current*Rtot 
            
                    
            it += 1
            e1 = np.abs(u-u_old).max() 
            e2 = np.abs(u_TEG_N_main-u_TEG_N_old).max() 
            e3 = np.abs(u_TEG_P_main-u_TEG_P_old).max() 
            e4 = (I_current_old - I_current)
            converged = e1 < tol and e2 < tol and e3 < tol and e4 < tol
            Erro_tot = e1 + e2 + e3 + e4 

            ## RELAX
            relax = 1
            u = relax*u+(1-relax)*u_old
            u_TEG_N_main = relax*u_TEG_N_main + (1-relax)*u_TEG_N_old
            u_TEG_P_main = relax*u_TEG_P_main + (1-relax)*u_TEG_P_old

            I_current_old = I_current    
            u_TEG_N_old = u_TEG_N_main.copy()
            u_TEG_P_old = u_TEG_P_main.copy()
            u_old = u.copy()
            
            if it > max_iter:
                print(t*dt, ' skiped with error: ', Erro_tot)
                break

        ## write cycle
        if ( t*dt == write or t == (N_time_steps/dt) ):
            
            ETA(start_time,t*dt,N_time_steps,Erro_tot,it)        

            print(round(t*dt, 2), ' Current: ', round(I_current, 2), ' U_matchedLoad: ', round(U_matchedLoad[0], 2), ' U0: ', round(U[0], 2), 
            ' Rtot: ', round(Rtot[0], 2),
            ' Pelect[W] ', round(U_matchedLoad[0]*I_current,2),
            ' Thot: ', round(u_TEG_N_main[0,0]- 273, 2) , ' Tcold: ', round(u_TEG_N_main[0,-1]- 273, 2)  )

            #print('\n', 't=',str(t*dt),' Energy balance= ', E_error, '\n', '\n')

            np.savetxt('results/{}_T2Dt.csv'.format(str(t*dt)), (u), delimiter=',')
            np.savetxt('results/{}_u_TEG_N.csv'.format(str(t*dt)), (u_TEG_N_main), delimiter=',')
            np.savetxt('results/{}_u_TEG_P.csv'.format(str(t*dt)), (u_TEG_P_main), delimiter=',')
            
            np.savetxt('results/{}_uLoo.csv'.format(str(t*dt)), (uoo), delimiter=',')
            np.savetxt('results/{}_uLcc.csv'.format(str(t*dt)), (ucc), delimiter=',')
            
            np.savetxt('results/{}_Q_exc.csv'.format(str(t*dt)), (Q_excess), delimiter=',')

            np.savetxt('results/{}_Qoo.csv'.format(str(t*dt)), (Qoo), delimiter=',')
            np.savetxt('results/{}_Qcc.csv'.format(str(t*dt)), (Qcc), delimiter=',')

            np.savetxt('results/{}_I.csv'.format(str(t*dt)), np.array(I_current, ndmin=1), delimiter=',')

                            
            if t == (N_time_steps/dt):
                print('Solver ran successfully')
            else:
                write = write + write_increment
                
                        
        
        # Update u_n before next step
        u_n = u.copy()
        
        u_TEG_N_n = u_TEG_N_main.copy()
        u_TEG_P_n = u_TEG_P_main.copy()


def SolveTEG():

    # initialization
    fun_TEG(0,0,0,0,0,0,0,0)
    Tamb = 20 + 273
    dt = 0.01
    t = 0
    
    tol = 1e-6
    I_current = 0
    global dx 
    global L_lat
    dx = 0.003
    L_lat = 0.003

    global write_increment
    write_increment = 0.1
    write = write_increment
    
    ## Declarations
    u_TEG_P_n = np.transpose( np.zeros([Nx_TEG+1]) )
    u_TEG_P_n[:] = Tamb

    u_TEG_N_n = u_TEG_P_n.copy()
    u_TEG_N_main = u_TEG_P_n.copy()
    u_TEG_P_main = u_TEG_P_n.copy()
    u_TEG_N_old = u_TEG_P_n.copy()
    u_TEG_P_old = u_TEG_P_n.copy()

    def condition():
        return True


    Q1 = 1
    Q2 = 0.5

    while condition():

        it = 0
        E_error = 1
        while E_error > tol:


            [u_TEG_N_main, u_TEG_P_main, R, vemf,delta_x,Area_TEG_N] = fun_TEG(Q1,Q2,
            u_TEG_P_n,u_TEG_N_n,u_TEG_N_main,u_TEG_P_main,dt,I_current)
            
            #u_TEG_N_main = u_TEG_n_N.copy()
            #u_TEG_P_main = u_TEG_n_P.copy()
            Rtot = R
            U = vemf
            
            if sum(Rtot) != 0:
                I_current = (2*sum(U))/(2*2*sum(Rtot)) # 2*Rtot -> Rload = Rinterno

            E_error = np.abs(u_TEG_N_main-u_TEG_N_old).max() +\
                 np.abs(u_TEG_P_main-u_TEG_P_old).max()


            it += 1


            ## RELAX
            relax = 0.3
            u_TEG_N_main = relax*u_TEG_N_main + (1-relax)*u_TEG_N_old
            u_TEG_P_main = relax*u_TEG_P_main + (1-relax)*u_TEG_P_old
                
            u_TEG_N_old = u_TEG_N_main.copy()
            u_TEG_P_old = u_TEG_P_main.copy()

            if str(E_error) == str(nan):
                sys.exit('Does not converge')

        plt.ion()
        PlotTEG(u_TEG_N_main, u_TEG_P_main)
        plt.pause(1)

        u_TEG_N_n = u_TEG_N_main.copy()
        u_TEG_P_n = u_TEG_P_main.copy()


        ## write cycle
        if ( t*dt == write ):
            
            ETA(start_time,t*dt,1000,E_error,it)        

            print('\n', 't=',str(t*dt),' Energy balance= ', E_error, '\n', '\n')

            np.savetxt('results/{}_u_TEG_N.csv'.format(str(t*dt)), (u_TEG_N_main), delimiter=',')
            np.savetxt('results/{}_u_TEG_P.csv'.format(str(t*dt)), (u_TEG_P_main), delimiter=',')
            np.savetxt('results/{}_I.csv'.format(str(t*dt)), np.array(I_current, ndmin=1), delimiter=',')
                       
            write = write + write_increment

        t += 1


''' PostProcess functions'''

def PlotTEG(uN, uP):

    fig = figure(1,clear=True)
    ax = fig.add_subplot(1, 1, 1)

    X = np.arange(0, len(uN), 1)

    ax.plot(X,uN, label='N-type', color='r')
    ax.plot(X,uP, label='P-type', color='b')
    ax.set_xlabel('TEG leg length')
    ax.set_ylabel('Temperature [K]')

    legend = ax.legend()


def Plots(u, uTEGn, uTEGp, uoo, ucc, Qexc, Pelect, Qoo, TEGeff, n):

    tStep = 0.5

    axis_font ={'fontname':'Arial', 'size':'14'}
    linewidth = 0.5

    fig = plt.figure(1, clear=True, figsize=(10,8))

    ''' 3D plot Temperature '''
    ax = fig.add_subplot(2, 2, 1, projection='3d')

    [yy, xx] = u.shape
    X = np.arange(0, xx, 1)
    Y = np.arange(0, yy, 1)
    X, Y = np.meshgrid(X, Y)

    surf = ax.plot_surface(X, Y, u, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    ax.set_title('Temperature 3D at {}'.format(n), axis_font)
    ax.set_xlabel('ny, k dT/dx = 0', axis_font)
    ax.set_ylabel('Exhaust direction', axis_font)

    # Customize the z axis.
    try:
        ax.set_zlim(u.min().min(), u.max().max())
    except:
        pass
    ax.zaxis.set_major_locator(LinearLocator(10))
    #ax.zaxis.set_major_formatter(FormatStrFormatter('.02f'))
    ax.zaxis.set_major_formatter(ticker.ScalarFormatter())

    ax.view_init(20,20)
    
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    ''' 2D Temperature levels '''
    ax = fig.add_subplot(2, 2, 2)

    [xx,yy] = uTEGn.shape
    X = np.arange(0, xx, 1)

    ax.plot(X,uoo, color='r', label='Exhaust')
    ax.plot(X,ucc, color='b', label='Cooler')
    ax.plot(X,uTEGn.iloc[:][0], '--', color='k', label='Hot Side TEG')
    ax.plot(X,uTEGn.iloc[:][yy-1], '--', color='gray', label='Cold side TEG')

    legend = ax.legend(loc='upper right', fontsize=axis_font['size'])
    ax.set_ylabel('Temperature [K]')

    ''' Total Electrical Power '''
    ax = fig.add_subplot(2, 2, 3)

    X = np.arange(tStep, n+tStep, tStep)

    pl1 = ax.plot(X,Pelect, color='g', label='Electrical Power', linewidth=linewidth)

    #ax.set_title('Total electrical Power', axis_font)
    ax.set_xlabel('Time [s]', axis_font)
    ax.set_ylabel('Electrical Power [w]', axis_font)

    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

    pl2 = ax2.plot(X, Qexc/1000, color='b', label='Absobed by HP', linewidth=linewidth)
    pl3 = ax2.plot(X, Qoo/1000, color='r', label='Availabel Exhaust', linewidth=linewidth)

    ax2.set_ylabel('Thermal Power [kw]', axis_font)

    lns = pl1 + pl2 + pl3 
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='upper left', fontsize=axis_font['size'])


    ''' Efficiency '''
    ax = fig.add_subplot(2, 2, 4)

    X = np.arange(tStep, n+tStep, tStep)

    labels = ['TEGeff', 'I_TEG', 'U_TEG']
    [x,y] = TEGeff.shape
    for i in range(y):
        ax.plot(X,TEGeff.loc[:][i], linewidth=linewidth, label=labels[i])

    ax.legend()
    #ax.set_title('TEG T', axis_font)
    ax.set_xlabel('Time [s]', axis_font)
    ax.set_ylabel('TEG ave. efficiency [%]', axis_font)

    fig.tight_layout() 


def PostProcess():

    fun_TEG(0,0,0,0,0,0,0,0)

    Tau_P = 0.00027 # [V/K] thomson coefficient
    Tau_N = -0.000156 # [V/K] thomson coefficient
    # n_leg divide by 6, not all legs fit in one slice 
    R_Exhaust (1,1,1,1, 60)        
    dx = 0.01
    Avoid = 0.2*( A_crosssec_N + A_crosssec_P) 
    n_leg = (L_lat*dx)/(Avoid + A_crosssec_N + A_crosssec_P)
    n_HiZ = 98




    Tamb = 20 + 273
    Temp_N = 20 + 273 # ref temp is the Tamb
    Temp_P = 20 + 273 # ref temp is the Tamb
    cp_Bi2Te3 = Prop_TEG_Cp (Temp_P)
    [alfa_Bi2Te3N, rho_Bi2Te3N, k_Bi2Te3N] = Prop_TEG_N(Temp_N)
    [alfa_Bi2Te3P, rho_Bi2Te3P, k_Bi2Te3P] = Prop_TEG_P(Temp_P)
    alfa_N_now_ref = alfa_Bi2Te3N
    alfa_P_now_ref = alfa_Bi2Te3P

    # empty declarations
    R_write = pd.DataFrame()
    U0 = pd.DataFrame()
    I_write = pd.DataFrame()
    J_N = pd.DataFrame()
    J_P = pd.DataFrame()

    Q_out_N = pd.DataFrame()
    Q_in_N = pd.DataFrame()
    Q_out_P = pd.DataFrame()
    Q_in_P = pd.DataFrame()
    Q_legs = pd.DataFrame()
    Q_eggcrate = pd.DataFrame()
    Q_in_module = pd.DataFrame()
    Q_out_module = pd.DataFrame()
    P_leg_load = pd.DataFrame()
    P_load = pd.DataFrame()
    eff_module = pd.DataFrame()
    U_TEG = pd.DataFrame()

    P_total = pd.DataFrame()
    Q_exc_acc = pd.DataFrame()
    Qin_absobed = pd.DataFrame()

    n = 0.5 # change with the write_increment variable || in the RESULTS folder check for the t times solved for, n = t0 - t1 
    t = 0
    t_step = 0.5

    tStamp = time.time()

    while True:
        print(n)
        
        try:
            u = pd.read_csv('results/{}_T2Dt.csv'.format(str(n)),',',index_col=False,header=None, sep=';')
            print(u)
            u_TEG_P = pd.read_csv('results/{}_u_TEG_P.csv'.format(str(n)),',',index_col=False,header=None, sep=';')
            u_TEG_N = pd.read_csv('results/{}_u_TEG_N.csv'.format(str(n)),',',index_col=False,header=None, sep=';')

            uLoo = pd.read_csv('results/{}_uLoo.csv'.format(str(n)),',',index_col=False,header=None)
            uLcc = pd.read_csv('results/{}_uLcc.csv'.format(str(n)),',',index_col=False,header=None)
            
            Q_exc = pd.read_csv('results/{}_Q_exc.csv'.format(str(n)),',',index_col=False,header=None)
            Qoo = pd.read_csv('results/{}_Qoo.csv'.format(str(n)),',',index_col=False,header=None)
            Qcc = pd.read_csv('results/{}_Qcc.csv'.format(str(n)),',',index_col=False,header=None)

            I = pd.read_csv('results/{}_I.csv'.format(str(n)),',',index_col=False,header=None)

            tStamp = time.time()

        except:

            try:
                u_TEG_P = pd.read_csv('results/{}_u_TEG_P.csv'.format(str(n)),',',index_col=False,header=None)
                u_TEG_N = pd.read_csv('results/{}_u_TEG_N.csv'.format(str(n)),',',index_col=False,header=None)
                I = pd.read_csv('results/{}_I.csv'.format(str(n)),',',index_col=False,header=None)

                u = pd.DataFrame(np.zeros([2, 2]))
                Q_exc = pd.DataFrame(np.zeros([2, 2]))
                Qoo = pd.DataFrame(np.zeros([2, 2]))
                uLoo = pd.DataFrame(np.zeros([2, 2]))
                uLcc = pd.DataFrame(np.zeros([2, 2]))
            
            except:
                
                dtStamp = time.time() - tStamp
                if dtStamp > 60:
                    print('No more data in 60s, ending ')
                    break
                continue

        k_Al = 200       
        cp_metal_HP = 910
        k_egg = 0.1
        width1 = 5.81/100 # TEG module width
        width2 = 5.81/100
        L_egg = 0.3/100
        A_eggcrate = Avoid #width1 * width2 - ((n_leg)*(leg_side_N**2) + (n_leg)*(leg_side_P**2)) # [m2]
        

        [xx, yy] = u_TEG_P.shape

        #Declarations
        rho_N_before = pd.DataFrame(np.zeros([yy]))
        rho_P_before = pd.DataFrame(np.zeros([yy]))
        density_before = pd.DataFrame(np.zeros([yy]))
        alfa_N_now = pd.DataFrame(np.zeros([yy]))
        alfa_P_now = pd.DataFrame(np.zeros([yy]))
        k_N_now = pd.DataFrame(np.zeros([yy]))
        k_P_now = pd.DataFrame(np.zeros([yy]))
        Cp_N_now = pd.DataFrame(np.zeros([yy]))
        Cp_P_now = pd.DataFrame(np.zeros([yy]))
        R_y = pd.DataFrame(np.zeros([yy]))
        Vemf = pd.DataFrame(np.zeros([yy]))


        Q_joulle_module = pd.DataFrame(np.zeros([xx]))
        Q_Peltier_N = pd.DataFrame(np.zeros([xx]))
        Q_Peltier_P = pd.DataFrame(np.zeros([xx]))


        Q_joulle_N_y = pd.DataFrame(np.zeros([xx,yy-1]))
        Q_joulle_P_y = pd.DataFrame(np.zeros([xx,yy-1]))
        Q_Peltier_N_Y = pd.DataFrame(np.zeros([xx,yy-1]))
        Q_Peltier_P_Y = pd.DataFrame(np.zeros([xx,yy-1]))        
        Z_N = pd.DataFrame(np.zeros([xx,yy-1]))
        Z_P = pd.DataFrame(np.zeros([xx,yy-1]))
        q_N = pd.DataFrame(np.zeros([xx,yy-1]))
        q_P = pd.DataFrame(np.zeros([xx,yy-1]))
        
        R_write = pd.concat( [R_write, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        U0 = pd.concat( [U0, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        J_N = pd.concat( [J_N, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        J_P = pd.concat( [J_P, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        U_TEG = pd.concat( [U_TEG, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )

        Q_out_N = pd.concat( [Q_out_N, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        Q_in_N = pd.concat( [Q_in_N, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        Q_out_P = pd.concat( [Q_out_P, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        Q_in_P = pd.concat( [Q_in_P, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        Q_legs = pd.concat( [Q_legs, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        Q_eggcrate = pd.concat( [Q_eggcrate, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        Q_in_module = pd.concat( [Q_in_module, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        Q_out_module = pd.concat( [Q_out_module, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        P_leg_load = pd.concat( [P_leg_load, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        P_load = pd.concat( [P_load, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )
        eff_module = pd.concat( [eff_module, pd.DataFrame(np.zeros([xx]))], axis=1, ignore_index=True )

        I_write = I_write .append( [I], ignore_index=True)

        for i in range(len(u)):
            
            for y in range(yy):

                if y >= N_connector and y <= N_connector + N_legs: 

                    Temp_N = u_TEG_N.iloc[i][y]
                    Temp_P = u_TEG_P.iloc[i][y]
                    Temp = (Temp_N + Temp_P)/2
                    
                    [alfa_Bi2Te3N, rho_Bi2Te3N, k_Bi2Te3N] = Prop_TEG_N(Temp_N)
                    [alfa_Bi2Te3P, rho_Bi2Te3P, k_Bi2Te3P] = Prop_TEG_P(Temp_P)
                    
                    rho_N_before[0][y] = rho_Bi2Te3N
                    rho_P_before[0][y] = rho_Bi2Te3P
                    
                    density_before[0][y] = 7700
                    
                else:
                    
                    rho_N_before[0][y] = 2.82e-08
                    rho_P_before[0][y] = 2.82e-08
                    
                    density_before[0][y] = 2712
                    
                # end material properties
                
                ## t+1 time
                if y >= N_connector and y <= N_connector + N_legs: # material properties
                    
                    Temp_N = u_TEG_N.iloc[i][y]
                    Temp_P = u_TEG_P.iloc[i][y]
                    Temp = (Temp_N + Temp_P)/2
                    
                    cp_Bi2Te3 = Prop_TEG_Cp (Temp)
                    [alfa_Bi2Te3N, rho_Bi2Te3N, k_Bi2Te3N] = Prop_TEG_N(Temp_N)
                    [alfa_Bi2Te3P, rho_Bi2Te3P, k_Bi2Te3P] = Prop_TEG_P(Temp_P)
                    
                    alfa_N_now[0][y] = alfa_Bi2Te3N
                    alfa_P_now[0][y] = alfa_Bi2Te3P
                    
                    k_N_now[0][y] = k_Bi2Te3N
                    k_P_now[0][y] = k_Bi2Te3P
                    
                    Cp_N_now[0][y] = cp_Bi2Te3
                    Cp_P_now[0][y] = cp_Bi2Te3
                    
                else:
                    alfa_N_now[0][y] = 3.5*10**-6
                    alfa_P_now[0][y] = 3.5*10**-6
                    
                    k_N_now[0][y] = k_Al
                    k_P_now[0][y] = k_Al
                    
                    Cp_N_now[0][y] = cp_metal_HP
                    Cp_P_now[0][y] = cp_metal_HP
                # end material properties


            for y in range(1,yy):

                Vemf.loc[y] = abs( alfa_N_now[0][y]*(u_TEG_N[y-1] - u_TEG_N[y]) ) + abs( alfa_P_now[0][y]*(u_TEG_P[y-1] - u_TEG_P[y]) ) # U0 pair
                R_y.loc[y] = (rho_P_before[0][y]*delta_x)/A_crosssec_P + (rho_N_before[0][y]*delta_x)/A_crosssec_N 

            #alpha_PN_0 = (alfa_N_now_ref - alfa_P_now_ref) + (Tau_P - Tau_N)*np.log( u_TEG_N.iloc[i][0]/Tamb )
            #alpha_PN_L = (alfa_N_now_ref - alfa_P_now_ref) + (Tau_P - Tau_N)*np.log( ((u_TEG_N.iloc[i][yy-1] + u_TEG_P.iloc[i][yy-1])/2)/Tamb )

            U0.loc[i][t] = ( sum(Vemf) * (L_lat * dx) ) / ( Avoid + A_crosssec_N + A_crosssec_P) /2

            R_write.loc[i][t] = ( sum(R_y) * (L_lat * dx) ) / ( Avoid + A_crosssec_N + A_crosssec_P) # Ri 
            #R0 = R_write.loc[i][t]

            U_TEG.loc[i][t] = U0.loc[i][t] - R_write.loc[i][t]*I_write.loc[t]
            
            #current density
            J_N.loc[i][t] = I_write.loc[t]/A_crosssec_N
            J_P.loc[i][t] = I_write.loc[t]/A_crosssec_P

            for y in range(yy):    
                                
                # [W] Joulle heat
                #Q_joulle_N_y.loc[i][y] = rho_N_before[0][y] * J_N.loc[i][t]**2
                #Q_joulle_P_y.loc[i][y] = rho_P_before[0][y] * J_P.loc[i][t]**2
                
                if y < yy-1:
                    # [W] Peltier heat
                    Q_Peltier_N_Y.loc[i][y] = I_write.loc[t] * ( abs(alfa_N_now[0][y]) ) * u_TEG_N.iloc[i][y]
                    Q_Peltier_P_Y.loc[i][y] = I_write.loc[t] * ( alfa_P_now[0][y] ) * u_TEG_P.iloc[i][y]
                
                
                Z_N.loc[i][y] = ( (1/rho_N_before[0][y])*alfa_N_now[0][y]**2 ) / k_N_now[0][y]
                Z_P.loc[i][y] = ( (1/rho_P_before[0][y])*alfa_P_now[0][y]**2 ) / k_P_now[0][y]
                
                if y < yy-1:

                    q_N.loc[i][y] = ( k_N_now[0][y]*(u_TEG_N.iloc[i][y] - u_TEG_N.iloc[i][y+1])/delta_x )*( alfa_N_now[0][y]*u_TEG_N.iloc[i][y]*( Z_N.loc[i][y] / (2*alfa_N_now[0][y]) ) + 1 )
                    q_P.loc[i][y] = ( k_P_now[0][y]*(u_TEG_P.iloc[i][y] - u_TEG_P.iloc[i][y+1])/delta_x )*( alfa_P_now[0][y]*u_TEG_P.iloc[i][y]*( Z_P.loc[i][y] / (2*alfa_P_now[0][y]) ) + 1 )
                
            
            
            #Q_joulle_N = np.sum(Q_joulle_N_y, axis=0) # sum(Q_joulle_N_y(i,:,t)) #  # [W]
            #Q_joulle_P = np.sum(Q_joulle_P_y, axis=0) # sum(Q_joulle_P_y(i,:,t)) # # [W]
            
            Q_joulle_module.loc[i] = ( R_write.loc[i][t]*I_write.loc[t]**2 )
            
            Q_Peltier_N.loc[i] = Q_Peltier_N_Y.loc[i][0]  # [W]
            Q_Peltier_P.loc[i] = Q_Peltier_P_Y.loc[i][0] # [W]
            
            Q_out_N.loc[i][t] = abs(q_N[size(q_N,axis=1) - 1][i]) * A_crosssec_N  # [W]
            Q_in_N.loc[i][t] = abs(q_N[0][i]) * A_crosssec_N  # [W]
            
            Q_out_P.loc[i][t] = abs(q_P[size(q_P,axis=1) - 1][i]) * A_crosssec_P  # [W]
            Q_in_P.loc[i][t] = abs(q_P[0][i]) * A_crosssec_P  # [W]
            
            #Q_legs.loc[i][t] = (n_leg/2)*(Q_in_N.loc[i][t]+Q_in_P.loc[i][t])  # [W] heat traversing the legs
            
            Q_eggcrate.loc[i][t] = k_egg * (u_TEG_N.iloc[i][0] -  u_TEG_N.iloc[i][yy-1] ) * A_eggcrate / L_egg  # [W] thermal power
            #Q_in_module.loc[i][t] = n_leg * (Q_in_N.loc[i][t] + Q_in_P.loc[i][t] + Q_Peltier_N.loc[i] + Q_Peltier_P.loc[i] )/2 + Q_eggcrate.loc[i][t] - Q_joulle_module.loc[i]  # [W] thermal power traversing thne module
            #Q_out_module.loc[i][t] = n_leg * (Q_out_N.loc[i][t] + Q_out_P.loc[i][t] + Q_Peltier_P_Y.loc[i][yy-2] + Q_Peltier_N_Y.loc[i][yy-2] )/2 + Q_eggcrate.loc[i][t] + Q_joulle_module.loc[i]  # [W] thermal power traversing thne module
            
            Q_in_module.loc[i][t] = n_leg * (Q_in_N.loc[i][t] + Q_in_P.loc[i][t] + Q_Peltier_N.loc[i] + Q_Peltier_P.loc[i] )/2 + Q_eggcrate.loc[i][t]   # [W] thermal power traversing thne module
            Q_out_module.loc[i][t] = n_leg * (Q_out_N.loc[i][t] + Q_out_P.loc[i][t]+ Q_Peltier_P_Y.loc[i][yy-2] + Q_Peltier_N_Y.loc[i][yy-2] )/2 + Q_eggcrate.loc[i][t]   # [W] thermal power traversing thne module

            # electrical power
            
            P_leg_load.loc[i][t] = I_write.loc[t]*U_TEG.loc[i][t]
            P_load.loc[i][t] = 2*P_leg_load.loc[i][t] # both sides of the TEG
            
            
            # efficiency
            
            eff_module.loc[i][t] = P_load.loc[i][t] / (abs(Q_in_module.loc[i][t]))
            
        
        #P_total = P_total.append( pd.DataFrame( [np.sum(P_load[t], axis=0)] ), ignore_index=True)
        

        Q_exc_acc = Q_exc_acc.append(sum(Q_exc), ignore_index=True)
        Qin_absobed = Qin_absobed.append(sum(Qoo), ignore_index=True)

        QinTEG = Q_in_module.mean()
        QoutTEG = Q_out_module.abs().mean()

        TEGeff = 100*eff_module.mean()
        U_TEGplot = 2*sum(U_TEG)
        P_total = P_total.append( pd.DataFrame( [2*I_write.loc[t]*U_TEGplot[t]] ), ignore_index=True)

        #plt.ion()

        #Plots(u, u_TEG_N, u_TEG_P, uLoo, uLcc, Q_exc_acc, P_total, Qin_absobed, pd.concat( [TEGeff, I_write] , ignore_index=True, axis=1), n)

        #plt.pause(0.0001)

        #print('Power per module: ', P_total.iloc[-1]/60/2)

        t = t + 1
        n = n + t_step

        ## MORE PLOTS
        
        '''axis_font ={'fontname':'Arial', 'size':'14'}
        linewidth = 0.5

        fig = plt.figure(2, clear=True)

        ax = fig.add_subplot(1, 1, 1)

        ax.set_title(QinTEG.iloc[-1]-QoutTEG.iloc[-1])

        ax.plot(QinTEG, color='r', label='QinTEG')
        ax.plot(QoutTEG, color='b', label='QoutTEG')

        ax.yaxis.tick_right()

        ax.legend()'''

    ''' save to Excel '''
    times = pd.DataFrame( np.arange(t_step, n, t_step) )
    Excel_save = pd.concat( [times, P_total, I_write, U_TEGplot, TEGeff, Qin_absobed, Q_exc_acc, QinTEG, QoutTEG], axis=1 )
    Excel_save.columns = ['Time [s]', 'Total Electrical Power [W]', 'Current [A]', 'Voltage [V]', 
    'TEG Module Efficiency [%]', 'Exhaust Power Absorbed [W]', 'Excess Power Absorbed by HP [W]',
    'Qin TEG [W]', 'Qout TEG [W]']
    Excel_save.to_excel('results_summary/Summary.xlsx', index=False)


def main():

    
    a_input = input('HP Solver [1], TEG Solver [2] or PostProcess [3]? ')

    if a_input == str(1):

        init_n_print()

        import multiprocessing 

        p1 = multiprocessing.Process(target=Solver)
        p2 = multiprocessing.Process(target=PostProcess)

        p1.start()
        p2.start()

        p1.join() # only moves on after p1 and p2 are complet
        p2.join()

    elif a_input == str(2):

        init_n_print()

        import multiprocessing 

        p1 = multiprocessing.Process(target=SolveTEG)
        p2 = multiprocessing.Process(target=PostProcess)

        p1.start()
        p2.start()

        p1.join() # only moves on after p1 and p2 are complet
        p2.join()

    else:
        PostProcess()


if __name__ == '__main__':

    main()
