# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 12:08:22 2022

@author: farhan
"""

"""
_____________________________________________________________
Versions:
    
python 3.9.13
numpy 1.21.5
panda 1.4.4
scipy 1.9.1
scikit-learn 1.0.2
matplotlib 3.5.2
CoolProp 6.4.1
_____________________________________________________________
"""
import numpy as np
import pandas as pd
from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI
import glob2 ,os, sys
import math
from scipy.optimize import fsolve
from uncertainties import ufloat,umath
try:
    from Partial_diff import air_props_un
except:
    from PythonFiles.Partial_diff import air_props_un
import copy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
from pathlib import Path

def air_flow(delta_P, un_delta_P, P_a, P_a_un, T_db, T_db_un, T_wb, T_wb_un,
             ND, un_ND, Test_name):
    '''
    Function to calculate the uncertainty in air flowrate
    delta_P: pressure drop across nozzle in Pa
    un_delta_P: uncertainty in delta P in Pa
    P_a: Air pressure at nozzle inlet in Pa
    un_P_a: Uncertainty in pressure at nozzle inlet in Pa
    T_db: Air dry bulb temperature at nozzle inlet in K
    un_T_db: Uncertainty in T_db in K
    T_wb: Air wet bulb temperature at nozzle inlet in K
    un_T_wb: Uncertainty in T_wb in K
    ND: Nozzle diamter , can be an int or array in m
    un_ND: uncertainty in ND in m
    '''

    # Air density and viscosity

    rho=1/(HAPropsSI('V', 'T', T_db,'P', P_a,
                        'Twb', T_wb))

    mu=HAPropsSI('mu', 'T', T_db, 'P', P_a,
                        'Twb', T_wb)

    # Uncertainty in Air density and viscosity
    
    #Calculate uncertainty in density calculation from CoolProp
    rho_un=air_props_un('Dens', T_db, T_db_un, T_wb,
                 T_wb_un, P_a, P_a_un, thres=1e-5)

    rho_un=rho_un[0]
    
    #Calculate uncertainty in viscosity calculation from CoolProp
    mu_un=air_props_un('mu', T_db, T_db_un, T_wb,
                 T_wb_un, P_a, P_a_un, thres=1e-5)

    mu_un=mu_un[0]

    # Air density and viscosity with uncertainties
    rho=ufloat(rho, rho_un, "rho")
    mu=ufloat(mu, mu_un, "mu")

    # Delta P
    delta_P=ufloat(delta_P, un_delta_P, "delta_P")

    # Alpha from ASHRAE standard 41.2-2018, page 18
    Alpha = P_a / (delta_P + P_a)

    # Expansion factor from ASHRAE standard 41.2-2018, page 18, eqn 9-22
    Y=1-(0.548*(1-Alpha))

    # Multiple nozzles or single nozzle
    if type(ND)==list or type(ND)==np.ndarray:
        NA=[] #Nozzle area
        Re=[]
        for i in range(len(ND)):
            NA.append((math.pi*ND[i]**2)/4)
    else:
        NA=[(math.pi*ND**2)/4]

    # Expansion coefficients
    CA=0

    def obj(Coeff):
        C=Coeff
        [Re,C1,Res]=[[],[],[]]
        
        for i in range(len(NA)):
            #ASHRAE 41.2-2018 page 18 eqn 9-20, Reynolds number
            Re.append(C[i]*(ND[i]*rho.nominal_value*Y.nominal_value*math.sqrt(2*(delta_P.nominal_value)/rho.nominal_value)/mu.nominal_value))

            #nozzle coefficient calculation, ASHRAE 41.2-2018 page 18 eqn 9-23
            C1.append(0.99855-(7.006/math.sqrt(Re[i]))+(134.6/Re[i]))
            Res.append(C[i]-C1[i])

        return(Res)

    C1=fsolve(obj,[1]*len(NA))

    Cn=[]

    ND_un=copy.deepcopy(ND)
    un_ND=copy.deepcopy(un_ND)

    [NDD, Re]=[[], []]
    
    for i in range(len(ND_un)):
        NDD.append(ufloat(ND_un[i], un_ND[i], "NDD"))
        NA[i]=(math.pi*NDD[i]**2)/4

        #ASHRAE 41.2-2018 page 18 eqn 9-20, Reynolds number
        Re.append(C1[i]*(NDD[i]*rho*umath.sqrt(2*delta_P/rho)/mu))
        
        #nozzle coefficient calculation,  ASHRAE 41.2-2018 page 18 eqn 9-23
        Cn.append(0.9986-7.006/umath.sqrt(Re[i])+134.6/Re[i])
        CA=CA+Cn[i]*NA[i]

    Q_SI=Y*umath.sqrt(2*delta_P/rho)*CA # Airflow rate SI

    
    Error=[]
    for (var, error) in Q_SI.error_components().items():
        percent_contrib=(error**2)*100/(Q_SI.s**2)
        Error.append(percent_contrib)
    
    fig1, ax1 = plt.subplots()
    labels=['$NDD$', '$mu$', '$rho$', '$Delta_P$', '$NDD$']
    sizes=Error[0:5]
    patches, texts = ax1.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, Error)]
    
    sort_legend = False
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.),
               fontsize=11)
    scr_loc = os.path.dirname( __file__ ) #Script Location
    prj_fold_loc = str(Path(scr_loc).parents[0]) #Project folder location
    plot_ind_dir = Test_name + "\\"
    plt_dir = prj_fold_loc + "\\06-plots\\" + plot_ind_dir
    plt_dir_exist = os.path.exists(plt_dir)
    if not plt_dir_exist:
        os.makedirs(plt_dir)
    plt.savefig(plt_dir + 'Uncertainty analysis of Air flow rate',
            dpi = 300,  bbox_inches = 'tight')
    #print (Error)
    return(Q_SI, Error)

# """
# TESTING THE FUNCTIONS
# """      
# if __name__=='__main__':
#     import matplotlib.pyplot as plt
#     import numpy as np
#     import matplotlib.ticker as mtick
    
#     Q=[]
#     results= []
# ##    delta_P=np.arange(200,300,5.25)
#     delta_P=[193.35]
#     fig, ax1 = plt.subplots()
#     for i in range(len(delta_P)):
#         un_delta_P=1.866
#         P_a=97080.00
#         T_db=273+17.2
#         T_wb=273+9.7
#         T_db_un=0.1
#         T_wb_un=0.1
#         P_a_un=30.94
#         ND=[0.1778, 0.1524]
#         un_ND=[0.002*ND[0], 0.002*ND[1]]
#         results = air_flow(delta_P[i],un_delta_P,P_a,P_a_un,T_db,T_db_un,T_wb,T_wb_un,ND,un_ND)
#         Q.append(results[0])
#         print(' the flowrate is ',Q[i].nominal_value*2118.88)
#         print('relative uncer is ', Q[i].s*100/Q[i].nominal_value)
#         ax1.scatter(Q[i].nominal_value*2118.88,Q[i].s*100/Q[i].nominal_value,color='blue')

#     fmt = '%.2f%%' # Format you want the ticks, e.g. '40%'
#     xticks = mtick.FormatStrFormatter(fmt)
#     ax1.yaxis.set_major_formatter(xticks)
#     ax1.set_xlabel('Air flow rate [CFM]',fontsize=14)
#     ax1.set_ylabel('Uncertainty [-]',fontsize=14)
#     plt.xticks(fontsize=14)
#     plt.yticks(fontsize=14)
#     plt.tight_layout()
#     plt.savefig('Uncer.png',dpi=450)
