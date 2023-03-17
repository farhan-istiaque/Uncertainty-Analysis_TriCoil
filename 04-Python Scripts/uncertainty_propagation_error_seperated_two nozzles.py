# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 12:08:22 2022
@Author Dr. Omar Sarfraz
@Modified by Farhan
"""

"""
_____________________________________________________________
Versions:
    
python 3.9.13
numpy 1.21.5
pandas 1.4.4
scipy 1.9.1
scikit-learn 1.0.2
matplotlib 3.5.2
CoolProp 6.4.1
____________________________________________________________
"""
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import numpy as np
from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI
from uncertainties import ufloat,umath,unumpy
import matplotlib.pyplot as plt
import math
import os
import glob
from pathlib import Path
from csv import writer

try:
    from Partial_diff import air_props_un
except:
    from PythonFiles.Partial_diff import air_props_un
    
try:
    from Partial_diff import ref_props_un
except:
    from PythonFiles.Partial_diff import ref_props_un
    
try:
    from Airflow_uncertainty import air_flow
except:
    from Airflow_uncertainty import air_flow

"""
Reading the file into numpy arrays
"""
scr_loc = os.path.dirname( __file__ ) #Script Location
prj_fold_loc = str(Path(scr_loc).parents[0]) #Project folder location
"""
Change the title of raw data file here
"""

data_file_loc = prj_fold_loc + "\\03-Trimmed Raw Data\\two_nozzles\\"
files = glob.glob(data_file_loc + "/*.csv")



for filename in files:
    df = pd.read_csv(filename)
    
    Test_name = filename[-14:]
    
    """
    Setting the plot directory
    """
    plot_ind_dir = Test_name + "\\"
    plt_dir = prj_fold_loc + "\\06-plots\\" + plot_ind_dir
    plt_dir_exist = os.path.exists(plt_dir)
    if not plt_dir_exist:
        os.makedirs(plt_dir)
    
    #Data headers
    P_atm = df[["Atm Pressure (PSI)"]].to_numpy()
    P_atm = P_atm * 6894.76 #PSI to Pa
    P_atm_avg = np.mean(P_atm)
    
    m_dot_w = df[["m_dot_w (kg/h)"]].to_numpy()
    m_dot_w = m_dot_w / 3600 #kg/h to kg/s
    
    T_w_in = df[["T_w_in (C)"]].to_numpy()
    T_w_in += 273.15 #C to K
    
    T_w_out = df[["T_w_out (C)"]].to_numpy()
    T_w_out += 273.15 #C to K
    
    P_w_in = df[["P_w_in (PSIG)"]].to_numpy()
    P_w_in = (P_w_in * 6894.76) + P_atm_avg #PSIG to Pa
    
    dp_w = df[["dp_w (in W.G.)"]].to_numpy()
    dp_w = dp_w * 248.84 #W.G. to Pa
    
    P_w_out = P_w_in + dp_w #Calculating outlet pressure of water from water deltaP
    
    P_a_in = df[["Pin_nozzle_TriCoil (in W.G)"]].to_numpy()
    P_a_in = P_a_in * 248.84 #in of W.G to Pa
    
    dp_a = df[["dP_nozzle_TriCoil (in W.G.)"]].to_numpy()
    dp_a = dp_a * 248.84 #in of W.G to Pa
    
    T_db_in = df[["Tdb_in (C)"]].to_numpy()
    T_db_in += 273.15 #C to K
    
    T_wb_in = df[["Twb_in (C)"]].to_numpy()
    T_wb_in += 273.15 #C to K
    
    T_db_out = df[["TriCoil Outlet (Avg)"]].to_numpy()
    T_db_out += 273.15 #C to K
    
    T_wb_out = df[["Twb_out(C)"]].to_numpy()
    T_wb_out += 273.15 #C to K
    
    Rh_in = df[["RH_in"]].to_numpy()
    
    T_r_evap_in = df[["T_evap_in (C)"]].to_numpy()
    T_r_evap_in += 273.15
    
    T_r_exv_in = df[["T_exv_in (C)"]].to_numpy()
    T_r_exv_in += 273.15 #C to K
    
    P_r_exv_in = df[["P_exv_in (PSIG)"]].to_numpy()
    P_r_exv_in = (P_r_exv_in * 6894.76) + P_atm_avg #PSIG to Pa
    
    T_r_evap_out = df[["T_evap_out (C)"]].to_numpy()
    T_r_evap_out += 273.15 #C to K
    
    P_r_evap_out = df[["P_evap_out (PSIG)"]].to_numpy()
    P_r_evap_out = (P_r_evap_out * 6894.76) + P_atm_avg #PSIG to Pa
    
    m_dot_r = df[["m_dot_r (kg/h)"]].to_numpy()
    m_dot_r = m_dot_r / 3600 #kg/h to kg/s
    
    """
    Mean values
    """
    m_dot_w_avg = np.mean(m_dot_w)
    T_w_in_avg = np.mean(T_w_in)
    T_w_out_avg = np.mean(T_w_out)
    P_w_in_avg = np.mean(P_w_in)
    P_w_out_avg = np.mean(P_w_out)
    dp_w_avg = np.mean(dp_w)
    T_db_in_avg = np.mean(T_db_in)
    T_wb_in_avg = np.mean(T_wb_in)
    T_db_out_avg = np.mean(T_db_out)
    T_wb_out_avg = np.mean(T_wb_out)
    Rh_in_avg = np.mean(Rh_in)
    T_r_exv_in_avg = np.mean(T_r_exv_in)
    P_r_exv_in_avg = np.mean(P_r_exv_in)
    T_r_evap_in_avg = np.mean(T_r_evap_in)
    T_r_evap_out_avg = np.mean(T_r_evap_out)
    P_r_evap_out_avg = np.mean(P_r_evap_out)
    m_dot_r_avg = np.mean(m_dot_r)
    dp_a_avg = np.mean(dp_a)
    P_a_in_avg = np.mean(P_a_in)
    
    """
    Random Uncertainties
    """
    rnd_un_T_db_in = 2 * np.std(T_db_in)
    rnd_un_T_wb_in = 2 * np.std(T_wb_in)
    rnd_un_T_db_out = 2 * np.std(T_db_out)
    rnd_un_T_wb_out = 2 * np.std(T_wb_out)
    rnd_un_dp_a = 2 * np.std(dp_a)
    rnd_un_Pa_in = 2 * np.std(P_a_in)
    rnd_un_Pa_atm = 2 * np.std(P_atm)
    rnd_un_T_r_exv_in = 2 * np.std(T_r_exv_in)
    rnd_un_P_r_exv_in = 2 * np.std(P_r_exv_in)
    rnd_un_T_r_evap_out = 2 * np.std(T_r_evap_out)
    rnd_un_P_r_evap_out = 2 * np.std(P_r_evap_out)
    rnd_un_m_dot_r = 2 * np.std(m_dot_r)
    rnd_un_T_w_in = 2 * np.std(T_w_in)
    rnd_un_T_w_out = 2 * np.std(T_w_out)
    rnd_un_P_w_in = 2 * np.std(P_w_in)
    rnd_un_dp_w = 2 * np.std(dp_w)
    rnd_un_m_dot_w = 2 * np.std(m_dot_w)
    
    """
    Systemetic Uncertainties
    """
    
    sys_un_T_db_in = 0.002 #K
    sys_un_T_wb_in = 0.006 #K
    sys_un_T_db_out = 0.600 #K
    sys_un_T_wb_out = 0.004 #K
    sys_un_dp_a = 1.8663 #Pa
    sys_un_Pa_in = 6.221 #Pa
    sys_un_Pa_atm = 0.001 #Pa
    sys_un_T_r_exv_in = 0.007 #K
    sys_un_P_r_exv_in =4481.594 #Pa
    sys_un_T_r_evap_out = 0.005 #K
    sys_un_P_r_evap_out = 4481.594 #Pa
    sys_un_m_dot_r = 0.001 * m_dot_r_avg #kg/s
    sys_un_T_w_in = 0.005
    sys_un_T_w_out = 0.004
    sys_un_P_w_in = 4481.594 #Pa
    sys_un_dp_w = 93.315 #Pa
    sys_un_m_dot_w =  0.001 * m_dot_w_avg #kg/s
    
    """
    Uncertaintes of measure variable
    """
    
    un_T_db_in = math.sqrt(sys_un_T_db_in**2 + rnd_un_T_db_in**2)#K
    un_T_wb_in = math.sqrt(sys_un_T_wb_in**2 + rnd_un_T_wb_in**2) #K
    un_T_db_out = math.sqrt(sys_un_T_db_out**2 + rnd_un_T_db_out**2)#K
    un_T_wb_out = math.sqrt(sys_un_T_wb_out**2 + rnd_un_T_wb_out**2)#K
    un_dp_a = math.sqrt(sys_un_dp_a**2 + rnd_un_dp_a**2)#Pa
    un_Pa_in = math.sqrt(sys_un_Pa_in**2 + rnd_un_Pa_in**2)#Pa
    un_Pa_atm = math.sqrt(sys_un_Pa_atm**2 + rnd_un_Pa_atm**2)#Pa
    un_T_r_exv_in = math.sqrt(sys_un_T_r_exv_in**2 + rnd_un_T_r_exv_in**2) #K
    un_P_r_exv_in = math.sqrt(sys_un_P_r_exv_in**2 + rnd_un_P_r_exv_in**2) #K
    un_T_r_evap_out = math.sqrt(sys_un_T_r_evap_out**2 + rnd_un_T_r_evap_out**2) #K
    un_P_r_evap_out = math.sqrt(sys_un_P_r_evap_out**2 + rnd_un_P_r_evap_out**2) #K
    un_m_dot_r = math.sqrt(sys_un_m_dot_r**2 + rnd_un_m_dot_r**2) #Kg/s
    un_T_w_in = math.sqrt(sys_un_T_w_in**2 + rnd_un_T_w_in**2)#K
    un_T_w_out = math.sqrt(sys_un_T_w_out**2 + rnd_un_T_w_out**2)#K
    un_P_w_in = math.sqrt(sys_un_P_w_in**2 + rnd_un_P_w_in**2)#Pa
    un_dp_w = math.sqrt(sys_un_dp_w**2 + rnd_un_dp_w**2)#Pa
    un_m_dot_w =  math.sqrt(sys_un_m_dot_w**2 + rnd_un_m_dot_w**2)#kg/s
    
    """
    Water
    """
    #Enthalpy
    h_w_in = PropsSI('H', 'P', P_w_in_avg, 'T', T_w_in_avg, 'water')
    h_w_out = PropsSI('H', 'P', P_w_out_avg, 'T', T_w_out_avg, 'water')
    
    #Calculating uncertainty of P_w_out and contribution of P_w_in and dP
    un_P_w_out = math.sqrt(un_P_w_in**2 + un_dp_w**2)
    
    #Calculating uncertainty of enthanlpy calculation
    un_h_w_in = ref_props_un('H', T_w_in_avg, un_T_w_in, P_w_in_avg, un_P_w_in, 'water', thres=1e-5)
    
    un_h_w_out = ref_props_un('H', T_w_out_avg, un_T_w_out, P_w_out_avg, un_P_w_out, 'water', thres=1e-5)
    
    h_w_in_comb = ufloat(h_w_in, un_h_w_in[0], tag="h_w_in")
    h_w_out_comb = ufloat(h_w_out, un_h_w_out[0], tag="h_w_out")
    m_dot_w_comb = ufloat(m_dot_w_avg, un_m_dot_w, tag="m_dot_w")
    
    #Calculating uncertainty propagation in heat balance calculation
    Q_water = m_dot_w_comb * (h_w_out_comb - h_w_in_comb)
    
    #Percentage contribution of variables in the uncertainty propagation on water heat balance
    Error_w=[]
    
    # Plot % contribution of uncertainty of factors to heat transfer uncertainty
    for (var, error) in Q_water.error_components().items():
        percent_contrib=(error**2)*100/(Q_water.s**2)
        Error_w.append(percent_contrib)
    
    # combine all the errors 
    Error_w=[Error_w[0],Error_w[1],Error_w[2]]
    
    #Plotting the contribution of enthalpy and mass flow rate in water capacity uncertaitnty
    fig1, ax1 = plt.subplots()
    labels=['$h_{w_{in}}$', '$h_{w_{out}}$', '$\dot{m}_w$']
    sizes=Error_w
    patches, texts = ax1.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, Error_w)]
    
    sort_legend = True
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.),
                fontsize=11)
    plt.savefig(plt_dir + 'Uncertainty analysis of Water Capacity.png',
                dpi = 300, bbox_inches = 'tight')
    
    #Plotting the contribution of temperature and pressure in water inlet enthalpy uncertaitnty
    fig4, ax4 = plt.subplots()
    labels=['$T_{w_{in}}$', '$P_{w_{in}}$']
    sizes=un_h_w_in[1:3]
    patches, texts = ax4.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, un_h_w_in[1:3])]
    
    sort_legend = True
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.),
                fontsize=11)
    plt.savefig(plt_dir + 'Uncertainty analysis of Water Inlet Enthalpy',
                dpi = 300,  bbox_inches = 'tight')
    
    #Plotting the contribution of temperature and pressure in water outlet enthalpy uncertaitnty
    fig5, ax5 = plt.subplots()
    labels=['$T_{w_{out}}$', '$P_{w_{out}}$']
    sizes=un_h_w_out[1:3]
    patches, texts = ax5.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, un_h_w_out[1:3])]
    
    sort_legend = True
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.),
                fontsize=11)
    plt.savefig(plt_dir + 'Uncertainty analysis of Water Outlet Enthalpy',
                dpi = 300,  bbox_inches = 'tight')
    """
    Air
    """
    #Q_air uncertainty propagation
    h_a_in = HAPropsSI('H', 'T', T_db_in_avg, 'B', T_wb_in_avg, 'P', P_atm_avg) #Need to check
    h_a_out = HAPropsSI('H', 'T', T_db_out_avg, 'B', T_wb_out_avg, 'P', P_atm_avg) #Need to check
    
    ND=[0.1778, 0.1524]
    un_ND=[0.002*ND[0], 0.002*ND[1]]
    
    #Determining uncertainty and variable contribution for calculating Air flow rate
    V_dot_a = air_flow(dp_a_avg, un_dp_a, P_atm_avg, un_Pa_atm, T_db_in_avg, un_T_db_in, T_wb_in_avg, un_T_wb_in, ND, un_ND, Test_name)
    
    ND_perc_Vdota = V_dot_a[1][0]
    mu_perc_Vdota = V_dot_a[1][1]
    rho_perc_Vdota = V_dot_a[1][2]
    dpa_perc_Vdota = V_dot_a[1][3]
    
    #Determining uncertainties and contribution for calculating enthalpy
    un_h_a_in = air_props_un('H',T_db_in_avg, un_T_db_in, T_wb_in_avg, un_T_wb_in, P_atm_avg, un_Pa_atm, thres=1e-5)
    un_h_a_out = air_props_un('H',T_db_out_avg, un_T_db_out, T_wb_out_avg, un_T_wb_out, P_atm_avg, un_Pa_atm, thres=1e-5)
    
    rho_a = 1 / HAPropsSI('V', 'T', T_db_in_avg, 'B', T_wb_in_avg, 'P', P_atm_avg)
    un_rho_a = air_props_un('Dens',T_db_out_avg, un_T_db_out, T_wb_out_avg, un_T_wb_out, P_a_in_avg, un_Pa_atm, thres=1e-5)
    
    h_a_in_comb = ufloat(h_a_in, un_h_a_in[0], tag="h_a_in")
    h_a_out_comb = ufloat(h_a_out, un_h_a_out[0], tag = "h_a_out")
    rho_a_comb = ufloat(rho_a, un_rho_a[0], tag="rho_a")
    
    T_db_in_avg_comb = ufloat(T_db_in_avg, un_T_db_in, tag="T_db_in_avg")
    T_db_out_avg_comb = ufloat(T_db_out_avg, un_T_db_out, tag="T_db_out_avg")
    
    Q_air = rho_a_comb * V_dot_a[0] * abs(h_a_out_comb - h_a_in_comb)
    Q_air_sens = rho_a_comb * V_dot_a[0] * abs(T_db_in_avg_comb - T_db_out_avg_comb) * 1007
    
    Error_a=[]
    
    # Plot % contribution of uncertainty of factors to heat transfer uncertainty
    for (var_2, error_2) in Q_air.error_components().items():
        percent_contrib_2=(error_2**2)*100/(Q_air.s**2)
        Error_a.append(percent_contrib_2)
    
    
    # #Finding contribution of each sensor
    Contr_Q_a = []
    Contr_Q_a.append(Error_a[0]) #h_a_in
    Contr_Q_a.append(Error_a[1]) #h_a_out
    Contr_Q_a.append(Error_a[2] + Error_a[6]) #NDD
    Contr_Q_a.append(Error_a[3]) #mu
    Contr_Q_a.append(Error_a[4] + Error_a[7]) #rho
    Contr_Q_a.append(Error_a[5]) #deltaP
    
    #Plotting the contribution of all variables in air capacity uncertainty
    fig2, ax2 = plt.subplots()
    labels=['$h_{a_{in}}$', '$h_{a_{out}}$', '$NDD$', '$mu$', '$rho$', '$delta_P$']
    sizes=Contr_Q_a
    patches, texts = ax2.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, Error_a)]
    
    sort_legend = True
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.), fontsize=11)
    plt.savefig(plt_dir + 'Uncertainty analysis of Air Capacity',
                dpi = 300, bbox_inches = 'tight')
    
    #Plotting the contribution of all variables in air inlet enthalpy uncertainty
    fig8, ax8 = plt.subplots()
    labels=['$T_{db_{in}}$', '$T_{wb_{in}}$', '$P_{air}$']
    sizes=un_h_a_in[1:4]
    patches, texts = ax8.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, un_h_a_in[1:4])]
    
    sort_legend = True
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.),
               fontsize=11)
    plt.savefig(plt_dir + 'Uncertainty analysis of Air Inlet Enthalpy',
                dpi = 300,  bbox_inches = 'tight')
    
    #Plotting the contribution of all variables in air outlet enthalpy uncertainty
    fig9, ax9 = plt.subplots()
    labels=['$T_{db_{out}}$', '$T_{wb_{out}}$', '$P_{air}$']
    sizes=un_h_a_out[1:4]
    patches, texts = ax9.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, un_h_a_out[1:4])]
    
    sort_legend = True
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.),
               fontsize=11)
    plt.savefig(plt_dir + 'Uncertainty analysis of Air Outlet Enthalpy',
                dpi = 300,  bbox_inches = 'tight')
    
    """
    Sensible air plotting
    """
    Error_a_sens=[]
    
    # Plot % contribution of uncertainty of factors to heat transfer uncertainty
    for (var_3, error_3) in Q_air_sens.error_components().items():
        percent_contrib_3=(error_3**2)*100/(Q_air_sens.s**2)
        Error_a_sens.append(percent_contrib_3)
    
    
    # #Finding contribution of each sensor
    Contr_Q_a_s = []
    Contr_Q_a_s.append(Error_a_sens[1]) #T_db_in_avg
    Contr_Q_a_s.append(Error_a_sens[0]) #T_db_out_avg
    Contr_Q_a_s.append(Error_a_sens[2] + Error_a_sens[6]) #NDD
    Contr_Q_a_s.append(Error_a[3]) #mu
    Contr_Q_a_s.append(Error_a[4] + Error_a[7]) #rho
    Contr_Q_a_s.append(Error_a[5]) #deltaP
    
    #Plotting the contribution of all variables in air capacity uncertainty
    fig10, ax10 = plt.subplots()
    labels=['$T_{db_{in}}$', '$T_{db_{out}}$', '$NDD$', '$mu$', '$rho$', '$delta_P$']
    sizes=Contr_Q_a_s
    patches, texts = ax10.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, Error_a_sens)]
    
    sort_legend = True
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.), fontsize=11)
    plt.savefig(plt_dir + 'Uncertainty analysis of Air Sensible Capacity',
                dpi = 300, bbox_inches = 'tight')
    
    """
    Refrigerant
    """
    
    h_r_in = PropsSI('H', 'P', P_r_exv_in_avg, 'T', T_r_exv_in_avg, 'R410A')
    h_r_out = PropsSI('H', 'P', P_r_evap_out_avg, 'T', T_r_evap_out_avg, 'R410A')
    
    #Calculating uncertainty of enthanlpy calculation
    un_h_r_in = ref_props_un('H', T_r_exv_in_avg, un_T_r_exv_in, P_r_exv_in_avg, un_P_r_exv_in, 'R410A', thres=1e-5)
    un_h_r_out = ref_props_un('H', T_r_evap_out_avg, un_T_r_evap_out, P_r_evap_out_avg, un_P_r_evap_out, 'R410A', thres=1e-5)
    
    h_r_in_comb = ufloat(h_r_in, un_h_r_in[0], tag="h_r_in")
    h_r_out_comb = ufloat(h_r_out, un_h_r_out[0], tag="h_r_out")
    m_dot_r_comb = ufloat(m_dot_r_avg, un_m_dot_r, tag="m_dot_r")
    
    #Calculating uncertainty propagation in heat balance calculation
    Q_refr = m_dot_r_comb * abs(h_r_out_comb - h_r_in_comb)
    #Percentage contribution of variables in the uncertainty propagation on water heat balance
    Error_r=[]
    # Plot % contribution of uncertainty of factors to heat transfer uncertainty
    for (var, error) in Q_refr.error_components().items():
        percent_contrib=(error**2)*100/(Q_refr.s**2)
        Error_r.append(percent_contrib)
    
    # combine all the errors 
    Error_r=[Error_r[0],Error_r[1],Error_r[2]]
    
    #Plotting the contribution of all variables in refrigerant capacity uncertainty
    fig3, ax3 = plt.subplots()
    labels=['$h_{r_{in}}$', '$h_{r_{out}}$', '$\dot{m}_r$']
    sizes=Error_r
    patches, texts = ax3.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, Error_r)]
    
    sort_legend = True
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.),
               fontsize=11)
    plt.savefig(plt_dir + 'Uncertainty analysis of Refrigerant Capacity',
                dpi = 300,  bbox_inches = 'tight')
    
    #Plotting the contribution of all variables in refr. inlet enthalpy uncertainty
    fig6, ax6 = plt.subplots()
    labels=['$T_{r_{exv,in}}$', '$P_{r_{exv,in}}$']
    sizes=un_h_r_in[1:3]
    patches, texts = ax6.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, un_h_r_in[1:3])]
    
    sort_legend = True
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.),
               fontsize=11)
    plt.savefig(plt_dir +'Uncertainty analysis of Refr Inlet Enthalpy',
                dpi = 300,  bbox_inches = 'tight')
    
    #Plotting the contribution of all variables in refr. outlet enthalpy uncertainty
    fig7, ax7 = plt.subplots()
    labels=['$T_{r_{evap,out}}$', '$P_{r_{evap,out}}$']
    sizes=un_h_r_out[1:3]
    patches, texts = ax7.pie(sizes, shadow=False, startangle=90)
    
    labels = ['{0} - {1:1.3f} %'.format(i,j) for i,j in zip(labels, un_h_r_out[1:3])]
    
    sort_legend = True
    
    #put the legend outside the pie
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, labels),
                                              key=lambda x: x[2],
                                              reverse=True))
    
    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.5, 1.),
               fontsize=11)
    plt.savefig(plt_dir + 'Uncertainty analysis of Refr Outlet Enthalpy',
                dpi = 300,  bbox_inches = 'tight')
    
    print("Test File Name: ", Test_name)
    print("Capacity of Air in W: ", Q_air)
    print("Sensible Capacity of Air in W: ", Q_air_sens)
    print("Capacity of Water in W: ", Q_water)
    print("Capacity of Refrigerant in W: ", Q_refr)
    a_w_heat_bal = (abs(Q_air_sens.nominal_value) - abs(Q_water.nominal_value))*100/(0.5*(abs(Q_air_sens.nominal_value) + abs(Q_water.nominal_value)))
    a_r_heat_bal = (abs(Q_air_sens.nominal_value) - abs(Q_refr.nominal_value))*100/(0.5*(abs(Q_air_sens.nominal_value) + abs(Q_refr.nominal_value)))
    w_r_heat_bal = (abs(Q_water.nominal_value) - abs(Q_refr.nominal_value))*100/(0.5*(abs(Q_water.nominal_value) + abs(Q_refr.nominal_value)))
    a_w_r_heat_bal = ((abs(Q_air_sens.nominal_value) + abs(Q_water.nominal_value)) - abs(Q_refr.nominal_value))*100/(0.5*((abs(Q_air_sens.nominal_value) + abs(Q_water.nominal_value)) + abs(Q_refr.nominal_value)))
    print("Air-Water heat balance: ", a_w_heat_bal)
    print("Air-Refr heat balance: ", a_r_heat_bal)
    print("Water-Refr heat balance: ", w_r_heat_bal)
    print("Air+Water_Refr heat balance: ", a_w_r_heat_bal)
    
    """
    Writing to CSV file
    """
    # Open our existing CSV file in append mode
    # Create a file object for this file
    #data_file_loc = prj_fold_loc + "\\03-Trimmed Raw Data\\" + Raw_file_name
    
    result_file_loc = prj_fold_loc + "\\09-Results\\"
    result_files = result_file_loc + "/results.csv"
    with open(result_files, 'a') as f_object:
     
        # Pass this file object to csv.writer()
        # and get a writer object
        writer_object = writer(f_object)
     
        # Pass the list as an argument into
        # the writerow()
        List = [Test_name, 
                T_db_in_avg-273.115, 
                T_wb_in_avg-273.15,
                Rh_in_avg, 
                T_db_out_avg-273.15, 
                V_dot_a[0].nominal_value*2118.87997, 
                dp_a_avg/248.84, 
                (P_a_in_avg-P_atm_avg)/248.84, 
                P_atm_avg/6894.76, 
                T_w_in_avg-273.15, 
                (P_w_in_avg-P_atm_avg)/6894.76, 
                T_w_out_avg-273.15, 
                dp_w_avg/248.84, 
                m_dot_w_avg*3600, 
                T_r_exv_in_avg-273.15, 
                (P_r_exv_in_avg-P_atm_avg)/6894.76,
                T_r_evap_in_avg-273.15, 
                T_r_evap_out_avg-273.15, 
                (P_r_evap_out_avg-P_atm_avg)/6894.76, 
                m_dot_r_avg*3600, 
                Q_air, 
                Q_air_sens, 
                Q_water, 
                Q_refr, 
                a_w_heat_bal, 
                a_r_heat_bal, 
                w_r_heat_bal,
                a_w_r_heat_bal]
        writer_object.writerow(List)
        # Close the file object
        f_object.close()
