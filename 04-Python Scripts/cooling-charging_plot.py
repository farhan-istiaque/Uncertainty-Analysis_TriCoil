# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 12:36:45 2023

@author: farhan
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import glob
from pathlib import Path
from csv import writer
import matplotlib as mpl
from matplotlib.lines import Line2D

"""
Reading the file into numpy arrays
"""
scr_loc = os.path.dirname( __file__ ) #Script Location
prj_fold_loc = str(Path(scr_loc).parents[0]) #Project folder location
result_file_loc = prj_fold_loc + "\\09-Results\\plot_results\\"
result_files = result_file_loc + "results_cc.csv"
plot_loc_1 = result_file_loc + "cc_water capacity.png"
plot_loc_2 = result_file_loc + "cc_refr capacity.png"
plot_loc_3 = result_file_loc + "cc_air sensible capacity.png"

df = pd.read_csv(result_files, encoding= 'unicode_escape')
    
"""
Editing matplotlibrc file
"""
mpl.rcParams['axes.labelsize'] = 'large' #Changing axes label size
mpl.rcParams['axes.labelpad'] = 12.0 #Changing axes label distance from axis
mpl.rcParams['axes.grid'] = 'true' #Showing grids
mpl.rcParams['grid.alpha'] = 0.5 #Changing grid transparency


x = df[["Tdb_in (°C)"]].to_numpy()
y_1 = abs(df[["Q_w (kW)"]].to_numpy())
y_1_err = df[["un_Q_w (kW)"]].to_numpy()
y_2 = abs(df[["Q_r (kW)"]].to_numpy())
y_2_err = df[["un_Q_r (kW)"]].to_numpy()
y_3 = abs(df[["Q_a_sens (kW)"]].to_numpy())
y_3_err = df[["un_Q_a_sens (kW)"]].to_numpy()

legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor = 'red', markersize=8, label='Inlet Water Temp = 30°C, Water mass flow = 1700 kg/h'),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor = 'green', markersize=8, label='Inlet Water Temp = 30°C, Water mass flow = 800 kg/h'),
                   Line2D([0], [0], marker='s', color='w', markerfacecolor = 'green', markersize=8, label='Inlet Water Temp = 27°C, Water mass flow = 800 kg/h'),
                   Line2D([0], [0], marker='s', color='w', markerfacecolor = 'red', markersize=8, label='Inlet Water Temp = 27°C, Water mass flow = 1700 kg/h')]

mass_flow = df[["m_w (kg/h)"]].to_numpy()
inlet_water = df[["T_w_in (°C)"]].to_numpy()
rows = np.size(x)

fig = plt.figure()
ax = fig.add_subplot(111)

for n in range(rows):
    if inlet_water[n]>29 and mass_flow[n]>1600:
       ax.scatter(x[n], y_1[n], marker = 'o', color = 'red', label='Inlet Water Temp = 30°C, Water mass flow = 1700 kg/h') #Added lable for the plot 
       
       ax.errorbar(x[n], y_1[n], yerr = y_1_err[n], capsize = 5.0, ecolor = 'red')
       
    elif inlet_water[n]>29 and mass_flow[n]<1600:
        ax.scatter(x[n], y_1[n], marker = 'o', color = 'green', label='Inlet Water Temp = 30°C, Water mass flow = 800 kg/h') #Added lable for the plot
        
        ax.errorbar(x[n], y_1[n], yerr = y_1_err[n], capsize = 5.0, ecolor = 'green')
        
    elif inlet_water[n]<29 and mass_flow[n]<1600:
        ax.scatter(x[n], y_1[n], marker = 's', color = 'green', label='Inlet Water Temp = 27°C, Water mass flow = 800 kg/h') #Added lable for the plot  
        
        ax.errorbar(x[n], y_1[n], yerr = y_1_err[n], capsize = 5.0, ecolor = 'green')
        
    elif inlet_water[n]<29 and mass_flow[n]>1600:
        ax.scatter(x[n], y_1[n], marker = 's', color = 'red', label='Inlet Water Temp = 27°C, Water mass flow = 1700 kg/h') #Added lable for the plot  
        
        ax.errorbar(x[n], y_1[n], yerr = y_1_err[n], capsize = 5.0, ecolor = 'red')
        
ax.set_xlabel('Inlet Air Dry Bulb Temperature (°C)')
ax.set_ylabel('Water Capacity (kW)')
ax.legend(handles = legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.25))
plt.savefig(plot_loc_1,
            dpi = 300, bbox_inches = 'tight')   


fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

for n in range(rows):
    if inlet_water[n]>29 and mass_flow[n]>1600:
       ax2.scatter(x[n], y_2[n], marker = 'o', color = 'red', label='Inlet Water Temp = 30°C, Water mass flow = 1720 kg/h') #Added lable for the plot 
       
       #ax2.errorbar(x[n], y_2[n], yerr = y_2_err[n], capsize = 5.0, ecolor = 'red')
       
    elif inlet_water[n]>29 and mass_flow[n]<1600:
        ax2.scatter(x[n], y_2[n], marker = 'o', color = 'green', label='Inlet Water Temp = 30°C, Water mass flow = 800 kg/h') #Added lable for the plot 
        
        #ax2.errorbar(x[n], y_2[n], yerr = y_2_err[n], capsize = 5.0, ecolor = 'green')
        
    elif inlet_water[n]<29 and mass_flow[n]<1600:
        ax2.scatter(x[n], y_2[n], marker = 's', color = 'green', label='Inlet Water Temp = 27°C, Water mass flow = 800 kg/h') #Added lable for the plot
        
        #ax2.errorbar(x[n], y_2[n], yerr = y_2_err[n], capsize = 5.0, ecolor = 'green')
        
    elif inlet_water[n]<29 and mass_flow[n]>1600:
        ax2.scatter(x[n], y_2[n], marker = 's', color = 'red', label='Inlet Water Temp = 27°C, Water mass flow = 1700 kg/h') #Added lable for the plot
        
        #ax2.errorbar(x[n], y_2[n], yerr = y_2_err[n], capsize = 5.0, ecolor = 'red')
        
ax2.set_xlabel('Inlet Air Dry Bulb Temperature (°C)')
ax2.set_ylabel('Refrigerant Capacity (kW)')
ax2.legend(handles = legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.25)) 
plt.savefig(plot_loc_2,
            dpi = 300, bbox_inches = 'tight')   

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)

for n in range(rows):
    if inlet_water[n]>29 and mass_flow[n]>1600:
       ax3.scatter(x[n], y_3[n], marker = 'o', color = 'red', label='Inlet Water Temp = 30°C, Water mass flow = 1700 kg/h') #Added lable for the plot 
       
       ax3.errorbar(x[n], y_3[n], yerr = y_3_err[n], capsize = 5.0, ecolor = 'red')
       
    elif inlet_water[n]>29 and mass_flow[n]<1600:
        ax3.scatter(x[n], y_3[n], marker = 'o', color = 'green', label='Inlet Water Temp = 30°C, Water mass flow = 800 kg/h') #Added lable for the plot 
        
        ax3.errorbar(x[n], y_3[n], yerr = y_3_err[n], capsize = 5.0, ecolor = 'green')
        
    elif inlet_water[n]<29 and mass_flow[n]<1600:
        ax3.scatter(x[n], y_3[n], marker = 's', color = 'green', label='Inlet Water Temp = 27°C, Water mass flow = 800 kg/h') #Added lable for the plot
        
        ax3.errorbar(x[n], y_3[n], yerr = y_3_err[n], capsize = 5.0, ecolor = 'green')
        
    elif inlet_water[n]<29 and mass_flow[n]>1600:
        ax3.scatter(x[n], y_3[n], marker = 's', color = 'red', label='Inlet Water Temp = 27°C, Water mass flow = 1700 kg/h') #Added lable for the plot
        
        ax3.errorbar(x[n], y_3[n], yerr = y_3_err[n], capsize = 5.0, ecolor = 'red')
        
ax3.set_xlabel('Inlet Air Dry Bulb Temperature (°C)')
ax3.set_ylabel('Air Sensible Capacity (kW)')
ax3.legend(handles = legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.25)) 
plt.savefig(plot_loc_3,
            dpi = 300, bbox_inches = 'tight')   







