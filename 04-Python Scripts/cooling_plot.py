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
result_files = result_file_loc + "results_cooling.csv"
plot_loc_1 = result_file_loc + "cooling_air_cap.png"
plot_loc_2 = result_file_loc + "cooling_ref_cap.png"


df = pd.read_csv(result_files, encoding= 'unicode_escape')
    
"""
Editing matplotlibrc file
"""
mpl.rcParams['axes.labelsize'] = 'large' #Changing axes label size
mpl.rcParams['axes.labelpad'] = 12.0 #Changing axes label distance from axis
mpl.rcParams['axes.grid'] = 'true' #Showing grids
mpl.rcParams['grid.alpha'] = 0.5 #Changing grid transparency


x = df[["Tdb_in (°C)"]].to_numpy()
y_1 = abs(df[["Q_a_sens (kW)"]].to_numpy())
y_1_err = df[["un_Q_a_sens (kW)"]].to_numpy()
y_2 = abs(df[["Q_r (kW)"]].to_numpy())
y_2_err = df[["un_Q_r (kW)"]].to_numpy()
comp_speed = df[["Compressor Speed (RPM)"]].to_numpy()

legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor = 'red', markersize=8, label='Comp Speed = 3725 RPM'),
                   Line2D([0], [0], marker='s', color='w', markerfacecolor = 'green', markersize=8, label='Comp Speed = 2050 RPM')]

rows = np.size(x)
fig = plt.figure()
ax = fig.add_subplot(111)

for n in range(rows):
    if comp_speed[n]>3000:
       ax.scatter(x[n], y_1[n], marker = 'o', color = 'red', label='Comp Speed = 3725 RPM') #Added lable for the plot 
       
       ax.errorbar(x[n], y_1[n], yerr = y_1_err[n], capsize = 5.0, ecolor = 'red')
        
    elif comp_speed[n]<3000:
        ax.scatter(x[n], y_1[n], marker = 's', color = 'green', label='Comp Speed = 3725 RPM, Water mass flow = 1720 kg/h') #Added lable for the plot  
        
        ax.errorbar(x[n], y_1[n], yerr = y_1_err[n], capsize = 5.0, ecolor = 'green')
        
ax.set_xlabel('Air Inlet Dry Bulb Temperature (°C)')
ax.set_ylabel('Air Sensible Capacity (kW)')
#ax.legend(loc = 'lower right') #Addeng legend in lower right location
ax.legend(handles = legend_elements, loc='center')
plt.savefig(plot_loc_1,
            dpi = 300, bbox_inches = 'tight')   

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

for n in range(rows):
    if comp_speed[n]>3000 :
       ax2.scatter(x[n], y_2[n], marker = 'o', color = 'red', label='Comp Speed = 3725 RPM') #Added lable for the plot 
       
       #ax2.errorbar(x[n], y_2[n], yerr = y_2_err[n], capsize = 5.0, ecolor = 'red')

    elif comp_speed[n]<3000:
        ax2.scatter(x[n], y_2[n], marker = 's', color = 'green', label='Comp Speed = 3725 RPM, Water mass flow = 1720 kg/h') #Added lable for the plot
        
        #ax2.errorbar(x[n], y_2[n], yerr = y_2_err[n], capsize = 5.0, ecolor = 'green')
        
ax2.set_xlabel('Air Inlet Dry Bulb Temperature (°C)')
ax2.set_ylabel('Refrigerant Capacity (kW)')
#ax.legend(loc = 'lower right') #Addeng legend in lower right location
ax2.legend(handles = legend_elements, loc='center')
plt.savefig(plot_loc_2,
            dpi = 300, bbox_inches = 'tight')   





