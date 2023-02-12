# -*- coding: utf-8 -*-
"""
Created on Sun Feb 27 21:43:39 2022

@author: farha
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
    
#Selecting csv file
root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
file_path = filedialog.askopenfilename(parent=root)

#Reading the file into numpy arrays
df = pd.read_csv(file_path)
m_dot_w = df[["Flow Meter (GL) (Kg/h)"]].to_numpy()
T_w_in = df[["Test Coil Inlet Temp (C)"]].to_numpy()
T_w_out = df[["Test Coil Outlet Temp (C)"]].to_numpy()
P_w_in = df[["Test Coil Inlet Pressure (PSIG)"]].to_numpy()
dp_w = df[["Test Coil Diff Pressure (in WC)"]].to_numpy()
T_db_in = df[["Test coil inlet Dry bulb (C)"]].to_numpy()
T_wb_in = df[["Test coil inlet Wet bulb (C)"]].to_numpy()
T_db_out = df[["Test coil outlet Dry bulb (C)"]].to_numpy()
T_wb_out = df[["Test coil outlet Wet bulb (C)"]].to_numpy()
P_a_in = df[["Nozzle Inlet Pressure (Pa)"]].to_numpy()
dp_a = df[["Nozzle Differential Pressure (Delta_P_a)"]].to_numpy()

#Mean values
m_dot_w_avg = np.mean(m_dot_w)
T_w_in_avg = np.mean(T_w_in)
T_w_out_avg = np.mean(T_w_out)
P_w_in_avg = np.mean(P_w_in)
T_db_in_avg = np.mean(T_db_in)
T_wb_in_avg = np.mean(T_wb_in)
T_db_out_avg = np.mean(T_db_out)
T_wb_out_avg = np.mean(T_wb_out)
P_a_in_avg = np.mean(P_a_in)
dp_a_avg = np.mean(dp_a)

#Variance
m_dot_w_var = np.var(m_dot_w)
T_w_in_var = np.var(T_w_in)
T_w_out_var = np.var(T_w_out)
P_w_in_var = np.var(P_w_in)
T_db_in_var = np.var(T_db_in)
T_wb_in_var = np.var(T_wb_in)
T_db_out_var = np.var(T_db_out)
T_wb_out_var = np.var(T_wb_out)
P_a_in_var = np.var(P_a_in)
dp_a_var = np.var(dp_a)

#Standard Deviation
m_dot_w_std = np.std(m_dot_w)
T_w_in_std = np.std(T_w_in)
T_w_out_std = np.std(T_w_out)
P_w_in_std = np.std(P_w_in)
T_db_in_std = np.std(T_db_in)
T_wb_in_std = np.std(T_wb_in)
T_db_out_std = np.std(T_db_out)
T_wb_out_std = np.std(T_wb_out)
P_a_in_std = np.std(P_a_in)
dp_a_std = np.std(dp_a)
