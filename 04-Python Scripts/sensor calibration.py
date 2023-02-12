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
_____________________________________________________________
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import os
from pathlib import Path

#df = pd.read_csv("RTD Calibration.csv") #reading csv file wtih rtd data and 
scr_loc = os.path.dirname( __file__ ) #Script Location
prj_fold_loc = str(Path(scr_loc).parents[0]) #Project folder location
"""
Change the title of calibration file here
"""
data_file_loc = prj_fold_loc + "\\07-Calibration Data\\RTD calibration.csv"

df = pd.read_csv(data_file_loc)
#reference rtd data
n = len(df) #number of rows

"""
_____________________________________________________________
Select sensor name
"""
rtd = "Tdb_in" 
"""
_____________________________________________________________
"""

x_fit = df[rtd].to_numpy() #selected sensor readings
y_fit = df["Refr. RTD"].to_numpy() #Reference RTD readings

plt.scatter(x_fit, y_fit) #plot actual output from sensor and Refr. sensor

def curvefit_func(x, a, b):
#linear equation for curvefitting
    return a*x + b

coeff, cov = curve_fit(curvefit_func, x_fit, y_fit) #curve fitting the actual sensor
#output to the reference RTD output

print("{} * {} + ({})".format(coeff[0], rtd, coeff[1])) #printing the curve 
#fit eqn

y_calib = curvefit_func(x_fit, coeff[0], coeff[1]) #curve fitted sensor output

plt.plot(y_calib, y_fit) #plotting the curve fitted output

print("r^2 error: ", r2_score(y_fit, y_calib)) #calculate and pring r^2 score

delta_y = y_fit - y_calib

error = np.mean(delta_y)