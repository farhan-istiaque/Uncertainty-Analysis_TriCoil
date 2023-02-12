# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:30:50 2022

@author: hamid
"""

import pandas as pd
import os
from pathlib import Path


#Standard Deviation and Mean of Dry Bulb Temperature
Recorded_Data_LabView = "TriCoil Test_Dec_06_2022_10_53_AM_trimmed.xlsx"


cols = ['m_dot_w (kg/h)', 'T_w_in (C)', 'T_w_out (C)', 'P_w_in (PSIG)',
            'dp_w (in W.G.)', 'Atm Pressure (PSI)', 'Pin_nozzle_TriCoil (PSI)', 
            'dP_nozzle_TriCoil (in W.G.)', 'Tdb_in (C)', 'Twb_in (C)', 'Tdb_out (C)',
            'Twb_out(C)', 'T_exv_in (C)', 'P_exv_in (PSIG)',
            'T_evap_out (C)', 'P_evap_out (PSIG)', 'm_dot_r (kg/h)']

scr_loc = os.path.dirname( __file__ ) #Script Location
prj_fold_loc = str(Path(scr_loc).parents[0]) #Project folder location
data_file_loc = prj_fold_loc + "\\03-Trimmed Raw Data\\" + Recorded_Data_LabView

df = pd.read_excel(data_file_loc, usecols=cols)
    
df_mean =[] # Creating empty lists to hold mean and random error values
df_rand_err = []
    
for cols in df:
            Rand_Err = 2*(df[cols].std())
            df_rand_err.append(Rand_Err) # Adding to the random error list
            Mean = df[cols].mean()
            df_mean.append(Mean) # Adding to the mean list
            print("Random Error of ", [cols], "is", Rand_Err, "\nMean value of ",
              [cols], "is =", Mean)
            syst = 0 # ADD VALUE FOR SYSTEMATIC ERROR HERE
            total_e = pow((pow(Rand_Err, 2) + pow(syst, 2)),0.5)
            print("Total Uncertainity of ", cols, "is", total_e)
            print("Resulting Total of ", cols, "is", Mean, u"\u00B1", total_e)
    
# Writing to Excel
df = pd.DataFrame({'Variable': df.columns,
                   'Mean': df_mean,
                   'Random Error' : df_rand_err})

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter("Random_Error_Hamid.xlsx", engine='xlsxwriter')

# Convert the dataframe to an XlsxWriter Excel object. Note that we turn off
# the default header and skip one row to allow us to insert a user defined
# header.
df.to_excel(writer, sheet_name='Sheet1', startrow=1, header=False, index=False)

# Get the xlsxwriter workbook and worksheet objects.
workbook  = writer.book
worksheet = writer.sheets['Sheet1']

# Add a header format.
header_format = workbook.add_format({
    'bold': True,
    'text_wrap': True,
    'valign': 'top',
    'fg_color': '#D7E4BC',
    'border': 1})

# Write the column headers with the defined format.
for col_num, value in enumerate(df.columns.values):
    worksheet.write(0, col_num, value, header_format)

# Close the Pandas Excel writer and output the Excel file.
writer.save()
