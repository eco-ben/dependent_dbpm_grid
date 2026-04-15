#!/usr/bin/env python3

# Loading libraries ----
import os
import numpy as np
import useful_functions as uf

# Define base variables ----
#Model resolution
res = '1deg'

#FAO region
reg = 'fao-47'

#Years to be processed - If not processing all files
yrs = np.arange(1741, 2011)

#Defining input and output folders
gridded_outputs = f'/g/data/vf71/fishmip_outputs/ISIMIP3a/fao_outputs/{reg}/fishing_runs-smoothed/{res}'

#Getting list of variables
varnames = np.unique([f.split(f'_{res}')[0] for f in os.listdir(gridded_outputs)])

#Applying merge function to all gridded outputs
for var in varnames:
    uf.merge_files(var, gridded_outputs, merge_by = 'year', years = yrs)
