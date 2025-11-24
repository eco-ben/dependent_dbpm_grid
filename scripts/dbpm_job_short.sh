#!/bin/bash

#PBS -P vf71
#PBS -q normalbw
#PBS -l ncpus=14
#PBS -l mem=126GB
#PBS -l walltime=00:10:00
#PBS -l storage=gdata/vf71+gdata/xp65
#PBS -M lilian.fierroarcos@utas.edu.au
#PBS -m abe
#PBS -l wd

module use /g/data/xp65/public/modules
module load conda/analysis3-25.11
python3 02_extracting_DBPM_regional_inputs.py
