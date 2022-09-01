#!/usr/bin/env python3

#####################
# A template script to lunch several suhmo runs 
# User should be at the root of the folder 
# Usage:
#   python multi_runs.py 
#####################
import sys
import os
import shutil
import csv
import numpy as np

# Paths to root and executable definition
root_dir = os.getcwd() 
for f in os.listdir(root_dir):
    if ( f.startswith("Suhmo2d") and f.endswith(".ex")):
        executable = f
# Go into root 
os.chdir(root_dir)

# Cases we want to run
patterns = ["run_C_starter", "run_C_1lev", "run_C_2lev", "run_C_3lev", "run_C_4lev", "run_C_5lev"]
# Process each runs -- run test case 
for i, pat in enumerate(patterns):
    print("Running case ", pat)
    exec_dir = root_dir + "/{}/".format(pat)
    os.chdir(exec_dir)
    os.system("mpirun -n 4 ../{} {}".format(executable, "input.hydro"))
    shutil.move(exec_dir+"pout.0", exec_dir+"{}".format("run1.output") )
    print(".. restart with bigger dt ")
    os.system("mpirun -n 4 ../{} {}".format(executable, "input.hydro_restart"))
    shutil.move(exec_dir+"pout.0", exec_dir+"{}".format("run2.output") )
    #os.system("python ../dump_PP.py {} > postproc.dat".format(pat)) 
