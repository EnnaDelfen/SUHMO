#!/usr/bin/env python3

#####################
# A template script to lunch several suhmo chanelizing runs with different space resolution
# in order to evaluate the convergence order.
# User should be in the CONV_ANA folder
# Usage:
#   python scripts/multi_runs.py 
#####################
import sys
import os
import shutil
import csv
import numpy as np

# Paths to root and executable definition
convergence_dir = os.getcwd()
root_dir = os.getcwd() + "/../" 
for f in os.listdir(root_dir):
    if ( f.startswith("Suhmo2d") and f.endswith(".ex")):
        executable = f
# Go into root 
os.chdir(root_dir)

# Cases we want to run
patterns = ["64x16","128x32","256x64","512x128","1024x256"]
# Process each runs -- run test case
for i, pat in enumerate(patterns):
    print("Running case ", pat)
    exec_dir = root_dir + "/{}/".format(pat)
    os.chdir(exec_dir)
    os.system("mpirun -n 4 ../{} {}".format(executable, "input.hydro"))
    shutil.move(exec_dir+"pout.0", exec_dir+"{}".format("run.output") )
    os.system("rm pout.1 pout.2 pout.3")
