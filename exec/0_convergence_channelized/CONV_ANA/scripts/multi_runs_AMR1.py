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
patterns = ["1lev_base", "2lev_base", "3lev_base", "4lev_base", "5lev_base"]
# Process each runs -- run test case AND PostProc
for i, pat in enumerate(patterns):
    print("Running case ", pat)
    exec_dir = root_dir + "/{}/".format(pat)
    os.chdir(exec_dir)
    os.system("mpirun -n 4 ../{} {}".format(executable, "input.hydro"))
    shutil.move(exec_dir+"pout.0", exec_dir+"{}".format("run.output") )
    print(".. running pp too ")
    os.system("../{} {}".format(executable, "input.hydro_pp"))
    shutil.move(exec_dir+"pout.0", exec_dir+"{}".format("postproc.output") )
    os.system("rm pout.1 pout.2 pout.3 plot003000.2d.hdf5 plot007000.2d.hdf5")
