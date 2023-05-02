#!/usr/bin/env python3

#####################
# A template script to lunch several suhmo Suite "A" runs 
# User should be at the root of the A_SHMIP folder 
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
patterns = ["A1", "A2", "A3", "A4", "A5", "A6"]
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
    print(".. analyze pp and produce datafile ")
    os.system("python ../dump_PP.py {} > postproc.dat".format(pat)) 
    os.system("rm pout.* plot000000.2d.hdf5") 
