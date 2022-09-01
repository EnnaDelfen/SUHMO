#!/usr/bin/env python3

#####################
# A template script to launch the postproc of hdf5 result plots from SUHMO
# with different space resolution in order to evaluate the convergence order.
# USES CHOMBO's chdiff tool !!
# User should be in the CONV_ANA folder
# Usage:
#   python scripts/launch_comparaison.py 
#####################
import sys
import os
import shutil
from tempfile import mkstemp
import csv
import numpy as np

# Path to CHOMBO's compare tool
CHDIR = os.getenv('CHOMBO_HOME')
chdiff_dir = CHDIR+"util/ChomboCompare/"
for f in os.listdir(chdiff_dir):
    if ( f.startswith("compare2d") and f.endswith(".ex")):
        executable = chdiff_dir+f 

# Paths to root and convergence dir
convergence_dir = os.getcwd() 
root_dir = os.getcwd() + "/../"
# Go into CONV_ANA folder
os.chdir(convergence_dir)   

# list of runs to post process with compare tool
patterns = ["1lev", "2lev", "3lev", "4lev", "5lev", "6lev", "7lev"]
LastIt   = ["12200", "12200", "12200", "12200", "12200", "12200", "12200"]
LastItSL = ["3200", "3200", "3200", "3200", "3200", "3200", "1600"]
# patterns to modify in the input file, inputs.compare
searchFor1 = "exactRoot"
searchFor2 = "computedRoot"
# process each inputs.compare
for i, pat in enumerate(patterns[0:-3]):
    print("Dealing with", pat+"_base2levs", "and", patterns[i+3])
    tmpFile, abs_path = mkstemp()
    with open(tmpFile,'w') as new_file:
        with open("inputs.compare", 'r') as old_file:
            for line in old_file:
                if (searchFor1 in line):
                    new_file.write("compare.exactRoot =  ../"+patterns[i+3]+"/plot00"+LastItSL[i+3]+".2d.hdf5")
                elif (searchFor2 in line): 
                    new_file.write("compare.computedRoot =  ../"+pat+"_base2levs/plot0"+LastIt[i]+".2d.hdf5")
                else:
                    new_file.write(line)
    tmp_file_loc = convergence_dir+"/tmp/inputs.compare_"+pat+"_base2levs-"+patterns[i+3]
    print("... moving input file to", tmp_file_loc)
    shutil.move(abs_path, tmp_file_loc)
    os.system("/{} {}".format(executable, tmp_file_loc))
    shutil.move("pout.0", "./tmp/ERR_"+pat+"_base2levs-"+patterns[i+3])
    print("... moving outpout file to", convergence_dir+"/tmp/ERR_"+pat+"_base2levs-"+patterns[i+3])

# process each inputs.compare_custom
for i, pat in enumerate(patterns[0:-3]):
    print("Dealing with", pat+"_base2levs", "and", patterns[i+3])
    tmpFile, abs_path = mkstemp()
    with open(tmpFile,'w') as new_file:
        with open("inputs.compare_custom", 'r') as old_file:
            for line in old_file:
                if (searchFor1 in line):
                    new_file.write("compare.exactRoot =  ../"+patterns[i+3]+"/plotCustom00"+LastItSL[i+3]+".2d.hdf5")
                elif (searchFor2 in line): 
                    new_file.write("compare.computedRoot =  ../"+pat+"_base2levs/plotCustom0"+LastIt[i]+".2d.hdf5")
                else:
                    new_file.write(line)
    tmp_file_loc = convergence_dir+"/tmp/inputs.compare_custom_"+pat+"_base2levs-"+patterns[i+3]
    print("... moving input file to", tmp_file_loc)
    shutil.move(abs_path, tmp_file_loc)
    os.system("/{} {}".format(executable, tmp_file_loc))
    shutil.move("pout.0", "tmp/ERR_"+pat+"_base2levs-"+patterns[i+3]+"_custom")
    print("... moving outpout file to", convergence_dir+"/tmp/ERR_"+pat+"_base2levs-"+patterns[i+3]+"_custom")

