import sys
import os
import shutil
from tempfile import mkstemp
import csv
import numpy as np

# Path to CHOMBO's compare tool
CHDIR = os.getenv('CHOMBO_HOME')
chdiff_dir = CHDIR+"/util/ChomboCompare/"
for f in os.listdir(chdiff_dir):
    if ( f.startswith("compare2d") and f.endswith(".ex")):
        chdiff = chdiff_dir+f 

# cv and root dir
convergence_dir = os.getcwd() 
root_dir = os.getcwd() + "/../"
os.chdir(root_dir)   

patterns = ["64x16","128x32","256x64","512x128","1024x256"]
LastIt   = ["5000", "5000", "5000", "5000", "5000"]
# patterns to modify
searchFor1 = "exactRoot"
searchFor2 = "computedRoot"

# process each inputs.compare
for i, pat in enumerate(patterns[0:-1]):
    print("Dealing with", pat, "and", patterns[i+1])
    tmpFile, abs_path = mkstemp()
    with open(tmpFile,'w') as new_file:
        with open(root_dir+"/CONV_ANA/inputs.compare", 'r') as old_file:
            for line in old_file:
                if (searchFor1 in line):
                    new_file.write("compare.exactRoot =  "+patterns[i+1]+"/plot00"+LastIt[i+1]+".2d.hdf5")
                elif (searchFor2 in line): 
                    new_file.write("compare.computedRoot =  "+pat+"/plot00"+LastIt[i]+".2d.hdf5")
                else:
                    new_file.write(line)
    tmp_file_loc = root_dir+"/tmp/inputs.compare_"+pat+"-"+patterns[i+1]
    print("... moving input file to", tmp_file_loc)
    shutil.move(abs_path, tmp_file_loc)
    os.system("/{} {}".format(chdiff, tmp_file_loc))
    shutil.move(root_dir+"/pout.0", root_dir+"/tmp/ERR_"+pat+"-"+patterns[i+1])
    print("... moving outpout file to", root_dir+"/tmp/ERR_"+pat+"-"+patterns[i+1])

