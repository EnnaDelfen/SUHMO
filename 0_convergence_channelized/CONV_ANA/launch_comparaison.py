import sys
import os
import shutil
from tempfile import mkstemp
import csv
import numpy as np

# exec 
chdiff="Users/ennadelfen/Documents/0_BERKELEY/0_CODES/Chombo_git/lib/util/ChomboCompare/compare2d.Darwin.64.mpic++.gfortran.DEBUG.MPI.ex"
# Files to modify
file_path ="/Users/ennadelfen/Documents/0_BERKELEY/0_CODES/SUHMO_git/0_convergence_channelized" 
# patterns to modify
searchFor1 = "exactRoot"
searchFor2 = "computedRoot"

patterns = ["1lev", "2lev", "3lev", "4lev", "5lev", "6lev", "7lev"]
LastIt   = ["3000", "3000", "3000", "3000", "3000", "3000", "3000"]

# put yourself in root 
os.chdir(file_path)
# process each inputs.compare
for i, pat in enumerate(patterns[0:-1]):
    print("Dealing with", pat, "and", patterns[i+1])
    tmpFile, abs_path = mkstemp()
    with open(tmpFile,'w') as new_file:
        with open(file_path+"/inputs.compare", 'r') as old_file:
            for line in old_file:
                if (searchFor1 in line):
                    new_file.write("compare.exactRoot =  "+patterns[i+1]+"/plot00"+LastIt[i+1]+".2d.hdf5")
                elif (searchFor2 in line): 
                    new_file.write("compare.computedRoot =  "+pat+"/plot00"+LastIt[i]+".2d.hdf5")
                else:
                    new_file.write(line)
    tmp_file_loc = file_path+"/tmp/inputs.compare_"+pat+"-"+patterns[i+1]
    print("... moving input file to", tmp_file_loc)
    shutil.move(abs_path, tmp_file_loc)
    os.system("/{} {}".format(chdiff, tmp_file_loc))
    shutil.move(file_path+"/pout.0", file_path+"/tmp/ERR_"+pat+"-"+patterns[i+1])
    print("... moving outpout file to", file_path+"/tmp/ERR_"+pat+"-"+patterns[i+1])


# process each inputs.compare_custom
for i, pat in enumerate(patterns[0:-1]):
    print("Dealing with", pat, "and", patterns[i+1])
    tmpFile, abs_path = mkstemp()
    with open(tmpFile,'w') as new_file:
        with open(file_path+"/inputs.compare_custom", 'r') as old_file:
            for line in old_file:
                if (searchFor1 in line):
                    new_file.write("compare.exactRoot =  "+patterns[i+1]+"/plotCustom00"+LastIt[i+1]+".2d.hdf5")
                elif (searchFor2 in line): 
                    new_file.write("compare.computedRoot =  "+pat+"/plotCustom00"+LastIt[i]+".2d.hdf5")
                else:
                    new_file.write(line)
    tmp_file_loc = file_path+"/tmp/inputs.compare_custom_"+pat+"-"+patterns[i+1]
    print("... moving input file to", tmp_file_loc)
    shutil.move(abs_path, tmp_file_loc)
    os.system("/{} {}".format(chdiff, tmp_file_loc))
    shutil.move(file_path+"/pout.0", file_path+"/tmp/ERR_"+pat+"-"+patterns[i+1]+"_custom")
    print("... moving outpout file to", file_path+"/tmp/ERR_"+pat+"-"+patterns[i+1]+"_custom")


                
