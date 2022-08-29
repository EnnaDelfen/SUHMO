import sys,os
import csv
import numpy as np

## STORAGE
nb_cells = []
level_nb = []
L2_head = []
L2_GH = []
L2_Pw = []
L2_Re = []
L2_moulin = []
L2_Diff = []
L2_Channel = []

listOfFiles = []

## LIST OF FILES 
#with open("list_files.dat", 'r') as infile:
#    lines_interp = infile.readlines()
lines_interp = ["tmp/ERR_1lev_base2levs-4lev", 
                "tmp/ERR_2lev_base2levs-5lev",
                "tmp/ERR_3lev_base2levs-6lev",
                "tmp/ERR_4lev_base2levs-7lev"]
for i, line in enumerate(lines_interp[0:]):
    raw_line = line.strip()
    listOfFiles.append(raw_line)
    listOfFiles.append(raw_line+"_custom")
    
## PROCESS LIST OF FILES 
for i, currfile in enumerate(listOfFiles):
    print("Processing", currfile)
    with open(currfile, 'r') as infile:
        lines_interp = infile.readlines()
    for i, line in enumerate(lines_interp[0:]):
        raw_line = line.strip()
        if 'head' in raw_line:
            tmp_line = raw_line.split(',')[1]
            L2_head.append(float(tmp_line))
        if 'gapHeight' in raw_line:
            tmp_line = raw_line.split(',')[1]
            L2_GH.append(float(tmp_line))
        if 'Pw:' in raw_line:
            tmp_line = raw_line.split(',')[1]
            L2_Pw.append(float(tmp_line))
        if 'Re:' in raw_line:
            tmp_line = raw_line.split(',')[1]
            L2_Re.append(float(tmp_line))
        if 'Channelization:' in raw_line:
            tmp_line = raw_line.split(',')[1]
            L2_Channel.append(float(tmp_line))
        if 'DiffusiveTerm:' in raw_line:
            tmp_line = raw_line.split(',')[1]
            L2_Diff.append(float(tmp_line))
        if 'RHS_hmoul:' in raw_line:
            tmp_line = raw_line.split(',')[1]
            L2_moulin.append(float(tmp_line))
        if 'computed Filename' in raw_line:
            tmp_line = raw_line.split('/')[1]
            #tmp_line2 = tmp_line.split('=')[1]
            tmp_line3 = tmp_line.split('lev')[0]
            nb_cells.append(2**(float(tmp_line3)-1)*32)
            level_nb.append(float(tmp_line3))

## CHECKS
len_of_file = len(listOfFiles) /2.

if (len_of_file != len(L2_head)):
    print("Not enough head data", len_of_file, len(L2_head))
    #stop
if (len_of_file != len(L2_GH)):
    print("Not enough GH data", len_of_file, len(L2_GH))
    #stop
if (len_of_file != len(L2_Pw)):
    print("Not enough Pw data", len_of_file, len(L2_Pw))
    stop
if (len_of_file != len(L2_Re)):
    print("Not enough Re data", len_of_file, len(L2_Re))
    stop
if (len_of_file != len(L2_moulin)):
    print("Not enough MS data", len_of_file, len(L2_moulin))
    stop
if (len_of_file != len(L2_Diff)):
    print("Not enough DT data", len_of_file, len(L2_Diff))
    stop
if (len_of_file != len(L2_Channel)):
    print("Not enough CD data", len_of_file, len(L2_Channel))
    stop
if (2*len_of_file != len(nb_cells)):
    print("Not enough cell number data", len_of_file, len(nb_cells))
    stop

## PRINT OUT
file_out    = "convergence_data_3Levels.dat"
csv_file = str(file_out)
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile,delimiter=' ',quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["#errorL2"])
    writer.writerow(["#case","head","gapHeight","Pw","Re","RHS_moulin","DT","CD"])
    for i in range(len(L2_Channel)):
        writer.writerow([nb_cells[2*i+1], L2_head[i], L2_GH[i], L2_Pw[i], L2_Re[i], L2_moulin[i], L2_Diff[i], L2_Channel[i]])
