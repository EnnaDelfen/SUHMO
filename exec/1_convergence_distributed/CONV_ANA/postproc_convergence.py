import sys,os
import csv
import numpy as np

## STORAGE
nb_cells = []
L2_head = []
L2_GH = []
L2_Pw = []
L2_Re = []

listOfFiles = []

## LIST OF FILES 
#with open("list_files.dat", 'r') as infile:
#    lines_interp = infile.readlines()
lines_interp = ["../tmp/ERR_128x32-256x64",
                "../tmp/ERR_256x64-512x128",
                "../tmp/ERR_512x128-1024x256",
                "../tmp/ERR_64x16-128x32"]
for i, line in enumerate(lines_interp[0:]):
    raw_line = line.strip()
    listOfFiles.append(raw_line)
    
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
        if 'computed Filename' in raw_line:
            tmp_line  = raw_line.split('/')[0]
            tmp_line2 = tmp_line.split('=')[1]
            tmp_line3 = tmp_line2.split('x')[0]
            nb_cells.append(tmp_line3)

## CHECKS
len_of_file = len(listOfFiles) 

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
if (len_of_file != len(nb_cells)):
    print("Not enough cell number data", len_of_file, len(nb_cells))
    stop

## PRINT OUT
file_out    = "convergence_data.dat"
csv_file = str(file_out)
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile,delimiter=' ',quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["#errorL2"])
    writer.writerow(["#case","head","gapHeight","Pw","Re"])
    for i in range(len(L2_head)):
        writer.writerow([nb_cells[i], L2_head[i], L2_GH[i], L2_Pw[i], L2_Re[i]])
