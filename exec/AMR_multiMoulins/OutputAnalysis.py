import numpy as np
import sys


cells_lev_0 = []
cells_lev_1 = []
cells_lev_2 = []
cells_lev_3 = []
cells_lev_4 = []

cells_advanced_per_timestep = []

lev0_count = 0
lev1_count = 0
lev2_count = 0
lev3_count = 0
lev4_count = 0

with open(sys.argv[1], 'r') as reader:
    for line in reader.readlines(): 
        line_split = line.split()
        if ("level 0 cells advanced" in line):
            cells_lev_0.append(float(line_split[-1])) 
            lev0_count+=1
        if ("level 1 cells advanced" in line):
            cells_lev_1.append(float(line_split[-1])) 
            lev1_count+=1
        if ("level 2 cells advanced" in line):
            cells_lev_2.append(float(line_split[-1])) 
            lev2_count+=1
        if ("level 3 cells advanced" in line):
            cells_lev_3.append(float(line_split[-1])) 
            lev3_count+=1
        if ("level 4 cells advanced" in line):
            cells_lev_4.append(float(line_split[-1])) 
            lev4_count+=1

for i in range(lev0_count):
    cells_advanced_per_timestep.append(cells_lev_0[i])

if (lev1_count > 0):
    if (lev0_count!=lev1_count):
        print("problem: level 1 has advanced less ?", lev0_count, lev1_count)
    else:
        for i in range(lev0_count):  
            cells_advanced_per_timestep[i] +=cells_lev_1[i]

if (lev2_count > 0):
    if (lev0_count!=lev2_count):
        print("problem: level 2 has advanced less ?", lev0_count, lev2_count)
    else:
        for i in range(lev0_count):  
            cells_advanced_per_timestep[i] +=cells_lev_2[i]

if (lev3_count > 0):
    if (lev0_count!=lev3_count):
        print("problem: level 3 has advanced less ?", lev0_count, lev3_count)
    else:
        for i in range(lev0_count):  
            cells_advanced_per_timestep[i] +=cells_lev_3[i]

if (lev4_count > 0):
    if (lev0_count!=lev4_count):
        print("problem: level 1 has advanced less ?", lev0_count, lev4_count)
    else:
        for i in range(lev0_count):  
            cells_advanced_per_timestep[i] +=cells_lev_4[i]

total_average_cells_advanced = 0
for i in range(lev0_count):
    total_average_cells_advanced += cells_advanced_per_timestep[i]
total_average_cells_advanced = total_average_cells_advanced/lev0_count

print("Averaged number of cells advanced in this run is: ", total_average_cells_advanced)
