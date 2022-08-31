import numpy as np
import sys

######## IMPORT STUFF
T_in_h = []
T_in_d = []
avgN = []
N_LB = []
N_MB = []
N_HB = []
rech = []
disch = []
case = sys.argv[1]
Header = ["#T_hrs_"+case, "T_days_"+case, "avgN_"+case, "N_LB_"+case, "N_MB_"+case, "N_HB_"+case, "rech_"+case, "dis_"+case ]
with open("run.output", 'r') as reader:
    lines=reader.readlines()
    for line in lines:
        if "avgN N_LB  N_M" in line:
            line_striped = (line.split("=")[1]).split("\r")[0]
            data = line_striped.split()
            T_in_h.append(float(data[0]))
            T_in_d.append(float(data[1]))
            avgN.append(float(data[2]))
            N_LB.append(float(data[3]))
            N_MB.append(float(data[4]))
            N_HB.append(float(data[5]))
        if "rech dis =" in line:
            line_striped = (line.split("=")[1]).split("\r")[0]
            data = line_striped.split()
            rech.append(float(data[2]))
            disch.append(float(data[3]))

    print(*Header)
    for i in range(len(T_in_h)):
        print(T_in_h[i], T_in_d[i], avgN[i], N_LB[i], N_MB[i], N_HB[i], rech[i], disch[i])
