import numpy as np
import sys

######## IMPORT STUFF
count_PP    = 0
count_lines = 0
line_start  = -1
line_end    = -1
case = sys.argv[1]
with open("postproc.output", 'r') as reader:
    lines=reader.readlines()
    for line in lines:
        if "POST PROC analysis" in line:
            line_end = int((line.split(")")[0]).split(" ")[-1])
        if "XaxisSUHMO_" in line:
            if (count_PP>0):
                line_start = count_lines
            else:
                count_PP += 1
        count_lines += 1

    Header = lines[line_start].split()
    for i in range(len(Header)):
        Header[i] = Header[i]+case
    print("#",*Header)
    for i in range(line_start+1,line_start+line_end+1):
        print(lines[i])
