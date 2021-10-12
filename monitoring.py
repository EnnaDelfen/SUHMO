import numpy as np
import sys

print("#Dealing with file", sys.argv[1])
print("#it max min")
######## IMPORT STUFF
Content = []
ItNum = []
Bmax  = []
Bmin  = []
i = int(sys.argv[2])
count = 0
with open(sys.argv[1], 'r') as reader:
    for line in reader.readlines():
        Content.append(line.split()[0])
        ItNum.append(i)
        max_arg = line.split("max =")[1].split(",")[0]
        Bmax.append(np.float(max_arg))
        min_arg = line.split("min =")[1].split(")")[0]
        Bmin.append(np.float(min_arg))
        i+=1
        count+=1

for j in range(count):
    print(ItNum[j], Bmax[j], Bmin[j]) 
        

