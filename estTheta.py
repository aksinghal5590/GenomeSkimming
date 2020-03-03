import os, re

#f = open("./ri_master.vcf")
#f = open(os.getcwd() + "/Desktop/CSE280A/smaller.txt")
f = open(os.getcwd() + "/Desktop/CSE280A/ri_master.vcf")
miDict = {}
for line in f:
    parsed = line.strip().replace('|', '\t').split('\t')
    if len(parsed) == 29+20 and parsed[0] != "#CHROM":
        values = list(map(int, parsed[9:]))
        entry1 = sum(values)
        if miDict.get(entry1) is not None:
            miDict[entry1] = miDict[entry1] + 1
        else:
            miDict[entry1] = 1


x = []
y = []         
for key in miDict.keys():
    val = miDict[key]
    x.append(key)
    y.append(key*val)
    
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.bar(x, y)
plt.show()

print(mean(y))

## calculate waterson's est for theta
#for i in miDict.keys():
#    for s in range(1, )