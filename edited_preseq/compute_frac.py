

import sys

fileName = sys.argv[1]
total_reads = sys.argv[2].split(",")

histFile = open(fileName, "r")

total_counts = 0
counts = []
for row in histFile:
	line = row.strip("\n").split("\t")
	counts.append([int(line[0])+1,float(line[1])])
	total_counts += float(line[1])*(int(line[0])+1)

print total_counts

fracFile = open(fileName[:fileName.find("_")] + "_fracs.txt", "w")
for i in range(0, len(counts), 1):
	fracFile.write(str(counts[i][0]) + "\t" + str(counts[i][1]/total_counts) + "\n")
	counts[i].append(counts[i][1]/total_counts)
fracFile.close()

for total in total_reads:
	new_histFile = open(fileName[:fileName.find("_")] + "_at_" + str(total) + "_reads.txt", "w")
	for i in range(0,len(counts),1):
		new_histFile.write(str(counts[i][0]) + "\t" + str(int(counts[i][2] * int(total))) + "\n")
	new_histFile.close()
