import matplotlib
import matplotlib.pyplot as plt
import pickle
import shutil
import math
import sys
import os
import commands
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

#name = sys.argv[1]
#pk_loc = sys.argv[2]
#total_reads = sys.argv[3]
#input arguments or static variables, for the name of the experiment, the names of the two samples being compared, the location
#of the majiq executable and gff3 file, number of threads, and the total number of reads of the sample being considered

exp_name= sys.argv[1]
subname1 = sys.argv[2]
subname2 = "Lung2.mm10.sorted"
majiq_loc = "/opt/majiq"
gff3_file = "/data/DB/mm10/ensembl.mm10.gff3"
nthreads = 8
preseq_loc = "/data/mrpaul/RED/preseq-master/"
fastq_loc = sys.argv[3]

##Calculating total amount of reads
print("Calculating total amount of reads")
#reads = commands.getoutput('wc -l ' + fastq_loc)
#total_reads = int(reads[:reads.find(" ")]) / 4.0
total_reads = 11743827.0
print("Total: " + str(total_reads))

print("Running majiq build")
name=exp_name
###Runs majiq build on the samples being considered
#os.system("export PYTHONPATH=" + majiq_loc)
#os.system(majiq_loc + "/majiq build " + gff3_file + " -conf " + name + ".config --nthreads " + str(nthreads) + " --output " + name + "_build")

print("Runnig majiq deltapsi")
###Runs majiq deltapsi
#os.system(majiq_loc + "/majiq deltapsi -grp1 " + name + "_build/" + subname1 + ".majiq -grp2 " + name + "_build/" \
#          + subname2 + ".majiq --minreads 2 --minpos 2 --nthreads " + str(nthreads) + " --output " + name + "_deltapsi --names Test1 Test2")
          
print("Parsing majiq deltapsi outout")          
###Loads the pickle file which contains the ID and number of reads for each junction found in that LSV
pp = pickle.load(open(name + "_deltapsi/clean_reads.Test1.pkl"))

###Parses pickle file to create a list of LSV IDs and their corresponding total read counts
lsvs = []
read_counts = []
distinct = 0
for x in pp:
    distinct += 1
    lsvs.append(x[0])
    read_counts.append(sum(x[1].tolist()))

lsv_reads = sum(read_counts)
print("Total LSV reads: " + str(lsv_reads))
print("Total Distinct LSV: " + str(distinct))
print("Calculating read count distribution")

###Calculates and stores histogram density values
###Uncomment following section to plot histogram
print("Max read count: " + str(max(read_counts)))
binsize = sys.argv[4]
counts, bins, bars = plt.hist(read_counts, bins=int(max(read_counts)/int(binsize)))
#plt.hist(read_counts, bins=4000, range=[0,1000])
#plt.xlabel("LSV Total Read Count")
#plt.ylabel("Frequency")
#plt.savefig(name + "_read_count_hist_zoomed.png")
hist_file_o = open(name + "_hist_values.txt","w")

for i in range(0,len(counts),1):
	hist_file_o.write(str(i+1) + "\t" + str(counts[i])+"\n")

hist_file_o.close()

###Extrapolates with preseq   
print("Running preseq")
preseq_command = preseq_loc + "preseq lc_extrap -o " + name + "_extra_bin" + "_" + str(binsize) + ".txt -H " + name + "_hist_values.txt"
os.system(preseq_command)

print("preseq completed")
###Opens the file that contains the extrapolated values from preseq 
scaled_c = open(name + "_extra_bin_" + str(binsize) + ".txt","r")

###Parses the parseq output 
ratio = total_reads/sum(read_counts)
print("LSV reads to Total Reads Ratio: " + str(ratio))

num_reads = []
num_distinct = [] 
low_bound = []
high_bound = []
scaled_c.next()
scaled = {}
for x in scaled_c:
    vals = x.strip("\n").split("\t") 
    num_reads.append(float(vals[0])*ratio)
    scaled[int(float(vals[0]))] = vals[1]
    num_distinct.append(vals[1])
    low_bound.append(vals[2])
    high_bound.append(float(vals[3]))

print int(lsv_reads/25000)*25000
print("Plotting extrapolated values")
###Plots the extrapolated values from preseq (only the first 150M reads)
plt.plot(num_reads, num_distinct)
plt.plot(num_reads, low_bound, "b--")
plt.plot(num_reads, high_bound, "b--")
plt.ylabel("# of Distinct LSVs (Unscaled)")
plt.xlabel("# of Reads")
plt.title("Extrapolated Complexity Plot")
plt.axis([0,150000000, 0, 1.10*max(high_bound)])
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(name +"_bin" + str(binsize) + "_extrapolated.png")
plt.close()



for population in [150000,200000,250000,300000,400000,500000]:

	ratio2 = float(scaled[population])/distinct
	print ratio2
        new_histFile = open(name + "_at_" + str(population) + "_lsv_reads.txt", "w")
	binBracket = 2
        for i in range(0,len(counts),1):
                new_histFile.write(str(binBracket) + "-" + str(binBracket+3)  + "\t" + str(counts[i] * ratio2) + "\n")
		binBracket += 4
        new_histFile.close()    








    
