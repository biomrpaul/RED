import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle as pkl
import shutil
import sys
import os
import commands
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import sys
sys.path.append("/opt/majiq")
from voila import vlsv

exp_name= sys.argv[1]
subname1 = sys.argv[2]
subname2 = sys.argv[3]
majiq_loc = "/opt/majiq"
gff3_file = "/data/DB/mm10/ensembl.mm10.gff3"
nthreads = 8
preseq_loc = "/data/mrpaul/RED/preseq-master/"
fastq_loc = sys.argv[4]
run_majiq = sys.argv[6]
##Calculating total amount of reads
count_reads = sys.argv[7]
run_build = sys.argv[8]

if count_reads == "count_reads":
	print("Calculating total amount of reads")
	reads = commands.getoutput('wc -l ' + fastq_loc)
	total_reads = int(reads[:reads.find(" ")]) / 4.0
else:	
	total_reads = 11743827.0
print("Total: " + str(total_reads))
name=exp_name

if run_majiq == "run_majiq":
	if run_build == "run_build":
		print("Running majiq build")
		###Runs majiq build on the samples being considered
		os.system("export PYTHONPATH=" + majiq_loc)
		os.system(majiq_loc + "/majiq build " + gff3_file + " -conf " + name + ".config --nthreads " + str(nthreads) + " --output " + name + "_build")

	print("Runnig majiq deltapsi")
###Runs majiq deltapsi
	os.system(majiq_loc + "/majiq deltapsi -grp1 " + name + "_build/" + subname1 + ".majiq -grp2 " + name + "_build/" \
          + subname2 + ".majiq --minreads 2 --minpos 2 --nthreads " + str(nthreads) + " --output " + name + "_deltapsi --names Test1 Test2")
 
         
print("Parsing majiq deltapsi outout")          
###Loads the pickle file which contains the ID and number of reads for each junction found in that LSV
pp = pkl.load(open(name + "_deltapsi/clean_reads.Test1.pkl"))

###Parses pickle file to create a list of LSV IDs and their corresponding total read counts
lsvs = {}
read_counts = []
for x in pp:
    lsvs[x[0]] = sum(x[1].tolist())
    read_counts.append(sum(x[1].tolist()))

lsv_reads = sum(read_counts)
lsv_amount = len(read_counts)
print("Total LSV reads: " + str(lsv_reads))
print("Total LSV amount: " + str(lsv_amount))
print("Calculating read count distribution")

###Calculates and stores histogram density values
###Uncomment following section to plot histogram
print("Max read count: " + str(max(read_counts)))
binsize = sys.argiv[5]
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
ratio = total_reads/lsv_reads
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

ll = pkl.load(open(exp_name + "_deltapsi/Test1_Test2.deltapsi.pickle"))

lsvsObj = ll.lsvs

bins = {}
for x in lsvsObj:
        bins[x.get_id()] = x.bins

sig_values = {}
for lsvObj in bins.keys():
        binsList = bins[lsvObj]
        sig = []
        for i in range(len(binsList)):
                sig.append(vlsv.get_expected(vlsv.collapse_matrix(binsList[i])))
        sig_values[lsvObj] = sig


sigReads = []

numSig = 0
for lsv in sig_values.keys():
        if sum([1 for x in sig_values[lsv] if float(x) > .8 or float(x) <.2 ]) > 0:
             sigReads.append(lsvs[lsv])   
	     print lsv
		
print len(sigReads)
print sigReads

counts, bins, bars = plt.hist(sigReads, bins=int(max(sigReads)/int(binsize)))
hist_file_o = open(name + "_delta_hist.txt","w")

for i in range(0,len(counts),1):
        hist_file_o.write(str(bins[i]) + "\t" + str(counts[i])+"\n")

hist_file_o.close()


print numSig
##Create new histograms based on the expected amount of captured distinct LSVs from preseq at <population>
for population in [50000,100000,150000,200000,250000,300000,400000,500000]:
	#Creates a ratio of the ls

	print(str(scaled[population]) + " expected LSVs from " + str(population) + "\n") 
	ratio2 = float(scaled[population])/lsv_amount
	print("Expected LSVs to actual LSVs " + str(ratio2) + "\n")
        new_histFile = open(name + "_at_" + str(population) + "_lsv_reads.txt", "w")
	binBracket = 2
        for i in range(0,len(counts),1):
                new_histFile.write(str(binBracket) + "-" + str(binBracket+3)  + "\t" + str(counts[i] * ratio2) + "\n")
		binBracket += 4
        new_histFile.close()    








    
