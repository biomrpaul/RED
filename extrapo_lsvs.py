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

def expected_dpsi(bins):
	return sum(np.array(bins) * np.arange(-1+1./len(bins), 1., 2./len(bins)))

exp_name = sys.argv[1]
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
	total_reads = 33088629.0
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
binsize = int(sys.argv[5])
counts, bins, bars = plt.hist(read_counts, bins=[min(read_counts)+i*binsize for i in range(int(np.ceil(max(read_counts)/binsize))) ])
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
    scaled[int(float(vals[0]))] = [float(vals[1]),float(vals[2]),float(vals[3])]
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
                sig.append(expected_dpsi(vlsv.collapse_matrix(binsList[i])))
        sig_values[lsvObj] = sig


sigReads = []

numSig = 0
for lsv in sig_values.keys():
        if sum([1 for x in sig_values[lsv] if float(abs(x)) >.2 ]) > 0:
             sigReads.append(lsvs[lsv])   
		

sigReadBins = {}
sigReadBinOrder = []
binBracket = 3
while binBracket < max(sigReads):
	sigReadBins[str(binBracket) + "-" + str(binBracket+3)] = 0
        sigReadBinOrder.append(  str(binBracket) + "-" + str(binBracket+binsize-1) )
	binBracket += binsize

numSig = 0
for i in range(len(sigReads)):
	binNum = float(sigReads[i] - 2) / float(binsize)
	if sigReadBins >= 10:
		numSig += 1
		sigReadBins[sigReadBinOrder[int(np.ceil(binNum)) - 1]] += 1

proports = []
for i in range(len(sigReadBinOrder)):
	proports.append(sigReadBins[sigReadBinOrder[i]] / counts[i])

print "Number of delta psi LSVs, with read cov. > 10 = " + str(numSig)

expected_deltaFile = open(name + "_expected_delta_lsvs_over_L_prime.txt", "w")
expected_deltaFile2 = open(name + "_expected_delta_lsvs_over_N_prime.txt", "w")
exp_delta_lsvs = []
##Create new histograms based on the expected amount of captured distinct LSVs from preseq at <population>
populations = [25000*i for i in range(106)]
for population in populations:
	#Creates a ratio of the ls
	ratio2 = float(scaled[population][0])/lsv_amount
	ratioLow = float(scaled[population][1])/lsv_amount
	ratioHigh = float(scaled[population][2])/lsv_amount
	total_amount = 0
	total_low = 0
	total_high = 0
	binBracket = 3
        for i in range(len(sigReadBinOrder)):
	
		if np.isnan(proports[i]) == False:
			total_amount += counts[i] * ratio2 * proports[i]
			total_low += counts[i] * ratioLow * proports[i]	
			total_high += counts[i] * ratioHigh * proports[i]


	exp_delta_lsvs.append([total_amount, total_low, total_high])
	expected_deltaFile.write(str(population) + "\t" + str(total_amount) + "\t" + str(total_low) + "\t" + str(total_high) +  "\n")
	expected_deltaFile.write(str(population * ratio) + "\t" + str(total_amount) + "\t" + str(total_low) + "\t" + str(total_high) + "\n")


plt.plot(populations, [amount[0] for amount in exp_delta_lsvs])
plt.plot(populations, [amount[1] for amount in exp_delta_lsvs], "b--")
plt.plot(populations, [amount[2] for amount in exp_delta_lsvs], "b--")
plt.plot(lsv_reads, numSig, "ro")
plt.ylabel("Expected # of Delta Psi LSVs")
plt.xlabel("# of Reads Supporting LSVs")
plt.title("Extrapolated Complexity Plot")
plt.axis([0,max(populations), 0, 1.10*max([amount[2] for amount in exp_delta_lsvs])])
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(name +"_exp_delta_LSVS_over_L_prime.png")
plt.close()

populations2 = [int(x)*ratio for x in populations]
plt.plot(populations2, [amount[0] for amount in exp_delta_lsvs])
plt.plot(populations2, [amount[1] for amount in exp_delta_lsvs], "b--")
plt.plot(populations2, [amount[2] for amount in exp_delta_lsvs], "b--")
plt.plot(total_reads, numSig, "ro")
plt.ylabel("Expected # of Delta Psi LSVs")
plt.xlabel("# of Sequenced Reads")
plt.title("Extrapolated Complexity Plot")
plt.axis([0,max(populations2), 0, 1.10*max([amount[2] for amount in exp_delta_lsvs])])
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(name +"_exp_delta_LSVS_over_N_prime.png")
plt.close()





    
