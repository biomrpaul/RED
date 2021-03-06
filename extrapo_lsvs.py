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
gff3_file = "~/annotations/ensembl.mm10.gff3"
nthreads = 8
preseq_loc = "~/Desktop/preseq-master/"
fastq_loc = sys.argv[4]
run_majiq = sys.argv[6]
##Calculating total amount of reads
count_reads = sys.argv[7]
run_build = sys.argv[8]
step_size = sys.argv[9]
max_extrap = sys.argv[10]
max_terms = sys.arv[11]
max_boot = sys.arv[12]

if count_reads == "count_reads":
	print("Calculating total amount of reads")
	reads = commands.getoutput('wc -l ' + fastq_loc)
	total_reads = int(reads.split(" ")[1]) / 2.0
else:	
	total_reads = 33088629.0 * 2
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
print sum(counts)
#plt.hist(read_counts, bins=4000, range=[0,1000])
#plt.xlabel("LSV Total Read Count")
#plt.ylabel("Frequency")
#plt.savefig(name + "_read_count_hist_zoomed.png")
hist_file_o = open(name + "_hist_values_" + str(binsize) + ".txt","w")
preseq_reads = 0
for i in range(0,len(counts),1):
	hist_file_o.write(str(i+1) + "\t" + str(counts[i])+"\n")
	preseq_reads += (i+1) * counts[i]
hist_file_o.close()

###Extrapolates with preseq   
print("Running preseq")
preseq_command = preseq_loc + "preseq lc_extrap -s " + str(max_step) + " -e " + max_extrap + " -x " + max_terms + " -n " + max_boot + " -o " + name + "_extra_bin" + "_" + str(binsize) + ".txt -H " + name + "_hist_values_" + str(binsize) + ".txt"
os.system(preseq_command)

preseq_command = preseq_loc + "preseq c_curve -o " + name + "_c_curve.txt -H " +  name + "_hist_values_" + str(binsize) + ".txt"
os.system(preseq_command)

print("preseq completed")
###Opens the file that contains the extrapolated values from preseq 
scaled_c = open(name + "_extra_bin_" + str(binsize) + ".txt","r")

###Parses the parseq output 
ratio = total_reads/lsv_reads
print("LSV reads to Total Reads Ratio: " + str(ratio))
preseq_to_lsv_reads_ratio = float(lsv_reads) / float(preseq_reads)
num_reads = []
num_distinct = [] 
low_bound = []
high_bound = []
scaled_c.next()
scaled = {}
for x in scaled_c:
    vals = x.strip("\n").split("\t") 
    num_reads.append(float(vals[0])*ratio*preseq_to_lsv_reads_ratio)
    scaled[int(float(vals[0]))] = [float(vals[1]),float(vals[2]),float(vals[3])]
    num_distinct.append(vals[1])
    low_bound.append(vals[2])
    high_bound.append(float(vals[3]))


##Parses c_curve values
curve = open(name + "_c_curve.txt", "r")
curve.next()
num_reads2 = []
num_distinct2 = []
scaled2 = {}
for x in curve:
	line = x.strip("\n").split("\t")
	num_reads2.append(float(line[0])*ratio*preseq_to_lsv_reads_ratio)
	num_distinct2.append(line[1])
	scaled2[float(line[0])] = float(line[1])

print "C_curve amount of reads : " + str(preseq_reads)

print("Plotting extrapolated values")
###Plots the extrapolated values from preseq (only the first 150M reads)
plt.plot(num_reads, num_distinct, label="Preseq	Extrapo.")
plt.plot(num_reads, low_bound, "b--", label="Upper CI")
plt.plot(num_reads, high_bound, "b--", label="Lower CI")
plt.plot(num_reads2, num_distinct2, "go", label="Preseq down samp.")
plt.plot(total_reads, lsv_amount, "ro", label="True Amount")
plt.ylabel("# of Distinct LSVs")
plt.xlabel("# of Reads")
plt.title("Extrapolated Complexity Plot")
plt.axis([0,250000000, 0, 25000])
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
binBracket = binsize - 1
while binBracket < max(sigReads) + binsize - 1:
	sigReadBins[str(binBracket) + "-" + str(binBracket+binsize-1)] = 0
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

expected_deltaFile = open(name + "_expected_delta_lsvs_over_L_prime" + str(binsize) + ".txt", "w")
expected_deltaFile2 = open(name + "_expected_delta_lsvs_over_N_prime" + str(binsize) + ".txt", "w")
exp_delta_lsvs = []
##Create new histograms based on the expected amount of captured distinct LSVs from preseq at <population>


populations = [25000*i for i in range(80)]
for population in populations:
	#Creates a ratio of the ls
	ratio2 = float(scaled[population][0])/lsv_amount
	ratioLow = float(scaled[population][1])/lsv_amount
	ratioHigh = float(scaled[population][2])/lsv_amount
	total_amount = 0
	total_low = 0
	total_high = 0
        for i in range(len(sigReadBinOrder)):
	
		if np.isnan(proports[i]) == False:
			total_amount += counts[i] * ratio2 * proports[i]
			total_low += counts[i] * ratioLow * proports[i]	
			total_high += counts[i] * ratioHigh * proports[i]


	exp_delta_lsvs.append([total_amount, total_low, total_high])
	expected_deltaFile.write(str(population*preseq_to_lsv_reads_ratio) + "\t" + str(total_amount) + "\t" + str(total_low) + "\t" + str(total_high) +  "\n")
	expected_deltaFile2.write(str(population * ratio*preseq_to_lsv_reads_ratio) + "\t" + str(total_amount) + "\t" + str(total_low) + "\t" + str(total_high) + "\n")


populations3 = [int(x)*preseq_to_lsv_reads_ratio for x in populations]
plt.plot(populations3, [amount[0] for amount in exp_delta_lsvs], label="Preseq extrapo.")
plt.plot(populations3, [amount[1] for amount in exp_delta_lsvs], "b--", label="Lower CI")
plt.plot(populations3, [amount[2] for amount in exp_delta_lsvs], "b--", label="Upper CI")
plt.plot(lsv_reads, numSig, "ro", label="True Data")
plt.ylabel("Expected # of Delta Psi LSVs")
plt.xlabel("# of Reads Supporting LSVs")
plt.title("Extrapolated Complexity Plot")
plt.axis([0,250000000/ ratio, 0, 1200])
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(name +"_exp_delta_LSVS_over_L_prime.png")
plt.close()

populations2 = [int(x)*ratio*preseq_to_lsv_reads_ratio for x in populations]
plt.plot(populations2, [amount[0] for amount in exp_delta_lsvs])
plt.plot(populations2, [amount[1] for amount in exp_delta_lsvs], "b--")
plt.plot(populations2, [amount[2] for amount in exp_delta_lsvs], "b--")
plt.plot(total_reads, numSig, "ro")
plt.ylabel("Expected # of Delta Psi LSVs")
plt.xlabel("# of Sequenced Reads")
plt.title("Extrapolated Complexity Plot")
plt.axis([0,250000000, 0, 1200])
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(name +"_exp_delta_LSVS_over_N_prime.png")
plt.close()

stats = open(name + "_stats" + str(binsize) + ".txt", "w")
stats.write("Total_reads\t" + str(total_reads) + "\n")
stats.write("lsv_reads\t" + str(lsv_reads) + "\n")
stats.write("preseq_bin_reads\t" + str(preseq_reads) + "\n")
stats.write("Distinct_lsvs\t" + str(lsv_amount) + "\n")
stats.write("Delta_psi_lsvs\t" + str(numSig) + "\n")
stats.write("Total_reads:lsv_reads\t" + str(float(total_reads)/float(lsv_reads)) + "\n")
stats.write("lsv_reads:bin_reads\t" + str(preseq_to_lsv_reads_ratio) + "\n")



    
