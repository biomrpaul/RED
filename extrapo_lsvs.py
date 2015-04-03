import matplotlib
import matplotlib.pyplot as plt
import pickle
import sys
import os
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

#name = sys.argv[1]
#pk_loc = sys.argv[2]
#total_reads = sys.argv[3]
#input arguments or static variables, for the name of the experiment, the names of the two samples being compared, the location
#of the majiq executable and gff3 file, number of threads, and the total number of reads of the sample being considered

exp_name= "test2"
subname1 = "Lung2.mm10.sorted"
subname2 = "Lung1.mm10.sorted"
majiq_loc = "/opt/majiq"
gff3_file = "/data/DB/mm10/ensembl.mm10.gff3"
nthreads = 8
preseq_loc = "/data/mrpaul/preseq-1.0.2.Linux_x86_64/"

#This could be calculating the number of aligned reads in corresponding bam file, 
#counting reads in fastq file, or some other metric
total_reads = 20000000


print("Running majiq build")
name=exp_name
###Runs majiq build on the samples being considered
os.system("export PYTHONPATH=" + majiq_loc)
os.system(majiq_loc + "/majiq build " + gff3_file + " -conf " + name + ".config --nthreads " + str(nthreads) + " --output " + name + "_build")

print("Runnig majiq deltapsi")
###Runs majiq deltapsi
os.system(majiq_loc + "/majiq deltapsi -grp1 " + name + "_build/" + subname1 + ".majiq -grp2 " + name + "_build/" \
          + subname2 + ".majiq --nthreads " + str(nthreads) + " --output " + name + "_deltapsi --names Test1 Test2")
          
print("Parsing majiq deltapsi outout")          
###Loads the pickle file which contains the ID and number of reads for each junction found in that LSV
pp = pickle.load(open(name + "_deltapsi/clean_reads.Test1.pkl"))

###Parses pickle file to create a list of LSV IDs and their corresponding total read counts
lsvs = []
read_counts = []
for x in pp:
    lsvs.append(x[0])
    read_counts.append(sum(x[1].tolist()))

print("Calculating read count distribution and scaling by ratio to total amount of reads")

###Calculates and stores histogram density values
###Uncomment following section to plot histogram
counts, bins, bars = plt.hist(read_counts, bins=max(read_counts))
#plt.hist(read_counts, bins=4000, range=[0,1000])
#plt.xlabel("LSV Total Read Count")
#plt.ylabel("Frequency")
#plt.savefig(name + "_read_count_hist_zoomed.png")

###Calculates ratio of the total reads of the sample to the total amount of reads supporting LSVs
###Proxy read counts are created by scaling the amount of reads counts supporting LSVs by the ratio
ratio = total_reads/sum(read_counts)
scaled_counts = [int(x*ratio) for x in counts]

###Because read counts of LSVs that have less than 10 supporting reads are not included, 
###an exponential is fit to [n_10,...n_30] to predict the values of [n_1,....,n_9]
fitted = np.polyfit(range(10,30,1), log(scaled_counts[:20]), 1)

scaled_counts_total = []
for i in range(1,10,1):
    scaled_counts_total.append(int(exp(float(fitted[1]))*exp(float(fitted[0])*i)))
    
scaled_counts_total += scaled_counts

###Outputs new histogram values from the scaled and extrapolated counts
histFile = open(name + "_hist_values.txt", "w")

for x in range(0, len(scaled_counts_total), 1):
    histFile.write(str(x+1) + "\t" + str(scaled_counts_total[x]) + "\n")

histFile.close()    

###Extrapolates with preseq   
print("Running preseq")
preseq_command = preseq_loc + "preseq lc_extrap -o " + name + "_extra.txt -H " + name + "_hist_values.txt"
os.system(preseq_command)

print("preseq completed")
###Opens the file that contains the extrapolated values from preseq 
scaled_c = open(name + "_extra.txt","r")

###Parses the parseq output 
num_reads = []
num_distinct = [] 
low_bound = []
high_bound = []
scaled_c.next()
for x in scaled_c:
    vals = x.strip("\n").split("\t")
    num_reads.append(vals[0])
    num_distinct.append(vals[1])
    low_bound.append(vals[2])
    high_bound.append(float(vals[3]))

pltdir=name+"_preseq_plots"
os.mkdir(pltdir)
print("Plotting extrapolated values")
###Plots the extrapolated line and confidence interval values as is
plt.plot(num_reads, num_distinct)
plt.plot(num_reads, low_bound, "b--")
plt.plot(num_reads, high_bound, "b--")
plt.xlabel("# of Distinct LSVs (Unscaled)")
plt.ylabel("# of Reads")
plt.title("Extrapolated Complexity Plot")
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(pltdir + "/"+name + "_unscaled.png")
plt.close()

###Plots the same values but zooms in so that it only looks at the effect of including up to 1B reads
plt.plot(num_reads, num_distinct)
plt.plot(num_reads, low_bound, "b--")
plt.plot(num_reads, high_bound, "b--")
plt.xlabel("# of Distinct LSVs (Unscaled)")
plt.ylabel("# of Reads")
plt.title("Extrapolated Complexity Plot")
plt.axis([0,1000000000, 0, 1.05*max(high_bound)])
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(pltdir + "/"+name + "_unscaled_zoomed_1B.png")
plt.close()

###Plots the same values but zooms in so that it only looks at the effect of including up to 150M reads
plt.plot(num_reads, num_distinct)
plt.plot(num_reads, low_bound, "b--")
plt.plot(num_reads, high_bound, "b--")
plt.xlabel("# of Distinct LSVs (Unscaled)")
plt.ylabel("# of Reads")
plt.title("Extrapolated Complexity Plot")
plt.axis([0,150000000, 0, 1.05*max(high_bound)])
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(pltdir + "/"+name + "_unscaled_zoomed_150M.png")
plt.close()

###Plots the extrapolated line and confidence interval values after being scaled back 
# to the actual amount of distinct LSVs, we should expect
num_distinct_scaled = [float(x)/ratio for x in num_distinct]
low_bound_scaled = [float(x)/ratio for x in low_bound]
high_bound_scaled = [float(x)/ratio for x in high_bound]
plt.plot(num_reads, num_distinct_scaled)
plt.plot(num_reads, low_bound_scaled, "b--")
plt.plot(num_reads, high_bound_scaled, "b--")
plt.xlabel("# of Distinct LSVs (Scaled)")
plt.ylabel("# of Reads")
plt.title("Extrapolated Complexity Plot")
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(pltdir + "/"+name + "_scaled.png")
plt.close()

###Plots the same values but zooms in so that it only looks at the effect of including up to 1B reads
num_distinct_scaled = [float(x)/ratio for x in num_distinct]
low_bound_scaled = [float(x)/ratio for x in low_bound]
high_bound_scaled = [float(x)/ratio for x in high_bound]
plt.plot(num_reads, num_distinct_scaled)
plt.plot(num_reads, low_bound_scaled, "b--")
plt.plot(num_reads, high_bound_scaled, "b--")
plt.xlabel("# of Distinct LSVs (Scaled)")
plt.ylabel("# of Reads")
plt.title("Extrapolated Complexity Plot")
plt.axis([0,1000000000, 0, 1.05*max(high_bound_scaled)])
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(pltdir + "/"+name + "_scaled_zoomed_1B.png")
plt.close()

###Plots the same values but zooms in so that it only looks at the effect of including up to 150M reads
num_distinct_scaled = [float(x)/ratio for x in num_distinct]
low_bound_scaled = [float(x)/ratio for x in low_bound]
high_bound_scaled = [float(x)/ratio for x in high_bound]
plt.plot(num_reads, num_distinct_scaled)
plt.plot(num_reads, low_bound_scaled, "b--")
plt.plot(num_reads, high_bound_scaled, "b--")
plt.xlabel("# of Distinct LSVs (Scaled)")
plt.ylabel("# of Reads")
plt.title("Extrapolated Complexity Plot")
plt.axis([0,150000000, 0, 1.05*max(high_bound_scaled)])
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig(pltdir + "/"+name + "_scaled_zoomed_150M.png")
plt.close()


    
    
