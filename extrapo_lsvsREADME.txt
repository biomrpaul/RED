Author: Matt R. Paul
Edited: 3 April 2015


extrapo_lsvs.py


0. Configure config file, <exp_name>.config. Set the two experimental conditions to:
Test1 = <consider_sample>
Test2 = <comparison_sample>
1. Before running, set variables in the beginning of script. (May change some/all to input arguments). These include:
   1. experiment name (exp_name)
   2. subname1 and subname2, which are equal to the considered sample and the compared sample names
   3. location of the majiq executabel (majiq_loc)
   4. location of relevant gff3 file (gff3_file)
   5. location of preseq executable (preseq_loc)
   6. total number of reads in considered sample (this could be calculated from bam or fastq file, but is now set to a fixed number)
1. Run the script, python extrapo_lscvs.py
The steps executed in the script are as follows:
     0. *Calculate total amount of reads? Right now given as input parameter.
1. Run majiq build
2. Run majiq deltapsi
3. Parse clean_reads file (deltapsi output)
4. Calculate histrogram density values
5. Scale based on ratio between total reads and total reads supporting lsvs
6. **Extrapolate first 9 read count frequencies from fitted exponential distribution (not included in clean_reads file)
7. Run preseq on new scaled values + extrapolated first 9 vaues
8. Plot lines of preseq output as is and scaled back to number of expected distinct LSVs


Notes and concerns:
*How should we quantify total amount of reads? Aligned reads in bam file? Total amount of reads in 1 fastq file, in both? Estimated by what to expect.?


**clean_reads.pkl does not contain read counts for lsvs with less than 10 counts, can we get these values?
n_x = number of times there were counts of x
n_2 cannot be greater than n_1 (which was the case with n_10 and n_11 in Lung1) becaus extrapolation saturates before the number of total reads doubles.
Had to extrapolate backwards to get n_1:9 using exponential fit.
