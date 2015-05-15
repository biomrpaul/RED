import matplotlib.pyplot as plt

from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import numpy as np

uno = open("Thymus1_read_counts.hist", "r")
counts = []
amounts = []
for row in uno:
    line = row.strip("\n").split("\t")
    if line[1] > 0:
        counts.append(float(line[0]))
        amounts.append(log(float(line[1])))
    
plt.bar(counts[:200], amounts[:200])
plt.ylabel("Log # of Distinct Reads")
plt.xlabel("Read Count")
plt.title("Thymus1 Count Histogram")
plt.show()
plt.savefig("Thymus1_count_hist_n200_log.png", format="png")
plt.close()

plt.bar(counts[:25], amounts[:25])
plt.ylabel("Log # of Distinct Reads")
plt.xlabel("Read Count")
plt.title("Thymus1 Count Histogram")
plt.show()
plt.savefig("Thymus1_count_hist_n25_log.png", format="png")
plt.close()

plt.bar(counts[1:25], amounts[1:25])
plt.ylabel("Log # of Distinct Reads")
plt.xlabel("Read Count")
plt.title("Thymus1 Count Histogram")
plt.show()
plt.savefig("Thymus1_count_hist_n25_log_sans1.png", format="png")
plt.close()

uno = open("Thymus1_read_counts.hist", "r")
counts = []
amounts = []
for row in uno:
    line = row.strip("\n").split("\t")
    if line[1] > 0:
        counts.append(float(line[0]))
        amounts.append(float(line[1]))
    
plt.bar(counts[:200], amounts[:200])
plt.ylabel("# of Distinct Reads")
plt.xlabel("Read Count")
plt.title("Thymus1 Count Histogram")
plt.show()
plt.savefig("Thymus1_count_hist_n200.png", format="png")
plt.close()

plt.bar(counts[:25], amounts[:25])
plt.ylabel("# of Distinct Reads")
plt.xlabel("Read Count")
plt.title("Thymus1 Count Histogram")
plt.show()
plt.savefig("Thymus1_count_hist_n25.png", format="png")
plt.close()

plt.bar(counts[1:25], amounts[1:25])
plt.ylabel("# of Distinct Reads")
plt.xlabel("Read Count")
plt.title("Thymus1 Count Histogram")
plt.show()
plt.savefig("Thymus1_count_hist_n25_sans1.png", format="png")
plt.close()

dos = open("ThymusVSpleen/ThymusVSpleen_hist_values_4.txt", "r")
dos_counts = []
dos_amounts = []
for row in dos:
    line = row.strip("\n").split("\t")
    if line[1] > 0:
        dos_counts.append(float(line[0]))
        dos_amounts.append(log(float(line[1])))
    
plt.bar(dos_counts[:200], dos_amounts[:200])
plt.ylabel("Log # of Distinct LSVs")
plt.xlabel("LSV Read Count")
plt.title("Thymus1 LSV Count Histogram")
plt.show()
plt.savefig("Thymus1_LSVs_count_hist_n200_log.png", format="png")
plt.close()

plt.bar(dos_counts[:25], dos_amounts[:25])
plt.ylabel("Log # of Distinct LSVs")
plt.xlabel("LSV Read Count")
plt.title("Thymus1 LSV Count Histogram")
plt.show()
plt.savefig("Thymus1_LSVs_count_hist_n25_log.png", format="png")
plt.close()

dos = open("ThymusVSpleen/ThymusVSpleen_hist_values_4.txt", "r")
dos_counts = []
dos_amounts = []
for row in dos:
    line = row.strip("\n").split("\t")
    if line[1] > 0:
        dos_counts.append(float(line[0]))
        dos_amounts.append(float(line[1]))
    
plt.bar(dos_counts[:200], dos_amounts[:200])
plt.ylabel("# of Distinct LSVs")
plt.xlabel("LSV Read Count")
plt.title("Thymus1 LSV Count Histogram")
plt.show()
plt.savefig("Thymus1_LSVs_count_hist_n200.png", format="png")
plt.close()

plt.bar(dos_counts[:25], dos_amounts[:25])
plt.ylabel("Log # of Distinct LSVs")
plt.xlabel("LSV Read Count")
plt.title("Thymus1 LSV Count Histogram")
plt.show()
plt.savefig("Thymus1_LSVs_count_hist_n25.png", format="png")
plt.close()
