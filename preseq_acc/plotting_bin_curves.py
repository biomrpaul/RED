import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter


twofive_2 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime2.txt"))
twofive_3 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime3.txt"))
twofive_4 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime4.txt"))
twofive_5 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime5.txt"))
twofive_6 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime6.txt"))
twofive_7 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime7.txt"))
twofive_8 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime8.txt"))

five_2 = np.transpose(np.loadtxt("ThymusVSpleen_.5_expected_delta_lsvs_over_N_prime2.txt"))
five_8 = np.transpose(np.loadtxt("ThymusVSpleen_.5_expected_delta_lsvs_over_N_prime8.txt"))

one_4 = np.transpose(np.loadtxt("ThymusVSpleen_expected_delta_lsvs_over_N_prime4.txt"))
one_6 = np.transpose(np.loadtxt("ThymusVSpleen_expected_delta_lsvs_over_N_prime6.txt"))
one_7 = np.transpose(np.loadtxt("ThymusVSpleen_expected_delta_lsvs_over_N_prime7.txt"))

plt.plot(twofive_2[0], twofive_2[1], "g-", label="0.25 bin_size=2")
plt.plot(twofive_3[0], twofive_3[1], "g--", label="0.25 bin_size=3")
plt.plot(twofive_4[0], twofive_4[1], "g-.",label="0.25 bin_size=4")
plt.plot(twofive_5[0], twofive_5[1], "g:",label="0.25 bin_size=5")
plt.plot(twofive_6[0], twofive_6[1], "g^",label="0.25 bin_size=6")
plt.plot(twofive_7[0], twofive_7[1], "g*",label="0.25 bin_size=7")
plt.plot(twofive_8[0], twofive_8[1], "gs",label="0.25 bin_size=8")
plt.plot(16544314, 206, "go", label="0.25 True Value")
plt.plot(five_2[0], five_2[1], "b-",label="0.50 bin_size=2")
plt.plot(five_8[0], five_8[1], "b--",label="0.50 bin_size=8")
plt.plot(33088628, 417, "bo", label="0.50 True Value")
plt.plot(one_4[0], one_4[1], "r-",label="1.00 bin_size=4")
plt.plot(one_6[0], one_6[1], "r--",label="1.00 bin_size=6")
plt.plot(one_7[0], one_7[1], "r-.",label="1.00 bin_size=7")
plt.plot(66177258, 738, "ro", label="1.00 True Value")
plt.ylabel("Expected # of Delta Psi LSVs")
plt.xlabel("# of Sequenced Reads")
plt.title("Bins versus Proportion of Sample")
plt.axis([0,250000000,0,900])
plt.legend(prop={'size':8})
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig("Bins_vs_size.png")
plt.close()
