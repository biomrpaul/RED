import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter


#twofive_2 = np.transpose(np.loadtxt("ThymusVSpleen_.25_extra_bin_2.txt"))
#twofive_3 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime3.txt"))
#twofive_4 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime4.txt"))
#twofive_5 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime5.txt"))
#twofive_6 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime6.txt"))
#twofive_7 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime7.txt"))
#twofive_8 = np.transpose(np.loadtxt("ThymusVSpleen_.25_expected_delta_lsvs_over_N_prime8.txt"))
#
#five_2 = np.transpose(np.loadtxt("ThymusVSpleen_.5_expected_delta_lsvs_over_N_prime2.txt"))
#five_8 = np.transpose(np.loadtxt("ThymusVSpleen_.5_expected_delta_lsvs_over_N_prime8.txt"))
#
#one_4 = np.transpose(np.loadtxt("ThymusVSpleen_expected_delta_lsvs_over_N_prime4.txt"))
#one_6 = np.transpose(np.loadtxt("ThymusVSpleen_expected_delta_lsvs_over_N_prime6.txt"))
#one_7 = np.transpose(np.loadtxt("ThymusVSpleen_expected_delta_lsvs_over_N_prime7.txt"))


twofive_2 = np.transpose(np.loadtxt("ThymusVSpleen_.25_extra_bin_2.txt", skiprows=1))
twofive_3 = np.transpose(np.loadtxt("ThymusVSpleen_.25_extra_bin_3.txt", skiprows=1))
twofive_4 = np.transpose(np.loadtxt("ThymusVSpleen_.25_extra_bin_4.txt", skiprows=1))
twofive_5 = np.transpose(np.loadtxt("ThymusVSpleen_.25_extra_bin_5.txt", skiprows=1))
twofive_6 = np.transpose(np.loadtxt("ThymusVSpleen_.25_extra_bin_6.txt", skiprows=1))
twofive_7 = np.transpose(np.loadtxt("ThymusVSpleen_.25_extra_bin_7.txt", skiprows=1))
twofive_8 = np.transpose(np.loadtxt("ThymusVSpleen_.25_extra_bin_8.txt", skiprows=1))

five_2 = np.transpose(np.loadtxt("ThymusVSpleen_.5_extra_bin_2.txt", skiprows=1))
five_8 = np.transpose(np.loadtxt("ThymusVSpleen_.5_extra_bin_8.txt", skiprows=1))

one_4 = np.transpose(np.loadtxt("ThymusVSpleen_extra_bin_4.txt", skiprows=1))
one_6 = np.transpose(np.loadtxt("ThymusVSpleen_extra_bin_6.txt", skiprows=1))
one_7 = np.transpose(np.loadtxt("ThymusVSpleen_extra_bin_7.txt", skiprows=1))



plt.plot(twofive_2[0]*2.17050132802, twofive_2[1], "g-", label="0.25 bin_size=2")
plt.plot(twofive_3[0]*3.15983953213, twofive_3[1], "g--", label="0.25 bin_size=3")
plt.plot(twofive_4[0]*4.08916340891, twofive_4[1], "g-.",label="0.25 bin_size=4")
plt.plot(twofive_5[0]*4.9671769935, twofive_5[1], "g:",label="0.25 bin_size=5")
plt.plot(twofive_6[0]*5.78774733301, twofive_6[1], "g^",label="0.25 bin_size=6")
plt.plot(twofive_7[0]*6.56215809285, twofive_7[1], "g*",label="0.25 bin_size=7")
plt.plot(twofive_8[0]*7.28336675579, twofive_8[1], "gs",label="0.25 bin_size=8")
plt.plot(130751, 6903, "go", label="0.25 True Value")
plt.plot(five_2[0] * 2.10687361696, five_2[1], "b-",label="0.50 bin_size=2")
plt.plot(five_8[0] * 7.65124514767, five_8[1], "b--",label="0.50 bin_size=8")
plt.plot(337045, 11414, "bo", label="0.50 True Value")
plt.plot(one_4[0] * 4.03831811334, one_4[1], "r-",label="1.00 bin_size=4")
plt.plot(one_6[0] * 5.93283057019, one_6[1], "r--",label="1.00 bin_size=6")
plt.plot(one_7[0] * 6.85163300726, one_7[1], "r-.",label="1.00 bin_size=7")
plt.plot(868513, 17197, "ro", label="1.00 True Value")
plt.ylabel("Expected # of LSVs")
plt.xlabel("# of LSV Reads")
plt.title("Bins versus Proportion of Sample")
plt.axis([0,1600000,0,22000])
plt.legend(prop={'size':8})
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig("Bins_vs_size_lsv_reads.png")
plt.close()


twofive_ratio = 126.532982539
five_ratio = 98.1727306443
one_ratio = 76.1960477276


plt.plot(twofive_2[0]*2.17050132802 * twofive_ratio, twofive_2[1], "g-", label="0.25 bin_size=2")
plt.plot(twofive_3[0]*3.15983953213* twofive_ratio, twofive_3[1], "g--", label="0.25 bin_size=3")
plt.plot(twofive_4[0]*4.08916340891* twofive_ratio, twofive_4[1], "g-.",label="0.25 bin_size=4")
plt.plot(twofive_5[0]*4.9671769935* twofive_ratio, twofive_5[1], "g:",label="0.25 bin_size=5")
plt.plot(twofive_6[0]*5.78774733301* twofive_ratio, twofive_6[1], "g^",label="0.25 bin_size=6")
plt.plot(twofive_7[0]*6.56215809285* twofive_ratio, twofive_7[1], "g*",label="0.25 bin_size=7")
plt.plot(twofive_8[0]*7.28336675579* twofive_ratio, twofive_8[1], "gs",label="0.25 bin_size=8")
plt.plot(16544314.0, 6903, "go", label="0.25 True Value")
plt.plot(five_2[0] * 2.10687361696*five_ratio, five_2[1], "b-",label="0.50 bin_size=2")
plt.plot(five_8[0] * 7.65124514767*five_ratio, five_8[1], "b--",label="0.50 bin_size=8")
plt.plot(33088628.0, 11414, "bo", label="0.50 True Value")
plt.plot(one_4[0] * 4.03831811334*one_ratio, one_4[1], "r-",label="1.00 bin_size=4")
plt.plot(one_6[0] * 5.93283057019*one_ratio, one_6[1], "r--",label="1.00 bin_size=6")
plt.plot(one_7[0] * 6.85163300726*one_ratio, one_7[1], "r-.",label="1.00 bin_size=7")
plt.plot(66177258.0, 17197, "ro", label="1.00 True Value")
plt.ylabel("Expected # of LSVs")
plt.xlabel("# of Sequenced Reads")
plt.title("Bins versus Proportion of Sample")
plt.axis([0,250000000,0,22000])
plt.legend(prop={'size':8})
ax=plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
plt.savefig("Bins_vs_size_total_reads.png")
plt.close()