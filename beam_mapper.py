from pynmeagps import NMEAReader
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

sns.set_theme()

# log_input = input("What file would you like to parse > ")
# if log_input == '':
# 	log_input = "../Jacksons_house_3/log_0000.nmea"
#
log_input = "../Jacksons_house_3/log_0000.nmea"
stream = open(log_input, 'rb')
nmr = NMEAReader(stream, nmeaonly=True)
alts = []
cnos = []
sv_corresponding = []
dict_sats = {}
for (raw_data, parsed_data) in nmr:
	if hasattr(parsed_data, "elv_01") and hasattr(parsed_data, "cno_01"):
		elv_01 = parsed_data.elv_01
		cno_01 = parsed_data.cno_01
		name_01 = parsed_data.svid_01
		if name_01 in dict_sats.keys():
			dict_sats[name_01]["elevation"].append(elv_01)
			dict_sats[name_01]["cno"].append(cno_01)
		else:
			dict_sats[name_01] ={}
			dict_sats[name_01]["elevation"] = [elv_01]
			dict_sats[name_01]["cno"] = [cno_01]
		if elv_01 != '' and cno_01 != '':
			alts.append(elv_01)
			cnos.append(cno_01)
	if hasattr(parsed_data, "elv_02") and hasattr(parsed_data, "cno_02"):
		elv_02 = parsed_data.elv_02
		cno_02 = parsed_data.cno_02
		name_02 = parsed_data.svid_02
		if name_02 in dict_sats.keys():
			dict_sats[name_02]["elevation"].append(elv_02)
			dict_sats[name_02]["cno"].append(cno_02)
		else:
			dict_sats[name_02] ={}
			dict_sats[name_02]["elevation"] = [elv_02]
			dict_sats[name_02]["cno"] = [cno_02]
		if elv_02 != '' and cno_02 != '':
			alts.append(elv_02)
			cnos.append(cno_02)
			sv_corresponding.append(name_02)
	if hasattr(parsed_data, "elv_03") and hasattr(parsed_data, "cno_03"):
		elv_03 = parsed_data.elv_03
		cno_03 = parsed_data.cno_03
		name_03 = parsed_data.svid_03
		if name_03 in dict_sats.keys():
			dict_sats[name_03]["elevation"].append(elv_03)
			dict_sats[name_03]["cno"].append(cno_03)
		else:
			dict_sats[name_03] ={}
			dict_sats[name_03]["elevation"] = [elv_03]
			dict_sats[name_03]["cno"] = [cno_03]
		if elv_03 != '' and cno_03 != '':
			alts.append(elv_03)
			cnos.append(cno_03)
			sv_corresponding.append(name_03)
	if hasattr(parsed_data, "elv_04") and hasattr(parsed_data, "cno_04"):
		elv_04 = parsed_data.elv_04
		cno_04 = parsed_data.cno_04
		name_04 = parsed_data.svid_04

		if name_04 in dict_sats.keys():
			dict_sats[name_04]["elevation"].append(elv_04)
			dict_sats[name_04]["cno"].append(cno_04)
		else:
			dict_sats[name_04] ={}
			dict_sats[name_04]["elevation"] = [elv_04]
			dict_sats[name_04]["cno"] = [cno_04]
		if elv_04 != '' and cno_04 != '':
			alts.append(elv_04)
			cnos.append(cno_04)
			sv_corresponding.append(name_04)

def convert_P(C_N):
	N_sys = 8.5  # dB
	Tant = 290  # K
	P = C_N + 10*np.log10(Tant + 290*(10**(N_sys/10) - 1)) - 228.6 # last i power
	return P
powers = convert_P(cnos)
direc = "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/plots_talk/"
plt.title("Receiver on Rooftop")
plt.scatter(alts, powers, c='k')
plt.xlabel("elevation [deg]")
plt.ylabel("P [dBW]")
plt.savefig(direc + "newest_beam_map_try.png", dpi=200)