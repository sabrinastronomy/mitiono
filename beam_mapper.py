"""
This file extracts data from NMEA files returned by the receiver to generate a beam map.

Written by Sabrina Berger
"""

from pynmeagps import NMEAReader
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import glob

data_direc_1 = "/Users/sabrinaberger/Library/CloudStorage/OneDrive-McGillUniversity/31622/"
data_direc_2 = "/Users/sabrinaberger/Library/CloudStorage/OneDrive-McGillUniversity/31722/"

sns.set_theme()

# Create a way for the user to input file they'd like to parse below
# log_input = input("What file would you like to parse > ")
# if log_input == '':
# 	log_input = "../Jacksons_house_3/log_0000.nmea"

def convert_P(C_N):
	"""
	Convert C/N_0 to power in dBw
	"""
	N_sys = 8.5  # dB
	Tant = 290  # K
	P = C_N + 10*np.log10(Tant + 290*(10**(N_sys/10) - 1)) - 228.6 # last i power
	return P

# Read NMEA data
nmr = []
for f in glob.glob(data_direc_1 + '/*.nmea')[:10]:
	stream = open(f, 'rb')
	nmr.append(NMEAReader(stream, nmeaonly=True))

# allocating numpy arrays for elevations (altitudes), CNOs, azimuths
alts = np.full(int(19e6), -1.)
cnos = np.full(int(19e6), -1.)
az = np.full(int(19e6), -1.)
sv_corresponding = np.full(int(19e6), -1.)
i = 0

for n in nmr:
	print("entered")
	for (raw_data, parsed_data) in n:
		for j in range(4): # there are max 4 satellites per NMEA message I believe
			if hasattr(parsed_data, f"elv_0{j}") and hasattr(parsed_data, f"cno_0{j}"):
				elv_01 = getattr(parsed_data, f"elv_0{j}")
				cno_01 = getattr(parsed_data, f"cno_0{j}")
				name_01 = getattr(parsed_data, f"svid_0{j}")
				az_01 = getattr(parsed_data, f"az_0{j}")
				if elv_01 != '' and cno_01 != '':
					alts[i] = elv_01
					cnos[i] = cno_01
					az[i] = az_01
					sv_corresponding[i] = name_01
					i += 1

sat_dict_direc = "../parsed_data/"

np.save(sat_dict_direc + "az.npy", az)
np.save(sat_dict_direc + "alts.npy", alts)
np.save(sat_dict_direc + "cnos.npy", cnos,)
np.save(sat_dict_direc + "sv_corresponding.npy", sv_corresponding)


