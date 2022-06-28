"""
This file extracts data from NMEA files returned by the receiver to generate a beam map. NMEA data does not
generate elevations or S/N accurate enough for meaningful beam maps.

Written by Sabrina Berger
Created in April, 2022
"""
# below are the packages needed to run this script
# you should be able to install most with pip, but
# I recommend making a Conda virtual environment to place them in.
from pynmeagps import NMEAReader # our NMEA parser
import seaborn as sns # a plot beautifier that's note being used yet
import numpy as np # we use numpy for our array creation/manipulation
import glob # this allows us to get all the files in a specific directory with a particular ending
import matplotlib.pyplot as plt

# directories are below
# data_direc_1 = "/Users/sabrinaberger/Library/CloudStorage/OneDrive-McGillUniversity/31622"
# data_direc_2 = "/Users/sabrinaberger/Library/CloudStorage/OneDrive-McGillUniversity/31722"
d3a_first_beam_map = "/Users/sabrinaberger/Desktop/June14_port_C2_GNSS/"

# setting the seaborn theme, ignore for now
sns.set_theme()

# Create a way for the user to input file they'd like to parse below
# log_input = input("What file would you like to parse > ")
# if log_input == '':
# 	log_input = "../Jacksons_house_3/log_0000.nmea"

# this function converts the receiver's C/N_0 measurement into a classical power measurement in dBw
def convert_P(C_N):
	"""
	Convert C/N_0 to power in dBw
	"""
	N_sys = 8.5  # dB
	Tant = 290  # K
	P = C_N + 10*np.log10(Tant + 290*(10**(N_sys/10) - 1)) - 228.6 # last i power
	return P

if __name__ == "__main__":
	### Below we start reading in NMEA data
	nmr = []

	# print(data_direc_1, '/*.nmea')
	for f in glob.glob(d3a_first_beam_map + '/*.nmea')[:10]: # gets files in data_direc_1 that end in .nmea
		stream = open(f, 'rb')
		print(stream)
		nmr.append(NMEAReader(stream, nmeaonly=True)) # appends the fully parsed file to nmr

	# allocating numpy arrays for elevations (altitudes), CNOs, azimuths
	alts = np.empty(int(19e6))
	cnos = np.empty(int(19e6))
	az = np.empty(int(19e6))
	sv_corresponding = np.empty(int(19e6))
	i = 0

	# iterating over each nmea message (n) and extracting the necessary information we need to make
	# beam maps
	for n in nmr:
		print(f"entered {n}")
		for (raw_data, parsed_data) in n:
			for j in range(4): # there are max 4 satellites per NMEA message I believe
				if hasattr(parsed_data, f"elv_0{j}") and hasattr(parsed_data, f"cno_0{j}"):
					# this line checks to see whether the NMEA message actually has an elevation of CN_0 value
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

	# directory where we'll save the dictionary
	sat_dict_direc = "../parsed_data/"

	# saving all the parsed data files as .npy files
	np.save(sat_dict_direc + "az.npy", az)
	np.save(sat_dict_direc + "alts.npy", alts)
	np.save(sat_dict_direc + "cnos.npy", cnos,)
	np.save(sat_dict_direc + "sv_corresponding.npy", sv_corresponding)
	slice = int(19e6)
	mask = cnos > 0  # get numeric version to use ask mask
	cnos = cnos[mask]
	power=convert_P(cnos[:slice])
	# normalized = power - max(power) # taking ratio of powers so subtracting
	az = az[mask]
	alts = alts[mask]
	plt.scatter(az[:slice], alts[:slice], c=power)
	plt.xlabel("az [deg]")
	plt.ylabel("alt [deg]")

	plt.colorbar(label="dBW", cmap="magma_r")
	plt.savefig("first_d3a_beam.png", dpi=300)