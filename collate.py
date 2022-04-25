"""
This file wraps fast_extract_sat_info.py to extract TEC, elevation, and time data from multiple directories.

Written by Sabrina Berger
"""

from fast_extract_sat_info import GetIonoMap
import numpy as np
from datetime import datetime, timedelta

# data directories with my .SBF files (output from Septentrio receiver)
data_direc_0 = "/Users/sabrinaberger/Library/CloudStorage/OneDrive-McGillUniversity/31522/"
data_direc_1 = "/Users/sabrinaberger/Library/CloudStorage/OneDrive-McGillUniversity/31622/"
data_direc_2 = "/Users/sabrinaberger/Library/CloudStorage/OneDrive-McGillUniversity/31722/"

# directory to store .npy files with parsed data
sat_dict_direc = "../parsed_data/"
test_iono = GetIonoMap(data_direc_1, "log_0000.sbf_measurements.txt", "log_0000.sbf_SBF_ChannelStatus.txt", plot=False, dict_name="../parsed_data/satellite_dict_all", include_elev=True)

WNc = 2201 # week number counter for GNSS time
test_iono = GetIonoMap("", plot=False, process=False) # empty instance of GetIonoMap
sat_dict = np.load(sat_dict_direc + "satellite_dict_all.npy", allow_pickle=True).item() # extracting saved dictionary
test_iono.replace_sat_dict(sat_dict) # replacing instance of GetIonoMap with saved dictionary
time_str = '16/3/2022 15:01' # NOT RIGHT, but should be initial time
date_format_str = '%d/%m/%Y %H:%M' # format of string for datetime
given_time = datetime.strptime(time_str, date_format_str)

test_iono.plot_process(WNc=WNc)


# Code below does parsing of data
# for elev in [0]: # min elevation is 0, no filtered satellites
#     for fn in range(1, 19):
#         if fn < 10:
#             fn_lab = "0" + str(fn)
#         else:
#             fn_lab = fn
#         filename_pr, filename_alt = 'log_00{}.sbf_measurements.txt'.format(fn_lab), 'log_00{}.sbf_SBF_SatVisibility1.txt'.format(fn_lab)
#         test_iono.update_files(data_direc_1, filename_pr, filename_alt)
#         test_iono.extract_tec_elev_times_sat_id()
#         print(f"Processed file {fn} from {data_direc_1}")
#
#
#     for fn in range(1, 16):
#         if fn < 10:
#             fn_lab = "0" + str(fn)
#         else:
#             fn_lab = fn
#         filename_pr, filename_alt = 'log_00{}.sbf_measurements.txt'.format(
#             fn_lab), 'log_00{}.sbf_SBF_SatVisibility1.txt'.format(fn_lab)
#         test_iono.update_files(data_direc_2, filename_pr, filename_alt)
#         test_iono.extract_tec_elev_times_sat_id()
#         print(f"Processed file {fn} from {data_direc_2}")
#
#     test_iono.convert_save_to_dict()
