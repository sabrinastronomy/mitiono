# mitiono
miti-iono (mitigate ionosphere): Code to extract Septentrio receiver data to create ionosphere and beam maps. Future work will include code to cross-correlate and decode GNSS baseband data. 

Files and their corresponding descriptions:
## beam_map_nmea.py
This file extracts data from NMEA files returned by the receiver to generate a beam map. NMEA data does not
generate elevations or S/N accurate enough for meaningful beam maps. However, this script is preserved for future
reference or when needing to parse NMEA data.
## beam_map_sbf.py
This file extracts data from SBF SatVisibility files returned by the receiver to generate a beam map.
## CSTUtils.py (written mostly by Vincent Mackay - vincent.mackay@mail.utoronto.ca)
Class containing beams as simulated by CST.
## extract_sbf_data.py
This file contains the ExtractSBF class to extract data from SBF files created by Septentrio receivers.
## iono_map_sbf.py
This script contains a class which takes in a data directory and parses the file for TEC, elevations, times, satellite IDs. It also has plotting capabilities.
## pull_unavco.py (written mostly Jonathan Sievers)
This script pulls GNSS data form UNAVCO specifically for satellites at DRAO.
## use_TLEs.py
This script contains helper functions for using satellite TLE data.
