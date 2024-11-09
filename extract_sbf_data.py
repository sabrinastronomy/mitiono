"""
This file contains the ExtractSBF class to extract data from SBF files created by Septentrio receivers.

Written by Sabrina Berger
"""

import numpy as np
import datetime
import os

### Code below does parsing of data for beam map


class ExtractSBF:
    """
    This class takes in a data directories and parses the file for all required information necessary for ionospheric
    mapping (GetIonoMap class) and beam mapping (GetBeamMap class).
    """
    def __init__(self, data_direcs, masking=True, min_elev=0, process=True, include_elev=True, mask_frequency="1", save_parsed_data_direc="parsed_data"):
        """
        :param data_direc - list of files, directory with receiver raw files
        """
        assert(type(data_direcs) == list)
        self.masking = masking
        self.data_direcs = data_direcs
        self.min_elev = min_elev
        self.include_elev = include_elev
        self.mask_frequency = mask_frequency
        self.save_parsed_data_direc = save_parsed_data_direc
        ## Set params
        ##########################################
        self.reset_params()
        ##########################################

        if process:
            self.parse_data_direcs() # actually does the extraction

    def reset_params(self):
        self.all_sat_dict = {}
        ### Keeping track of where we are in huge list
        self.j = 0
        self.elev_j = 0

        ### Variables to be extracted from Measurement Files ###

        self.pr = np.full(int(19e7), -1.)
        self.all_times = np.full(int(19e7), -1.)
        self.wnc = np.full(int(19e7), -1.)
        self.c_n_0 = np.full(int(19e7), -1.)
        self.sat_id = np.empty(int(19e7), dtype=object)
        self.sig_type = np.empty(int(19e7), dtype=object)

        ### Variables to be from SatVisibility Files ###
        self.az = np.full(int(19e7), -1.)
        self.elev_new_times = np.full(int(19e7), -1.)  # allocating numpy array for times values
        self.elev_new_times_WNC = np.full(int(19e7), -1.)
        self.elev = np.full(int(19e7), -1.)  # allocating numpy array for elevations
        self.elev_new_sat_id = np.empty(int(19e7), dtype=object)  # allocating numpy array for satellite name values

    def parse_data_direcs(self, convert_to_dict=True):
        ### grab all measurement and SatVisibility files from data_direcs inputted ###
        for direc in self.data_direcs:
            self.reset_params()
            measures = []
            sat_vises = []
            for file in sorted(os.listdir(direc)):
                if file.endswith("_measurements.txt"):
                    measures.append(file)
                if file.endswith("_SBF_SatVisibility1.txt"):
                    sat_vises.append(file)
            dict_name_from_direc = direc.split("/")[-2]
            print(f"direc {direc}")
            if self.masking:
                self.dict_name = f"{self.save_parsed_data_direc}/{dict_name_from_direc}_satellite_dict_all_{self.mask_frequency}"
            else:
                self.dict_name = f"{self.save_parsed_data_direc}/{dict_name_from_direc}_satellite_dict_all"
            print(f"Saving all parsed data to {self.dict_name}.")
            for measure, sat_vis in zip(measures, sat_vises):
                self.filename_sats_pr, self.filename_sats_elevs = direc + measure, direc + sat_vis
                self.update_files(direc, measure, sat_vis)
                self.update_map()

            if convert_to_dict:
                self.convert_save_to_dict_all()

    def update_files(self, data_direc, filename_sats_pr, filename_sats_elevs):
        """
        Update files such that you can continue to parse and add to dictionary.
        """
        self.filename_sats_pr, self.filename_sats_elevs = data_direc + filename_sats_pr, data_direc + filename_sats_elevs

    def update_map(self):
        self.extract_measurements_file()
        self.extract_satvisibility_file()
        print(f"Updated instance to include {self.filename_sats_pr} and {self.filename_sats_elevs}.")


    def convert_save_to_dict_all(self):
        """
        This method converts all the TEC, elevation, sat_id, and times array into a dictionary where each key is a
        satellite and the value is a dictionary containing properties of that satellite.
        """
        # Below get dictionary with dictionary (id_dict_pr) containing sat_ids and corresponding numeric values
        # and array of ids in numeric version (id_vec_numeric_sat)
        id_dict_pr, id_vec_pr_numeric_sat = self.helper_convert_int_keys(self.sat_id)
        id_dict_elev, id_vec_elev_numeric_sat = self.helper_convert_int_keys(self.elev_new_sat_id)

        for num, key in enumerate(id_dict_pr):
            # getting masks for pseudoranges
            targ_pr = id_dict_pr[key]  # pick out a satellite
            mask = id_vec_pr_numeric_sat == targ_pr # get numeric version to use ask mask

            # getting masks for measurements
            curr_times = self.all_times[mask]
            curr_wnc = self.wnc[mask]
            curr_cnos = self.c_n_0[mask] # grabbing cnos and times for a specific satellite
            curr_pr = self.pr[mask] # grabbing cnos and times for a specific satellite
            curr_sig_type = self.sig_type[mask] # grabbing signal types for each observation
            # getting masks for SatVisibility
            if key not in id_dict_elev.keys():
                continue
            targ_elev = id_dict_elev[key]  # pick out a satellite
            mask_elev = id_vec_elev_numeric_sat == targ_elev
            curr_elev = self.elev[mask_elev]
            curr_az = self.az[mask_elev]
            time_elev = self.elev_new_times[mask_elev]
            wnc_time_elev = self.elev_new_times_WNC[mask_elev]

            if key not in self.all_sat_dict:
                self.all_sat_dict[key] = {}
                self.all_sat_dict[key]["cno"] = []
                self.all_sat_dict[key]["times"] = []
                self.all_sat_dict[key]["wnc"] = []
                self.all_sat_dict[key]["pr"] = []
                self.all_sat_dict[key]["sig_type"] = []

                if self.include_elev:
                    self.all_sat_dict[key]["elevations"] = []
                    self.all_sat_dict[key]["azimuths"] = []
                    self.all_sat_dict[key]["elev_time"] = []
                    self.all_sat_dict[key]["wnc_elev_time"] = []

            # adding to current satellite key
            self.all_sat_dict[key]["cno"].append(curr_cnos)
            self.all_sat_dict[key]["times"].append(curr_times)
            self.all_sat_dict[key]["wnc"].append(curr_wnc)
            self.all_sat_dict[key]["pr"].append(curr_pr)
            self.all_sat_dict[key]["sig_type"].append(curr_sig_type)

            if self.include_elev:
                self.all_sat_dict[key]["elevations"].append(curr_elev)
                self.all_sat_dict[key]["azimuths"].append(curr_az)
                self.all_sat_dict[key]["elev_time"].append(time_elev)
                self.all_sat_dict[key]["wnc_elev_time"].append(wnc_time_elev)
        self.save_sat_dict()

    def save_sat_dict(self):
        np.save("{}.npy".format(self.dict_name), self.all_sat_dict)


    def get_sat_dict(self):
        """
        Return current satellite dict for instance
        """
        return self.all_sat_dict

    def replace_sat_dict(self, sat_dict):
        """
        Replace the satellite dictionary. Useful to do plotting on preprocessed data held in a .npy saved dictionary
        """
        self.all_sat_dict = sat_dict

    def extract_measurements_file(self):
        """
        Reads a pseudorange file (filenam_alt)
        mask_frequency: single frequency number as string
        """
        filename_pr = self.filename_sats_pr
        with open(filename_pr, 'r', encoding="ISO-8859-1") as f:
            all_lines = f.readlines()
        f.close()
        lines = all_lines[2:]
        n = len(lines)
        print(f"Number of lines: {n}")

        for i in range(n):
            tags = lines[i].split(',')
            if self.masking and (tags[5] == '' or tags[5] == 0.0 or self.mask_frequency not in tags[3]):
                continue

            self.pr[self.j] = float(tags[5])
            self.sat_id[self.j] = tags[2]
            self.sig_type[self.j] = tags[3]
            if self.sat_id[self.j] == "G14":
                print("A G14 signal type: ")
                print(self.sig_type[self.j])
            self.all_times[self.j] = float(tags[0])
            self.wnc[self.j] = float(tags[1])
            self.c_n_0[self.j] = float(tags[-2])
            self.j += 1

    def extract_satvisibility_file(self):
        """
        Reads a SatVisibility1 file (filename_alt)
        """
        filename_alt = self.filename_sats_elevs
        with open(filename_alt, 'r', encoding="ISO-8859-1") as f:
            all_lines = f.readlines()
        f.close()
        lines = all_lines[6:]
        n = len(lines)
        for i in range(n):
            tags = lines[i].split(',')
            if tags[3] == '':
                continue
            self.elev_new_times[self.elev_j] = float(tags[7])/1000 # time values in SatVisibility are reported in milliseconds
            self.elev_new_times_WNC[self.elev_j] = float(tags[8])
            self.elev_new_sat_id[self.elev_j] = self.convert_prn2rinexsatcode(tags[12])
            self.elev[self.elev_j] = tags[-3]
            self.az[self.elev_j] = tags[-5]
            self.elev_j += 1
            # print(self.elev_j)


    @staticmethod
    def weeksecondstoutc(gpsweek, gpsseconds, leapseconds=0):
        # this is stolen from a random stack overflow post
        datetimeformat = "%Y-%m-%d %H:%M:%S"
        epoch = datetime.datetime.strptime("1980-01-06 00:00:00", datetimeformat)
        elapsed = datetime.timedelta(days=(gpsweek * 7), seconds=(gpsseconds - leapseconds))
        utc_time = datetime.datetime.strftime(epoch + elapsed, datetimeformat)
        # print(f"converted utc_time {utc_time}")
        return utc_time

    @staticmethod
    def convert_GPStime_wrapper(times_GPS, Wnc):
        """
        Convert list of GPS times to a list of UTC times given a WNc (Week number counter)
        """
        times_str = np.empty(len(times_GPS), dtype=object)
        for i in range(len(times_GPS)):
            times_str[i] = ExtractSBF.weeksecondstoutc(gpsweek=Wnc[i], gpsseconds=times_GPS[i])
        return times_str

    @staticmethod
    def helper_convert_int_keys(sat_id):
        """
        Convert satellite IDs into numeric values stored in dictionary. Some of this written by Jon Sievers.
        :param sat_id - satellite IDs

        :return
        id_dict - dictionary where keys are satellite names and values are their corresponding IDs
        """
        ids = list(set(sat_id)) # get unique satellite values
        id_dict = {}
        for i, id in enumerate(ids): # create a dictionary where each value corresponds to an integer to represent a satellite
            id_dict[id] = i
        n = len(sat_id) # get length of saellite

        # now that we have the dictionary, we'll go back to our list and
        # convert the text names into the int keys
        id_vec = np.empty(len(sat_id), dtype='int')

        # convert sat names (sat_id) to integers (in id_vec)
        for i in range(n):
            if sat_id[i] == "" and sat_id[i] == None:
                continue
            id_vec[i] = id_dict[sat_id[i]] # changing sat_id into an array where each element corresponds to numeric representation of that satellite
        return id_dict, id_vec

    @staticmethod
    def convert_prn2rinexsatcode(prn):
        """
        This function converts PRNs to RINEX satellite codes.
        :param prn - pseudorandom number from a satellite
        :returns rinex satellite code needed to match with other parameters
        """
        prn = int(prn)
        if prn >= 1 and prn <= 37: # GPS
            if prn < 10:
                return f"G0{prn}"
            else:
                return f"G{prn}"
        if prn >= 38 and prn <= 61:  # GLONASS
            new_prn = prn-37
            if new_prn < 10:
                return f"R0{new_prn}"
            else:
                return f"R{new_prn}"

        if prn >= 63 and prn <= 68:  # GLONASS
            new_prn = prn-38
            if new_prn < 10:
                return f"R0{new_prn}"
            else:
                return f"R{new_prn}"

        if prn >= 71 and prn <= 106:  # Galileo
            new_prn = prn-70
            if new_prn < 10:
                return f"E0{new_prn}"
            else:
                return f"E{new_prn}"
        if prn >= 120 and prn <= 140:  # SBAS
            new_prn = prn-100
            return f"S{prn-100}"
        if prn >= 141 and prn <= 180:  # Beidou
            new_prn = prn-140
            if new_prn < 10:
                return f"C0{new_prn}"
            else:
                return f"C{new_prn}"
        if prn >= 181 and prn <= 187:  # QZSS
            new_prn = prn-180
            if new_prn < 10:
                return f"J0{new_prn}"
            else:
                return f"J{new_prn}"

        if prn >= 223 and prn <= 245:  # Beidou
            new_prn = prn - 182
            if new_prn < 10:
                return f"C0{new_prn}"
            else:
                return f"C{new_prn}"
        else:
            return None