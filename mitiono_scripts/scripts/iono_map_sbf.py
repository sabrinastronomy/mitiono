import matplotlib.pyplot as plt
import numpy as np
from extract_sbf_data import ExtractSBF
import datetime
import matplotlib.dates as dates

# IN PROGRESS

class GetIonoMap(ExtractSBF):
    """
    This class takes in a data directory and parses the file for TEC, elevations, times, satellite IDs. It also has
    plotting capabilities/
    """
    def __init__(self, data_direc, filename_sats_pr='', filename_sats_alts='', dict_name="satellite_dict", min_elev=0, process=False, plot=False, include_elev=False, plot_name=""):
        """
        :param data_direc - directory with receiver raw files
        :param filename_sats_pr - file containing pseudorange data
        :param filename_sats_alts -  file containing elevation/altitude data
        :param dict_name - dictionary name to save parsed data to
        :param min_elev - minimum elevation of included satellites
        :param process - whether or not to parse data
        :param plot - whether or not to plot data
        :param include_elev - whether or not to extract elevations
        """

        self.all_sat_dict = {}
        self.filename_sats_pr, self.filename_sats_alts = data_direc + filename_sats_pr, data_direc + filename_sats_alts
        self.min_elev = min_elev
        self.include_elev = include_elev
        self.dict_name = dict_name
        self.tec = np.full(int(19e6), -1.)  # allocating numpy array for possible TEC values
        self.new_times = np.empty(int(19e6), dtype=object)  # allocating numpy array for times values

        self.new_sat_id = np.empty(int(19e6), dtype=object)  # allocating numpy array for satellite name values

        ### Extracted from SatVisibility Files ###
        self.elev_new_times = np.full(int(19e6), -1.)  # allocating numpy array for times values
        self.elev_new_times_WNC = np.full(int(19e6), -1.)
        self.elev = np.full(int(19e6), -1.)  # allocating numpy array for elevations
        self.elev_new_sat_id = np.empty(int(19e6), dtype=object)  # allocating numpy array for satellite name values
        self.extract_elev_time()
        ##########################################

        ### Extracted from Measurement Files ###
        self.pr = np.full(int(19e6), -1.)
        self.all_times = np.full(int(19e6), -1.)
        self.wnc = np.full(int(19e6), -1.)
        self.c_n_0 = np.full(int(19e6), -1.)
        self.sat_id = np.empty(int(19e6), dtype=object)
        self.sig_type = np.empty(int(19e6), dtype=object)
        self.extract_measurements_file()

        self.j = 0
        self.elev_j = 0
        self.plot_name = plot_name

        if process: # process file to extract data
            self.extract_tec_elev_times_sat_id()
            self.convert_save_to_dict_iono_map()
        if plot: # plot data
            self.plot_process()

    # def create_time2elevation_func(self, time, elevations):
    #     """
    #     This function returns a function that takes in an input time (in the following form: "%Y-%m-%d %H:%M:%S")
    #     and outputs an elevation. Time range will be limited between min and max of interpolated times.
    #     TODO: Make sure this is working correctly
    #     :param elevations - to be interpolated elevations
    #     :param times - corresponding times to be independent variable in interpolation
    #     :return func_intermediate (function that takes in input time and outputs corresponding time,
    #     origin is minimum time for this sat)
    #     min_time,
    #     min_time
    #     """
    #     form = "%Y-%m-%d %H:%M:%S"
    #     min_time = datetime.datetime.strptime(time[0], form)
    #     max_time = datetime.datetime.strptime(time[-1], form)
    #     new_time_deltas = np.empty(len(time))
    #     for i, t in enumerate(time):
    #         t = datetime.datetime.strptime(t, form)
    #         new_time_deltas[i] = (t - min_time).seconds # origin is first time
    #     u, inds = np.unique(new_time_deltas, return_index=True) # getting unique time indices
    #     new_time_deltas = new_time_deltas[inds]
    #     elevations = elevations[inds]
    #     # func_intermediate = interpolate.interp1d(new_time_deltas, elevations, kind='linear', bounds_error=False)
    #     return func_intermediate, min_time, max_time

    def time2elevation_func(self, time_formatted, func_intermediate, min_time, max_time):
        """
        Function to get elevation at closest time.
        """
        form = "%Y-%m-%d %H:%M:%S"
        time_formatted = datetime.datetime.strptime(time_formatted, form)
        tolerance = datetime.timedelta(minutes=10)
        assert time_formatted > (min_time - tolerance)
        assert time_formatted < (max_time + tolerance)

        curr_time_secs_from_min = (time_formatted - min_time).seconds
        elev_curr_time = func_intermediate(curr_time_secs_from_min)
        return elev_curr_time

    def update_files(self, data_direc, filename_sats_pr, filename_sats_alts):
        """
        Update files such that you can continue to parse and add to dictionary.
        """
        self.filename_sats_pr, self.filename_sats_alts = data_direc + filename_sats_pr, data_direc + filename_sats_alts

    def convert_GPStime_wrapper(self, times_GPS, Wnc):
        """
        Convert list of GPS times to a list of strings given a WNc (Week number counter)
        """
        times_str = np.empty(len(times_GPS), dtype=object)
        for i in range(len(times_GPS)):
            times_str[i] = self.weeksecondstoutc(gpsweek=Wnc[i], gpsseconds=times_GPS[i])
        return times_str


    def extract_tec_elev_times_sat_id(self):
        """
        This method extracts information from the input file and places them into the already allocated
        TEC, elevations, times, and satellite_ids array. It creates vectorized data that allows for efficient
        processing and manipulation.
        """

        wnc, all_times, sig_type, sat_id, pr = self.wnc, self.all_times, self.sig_type, self. sat_id, self.pr # note this contains extracted information for all files parsed so far
        utc_time = self.convert_GPStime_wrapper(all_times, wnc)
        i = 1
        while i < len(pr): # TODO: vectorize this. get rid of while loop
            sat_name = sat_id[i]
            if sig_type[i] != None and sig_type[i] != None and sat_id[i] == sat_id[i - 1] and "1" in sig_type[i - 1] and ("2" in sig_type[i] or "5" in sig_type[i]):
                # checks if previous sat ID equals current one (i) and that L1 and L2/L5 are in the right places
                delta_pr = pr[i] - pr[i - 1]  # subtracting lower frequency PR from higher frequency PR
                next_tec = self.get_dTEC_first_order(delta_pr, signal_types=[sig_type[i - 1], sig_type[i]])
                # using class variable self.j to keep track of where you are in huge array of tecs, times, and sat_IDs
                self.tec[self.j] = next_tec
                self.new_times[self.j] = utc_time[i]
                self.new_sat_id[self.j] = sat_name
                i += 2
                self.j += 1
            i += 1

    def convert_save_to_dict_iono_map(self):
        """
        This method converts all the TEC, elevation, sat_id, and times array into a dictionary where each key is a
        satellite and the value is a dictionary containing properties of that satellite.
        """
        # Below get dictionary with dictionary (id_dict_pr) containing sat_ids and corresponding numeric values
        # and array of ids in numeric version (id_vec_numeric_sat)
        id_dict_pr, id_vec_pr_numeric_sat = self.helper_convert_int_keys(self.new_sat_id)
        if self.include_elev:
            id_dict_elev, id_vec_elev_numeric_sat = self.helper_convert_int_keys(self.new_sat_id)

        for num, key in enumerate(id_dict_pr):
            print(key)
            # getting masks for pseudoranges
            targ_pr = id_dict_pr[key]  # pick out a satellite
            mask = id_vec_pr_numeric_sat == targ_pr # get numeric version to use ask mask

            # getting masks for elevations
            if self.include_elev:
                targ_elev = id_dict_elev[key]  # pick out a satellite
                mask_elev = id_vec_elev_numeric_sat == targ_elev
                curr_elev = self.elev[mask_elev]
                time_elev = self.elev_new_times[mask_elev]
                wnc_time_elev = self.elev_new_times_WNC[mask_elev]
                elev_approx = np.max(curr_elev)
            else:
                elev_approx = 100 # random elev so it will enter if statement below when not including elevations

            if elev_approx > self.min_elev:
                if key not in self.all_sat_dict:
                    self.all_sat_dict[key] = {}
                    self.all_sat_dict[key]["tec"] = []
                    self.all_sat_dict[key]["times"] = []
                    if self.include_elev:
                        self.all_sat_dict[key]["elev_approx"] = []
                        self.all_sat_dict[key]["elev_time_approx"] = []
                        self.all_sat_dict[key]["wnc_elev_time_approx"] = []

                curr_tec = self.tec[mask]
                curr_new_times = self.new_times[mask]
                self.all_sat_dict[key]["tec"].append(curr_tec)
                self.all_sat_dict[key]["times"].append(curr_new_times)
                if self.include_elev:
                    self.all_sat_dict[key]["elev_approx"].append(elev_approx)
                    self.all_sat_dict[key]["elev_time_approx"].append(time_elev)

        np.save("{}.npy".format(self.dict_name), self.all_sat_dict)
        return

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

    def extract_measurements_file(self, filename_pr):
        """
        Reads a pseudorange file (filenam_alt)
        """
        with open(filename_pr, 'r', encoding="ISO-8859-1") as f:
            all_lines = f.readlines()
        f.close()
        lines = all_lines[2:]
        n = len(lines)
        print(f"Number of lines: {n}")

        for i in range(n):
            tags = lines[i].split(',')
            if tags[5] == '' or tags[5] == 0.0:
                continue
            self.pr[self.j] = float(tags[5])
            self.sat_id[self.j] = tags[2]
            self.sig_type[self.j] = tags[3]
            self.all_times[self.j] = float(tags[0])
            self.wnc[self.j] = float(tags[1])
            self.c_n_0[self.j] = float(tags[-2])
            print(f"c_n_0 is { self.c_n_0[self.j]}")
            self.j += 1
        return

    def extract_elev_time(self, filenam_alt):
        """
        Reads a SatVisibility1 file (filename_alt)
        """
        with open(filenam_alt, 'r', encoding="ISO-8859-1") as f:
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
            self.elev_new_sat_id[self.elev_j] = tags[12]
            self.elev[self.elev_j] = tags[-3]
            self.az[self.elev_j] = tags[-6]
            self.elev_j += 1
        return

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
            if sat_id[i] == "":
                continue
            id_vec[i] = id_dict[sat_id[i]] # changing sat_id into an array where each element corresponds to numeric representation of that satellite
        return id_dict, id_vec

    @staticmethod
    def get_dTEC_first_order(delta_rho, signal_types=None):
        """
        Return a sTEC measurement from a given pseudorange difference (delta_rho)

        :param delta_rho - difference in pseudoranges
        :signal_types - types of signals we are using

        :return sTEC measurement
        """
        if signal_types is None:
            signal_types = []
        r_e = 2.8179e-15  # [m]
        c = 3e8  # [m/s]
        kappa = c ** 2 * r_e / (2 * np.pi)
        L1 = 1575e6  # Hz
        if "5" in signal_types[1]:
            L5 = 1176e6  # Hz
            f_sq = np.abs(1 / L5 ** 2 - 1 / L1 ** 2)
        elif "2" in signal_types[1]:
            L2 = 1227e6  # Hz
            f_sq = np.abs(1 / L2 ** 2 - 1 / L1 ** 2)

        rho_iono = np.asarray(delta_rho) / f_sq
        diff_TEC_TECu = rho_iono / (kappa * 10 ** 16)  # this returns the TEC in TECu
        return diff_TEC_TECu  # returns differential TEC measurement


    def plot_process(self, TEC=False, no_negative_stec=True):
        """
        :param WNc - week number count (GNSS time thing)
        :param TEC - whether to plot sTEC or TEC, not sure if this is working
        Plot satellites sTECs versus time.
        """
        num_curr = 0
        print(f"total number of satellites: {len(self.all_sat_dict.keys())}")
        # for key in self.all_sat_dict.keys():
        for key in ["G12"]:

            if key == None:
                continue

            num_curr += 1
            stecs = np.concatenate(self.all_sat_dict[key]["tec"]).ravel()
            times = np.concatenate(self.all_sat_dict[key]["times"]).ravel()
            stec_length = len(stecs)

            # if stec_length > 1e8:
            #     mod_stecs = stec_length % 100
            #     stecs = stecs[:mod_stecs*100]
            #     stecs = np.max(stecs.reshape(-1, 100), axis=1) # use this to average every 100 elements
            #
            #     times = times[:mod_stecs*100]
            #     times = np.max(times.reshape(-1, 100), axis=1) # use this to average every 100 elements

            if min(stecs) < 0 and no_negative_stec:
                continue

            # times_str = self.convert_GPStime_wrapper(times, self.WNc)

            # elevation = np.concatenate(self.all_sat_dict[key]["elev_approx"]).ravel()
            # time_elevation = np.concatenate(self.all_sat_dict[key]["elev_time_approx"]).ravel()/1000 # times reported in 0.001s for SatVis files for some reason ugh
            # times_elev_str = self.convert_GPStime_wrapper(time_elevation, self.WNc)

            # interpolate elevations
            # make function (func_intermediate), TODO below interpolation of elevations may not be working correctly

            # func_intermediate, min_time, max_time = self.create_time2elevation_func(times_elev_str, elevation)
            # new_elevations = np.empty(len(times))
            # for k, tim in enumerate(times_str):
            #     form = "%Y-%m-%d %H:%M:%S"
            #     t_test = datetime.datetime.strptime(tim, form)
            #     if t_test > min_time and t_test < max_time:
            #         new_elevations[k] = self.time2elevation_func(tim, func_intermediate, min_time, max_time)
            #     else:
            #         new_elevations[k] = None

            plt_dates = dates.datestr2num(times_str)
            plt.xticks(rotation=90)
            print(f"plotting {num_curr}")

            if TEC:
                tecs = stecs *  1/np.cos(np.deg2rad(90-new_elevations))
                plt.scatter(plt_dates, tecs, alpha=0.5, s=1, label=key + " TEC")
                plt.scatter(plt_dates, stecs, alpha=0.5, s=1, label=key + " sTEC")
                # plt.xlim(plt_dates[1000])

            else:
                plt.scatter(plt_dates, stecs, alpha=1, s=1, label=key)
            plt.gca().xaxis.set_major_formatter(dates.DateFormatter('%m-%d %H:%M'))
            plt.gca().xaxis.set_major_locator(dates.HourLocator(interval=5))

        if TEC:
            plt.ylabel("TEC [TECu]")
        else:
            plt.ylabel("sTEC [TECu]")
        plt.tight_layout()
        plt.legend(prop={'size': 6})
        # plt.legend(ncol=10, prop={'size': 3})
        plt.savefig(f"../plots/{self.plot_name}_{self.min_elev}.png", dpi=300)
        plt.close()

