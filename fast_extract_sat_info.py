import matplotlib.pyplot as plt
import numpy as np
import datetime
import datetime
import matplotlib.dates as md


class GetIonoMap:
    def __init__(self, data_direc, filename_sats_pr='', filename_sats_alts='', dict_name="satellite_dict", min_elev=0, process=False, plot=False, include_elev=False):
        self.all_sat_dict = {}
        self.filename_sats_pr, self.filename_sats_alts = data_direc + filename_sats_pr, data_direc + filename_sats_alts
        self.min_elev = min_elev
        self.include_elev = include_elev
        self.dict_name = dict_name
        self.tec = np.full(int(19e6), -1.)  # allocating numpy array for possible TEC values
        self.new_times = np.full(int(19e6), -1.)  # allocating numpy array for times values
        self.new_sat_id = np.empty(int(19e6), dtype=object)  # allocating numpy array for satellite name values

        self.elev_new_times = np.full(int(19e6), -1.)  # allocating numpy array for times values
        self.elev_new_sat_id = np.empty(int(19e6), dtype=object)  # allocating numpy array for satellite name values
        self.elev = np.full(int(19e6), -1.)  # allocating numpy array for elevations
        self.j = 0
        self.elev_j = 0
        if process:
            self.extract_tec_elev_times_sat_id()
            self.convert_to_dict()
        if plot:
            self.plot_process()

    def reset_initial_time(self):
        self.begin_time_arr = 0

    def update_files(self, data_direc, filename_sats_pr, filename_sats_alts):
        self.filename_sats_pr, self.filename_sats_alts = data_direc + filename_sats_pr, data_direc + filename_sats_alts

    def plot_process(self, WNc):
        i = 0
        print(f"total number of satellites: {len(self.all_sat_dict)}")
        for key in self.all_sat_dict.keys():
            print(f"elevations {self.all_sat_dict[key]['elev_approx']}")
            i += 1

            tecs = np.concatenate(self.all_sat_dict[key]["tec"]).ravel()
            tecs = tecs[::100]
            if min(tecs) < 0:
                continue

            times = np.concatenate(self.all_sat_dict[key]["times"]).ravel()
            times = times[::100]

            times_str = np.empty(len(times), dtype=object)
            for i in range(len(times)):
                # print(f"times[i] {times[i]}")
                times_str[i] = self.weeksecondstoutc(gpsweek=WNc, gpsseconds=times[i])

            # xfmt = md.DateFormatter('%H:%M')
            # ax = plt.gca()
            plt.xticks(rotation=90)
            # ax.xaxis.set_major_formatter(xfmt)
            print(f"plotting {i}")
            print(times_str)
            import matplotlib.dates as dates

            plt_dates = dates.datestr2num(times_str)
            plt.scatter(plt_dates, tecs, label=key, alpha=0.5)
            plt.gca().xaxis.set_major_formatter(md.DateFormatter('%m-%d %H:%M'))
            plt.gca().xaxis.set_major_locator(md.HourLocator(interval=5))
            # plt.setp(plt.gca().get_xticklabels(), rotation=60, ha="right")

        plt.tight_layout()
        # plt.xlabel("time of day [hours]")
        plt.ylabel("sTEC [TECu]")
        plt.legend()
        print("if paused below here, showing figure is taking forever.")
        plt.show()
        print("if paused below here, saving figure is taking forever.")
        plt.savefig(f"stecs_time_elev_min_{self.min_elev}.png")
        plt.close()

    def extract_tec_elev_times_sat_id(self):
        """
        This method extracts information from the input file and places them into the already allocated
        TEC, elevations, times, and satellite_ids array. It creates vectorized data that allows for efficient
        processing and manipulation.
        """
        if self.include_elev:
            self.extract_elev_time(self.filename_sats_alts)

        all_times, sig_type, sat_id, pr = self.extract_pr_time(self.filename_sats_pr)

        i = 1
        while i < len(pr): # TODO: vectorize this. get rid of while loop
            # print(f"currently working on pr {i}")
            sat_name = sat_id[i]
            if sat_id[i] == sat_id[i - 1] and "1" in sig_type[i - 1] and ("2" in sig_type[i] or "5" in sig_type[i]):
                # checks if previous sat ID equals current one (i) and that L1 and L2/L5 are in the right places
                delta_pr = pr[i] - pr[i - 1]  # subtracting lower frequency PR from higher frequency PR
                next_tec = self.get_dTEC_first_order(delta_pr, signal_types=[sig_type[i - 1], sig_type[i]])
                # using class variable self.j to keep track of where you are in huge array of tecs, times, and sat_IDs
                self.tec[self.j] = next_tec
                self.new_times[self.j] = all_times[i]
                self.new_sat_id[self.j] = sat_name
                i += 2
                self.j += 1
            i += 1

    def convert_to_dict(self):
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

            # getting masks for pseudoranges
            targ_pr = id_dict_pr[key]  # pick out a satellite
            mask = id_vec_pr_numeric_sat == targ_pr # get numeric version to use ask mask

            # getting masks for elevations
            if self.include_elev:
                targ_elev = id_dict_elev[key]  # pick out a satellite
                mask_elev = id_vec_elev_numeric_sat == targ_elev
                curr_elev = self.elev[mask_elev]
                elev_approx = curr_elev[0]
            else:
                elev_approx = 100 # random elev so it will enter if statement below when not including elevations

            if elev_approx > self.min_elev:
                if key not in self.all_sat_dict:
                    self.all_sat_dict[key] = {}
                    self.all_sat_dict[key]["tec"] = []
                    self.all_sat_dict[key]["times"] = []
                    if self.include_elev:
                        self.all_sat_dict[key]["elev_approx"] = []

                curr_tec = self.tec[mask]
                curr_new_times = self.new_times[mask]
                self.all_sat_dict[key]["tec"].append(curr_tec)
                self.all_sat_dict[key]["times"].append(curr_new_times)
                if self.include_elev:
                    self.all_sat_dict[key]["elev_approx"].append(elev_approx)

        np.save("../data/{}.npy".format(self.dict_name), self.all_sat_dict)
        return


    def get_sat_dict(self):
        return self.all_sat_dict

    def replace_sat_dict(self, sat_dict):
        self.all_sat_dict = sat_dict

    def extract_pr_time(self, filename_pr):
        with open(filename_pr, 'r', encoding="ISO-8859-1") as f:
            all_lines = f.readlines()
        f.close()
        lines = all_lines[2:]
        n = len(lines)
        pr = np.empty(n)
        all_times = np.empty(n)
        sat_id = [None]*n
        sig_type = [None]*n
        # init_time = float(lines[0].split(',')[0])
        for i in range(n):
            tags = lines[i].split(',')
            pr[i] = float(tags[5])
            sat_id[i] = tags[2]
            sig_type[i] = tags[3]
            all_times[i] = float(tags[0])
            # print(f"should be TOW: {all_times[i]}")
        return all_times, sig_type, sat_id, pr

    def extract_elev_time(self, filenam_alt):
        # READS SATVISIBILITY1 FILE
        with open(filenam_alt, 'r', encoding="ISO-8859-1") as f:
            all_lines = f.readlines()
        f.close()
        lines = all_lines[6:]
        n = len(lines)
        for i in range(n):
            tags = lines[i].split(',')
            if tags[3] == '':
                continue
            self.elev_new_times[self.elev_j] = float(tags[7])
            self.elev_new_sat_id[self.elev_j] = tags[12]
            self.elev[self.elev_j] = tags[-3]
            self.elev_j += 1

        return

    @staticmethod
    def weeksecondstoutc(gpsweek, gpsseconds, leapseconds=0):
        # random stack overflow post
        datetimeformat = "%Y-%m-%d %H:%M:%S"
        epoch = datetime.datetime.strptime("1980-01-06 00:00:00", datetimeformat)
        elapsed = datetime.timedelta(days=(gpsweek * 7), seconds=(gpsseconds - leapseconds))
        utc_time = datetime.datetime.strftime(epoch + elapsed, datetimeformat)
        # print(f"converted utc_time {utc_time}")
        return utc_time

    @staticmethod
    def helper_convert_int_keys(sat_id):
        # ELEVS
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



