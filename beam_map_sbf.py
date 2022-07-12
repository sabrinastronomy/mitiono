"""
This file extracts data from SBF SatVisibility files returned by the receiver to generate a beam map.

Written by Sabrina Berger
"""

from extract_sbf_data import ExtractSBF
import numpy as np
import matplotlib.pyplot as plt


class GetBeamMap(ExtractSBF):
    """
    This class takes in data directories and created a beam map.
    """
    def __init__(self, data_direcs):
        super().__init__(data_direcs=data_direcs, include_elev=True, process=True)

    def print_dictionary(self):
        print(self.all_sat_dict)
        return

    def replace_dictionary(self, dict):
        self.all_sat_dict = dict

    def make_plot(self, all_sat, plot_title, sat_list=[], direc="/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/plots/", filename="all_sat.png"):
        """
        Plot all the satellites or a subset
        param: all_sat - boolean when true all satellites are plotted, when false, looks for satellite in sat_list
        param: plot_title - title above plot to be generated
        param: direc - directory where to save the plot
        param: filename - name of plot you're making
        """
        plt.close()
        max_power = -200
        min_power = 0
        if all_sat:
            sats_to_plot = self.all_sat_dict
        else:
            sats_to_plot = sat_list

        for key in sats_to_plot:
            satellite = sats_to_plot[key]
            times_cno = satellite["times"][0]
            times_elevs = satellite["elev_time"][0]
            elevs = satellite["elevations"][0]
            az = satellite["azimuths"][0]
            cnos = satellite["cno"][0]
            times_elevs, cnos, elevs, az = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
            power = self.convert_P(cnos)  # dBw
            # close_to_beam_bool = np.logical_and(elevs > 80, elevs < 82.5).any()


            elev_beam = np.full(len(elevs), 81.41)
            diff = np.abs(elevs - elev_beam)
            print(diff)
            elev_bool_tups = diff < 2


            az_beam = np.full(len(az), 180)
            diff = np.abs(az - az_beam)
            az_bool_tups = diff < 2
            bool_tups = np.logical_and(elev_bool_tups, az_bool_tups)

            print(bool_tups.any())
            # if not bool_tups.any():
            #     continue
            # difference_array = np.absolute(elevs - 81.41)
            # TO DO add azimuthal dependence
            # index = difference_array.argmin()
            if len(power) < 1:
                continue
            max_power_current = max(power)
            if max_power_current > max_power:
                max_power = max_power_current

            min_power_current = min(power)
            if min_power_current < min_power:
                min_power = min_power_current

            # normalized_power = power - max_cno
            plt.scatter(az, elevs, c=power, cmap='inferno_r', label=key)

        print(f"max_power {max_power}")
        print(f"min_power {min_power}")

        plt.colorbar(label="Power [dBW]")
        plt.xlim(160, 190)
        plt.ylim(70, 90)
        plt.clim(-130, -180)
        # plt.legend()
        plt.scatter(180, 81.41, s=100, facecolors='none', edgecolors='y', lw=2, zorder=100)
        plt.xlabel("Azimuth [deg]")
        plt.ylabel("Elevation [deg]")
        plt.title(plot_title)

        plt.savefig(direc + filename, dpi=400)
        plt.close()

    @staticmethod
    def convert_P(C_N):
        """
        This function converts the receiver's C/N_0 measurement into a classical power measurement in dBw
        Convert C/N_0 to power in dBw.
        """
        N_sys = 8.5  # dB
        Tant = 290  # K
        P = C_N + 10 * np.log10(Tant + 290 * (10 ** (N_sys / 10) - 1)) - 228.6  # last i power
        return P

    @staticmethod
    def match_elevs(times_cno, times_elevs, cnos, elevs, az):

        indices = np.intersect1d(times_elevs, times_cno, return_indices=True)
        elev_indices = indices[1]
        cnos_indices = indices[2]

        times_elevs = times_elevs[elev_indices]
        elevs = elevs[elev_indices]
        az = az[elev_indices]
        cnos = cnos[cnos_indices]
        return times_elevs, cnos, elevs, az



if __name__ == "__main__":
    ## Sample code to grab parsed dictionary file ######
    directories = ["/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/parsed_data/"]
    # days = ["June14_port_C2_GNSS_satellite_dict_all", "June_16_port_C2_GNSS_satellite_dict_all", "June_16_port_C2_GNSS_part2_satellite_dict_all"]
    days = ["June_16_port_C2_GNSS_part2_satellite_dict_all"]
    # days = []
    a = []
    for day in days:
        sat_dict = np.load(directories[0] + day + ".npy", allow_pickle=True).item()  # extracting saved dictionary
        beam_map = GetBeamMap(data_direcs=[])
        beam_map.replace_dictionary(sat_dict)
        beam_map.make_plot(all_sat=True, plot_title="Zoomed In Around Beam", filename="all_zoomed_" + day)

    # directories = ["/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/D3A(6)_Data/June14_port_C2_GNSS/", "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/D3A(6)_Data/June_16_port_C2_GNSS/", "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/D3A(6)_Data/June_16_port_C2_GNSS_part2/"]
    # beam_map = GetBeamMap(data_direcs=directories)
