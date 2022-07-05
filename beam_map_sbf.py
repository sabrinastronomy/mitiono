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

            # difference_array = np.absolute(elevs - 81.41)
            # TO DO add azimuthal dependence
            # index = difference_array.argmin()

            max_cno = max(cnos)
            normalized_power = power - max_cno
            plt.scatter(az, elevs, c=normalized_power)

        plt.scatter(180, 81.41, s=10)
        plt.xlabel("Azimuth [deg]")
        plt.ylabel("Elevation [deg]")
        plt.colorbar(label="Power [dBW]")
        plt.title(plot_title)

        plt.savefig(direc + filename, dpi=400)

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

#### Sample code to grab parsed dictionary file ######
# directories = ["/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/mitiono_scripts/parsed_data/"]
# sat_dict = np.load(directories[0] + "june14_satellite_dict_all.npy", allow_pickle=True).item()  # extracting saved dictionary
#####

if __name__ == "__main__":
    directories = ["/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/D3A(6)_Data/June16_port_C2_GNSS_part2/"]
    beam_map_14 = GetBeamMap(data_direcs=directories)
    beam_map_14.make_plot(all_sat=True, plot_title="All Satellite Beam Plot")
