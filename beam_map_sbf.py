"""
This file extracts data from SBF SatVisibility files returned by the receiver to generate a beam map.

Written by Sabrina Berger
"""


# import matplotlib.font_manager
# from IPython.core.display import HTML
# font = "Gill Sans MT"
#
# def make_html(fontname):
#     return "<p>{font}: <span style='font-family:{font}; font-size: 24px;'>{font}</p>".format(font=fontname)
#
# code = "\n".join([make_html(font) for font in sorted(set([f.name for f in matplotlib.font_manager.fontManager.ttflist]))])
#
# # HTML("<div style='column-count: 2;'>{}</div>".format(code))


from extract_sbf_data import ExtractSBF
import numpy as np
import matplotlib.pyplot as plt
from scipy import special

D3A_ALT_deg = 81.41
D3A_AZ_deg = 180

D3A_ALT = 81.41 * np.pi/180
D3A_AZ = 180 * np.pi/180

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
        plt.close()
        for key in sats_to_plot:
            satellite = sats_to_plot[key]
            times_cno = satellite["times"][0]
            times_elevs = satellite["elev_time"][0]
            elevs = satellite["elevations"][0]
            az = satellite["azimuths"][0]
            cnos = satellite["cno"][0]
            times_elevs, cnos, elevs, az = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
            power = self.convert_P(cnos)  # dBw
            #
            # mask = power > 20
            # az = az[mask]
            # elevs = elevs[mask]
            # power = power[mask]

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
            power -= 59.875
            # normalized_power = power - max_cno
            plt.scatter(az, elevs, s=0.1, c=power)

        print(f"max_power {max_power}")
        print(f"min_power {min_power}")

        fig = plt.gcf()
        ax = fig.gca()


        az_mask = az
        plt.colorbar(label="Log Power [uncalibrated]")
        plt.xlim(160, 190)
        plt.ylim(70, 85)
        # _min_power =
        # _max_power =
        plt.clim(-60, 0)
        # plt.legend()
        circle1 = plt.Circle((D3A_AZ_deg, D3A_ALT_deg), radius=1.25, color='k', fill=None)
        ax.add_patch(circle1)
        plt.scatter(None, None, c="k", label="D3A(6) Beam")

        # plt.scatter(D3A_AZ_deg, D3A_ALT_deg, s=100, , facecolors='none', edgecolors='k', lw=2)
        plt.xlabel("Azimuth [deg]")
        plt.ylabel("Elevation [deg]")
        plt.legend()
        plt.title("Satellites Near D3A(6) Beam over 24 Hours")
        plt.savefig(direc + filename, dpi=400, bbox_inches='tight')
        plt.close()

    @staticmethod
    def get_bessel_0(ka_sintheta):
        return special.j0(ka_sintheta)

    @staticmethod
    def get_bessel_1(ka_sintheta):
        return special.j1(ka_sintheta)

    def get_airy_total_intensity(self, theta, max_inten, k=2*np.pi/0.2, a=3):
        ka_sintheta = k * a * np.sin(theta) # theta is in rad
        bessel_1 = self.get_bessel_1(ka_sintheta)
        return max_inten * (2*bessel_1 / ka_sintheta)**2

    def get_airy_total_pow(self, theta, max_power, k=2*np.pi/0.2, a=3):
        ka_sintheta = k * a * np.sin(theta) # theta is in rad
        bessel_0 = self.get_bessel_0(ka_sintheta)**2
        bessel_1 = self.get_bessel_1(ka_sintheta)**2
        return max_power * (1 - bessel_0 - bessel_1)

    def hist_single_sat(self):
        """
        1D profile where y = power, x = distance from pointing center
        Really quickly - the el vs. az is because the actual angular distance is az*cos(el).
        You won’t get more likelihood of direct beam center hit with more dishes because they are all pointed in the same direction.
        However, for beam modelling, it’s not that important to get a bulls-eye, because the beams ought to be smooth.
        The single-satellite pass makes a lot more sense to me.
        I would take the data from that and turn it into a 1-d profile where you actually calculate the distance from the pointing center.
        Because I suck at spherical trig, the way I do it is to turn the pointing center into an x,y,z vector,
        the satellite position into an x,y,z vector, and take the dot product to get cos(angle).
        """

        satellite = self.all_sat_dict["E03"]

        elevs = satellite["elevations"][0]
        az = satellite["azimuths"][0]
        cnos = satellite["cno"][0]
        times_cno = satellite["times"][0]
        times_elevs = satellite["elev_time"][0]

        times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
        elevs = elevs_deg * np.pi/180
        az = az_deg * np.pi/180

        # masking
        mask_az = az_deg < D3A_AZ_deg # left
        mask = mask_az

        left_az = az[mask]
        left_elevs = elevs[mask]
        left_cnos = cnos[mask]

        left_angles_rad = self.convert_angular_distance_from_center_beam(left_elevs, left_az)
        left_angles_deg = left_angles_rad * 180/np.pi
        left_angles_deg = -left_angles_deg
        left_powers = self.convert_P(left_cnos)

        mask_az = az_deg > D3A_AZ_deg # right
        mask = mask_az

        right_az = az[mask]
        right_elevs = elevs[mask]
        right_cnos = cnos[mask]

        right_angles_rad = self.convert_angular_distance_from_center_beam(right_elevs, right_az)
        right_angles_deg = right_angles_rad * 180/np.pi
        right_powers = self.convert_P(right_cnos)

        # getting perfect beam
        # angles = np.concatenate((left_angles_rad, right_angles_rad))

        # _right_angles_perfect = np.linspace(1e-3, np.pi/2, 1000)
        # _left_angles_perfect = np.linspace(1e-3, np.pi/2, 1000)
        #
        # pow_W = 10 ** (max(left_powers)/ 10)
        # p_perfect = self.get_airy_total_intensity(_left_angles_perfect, pow_W)
        # pow_dbW = 10 * np.log10(p_perfect)
        # _left_angles_perfect *= (180/np.pi)
        # plt.plot(_left_angles_perfect, pow_dbW, c="red")
        #
        # print(max(right_powers))
        # pow_W = 10 ** (max(right_powers)/ 10)
        # p_perfect = self.get_airy_total_intensity(_right_angles_perfect, pow_W)
        # pow_dbW = 10 * np.log10(p_perfect)
        # _right_angles_perfect *= -(180/np.pi)

        # plt.plot(_right_angles_perfect, pow_dbW, c="red", label="Analytic Beam Pattern")
        # plt.show()
        # exit()

        direc = "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/notebook_and_beams/"

        # print(np.shape(powers))
        right_angles_deg_sim = np.load(direc + "theta_right.npy")
        left_angles_deg_sim = np.load(direc + "theta_left.npy")
        left_powers_sim = np.load(direc + "power_left.npy")
        right_powers_sim = np.load(direc + "power_right.npy")

        plt.plot(left_angles_deg_sim, left_powers_sim - 59.875, c="r")
        plt.plot(right_angles_deg_sim, right_powers_sim - 59.875, c="r", label="D3A(6) Simulated Beam Shape at 1.5 GHz")

        plt.scatter(left_angles_deg, left_powers - 59.875, c="k", s=0.5)
        plt.scatter(right_angles_deg, right_powers - 59.875, c="k", s=0.5, label="GNSS Beam Measurement at 1.575 GHz")
        plt.xlim(-20, 20)

        print(max(right_powers) - max(right_powers_sim))

        # plt.scatter(angles, p_perfect, label="Airy Pattern Beam Model")
        plt.xlabel(r"$\theta$, the angle from beam center [deg]")
        plt.ylabel("Log Power [uncalibrated]")
        plt.title("Galileo Satellte: E03")
        plt.legend()
        direc = "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/plots/"
        plt.savefig(direc + "poster_ang_dist.png", dpi=400, bbox_inches='tight')
        plt.show()



    @staticmethod
    def convert_angular_distance_from_center_beam(alt, az):

        x_beam = np.cos(D3A_AZ) * np.cos(D3A_ALT)
        y_beam = np.sin(D3A_AZ) * np.cos(D3A_ALT)
        z_beam = np.sin(D3A_ALT)

        x_sat = np.cos(az) * np.cos(alt)
        y_sat = np.sin(az) * np.cos(alt)
        z_sat = np.sin(alt)
        # sat_loc_xyz = np.transpose((x_sat, y_sat, z_sat))
        # beam_loc_xyz = np.transpose((np.full(len(x_sat), x_beam), np.full(len(y_sat), y_beam), np.full(len(z_sat), z_beam)))

        cos_ang = x_beam*x_sat + y_beam*y_sat + z_beam*z_sat
        ang = np.arccos(cos_ang)
        print(np.shape(ang))
        return ang

    @staticmethod
    def convert_P(C_N, nothing=True):
        """
        This function converts the receiver's C/N_0 measurement into a classical power measurement in dBw
        Convert C/N_0 to power in dBw.
        """
        if nothing:
            return C_N
        N_sys = 8.5  # dB
        Tant = 30  # K
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
        beam_map.hist_single_sat()
        # print(beam_map.con0vert_GPStime_wrapper([423450.000],[2214]))
        # print(beam_map.convert_GPStime_wrapper([509250.000],[2214]))
        beam_map.make_plot(True, "All Satellites Visible by D3A(6) Over 24 Hours")

    # directories = ["/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/D3A(6)_Data/June14_port_C2_GNSS/", "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/D3A(6)_Data/June_16_port_C2_GNSS/", "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/D3A(6)_Data/June_16_port_C2_GNSS_part2/"]
    # beam_map = GetBeamMap(data_direcs=directories)
