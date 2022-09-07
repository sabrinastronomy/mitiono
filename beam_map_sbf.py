"""
This file extracts data from SBF SatVisibility files returned by the receiver to generate a beam map.

Written by Sabrina Berger and Vincent MacKay
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
from CSTUtils import *


D3A_ALT_deg = 80.5 # This is the elevation of the dish from horizon
D3A_AZ_deg = 0 # old
D3A_ALT = D3A_ALT_deg * np.pi/180
D3A_AZ = D3A_AZ_deg * np.pi/180

class GetBeamMap(ExtractSBF):
    """
    This class takes in data directories and creates a beam map.
    """
    def __init__(self, data_direcs):
        super().__init__(data_direcs=data_direcs, include_elev=True, process=True)
        self.close_sat_beams = []

    def print_dictionary(self):
        print(self.all_sat_dict)
        return

    def replace_dictionary(self, dict):
        self.all_sat_dict = dict

    def make_plot(self, all_sat, plot_title, sat_list=[], shift=True, direc="/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/plots/", filename="all_sat.png"):
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
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for key in sats_to_plot:
            if key is not None:
                satellite = self.all_sat_dict[key]
                times_cno = satellite["times"][0]
                times_elevs = satellite["elev_time"][0]
                elevs = satellite["elevations"][0]
                az = satellite["azimuths"][0]
                if shift:  # shift > 180 to negative values
                    az = np.asarray([-(360 - x) if x > 180 else x for x in az])

                cnos = satellite["cno"][0]
                times_elevs, cnos, elevs, az = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)

                power = self.convert_P(cnos)  # dBw

                elev_beam = np.full(len(elevs), D3A_ALT_deg)
                diff_elev = np.abs(elevs - elev_beam)
                az_beam = np.full(len(az), D3A_AZ_deg)
                diff_az = np.abs(az - az_beam)

                # Checking if current satellite passes close to center of beam
                if (diff_az < 5).any() and (diff_elev < 5).any():
                    print(f"{key} passes close to center of beam.")
                    self.close_sat_beams.append(key)
                else: # not including satellites that do NOT pass close to beam center
                    continue

                if len(power) < 1:
                    continue
                max_power_current = max(power)
                if max_power_current > max_power:
                    max_power = max_power_current

                min_power_current = min(power)
                if min_power_current < min_power:
                    min_power = min_power_current
                power -= 63.0625 # max power
                cmap_reversed = plt.cm.get_cmap('viridis_r')
                scat = ax.scatter(az, elevs, s=0.1, c=power, cmap=cmap_reversed)
        print(max_power)
        fig.subplots_adjust(right=0.8)

        fig.colorbar(scat)
        circle1 = plt.Circle((D3A_AZ_deg, D3A_ALT_deg), radius=1.25, color='k', fill=None)
        ax.add_patch(circle1)

        circle1 = plt.Circle((21.289999, 81.349998), radius=4, color='r', fill=None)
        ax.add_patch(circle1)
        scat.set_clim(-10, -30)
        ax.scatter(None, None, c="k", label="D3A(6) Beam")
        ax.scatter(None, None, c="r", label="weird glitch")
        ax.legend()

        ax.set_xlabel("Azimuth [deg]")
        ax.set_ylabel("Elevation [deg]")
        # ax.set_xlim(0, 50)
        # ax.set_ylim(75, 85)

        ax.set_title(plot_title)
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

    def compare_sats(self, list_sats):
        for sat in list_sats:
            satellite = self.all_sat_dict[sat]
            elevs = satellite["elevations"][0]
            az = satellite["azimuths"][0]
            cnos = satellite["cno"][0]
            times_cno = satellite["times"][0]
            times_elevs = satellite["elev_time"][0]
            times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
            plt.scatter(times_elevs, cnos, label=sat, s=0.1)
            plt.title(f"L{self.mask_frequency}")
        # plt.xlim(490000, 495000)
        plt.legend()
        plt.savefig(f"../plots/compare_sats_{self.mask_frequency}.png", dpi=400)

    def hist_single_sat(self, sat_name, i_pol, i_freq, dB_offset_sim=0, dB_offset_data=0, zenith_rot=0, ns_rot=0, ew_rot=0,
                        norm_power=True, shift=True):
        """
        1D profile where y = power, x = distance from pointing center
        """

        satellite = self.all_sat_dict[sat_name]
        elevs = satellite["elevations"][0]
        az = satellite["azimuths"][0]
        cnos = satellite["cno"][0]
        times_cno = satellite["times"][0]
        times_elevs = satellite["elev_time"][0]
        times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)

        if shift: # shift > 180 to negative values
            print(az_deg)
            az_deg = np.asarray([-(360 - x) if x > 180 else x for x in az_deg])
            print(az_deg)

        elevs = elevs_deg * np.pi / 180
        az = az_deg * np.pi / 180

        # masking
        mask_az = az_deg < D3A_AZ_deg  # left
        mask = mask_az
        left_az = az[mask]
        left_elevs = elevs[mask]
        left_cnos = cnos[mask]
        left_angles_rad = self.convert_angular_distance_from_center_beam(left_elevs, left_az)
        left_angles_deg = left_angles_rad * 180 / np.pi
        left_angles_deg = -left_angles_deg
        left_powers = self.convert_P(left_cnos)
        mask_az = az_deg > D3A_AZ_deg  # right
        mask = mask_az
        right_az = az[mask]
        right_elevs = elevs[mask]
        right_cnos = cnos[mask]
        right_angles_rad = self.convert_angular_distance_from_center_beam(right_elevs, right_az)
        right_angles_deg = right_angles_rad * 180 / np.pi
        right_powers = self.convert_P(right_cnos)

        # hot check for min
        print(right_angles_deg[np.argmin(right_powers)])
        print(180 * right_elevs[np.argmin(right_powers)] / np.pi)
        print(180 * right_az[np.argmin(right_powers)] / np.pi)

        # The beam is currently defined in theta and phi,
        # with the beam's boresight at beam.theta = 0 degrees,
        # and its horizon is at beam.theta = 90 degrees

        # There are three issues to fix:
        # (1) The satellite positions are defined in azimuth and elevation (azel),
        #     so we would like to either re-express beam's theta/phi coordinates in
        #     terms of azel, or re-express the satellite positions in terms of
        #     theta and phi.
        # (2) The beam boresight is at an elevation of D3A_ALT_deg = 81.41,
        #     which, in phi/theta coordinates, is 90-81.41 = 8.59 deg.
        #     We would like to rotate the beam so that it is pointed at 81.41
        #     degrees of elevation, or at theta = 8.59 deg. Note that since
        #     the beam's resolution is 0.5 degrees, this will need to be rounded
        #     to theta = 8.5 deg.
        # (3) The beam may need to be rotated around its boresight (likely 45, 135,
        #     225, or 315 degrees, to account for the feed's orientation).

        # We can readily fix issues (2) and (3) when calling the beam, with
        # the optional variables "zenith_rot," "ns_rot," and "ew_rot," e.g.
        # beam = CSTBeam(beams_folder,zenith_rot = 45, ns_rot = 9.5)
        # The variable zenith_rot is the rotation around boresight,
        # which is performed first. In other words, before tilting the dish
        # away from zenith, we look at it from above, and rotate it by
        # the amount given by zenith_rot (in degrees), to account for feed rotation.
        # The variables ns_rot and ew_rot are the rotation around a tilting axis.
        # In other words, we tilt the dish away from zenith. A positive
        # ns_rot tilts it towards the south, and negative, towards
        # the north. Similarly for ew_rot, although D3A/6 dishes don't tilt
        # in that direction.

        # To load the beams, you only need to specify the folder.
        # The zenith_rot, ns_rot, and ew_rot variables are optional.
        # In our case, I added those variables to the hist_single_sat() function,
        # so that we can easily change the orientation of the beam when
        # calling that function.
        beams_folder = '/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/notebook_and_beams/new_beams/'
        beam = CSTBeam(beams_folder, zenith_rot=zenith_rot, ns_rot=ns_rot, ew_rot=ew_rot)

        # Now, let's solve problem (1).
        # We would have two options: converting the beam to azel,
        # or converting the satellite to theta/phi.
        # I chose to try the latter, because once we got the theta/phi
        # coordinates, it's going to be easy to retrieve the corresponding
        # beam indices.

        # Satellite positions in theta/phi:
        # Thankfully, we're using the "physics" theta/phi convention (not the
        # "maths" one), so the conversion is (almost) trivial.
        left_sat_theta = np.pi / 2 - left_elevs
        left_sat_phi = left_az
        right_sat_theta = np.pi / 2 - right_elevs
        right_sat_phi = right_az

        # Convert to degrees, and don't forget to
        # round to the nearest beam.theta_step & beam.phi_step
        left_sat_theta_deg = beam.theta_step * np.round((180 / np.pi) * left_sat_theta / beam.theta_step)
        left_sat_phi_deg = beam.phi_step * np.round((180 / np.pi) * left_sat_phi / beam.phi_step)
        right_sat_theta_deg = beam.theta_step * np.round((180 / np.pi) * right_sat_theta / beam.theta_step)
        right_sat_phi_deg = beam.phi_step * np.round((180 / np.pi) * right_sat_phi / beam.phi_step)

        # Now we want to find what points on the beam correspond to those
        # coordinates we just found (sat_theta_deg and sat_phi_deg).
        # The way I do it is to use beam.phi_step and beam.theta_step.
        # Say a point is at theta = 9, the corresponding index in the
        # theta axis will be theta / beam.theta_step = 9 / 0.5 = 18.
        i_sat_theta_left = (left_sat_theta_deg // beam.theta_step).astype('int')
        i_sat_phi_left = (left_sat_phi_deg // beam.phi_step).astype('int')
        i_sat_left = (i_sat_phi_left, i_sat_theta_left)
        i_sat_theta_right = (right_sat_theta_deg // beam.theta_step).astype('int')
        i_sat_phi_right = (right_sat_phi_deg // beam.phi_step).astype('int')
        i_sat_right = (i_sat_phi_right, i_sat_theta_right)
        # Then, the corresponding beam cut will be beam.beams_dir[i_pol,i_freq,i_sat_phi,i_sat_theta,:]

        # Plotting time!
        fig, ax = plt.subplots(1, 1, figsize=[10, 8])

        # Plot the satellite
        max_sat_power = max(left_powers.max(), right_powers.max())
        # I added a variable called norm_power, that is at True by default.
        # If True, the norm_power variable makes it so that the plotted power
        # (both the data and simulated beam) is max subtracted in dB
        # (i.e. normalized by max, in linear scale).
        # If False, then the dB_offset_data and dB_offset_beam variables can be
        # used to adjust the height of either curve.
        if norm_power:
            ax.scatter(left_angles_deg, left_powers - max_sat_power, c="k", s=0.5)
            ax.scatter(right_angles_deg, right_powers - max_sat_power, c="k", s=0.5,
                       label="GNSS Beam Measurement at ~1.175 GHz")
        else:
            ax.scatter(left_angles_deg, left_powers - 59.875 + dB_offset_data, c="k", s=0.5)
            ax.scatter(right_angles_deg, right_powers - 59.875 + dB_offset_data, c="k", s=0.5,
                       label="GNSS Beam Measurement at ~1.175 GHz")

        # Plot the simulated beam
        line_props = {
            'linestyle': '-',
            'color': 'r',
            'linewidth': 2,
            'alpha': 1,
            'color': 'b'
        }
        # max_beam_power = 10 * np.log10(
        #     max(beam.beams_dir[i_pol, i_freq][i_sat_left].max(), beam.beams_dir[i_pol, i_freq][i_sat_right].max()))
        # if norm_power:
        #     ax.plot(left_angles_deg, 10 * np.log10(beam.beams_dir[i_pol, i_freq][i_sat_left]) - max_beam_power,
        #             **line_props, label='Simulated beam at {:.3g} GHz'.format(beam.freqs[i_freq]))
        #     ax.plot(right_angles_deg, 10 * np.log10(beam.beams_dir[i_pol, i_freq][i_sat_right]) - max_beam_power,
        #             **line_props)
        #
        # else:
        #     ax.plot(left_angles_deg, 10 * np.log10(beam.beams_dir[i_pol, i_freq][i_sat_left]) + dB_offset_sim,
        #             **line_props, label='Simulated beam at {:.3g} GHz'.format(beam.freqs[i_freq]))
        #     ax.plot(right_angles_deg, 10 * np.log10(beam.beams_dir[i_pol, i_freq][i_sat_right]) + dB_offset_sim,
        #             **line_props)

        ax.set_xlim([-20, 20])

        ax.set_xlabel(r"$\theta$, the angle from beam center [deg]")
        ax.set_ylabel("Log Power [uncalibrated]")
        ax.set_title(f"Galileo Satellte: {sat_name}")
        ax.legend(loc='lower center')
        print("Saving single satellite plot at " + "/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/plots/sim_beam_{:.3g}.png".format(beam.freqs[i_freq]))
        fig.savefig(f"../plots/sim_beam_{beam.freqs[i_freq]}_{sat_name}.png", dpi=300, bbox_inches='tight')

    @staticmethod
    def convert_az_el_to_beam_theta_phi(az, el, pitch_offset, roll_offset):
        """
        This function converts the azimuth and elevation coordinates of the
        satellite's position to theta and phi coordinates centered around
        the beam boresight.

        The function has three steps.
        First, it converts to cartesian coordinates, where doing a rotation
        is easier.
        Second, it performs the rotation.
        Third, it converts to theta/phi.
        """

        # First, convert (az,el) to cartesian, just copy-pasting
        # the convert_angular_distance_from_center_beam function.
        x_sat = np.cos(az) * np.cos(el)
        y_sat = np.sin(az) * np.cos(el)
        z_sat = np.sin(el)

        # Then, rotate those coordinates so that they are aligned with
        # the beam.
        # It is sufficient to just rotate around two axes in this context,
        # so I do z-axis, followed by y-axis.
        theta_z = 0 + roll_offset * (np.pi / 180)
        theta_y = D3A_ALT + pitch_offset * (np.pi / 180)
        Rz = np.matrix([[np.cos(theta_z), -np.sin(theta_z), 0],
                        [np.sin(theta_z), np.cos(theta_z), 0],
                        [0, 0, 1]])
        Ry = np.matrix([[np.cos(theta_y), 0, np.sin(theta_y)],
                        [0, 1, 0],
                        [-np.sin(theta_y), 0, np.cos(theta_y)]])
        # cart_rot = np.matmul(Ry,np.matmul(Rz,np.array([x_sat,y_sat,z_sat])))
        cart_rot = np.array(Ry.dot(Rz.dot(np.array([x_sat, y_sat, z_sat]))))

        # Convert back to theta/phi
        theta_sat = np.arctan2((cart_rot[0] ** 2. + cart_rot[1] ** 2.) ** 0.5, cart_rot[2])
        phi_sat = np.arctan2(cart_rot[1], cart_rot[0])

        return theta_sat, phi_sat

    @staticmethod
    def convert_angular_distance_from_center_beam(alt, az):

        x_beam = np.cos(D3A_AZ) * np.cos(D3A_ALT)
        y_beam = np.sin(D3A_AZ) * np.cos(D3A_ALT)
        z_beam = np.sin(D3A_ALT)

        x_sat = np.cos(az) * np.cos(alt)
        y_sat = np.sin(az) * np.cos(alt)
        z_sat = np.sin(alt)

        cos_ang = x_beam*x_sat + y_beam*y_sat + z_beam*z_sat
        ang = np.arccos(cos_ang)
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
    direc = "/Volumes/Photos Mac/all_one_drive_MCGILL/D3A(6)_Data_Long_2022/"
    days = ["August_25_port_2B_part1_Attenuator_GNSS_different_settings"]
    new_days = []
    for i, day in enumerate(days):
        new_days.append(direc + day + "/")
    # beam_map = GetBeamMap(data_direcs=new_days)


    directories = ["/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/parsed_data/"]
    # days = ["June14_port_C2_GNSS_satellite_dict_all", "June_16_port_C2_GNSS_satellite_dict_all", "June_16_port_C2_GNSS_part2_satellite_dict_all"]
    # days = ["June_16_port_C2_GNSS_part2_satellite_dict_all"]

    ## Sample code to grab parsed dictionary file ######
    for day in days:
        sat_dict = np.load(directories[0] + day + "_satellite_dict_all" + ".npy", allow_pickle=True).item()  # extracting saved dictionary

        beam_map = GetBeamMap(data_direcs=[])
        beam_map.replace_dictionary(sat_dict)
        beam_map.make_plot(all_sat=True, plot_title="24 hours of satellites tracks close to the center of the beam", filename=f"{day}_all_sat.png")
        i_freq = 8
        dB_offset_sim = 0
        dB_offset_data = 35
        zenith_rot = 135
        ns_rot = 8.5
        ew_rot = 0
        i_pol = 1
        x = np.arange(0, 10, 1)
        y = np.arange(85, 95, 1)
        xv, yv = np.meshgrid(x, y, indexing='ij')
        # for i in range(10):
        #     for j in range(10):
        #         which_5_deg = beam_map.convert_angular_distance_from_center_beam(xv[i][j], yv[i][j])
        #         if (np.abs(which_5_deg - 3.50) < 2).any():
        #             print(which_5_deg)
        #             print(xv[i][j], yv[i][j])
        # for sat in beam_map.close_sat_beams:
        for sat in ["C36"]:
            beam_map.hist_single_sat(sat, i_pol,
                                     i_freq,
                                     norm_power=True,
                                     dB_offset_data=dB_offset_data,
                                     zenith_rot=zenith_rot,
                                     ns_rot=ns_rot,
                                     ew_rot=ew_rot)


