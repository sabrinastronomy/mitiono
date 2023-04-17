"""
This file extracts data from SBF SatVisibility files returned by the receiver to generate a beam map.

Written by Sabrina Berger and Vincent MacKay
"""
import matplotlib.pyplot as plt
import numpy as np
from extract_sbf_data import ExtractSBF
import matplotlib
from scipy import special
from CSTUtils import *
from math import radians

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# June, 2022
# D3A_ALT_deg = 81.41
# D3A_AZ_deg = 180
#
# D3A_ALT = 81.41 * np.pi/180
# D3A_AZ = 180 * np.pi/180

# Post August, 2022
D3A_ALT_deg = 80.5 # This is the elevation of the dish from horizon
D3A_AZ_deg = 0 # new
D3A_ALT = D3A_ALT_deg * np.pi/180
D3A_AZ = D3A_AZ_deg * np.pi/180

class GetBeamMap(ExtractSBF):
    """
    This class takes in data directories and creates a beam map.
    """
    def __init__(self, data_direcs, mask_frequency="1", save_parsed_data_direc="parsed_data", plot_dir="/Users/sabrinaberger/Desktop/ursi/", masking=False):
        super().__init__(data_direcs=data_direcs, include_elev=True, process=True, mask_frequency=mask_frequency, save_parsed_data_direc=save_parsed_data_direc, masking=masking)
        self.plot_dir = plot_dir
        self.panel_plot_num = 0 # number of panel plot so we can plot multiple satellites
        if self.all_sat_dict != None:
            self.get_close_sats()

    def print_dictionary(self):
        print(self.all_sat_dict)
        return

    def replace_dictionary(self, dict):
        self.all_sat_dict = dict
        self.get_close_sats()

    def get_close_sats(self, tol_beam=10):
        self.close_sat_beams = []

        min_diff_az = 100
        min_diff_el = 100
        min_nam = ""
        for key in self.all_sat_dict:
            if key is not None:
                satellite = self.all_sat_dict[key]
                times_cno = satellite["times"][0]
                times_elevs = satellite["elev_time"][0]
                cnos = satellite["cno"][0]
                elevs = satellite["elevations"][0]
                az = satellite["azimuths"][0]

                times_elevs, cnos, elevs, az = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
                elev_beam = np.full(len(elevs), D3A_ALT_deg)
                diff_elev = np.abs(elevs - elev_beam)
                az_beam = np.full(len(az), D3A_AZ_deg)
                diff_az = np.abs(az - az_beam)

                # if min_diff_az > min(diff_az) and min_diff_el > min(diff_elev):
                #     min_diff_az = min(diff_az)
                #     min_diff_el = min(diff_elev)
                #     min_nam = key
                # Checking if current satellite passes close to center of beam
                if (diff_az < tol_beam).any() and (diff_elev < tol_beam).any():
                    print(f"{key} passes close to center of beam.")
                    min_diff_az = min(diff_az)
                    min_diff_el = min(diff_elev)
                    print(f"min diff az: {min_diff_az}")
                    print(f"min diff elev: {min_diff_el}")

                    self.close_sat_beams.append(key)
                else:  # not including satellites that do NOT pass close to beam center
                    continue
        print(f"Total number of near beam satellites is {len(self.close_sat_beams)}.")

    def get_min_max(self, arrs):
        arr_concatenated = np.concatenate(arrs, axis=0)
        min_arr, max_arr = arr_concatenated.min(), arr_concatenated.max()
        return min_arr, max_arr

    def make_plot(self, all_sat, plot_title, sat_list=[], shift=True, show_power=False, filename="all_sat.png"):
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
        min_time = 1e10
        max_time = 0

        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        powers_arr = []

        for i, key in enumerate(sats_to_plot):
            satellite = self.all_sat_dict[key]
            cnos = satellite["cno"][0]
            times_cno = satellite["times"][0]
            times_elevs = satellite["elev_time"][0]
            elevs = satellite["elevations"][0]
            az = satellite["azimuths"][0]

            _, cnos, _, _ = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
            powers_arr.append(cnos)
        min_power, max_power = self.get_min_max(powers_arr)

        for i, key in enumerate(sats_to_plot):

            if key is not None:
                satellite = self.all_sat_dict[key]
                times_cno = satellite["times"][0]
                times_elevs = satellite["elev_time"][0]
                wnc = satellite["wnc"][0]

                elevs = satellite["elevations"][0]
                x = ExtractSBF.weeksecondstoutc(gpsweek=wnc[0], gpsseconds=times_cno[0])
                # getting time over which observations took place
                if min_time > min(times_elevs):
                    min_time = min(times_elevs)
                if max_time < max(times_elevs):
                    max_time = max(times_elevs)

                az = satellite["azimuths"][0]
                if shift:  # shift > 180 to negative values
                    print("Shifted azimuths...")
                    az = self.shift_az(az)

                cnos = satellite["cno"][0]
                times_elevs, cnos, elevs, az = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
                power = self.convert_P(cnos)
                powers_arr.append(power)

                if len(power) < 1:
                    continue
                if show_power:
                    cmap_reversed = plt.cm.get_cmap('viridis_r')
                    scat = ax.scatter(np.deg2rad(az), 90-elevs, s=5, c=power, cmap=cmap_reversed, vmin=min_power, vmax=max_power)
                else:
                    ax.scatter(np.deg2rad(az), 90-elevs, s=0.01, c="k")
                    print(np.deg2rad(az[:2]), elevs[:2])

        # fig.subplots_adjust(right=0.8)
        # N = 100
        # theta = np.random.rand(N) * np.pi * 2
        # r = np.cos(theta * 2) + np.random.randn(N) * 0.1
        # ax.scatter(theta, r)
        s = 100*1.22 * 0.2 / 6
        circle = plt.Circle((np.deg2rad(D3A_AZ), 90-D3A_ALT_deg), s, fill=False, transform=ax.transData._b, color="cyan")
        ax.add_artist(circle)
        # cax = plt.axes([0.85, 0.1, 0.075, 0.8])
        # plt.colorbar(cax=cax)
        # # print(f"beginning of observation time (GPS seconds): {min_time}")
        # print(f"beginning of observation time (GPS seconds): {min_time}")
        # print(f"end of observation time (GPS seconds): {max_time}")

        # plt.scatter(None, None, c="green",  label="Center of D3A Beam")
        # ax.set_xlabel(r"${\rm Azimuth ~ [deg]}$")
        # ax.set_ylabel(r"${\rm Elevation ~ [deg]}$")
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_yticks(range(0, 90 + 10, 10))  # Define the yticks
        yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
        ax.set_yticklabels(yLabel)
        ax.set_title(plot_title)
        if show_power:
            fig.colorbar(scat)
        # plt.legend()
        plt.savefig(self.plot_dir + filename, bbox_inches='tight', dpi=300)
        plt.close()

    def shift_az(self, az_deg):
        az_deg = np.asarray([-(360 - x) if x > 180 else x for x in az_deg])
        return az_deg

    def mask_array_for_one_d(self, az, elevs, cnos, mask, side):
        assert(side == "right" or side == "left")
        az = az[mask]
        elevs = elevs[mask]
        cnos = cnos[mask]
        angles_rad = self.convert_angular_distance_from_center_beam(elevs, az)
        angles_deg = angles_rad * 180 / np.pi
        if "left":
            angles_deg = -angles_deg # getting sign right for left
        powers = self.convert_P(cnos)
        return angles_deg, powers

    def save_particular_sat(self, sat_name):
        satellite = self.all_sat_dict[sat_name]
        elevs = satellite["elevations"][0]
        az = satellite["azimuths"][0]
        cnos = satellite["cno"][0]
        times_cno = satellite["times"][0]
        times_elevs = satellite["elev_time"][0]
        times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)

        np.save(f"{sat_name}_times.npy", times_elevs)
        np.save(f"{sat_name}_cnos.npy", cnos)
        np.save(f"{sat_name}_elev.npy", elevs_deg)
        np.save(f"{sat_name}_az.npy", az_deg)

    def get_angles_for_one_d_prof(self, sat_name, shift=True, select_pass=False):
        """
        Gives you a 1D profile for a satellite as a function of angular distance from beam
        :param sat_name - name of satellite
        :param shift - whether or not to shift azimuths with shift_az
        :return angles_deg, powers - theta and power of a one d beam profile for a given satellite
        """
        satellite = self.all_sat_dict[sat_name]
        elevs = satellite["elevations"][0]
        az = satellite["azimuths"][0]
        cnos = satellite["cno"][0]
        times_cno = satellite["times"][0]
        times_elevs = satellite["elev_time"][0]
        times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
        i = 0
        if select_pass:
            diff_times = np.diff(times_elevs)
            indices = np.argwhere(diff_times > 10000).flatten()
            print(indices[i])
            elevs = elevs_deg[:indices[i]] * np.pi / 180
            az = az_deg[:indices[i]] * np.pi / 180
            az_deg = az_deg[:indices[i]]
            cnos = cnos[:indices[i]]
        else:
            elevs = elevs_deg * np.pi / 180
            az = az_deg * np.pi / 180
        if shift: # shift > 180 to negative values
            az_deg = self.shift_az(az_deg)

        # masking
        right_mask = az_deg > D3A_AZ_deg  # right
        left_mask = az_deg < D3A_AZ_deg  # left
        # this takes in radians
        right_angles_deg, right_powers = self.mask_array_for_one_d(az, elevs, cnos, right_mask, "right")
        left_angles_deg, left_powers = self.mask_array_for_one_d(az, elevs, cnos, left_mask, "left")

        # powers = np.concatenate((left_powers, right_powers))
        # angles_deg = np.concatenate((left_angles_deg, right_angles_deg))
        return left_angles_deg, left_powers, right_angles_deg, right_powers

    def plot_panel_sats(self, rows=6, cols=3, start=0, end=19, chosen=True, figsize=(8, 12), offset=True, want_theta=False, want_time=False):
        fig, axes = plt.subplots(rows, cols, figsize=figsize)
        var = 0
        print("num close_sats", len(self.close_sat_beams))
        if chosen:
            close_sat_beams = self.chosen_list
        else:
            close_sat_beams = self.close_sat_beams[start:end] # local to function version, NO SELF
        print(len(close_sat_beams))
        for i in range(rows):
            for j in range(cols):
                sat_name = close_sat_beams[var]
                if want_theta:
                    left_angles_deg, left_powers, right_angles_deg, right_powers = self.get_angles_for_one_d_prof(
                        sat_name, shift=True)
                    axes[i][j].scatter(left_angles_deg, left_powers, s=0.5, label=sat_name, c="k")
                    axes[i][j].scatter(-right_angles_deg, right_powers, s=0.5, c="k")
                    axes[i][j].set_xlabel(r"${\theta \rm ~ [deg]}$")
                    axes[i][j].set_ylabel(r'${C/N_0 \rm ~ [dB-Hz]}$')
                    axes[i][j].set_title(rf"${sat_name}$")
                    axes[i][j].set_xlim((-90, 90))

                if want_time:
                    satellite = self.all_sat_dict[sat_name]
                    times_cno = satellite["times"][0]
                    times_elevs = satellite["elev_time"][0]
                    cnos = satellite["cno"][0]
                    elevs = satellite["elevations"][0]
                    az = satellite["azimuths"][0]
                    times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
                    diff_times = np.diff(times_elevs)
                    indices = np.argwhere(diff_times > 10000).flatten()
                    past_index = 0
                    k = 0
                    while k < len(indices):
                        next_index = indices[k]
                        times = times_elevs[past_index:next_index]
                        time_print = ExtractSBF.weeksecondstoutc(2225, times_elevs[next_index])
                        print(time_print)
                        times = times - times_elevs[past_index + 1]
                        if len(times) < 250:
                            past_index = indices[k]
                            k += 1
                            continue
                        if offset:
                            axes[i][j].scatter(times, cnos[past_index:next_index] + 10 * k, c="k", s=0.05)
                        else:
                            axes[i][j].scatter(times, cnos[past_index:next_index], c="k", s=0.05)

                        past_index = indices[k]
                        times_max = max(times) # need to take max of times here to ensure it got through the if loop
                        k += 1
                    axes[i][j].set_xlabel(r'$Time$')
                    axes[i][j].set_ylabel(r'${C/N_0 \rm ~ [dB-Hz]}$')
                    axes[i][j].set_title(rf"${sat_name}$")
                    axes[i][j].set_xlim((0, times_max))
                var += 1
        plt.tight_layout()
        fig.savefig(self.plot_dir + f"panel_{self.panel_plot_num}.png", dpi=900)
        self.panel_plot_num += 1


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
        plt.legend()
        plt.savefig(f"../plots/compare_sats_{self.mask_frequency}.png", dpi=400)

    def fix_hist_single_sat(self, sat_name, i_pol, i_freq, dB_offset_sim=0, dB_offset_data=0, zenith_rot=0, ns_rot=0, ew_rot=0,
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
            az_deg = np.asarray([-(360 - x) if x > 180 else x for x in az_deg])

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

        beams_folder = '/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/Current Research/mitiono/notebook_and_beams/new_new_beams/'
        new_beam = CSTBeam(beams_folder)
        beam = new_beam.rotate(zenith_rot=D3A_ALT_deg, ns_rot=90-D3A_AZ_deg)

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

        # # Convert to degrees, and don't forget to
        # # round to the nearest beam.theta_step & beam.phi_step
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
        # Then, the corresponding beam cut will be beam.directivity[i_pol,i_freq,i_sat_phi,i_sat_theta,:]

        # # Plotting time!
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
            ax.scatter(left_angles_deg, left_powers - max_sat_power, c="k", s=0.1)
            ax.scatter(right_angles_deg, right_powers - max_sat_power, c="k", s=0.1,)
        else:
            ax.scatter(left_angles_deg, left_powers - 59.875 + dB_offset_data, c="k", s=0.1)
            ax.scatter(right_angles_deg, right_powers - 59.875 + dB_offset_data, c="k", s=0.1)
        ax.scatter(left_angles_deg, left_powers, c="k", s=0.1)
        ax.scatter(right_angles_deg, right_powers, c="k", s=0.1)

        # Plot the simulated beam
        line_props = {
            'linestyle': '-',
            'color': 'r',
            'linewidth': 2,
            'alpha': 1,
            'color': 'b'
        }
        max_beam_power = 10 * np.log10(
            max(beam.directivity[i_pol, i_freq][i_sat_left].max(), beam.directivity[i_pol, i_freq][i_sat_right].max()))
        if norm_power:
            ax.plot(left_angles_deg, 10 * np.log10(beam.directivity[i_pol, i_freq][i_sat_left]) - max_beam_power,
                    **line_props, label='Simulated beam at {:.3g} GHz'.format(beam.freqs[i_freq]))
            ax.plot(right_angles_deg, 10 * np.log10(beam.directivity[i_pol, i_freq][i_sat_right]) - max_beam_power,
                    **line_props)

        else:
            ax.plot(left_angles_deg, 10 * np.log10(beam.directivity[i_pol, i_freq][i_sat_left]) + dB_offset_sim,
                    **line_props, label='Simulated beam at {:.3g} GHz'.format(beam.freqs[i_freq]))
            ax.plot(right_angles_deg, 10 * np.log10(beam.directivity[i_pol, i_freq][i_sat_right]) + dB_offset_sim,
                    **line_props)

        ax.set_xlim([-20, 20])
        print(beam.freqs)
        ax.set_xlabel(r"$\theta$, the angle from beam center [deg]")
        ax.set_ylabel("Log Power [uncalibrated]")
        ax.set_title(f"Galileo Satellte: {sat_name}")
        ax.legend(loc='lower center')
        # print("Saving single satellite plot at " +
        #       f"/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/Current Research/mitiono/plots/" + f"sim_beam_{beam.freqs[i_freq]}_{sat_name}.png".format(beam.freqs[i_freq], sat_name))
        fig.savefig(f"test_sim_beam_{sat_name}.png", dpi=300, bbox_inches='tight')

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

    def overplot_days(self, sat_name):
        satellite = self.all_sat_dict[sat_name]
        times_cno = satellite["times"][0]
        times_elevs = satellite["elev_time"][0]
        cnos = satellite["cno"][0]
        elevs = satellite["elevations"][0]
        az = satellite["azimuths"][0]
        times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)

        # below I'm splicing the array when the difference between times in subsequent elements is greater than 10k
        # this splicing tolerance might have to be changed
        diff_times = np.diff(times_elevs)
        indices = np.argwhere(diff_times > 10000).flatten()
        print(indices)
        past_index = 0
        i = 0
        while i < len(indices):
            next_index = indices[i]
            print(next_index)
            times = times_elevs[past_index:next_index]
            times = times - times_elevs[past_index+1]
            if len(times) < 250:
                past_index = indices[i]
                i += 1
                continue
            plt.scatter(times, cnos[past_index:next_index] + 10*i, label=f"Day {i}", s=0.5)
            past_index = indices[i]
            i += 1

        plt.ylabel(r'${C/N_0 \rm ~ [dB-Hz]}$')
        plt.legend()
        plt.title(sat_name)
        plt.xlim(0, 25000)
        plt.savefig(self.plot_dir + "overplotted_offset.png")

    def plot_one_d_prof(self, sat_names):
        """
        Plots 1D beam profile in time and angle
        """
        fig, axes = plt.subplots()

        for sat_name in sat_names:
            left_angles_deg, left_powers, right_angles_deg, right_powers = self.get_angles_for_one_d_prof(sat_name, shift=True)
            axes.scatter(left_angles_deg, left_powers, s=0.5, label=sat_name, c="k")
            axes.scatter(-right_angles_deg, right_powers, s=0.5, c="k")


            # plt.close()
            # cnos = satellite["cno"][0]
            # satellite = self.all_sat_dict[sat_name]
            # times_cno = satellite["times"][0]
            # plt.scatter(times_cno, cnos, c="k", s=0.1)
            # plt.xlabel(r"${\rm Time ~ [s]}$")
            # plt.ylabel(r'${C/N_0 \rm ~ [dB-Hz]}$')
            # plt.title(sat_name)
            # plt.savefig(self.plot_dir + sat_name + "_time_indiv.png", bbox_inches='tight')
        plt.xlabel(r"${\theta \rm ~ [deg]}$")
        plt.ylabel(r'${C/N_0 \rm ~ [dB-Hz]}$')
        plt.legend()
        plt.title(sat_name[:3])
        plt.savefig(self.plot_dir + sat_name + "_power_indiv.png", bbox_inches='tight', dpi=300)

    @staticmethod
    def get_bessel_0(ka_sintheta):
        return special.j0(ka_sintheta)

    @staticmethod
    def get_bessel_1(ka_sintheta):
        return special.j1(ka_sintheta)

    @staticmethod
    def convert_angular_distance_from_center_beam(alt, az):
        """
        Finds the angular distance from the beam center for a given alt and az using the dot product (lazy trig).
        """
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
        """
        Septentrio returns cnos and postiions of different lengths. This function just finds their intersection.
        :returns intersected times_elevs, cnos, elevs, az
        """
        # TODO fix duplicates
        indices = np.intersect1d(times_elevs, times_cno, return_indices=True)
        elev_indices = indices[1]
        cnos_indices = indices[2]

        times_elevs = times_elevs[elev_indices]
        elevs = elevs[elev_indices]
        az = az[elev_indices]
        cnos = cnos[cnos_indices]
        return times_elevs, cnos, elevs, az



# initial data run
# days = ["June14_port_C2_GNSS_satellite_dict_all", "June_16_port_C2_GNSS_satellite_dict_all", "June_16_port_C2_GNSS_part2_satellite_dict_all"]

if __name__ == "__main__":
    ### example to feed in directories

    # data directories with ability to parse multiple folders and days (currently just one)
    data_directories = ["/Users/sabrinaberger/Desktop/"]
    parse_days = ["August_29_port2B/"]

    parsed_data_directory = "../parsed_data"
    mask_frequencies = ["L1"]
    parse = False
    ## Sample code to parse data files and select frequencies#####

    for freq in mask_frequencies:
        if parse:
            beam_map_1 = GetBeamMap(data_direcs=[data_directories[0] + parse_days[0]], masking=True, mask_frequency=freq, save_parsed_data_direc=parsed_data_directory)
    days = ["August_29_port2B_satellite_dict_all_L1.npy"]
    # x = np.load("parsed_data/August_29_port2B_satellite_dict_all_E5a.npy")
    ## Sample code to grab parsed dictionary file ######
    # days = ["August_29_port2B"]
    parsed_data_directory = "../parsed_data/"

    for day in days:
        sat_dict = np.load(parsed_data_directory + day , allow_pickle=True).item()  # extracting saved dictionary

        beam_map = GetBeamMap(data_direcs=[])
        beam_map.replace_dictionary(sat_dict)
        # beam_map.plot_one_d_prof(["G07"])
        # beam_map.plot_panel_sats(start=0, end=19)
        # beam_map.plot_panel_sats(start=19, end=38)
        # beam_map.plot_panel_sats(start=37, end=56)
        # beam_map.plot_panel_sats(start=len(beam_map.close_sat_beams)-19, end=len(beam_map.close_sat_beams))

        # good_sats = ["G07", "G30", "G10", "G29"]
        # beam_map.chosen_list = good_sats
        # beam_map.plot_panel_sats(rows=2, cols=2, chosen=True, want_theta=True, figsize=(12,8))
        beam_map.make_plot(show_power=False, shift=True, all_sat=True, plot_title=r"${\rm Approx. 72 ~ hours ~ of ~ satellites' ~ tracks ~ at ~ D3A(6)}$", filename=f"{day}_all_sat.png")
