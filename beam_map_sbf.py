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
from pyuvdata import UVBeam



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
    def __init__(self, data_direcs, mask_frequency="1", save_parsed_data_direc="parsed_data", plot_dir="/Users/sabrinaberger/Desktop/beam_paper_plots/", masking=False):
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
                times_cno = satellite["times"][0] # SECONDS
                times_elevs = satellite["elev_time"][0] # SECONDS
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

    def make_plot(self, all_sat, plot_title, sat_list=[], shift=True, show_power=False, filename="all_sat.png", zoom=True):
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
                    cmap = plt.cm.get_cmap('viridis')
                    scat = ax.scatter(np.deg2rad(az), 90-elevs, s=5, c=power, cmap=cmap, vmin=min_power, vmax=max_power)
                else:
                    ax.scatter(np.deg2rad(az), 90-elevs, s=0.01, c="k")
                    print(np.deg2rad(az[:2]), elevs[:2])

        # fig.subplots_adjust(right=0.8)
        # N = 100
        # theta = np.random.rand(N) * np.pi * 2
        # r = np.cos(theta * 2) + np.random.randn(N) * 0.1
        # ax.scatter(theta, r)
        # s = 100*1.22 * 0.2 / 6
        # circle = plt.Circle((0, np.deg2rad(90-D3A_ALT_deg)), s, fill=False, color="cyan")
        # ax.add_artist(circle)

        # Define the center of the circle
        r_center = 90 - 80.5  # altitude (radius)
        theta_center = 0  # azimuth (angle in radians, converted from degrees)

        # Convert polar coordinates (r_center, theta_center) to Cartesian coordinates
        x_center = r_center * np.cos(theta_center)
        y_center = r_center * np.sin(theta_center)

        # Define the radius of the circle
        circle_radius = 100*1.22 * 0.2 / 6

        # Create the circle in Cartesian coordinates
        # circle = plt.Circle((0, 1), circle_radius, fill=False, color='cyan')
        circle = plt.Circle((x_center, y_center), circle_radius, fill=False, color='black',
                            transform=ax.transData._b, linewidth=1)
        # Add the circle to the plot
        ax.add_artist(circle)
        # cax = plt.axes([0.85, 0.1, 0.075, 0.8])
        # plt.colorbar(cax=cax)
        # # print(f"beginning of observation time (GPS seconds): {min_time}")
        # print(f"beginning of observation time (GPS seconds): {min_time}")
        # print(f"end of observation time (GPS seconds): {max_time}")

        # plt.scatter(None, None, c="green",  label="Center of D3A Beam")
        # ax.set_xlabel(r"${\rm Azimuth ~ [deg]}$")
        # ax.set_ylabel(r"${\rm Elevation ~ [deg]}$")



        ax.set_theta_zero_location('E')
        # ax.set_theta_direction(-1)

        ax.set_title(plot_title)
        if show_power:
            colorbar = fig.colorbar(scat)
            colorbar.set_label(r'$C/N_0$ [db-Hz]', rotation=270, labelpad=15)

        # plt.legend()
        # if zoom:
        #     ax.set_thetalim(0, np.pi)
        circle_radius = 3 * 100 * 1.22 * 0.2 / 6
        # convert circle_radius to polar plot limits
        zoom_radius = circle_radius  # Add some margin around the circle
        ax.set_yticks(range(0, 90 + 10, 10))  # Define the yticks
        # yLabel = ['90', '80', '70', '60', '50', '', '30', '', '', '']
        yLabel = ['90', '80', '70', '', '', '', '', '', '', '']
        ax.set_ylim(90 - (D3A_ALT_deg + zoom_radius), 90 - (D3A_ALT_deg - zoom_radius))


        ax.set_yticklabels(yLabel)
        plt.savefig(self.plot_dir + filename, bbox_inches='tight', dpi=300)
        plt.close()

    def shift_az(self, az_deg):
        az_deg = np.asarray([-(360 - x) if x > 180 else x for x in az_deg])
        return az_deg

    def mask_array_for_one_d(self, az, elevs, cnos, mask, side, times):
        assert(side == "right" or side == "left")
        az = az[mask]
        elevs = elevs[mask]
        cnos = cnos[mask]
        times = times[mask]
        angles_rad = self.convert_angular_distance_from_center_beam(elevs, az)
        angles_deg = angles_rad * 180 / np.pi
        if "left":
            angles_deg = -angles_deg # getting sign right for left
        powers = self.convert_P(cnos)
        return angles_deg, powers, times

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
        right_angles_deg, right_powers, right_times = self.mask_array_for_one_d(az, elevs, cnos, right_mask, "right", times_elevs)
        left_angles_deg, left_powers, left_times = self.mask_array_for_one_d(az, elevs, cnos, left_mask, "left", times_elevs)

        plt.close()
        plt.scatter(right_times, right_angles_deg)
        plt.show()
        plt.close()

        # powers = np.concatenate((left_powers, right_powers))
        # angles_deg = np.concatenate((left_angles_deg, right_angles_deg))
        return left_angles_deg, left_powers, right_angles_deg, right_powers

    def convert_C_n_to_SNR_linear(self, c_n, bw=2e6):
        # assuming nominal L1 bandwidth of 2Mhz
        c_n = np.asarray(c_n)
        bw = np.log10(bw)
        s_n_db = c_n - bw
        s_n_linear = 10**(s_n_db/10)
        return s_n_linear

    def convert_SNR_linear_C_n(self, SNR_linear, bw=2e6):
        # assuming nominal L1 bandwidth of 2Mhz
        s_n_db = 10 * np.log10(SNR_linear)
        bw = np.log10(bw)
        c_n = s_n_db + bw
        return c_n

    def average_C_N_0(self, c_n, angles, bw=2e6, linear_conversion=False):
        unique_angles_deg = np.unique(angles)
        if linear_conversion:
            s_n_linear = self.convert_C_n_to_SNR_linear(c_n)
            s_n_linear_average = np.array([np.mean(s_n_linear[angles == value]) for value in unique_angles_deg])
            s_n_linear_rms = np.array([np.std(s_n_linear[angles == value], ddof=0) for value in unique_angles_deg])
            c_n_average = self.convert_SNR_linear_C_n(s_n_linear_average)
            s_n_linear_rms_mask = s_n_linear_rms == 0
            s_n_linear_rms[s_n_linear_rms_mask] = 1
            c_n_rms = 10 * np.log10(s_n_linear_rms)
        else:
            s_n_linear = self.convert_C_n_to_SNR_linear(c_n)
            s_n_linear_average = np.array([np.mean(s_n_linear[angles == value]) for value in unique_angles_deg])
            c_n_average = self.convert_SNR_linear_C_n(s_n_linear_average)
            c_n_rms = np.array([np.std(c_n[angles == value], ddof=0) for value in unique_angles_deg])

        return unique_angles_deg, c_n_average, c_n_rms


    def plot_panel_sats(self, rows=6, cols=3, start=0, end=19, chosen=True, figsize=(8, 12), offset=True, want_theta=False, want_time=False, rms=False):
        fig, axes = plt.subplots(rows, cols, figsize=figsize)
        var = 0
        print("num close_sats", len(self.close_sat_beams))
        if chosen:
            close_sat_beams = self.chosen_list
        else:
            close_sat_beams = self.close_sat_beams[start:end] # local to function version, NO SELF

        sat_names = []
        peak_sigmas = []
        peak_powers = []

        for i in range(rows):
            for j in range(cols):
                print("row", i, "col", j)
                sat_name = close_sat_beams[var]
                if want_theta:
                    left_angles_deg, left_powers, right_angles_deg, right_powers = self.get_angles_for_one_d_prof(
                        sat_name, shift=True)
                    unique_angles_deg = np.unique(left_angles_deg)
                    counts_left = np.array([np.sum(left_angles_deg == value) for value in unique_angles_deg])

                    # get counts left = 1
                    counts_left_mask = counts_left < 2

                    unique_angles_deg = np.unique(right_angles_deg)
                    counts_right = np.array([np.sum(right_angles_deg == value) for value in unique_angles_deg])
                    counts_right_mask = counts_right < 2


                    print("average counts right")
                    print(np.average(counts_right))
                    unique_left_angles_deg, avg_left_powers, avg_left_sigma = self.average_C_N_0(left_powers, left_angles_deg)
                    unique_right_angles_deg, avg_right_powers, avg_right_sigma = self.average_C_N_0(right_powers, right_angles_deg)

                    print("average sigma value")
                    print(np.average(avg_right_sigma))
                    all_sigma = np.concatenate((avg_left_sigma, avg_right_sigma))
                    all_power = np.concatenate((avg_left_powers, avg_right_powers))

                    sat_names.append(sat_name)
                    peak_sigmas.append(np.max(all_sigma))
                    peak_powers.append(np.max(all_power))


                    axes[i][j].scatter(unique_left_angles_deg, avg_left_powers, s=0.1, label=sat_name, c="k")
                    axes[i][j].scatter(-unique_right_angles_deg, avg_right_powers, s=0.1, c="k")

                    axes[i][j].scatter(unique_left_angles_deg[counts_left_mask], avg_left_powers[counts_left_mask], s=1, label=sat_name, c="red")
                    axes[i][j].scatter(-unique_right_angles_deg[counts_right_mask], avg_right_powers[counts_right_mask], s=1, c="red")

                    axes[i][j].fill_between(unique_left_angles_deg, avg_left_powers - avg_left_sigma, avg_left_powers + avg_left_sigma, color='blue', alpha=0.5)
                    axes[i][j].fill_between(-unique_right_angles_deg, avg_right_powers - avg_right_sigma, avg_right_powers + avg_right_sigma, color='blue', alpha=0.5)

                    # axes[i][j].scatter(left_angles_deg, np.average(left_powers), s=0.5, label=sat_name, c="k")
                    # axes[i][j].scatter(-right_angles_deg, np.average(right_powers), s=0.5, c="k")
                    axes[i][j].set_xlabel(r"${\theta \rm ~ [deg]}$")
                    axes[i][j].set_ylabel(r'${C/N_0 \rm ~ [dB-Hz]}$')
                    axes[i][j].set_title(rf"${sat_name}$")
                    axes[i][j].set_xlim((-90, 90))

                elif want_time:
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
                    if rms:
                        cnos_rms = []
                        times_rms = []
                    pass_num = 0
                    while k < len(indices):
                        next_index = indices[k]
                        times = times_elevs[past_index:next_index]
                        time_print = ExtractSBF.weeksecondstoutc(2225, times_elevs[next_index])
                        # print(time_print)
                        times = times - times_elevs[past_index + 1]  # CONVERTING TO MINUTES
                        times /= 360
                        if len(times) < 500:
                            past_index = indices[k]
                            k += 1
                            continue
                        # print("passes must be greater than this length")
                        # print(times[250] - times[0])
                        # print(times_elevs[250] - times_elevs[0])
                        if rms:
                            length = next_index-past_index
                            cnos_rms.append(cnos[past_index:next_index])
                            times_rms.append(times)

                        else:
                            if offset:
                                axes[i][j].scatter(times, cnos[past_index:next_index] + 20 * k, s=1, label = f"Pass {pass_num}")
                            else:
                                axes[i][j].scatter(times, cnos[past_index:next_index], c="k", s=0.05)

                        past_index = indices[k]
                        times_max = max(times) # need to take max of times here to ensure it got through the if loop
                        k += 1
                        pass_num += 1

                    if rms:
                        # Determine the maximum length of the lists
                        min_length = min(len(lst) for lst in cnos_rms)

                        trimmed_lists = [lst[:min_length] for lst in cnos_rms]

                        # cnos_rms = np.transpose(cnos_rms) # switching to get in right order for stds
                        cnos_std = np.std(trimmed_lists, axis=1)
                        axes[i][j].plot(cnos_std, c="k")
                        axes[i][j].set_xlabel(r'Time [hr]')
                        axes[i][j].set_ylabel(r'$\sigma_{C/N_0}$')
                        axes[i][j].set_title(rf"${sat_name}$")
                    else:
                        axes[i][j].legend(loc='upper right')
                        axes[i][j].set_xlabel(r'Time [hr]')
                        axes[i][j].set_ylabel(r'Offset ${C/N_0 \rm ~ [dB-Hz]}$')
                        axes[i][j].set_title(rf"${sat_name}$")
                        axes[i][j].set_xlim((0, times_max))
                var += 1
        plt.tight_layout()
        if rms:
            fig.savefig(self.plot_dir + f"RMS_panel_{self.panel_plot_num}.png", dpi=300)
        else:
            if want_time:
                fig.savefig(self.plot_dir + f"time_panel_{self.panel_plot_num}.png", dpi=300)
            elif want_theta:
                fig.savefig(self.plot_dir + f"theta_panel_{self.panel_plot_num}.png", dpi=300)
            else:
                print("You didn't specify theta, time, or RMS.")
        self.panel_plot_num += 1

        # Start constructing the LaTeX table string
        latex_table = r"\begin{table}[ht]\n\centering\n\begin{tabular}{|c|c|c|c|}\n\hline\n"
        latex_table += "Satellite & Peak $\sigma$ & Peak $C/N_0$ \\\\ \\hline\n"

        # Single loop to fill in table rows
        for idx, (sat, sigma, peak_power) in enumerate(zip(sat_names, peak_sigmas, peak_powers)):
            latex_table += f"{sat} & {sigma:.2f} & {peak_power:.2f} \\\\ \\hline\n"

        # Finish the LaTeX table
        latex_table += r"\end{tabular}\n\caption{Satellite Data Table}\n\end{table}"

        # Output LaTeX table string
        print(latex_table)


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

    def plot_one_d_prof(self, sat_names, with_sim=False):
        """
        Plots 1D beam profile in time and angle
        """
        fig, axes = plt.subplots()

        for sat_name in sat_names:
            left_angles_deg, left_powers, right_angles_deg, right_powers = self.get_angles_for_one_d_prof(sat_name, shift=False)
            axes.scatter(-left_angles_deg, left_powers-40, s=0.5, label=sat_name, c="k")
            axes.scatter(right_angles_deg, right_powers-40, s=0.5, c="k")
            power_beam = beam.efield_to_power(inplace=False)

            if with_sim:
                axes.plot(power_beam.axis2_array * 180 / np.pi, power_beam.data_array[0, 0, -1, :, 0])

                # axes.scatter(beam.axis2_array * 180 / np.pi, beam.data_array[1, 0, 48, :, 0], s=0.5)
                # axes.scatter(-beam.axis2_array * 180 / np.pi, beam.data_array[1, 0, 48, :, 0], s=0.5)
            plt.xlim(-90, 90)
            x = beam.data_array[0, 0, :, 0, 0]
            print("x: ", x)
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
        plt.savefig(self.plot_dir + sat_name + "_power_indiv_new.png", bbox_inches='tight', dpi=300)

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
    data_directories = ["/Users/sabrinaberger/d3a_data/"]
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
        # beam_map.plot_one_d_prof(["G14"])
        # beam_map.plot_panel_sats(start=len(beam_map.close_sat_beams)-19, end=len(beam_map.close_sat_beams), chosen=False)

        good_sats = ["G04", "G07", "G10", "G23", "G29", "G30", "E36", "E15"]
        beam_map.chosen_list = good_sats
        beam_map.plot_panel_sats(start=0, end=8, rows=4, cols=2, chosen=True, want_time=True, figsize=(8,12))
        # beam_map.plot_panel_sats(start=0, end=8, rows=4, cols=2, chosen=True, want_time=True, figsize=(8,12), rms=True, offset=False)
        beam_map.plot_panel_sats(start=0, end=8, rows=4, cols=2, chosen=True, want_theta=True, figsize=(8,12), rms=False, offset=False)

        # beam_map.plot_panel_sats(start=8, end=16, rows=4, cols=2, chosen=False, want_time=True, figsize=(8,12))
        # beam_map.plot_panel_sats(start=16, end=24, rows=4, cols=2, chosen=False, want_time=True, figsize=(8,12))
        # beam_map.plot_panel_sats(start=24, end=32, rows=4, cols=2, chosen=False, want_time=True, figsize=(8,12))
        # beam_map.plot_panel_sats(start=32, end=40, rows=4, cols=2, chosen=False, want_time=True, figsize=(8,12))
        # beam_map.plot_panel_sats(start=40, end=46, rows=3, cols=2, chosen=False, want_time=True, figsize=(8,12))

        beam_map.make_plot(show_power=True, shift=True, all_sat=True, plot_title=r"${\rm Approx. 72 ~ hours ~ of ~ GNSS ~ satellite ~ tracks ~ at ~ D3A}$", filename=f"{day}_all_sat.png")
