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
import healpy
from functools import reduce
from scipy.stats import binned_statistic_2d
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import scipy
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

# class GetBeamMap(ExtractSBF):
#     """
#     This class takes in data directories and creates a beam map.
#     """
#     def __init__(self, data_direcs, mask_frequency="1", save_parsed_data_direc="parsed_data", plot_dir="/Users/sabrinaberger/Desktop/beam_paper_plots/", masking=False):
#         super().__init__(data_direcs=data_direcs, include_elev=True, process=True, mask_frequency=mask_frequency, save_parsed_data_direc=save_parsed_data_direc, masking=masking)
#         self.plot_dir = plot_dir
#         self.panel_plot_num = 0 # number of panel plot so we can plot multiple satellites
#         if self.all_sat_dict != None:
#             self.get_close_sats()
#
#     def print_dictionary(self):
#         print(self.all_sat_dict)
#         return
#
#     def replace_dictionary(self, dict):
#         self.all_sat_dict = dict
#         self.get_close_sats()
#
#     def get_close_sats(self, tol_beam=20):
#         self.close_sat_beams = []
#
#         min_diff_az = 100
#         min_diff_el = 100
#         min_nam = ""
#         for key in self.all_sat_dict:
#             if key is not None:
#                 satellite = self.all_sat_dict[key]
#                 times_cno = satellite["times"][0] # SECONDS
#                 times_elevs = satellite["elev_time"][0] # SECONDS
#                 cnos = satellite["cno"][0]
#                 elevs = satellite["elevations"][0]
#                 az = satellite["azimuths"][0]
#                 az = self.shift_az(az) # shifting azimuth
#                 times_elevs, cnos, elevs, az = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
#                 elev_beam = np.full(len(elevs), D3A_ALT_deg)
#                 diff_elev = np.abs(elevs - elev_beam)
#                 az_beam = np.full(len(az), D3A_AZ_deg)
#                 diff_az = np.abs(az - az_beam)
#
#                 # if min_diff_az > min(diff_az) and min_diff_el > min(diff_elev):
#                 #     min_diff_az = min(diff_az)
#                 #     min_diff_el = min(diff_elev)
#                 #     min_nam = key
#                 # Checking if current satellite passes close to center of beam
#                 if (diff_az < tol_beam).any() and (diff_elev < tol_beam).any():
#                     print(f"{key} passes close to center of beam.")
#                     # min_diff_az = min(diff_az)
#                     # min_diff_el = min(diff_elev)
#                     # print(f"min diff az: {min_diff_az}")
#                     # print(f"min diff elev: {min_diff_el}")
#
#                     self.close_sat_beams.append(key)
#                 else:  # not including satellites that do NOT pass close to beam center
#                     continue
#         print(f"Total number of near beam satellites is {len(self.close_sat_beams)}.")
#
#     def get_min_max(self, arrs):
#         arr_concatenated = np.concatenate(arrs, axis=0)
#         min_arr, max_arr = arr_concatenated.min(), arr_concatenated.max()
#         return min_arr, max_arr
#
#     def make_healpix_plot(self, all_sat, plot_title, sat_list=[], filename="heal_pix_all_sat.png"):
#         # Create synthetic HEALPix data
#         if all_sat:
#             sats_to_plot = self.all_sat_dict
#         else:
#             sats_to_plot = sat_list
#         # sorry to not vectorize the below yikes
#         powers_arr = []
#         az_arr = []
#         alt_arr = []
#         times_arr = []
#
#         # below in this for loop I'm binning each az,alt power
#         for i, key in enumerate(sats_to_plot):
#             satellite = self.all_sat_dict[key]
#             for i, sat_typ in enumerate(satellite['sig_type'][0]): # checking to make sure only using L1
#                 if sat_typ is not None:
#                     if "L1" not in sat_typ:
#                         print(sat_typ)
#                         print("Not using L1")
#                         exit()
#             cnos = satellite["cno"][0]
#             times_cno = satellite["times"][0]
#             times_elevs = satellite["elev_time"][0]
#             elevs = satellite["elevations"][0]
#             az = satellite["azimuths"][0]
#
#             times_elevs, cnos_matched, elevs_matched, az_matched = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
#             powers_arr.append(cnos_matched)
#             alt_arr.append(elevs_matched)
#             az_arr.append(az_matched)
#             times_elevs
#
#         # flattening ragged numpy arrays
#         powers_arr = np.concatenate(powers_arr)
#         alt_arr = np.concatenate(alt_arr)
#         # shifting elevs
#         alt_arr = 90 - alt_arr
#         az_arr = np.concatenate(az_arr)
#
#         np.save("az_arr.npy", az_arr)
#         np.save("alt_arr.npy", alt_arr)
#         np.save("powers_arr.npy", powers_arr)
#         # Define HEALPix parameters
#         nside = 5  # HEALPix resolution parameter
#         npix = healpy.nside2npix(nside)
#
#         print("Shifted azimuths 180 degrees...")
#         az_arr = self.shift_az(az_arr)
#
#         # Convert RA/Dec to radians
#         az_rad = np.radians(az_arr)
#         alt_rad = np.radians(alt_arr)  # Convert altitude to declination
#
#         # Convert RA/Dec to HEALPix pixel indices (vectorized)
#         pix_indices = healpy.ang2pix(nside, alt_rad, az_rad)
#
#         # Initialize a HEALPix map to hold the accumulated power values
#         healpix_map = np.zeros(npix)
#
#         # Accumulate power values in the HEALPix map using NumPy's accumarray
#         np.add.at(healpix_map, pix_indices, powers_arr)
#
#         # # Convert RA/Dec to HEALPix pixel indices and accumulate power values
#         # for ra, dec, power in zip(az_arr, alt_arr, powers_arr):
#         #     pix = healpy.ang2pix(nside, np.radians(alt_arr), np.radians(az_arr))  # Convert to radians
#         #     healpix_map[pix] += power  # Accumulate power values
#
#         # Healpix map to Mollweide projection
#         healpy.orthview(healpix_map,  coord='C', cmap='viridis', norm='hist', half_sky=True) #"C" is celestial coordinates
#         # Display the colorbar
#         # plt.colorbar(label='Accumulated Power')
#         # Show the plot
#         plt.savefig(self.plot_dir + "healpy.png", dpi=300)
#
#     def make_plot(self, all_sat, plot_title, sat_list=[], show_power=False, filename="all_sat.png"):
#         """
#         Plot all the satellites or a subset
#         param: all_sat - boolean when true all satellites are plotted, when false, looks for satellite in sat_list
#         param: plot_title - title above plot to be generated
#         param: direc - directory where to save the plot
#         param: filename - name of plot you're making
#         """
#         if all_sat:
#             sats_to_plot = self.all_sat_dict
#         else:
#             sats_to_plot = sat_list
#
#         fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
#
#         # sorry to not vectorize the below yikes
#         powers_arr = []
#         az_arr = []
#         alt_arr = []
#
#         # below in this for loop I'm binning each az,alt power
#         for i, key in enumerate(sats_to_plot):
#             satellite = self.all_sat_dict[key]
#             for i, sat_typ in enumerate(satellite['sig_type'][0]): # checking to make sure only using L1
#                 if sat_typ is not None:
#                     if "L1" not in sat_typ:
#                         print(sat_typ)
#                         print("Not using L1")
#                         exit()
#             cnos = satellite["cno"][0]
#             times_cno = satellite["times"][0]
#             times_elevs = satellite["elev_time"][0]
#             elevs = satellite["elevations"][0]
#             az = satellite["azimuths"][0]
#
#             if show_power:
#                 times_elevs, cnos_matched, elevs_matched, az_matched = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
#                 powers_arr.append(cnos_matched)
#                 alt_arr.append(elevs_matched)
#                 az_arr.append(az_matched)
#             else: # just plotting tracks in black
#                 ax.scatter(np.deg2rad(az), 90 - elevs, s=0.01, c="k")
#                 print(np.deg2rad(az[:2]), elevs[:2])
#
#         if show_power:
#             # flattening ragged numpy arrays
#             powers_arr = np.concatenate(powers_arr)
#             alt_arr = np.concatenate(alt_arr)
#             # shifting elevs
#             alt_arr = 90 - alt_arr
#             az_arr = np.concatenate(az_arr)
#             # 2D binning getting maximums, using ChatGPT to help accelerate generating this plot
#             # Step 1: Define bin edges for azimuth and altitude (in degrees)
#             alt_bins = np.linspace(70, 90, int(100))  # Altitude bins (0 to 90 degrees), binwidth = 0.1 deg
#             az_bins = np.linspace(0, 360, int(100))  # Azimuth bins (-180 to 180 degrees), binwidth = 0.1 deg
#
#             # Step 2: Perform 2D binning using the average of `powers_arr` for each bin
#             peaks, az_edges, alt_edges, binnumber = binned_statistic_2d(
#                 az_arr, alt_arr, powers_arr, statistic='max', bins=[az_bins, alt_bins]
#             )
#
#             # Step 3: Convert bin edges to centers for plotting
#             az_centers = (az_edges[:-1] + az_edges[1:]) / 2  # Convert azimuth edges to bin centers
#             alt_centers = (alt_edges[:-1] + alt_edges[1:]) / 2  # Convert altitude edges to bin centers
#
#             # Step 4: Convert azimuth and altitude to radians for polar plotting
#             azimuth_bin_centers, altitude_bin_centers = np.meshgrid(np.radians(az_centers), alt_centers)
#
#             # Step 5: Create the polar plot with pcolormesh
#             # pcolormesh using binned azimuth and altitude
#             c = ax.pcolormesh(azimuth_bin_centers, altitude_bin_centers, peaks.T[:-1, :-1], cmap='viridis', shading="flat")
#             cbar = fig.colorbar(c, ax=ax)
#             cbar.set_label(r'$C/N_0$ [db-Hz]', rotation=270, labelpad=15)
#
#
#         # Define the center of the circle
#         r_center = np.radians(90 - D3A_ALT_deg)  # altitude (radius)
#         theta_center = np.radians(D3A_AZ_deg)  # azimuth (angle in radians, converted from degrees)
#
#         # Convert polar coordinates (r_center, theta_center) to Cartesian coordinates
#         x_center = r_center * np.cos(theta_center)
#         y_center = r_center * np.sin(theta_center)
#
#         # Define the radius of the circle
#         circle_radius = 100*1.22 * 0.2 / 6
#
#         # Create the circle in Cartesian coordinates
#         circle = plt.Circle((x_center, y_center), circle_radius, fill=False, color='black',
#                             transform=ax.transData._b, linewidth=1)
#         # Add the circle to the plot
#         ax.add_artist(circle)
#
#         ax.set_theta_zero_location('E')
#         ax.set_title(plot_title)
#         # convert circle_radius to polar plot limits
#         ax.set_yticks(range(0, 90 + 10, 10))  # Define the yticks
#         yLabel = ['90', '80', '70', '', '', '', '', '', '', '']
#         ax.set_ylim(70, 90)  # This would make the center at 90 degrees
#
#         ax.set_yticklabels(yLabel)
#         plt.savefig(self.plot_dir + filename, bbox_inches='tight', dpi=300)
#         plt.close()
#
#     def shift_az(self, az_deg):
#         """
#         :param az_deg:  azimuth in degrees between 0 and 360 degrees
#         :return: shifted azimuth in degrees between -180 and 180 degrees
#         """
#         az_deg = np.asarray([-(360 - x) if x > 180 else x for x in az_deg])
#         return az_deg
#
#
#     def save_particular_sat(self, sat_name):
#         satellite = self.all_sat_dict[sat_name]
#         elevs = satellite["elevations"][0]
#         az = satellite["azimuths"][0]
#         cnos = satellite["cno"][0]
#         times_cno = satellite["times"][0]
#         times_elevs = satellite["elev_time"][0]
#         times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
#
#         np.save(f"{sat_name}_times.npy", times_elevs)
#         np.save(f"{sat_name}_cnos.npy", cnos)
#         np.save(f"{sat_name}_elev.npy", elevs_deg)
#         np.save(f"{sat_name}_az.npy", az_deg)
#
#     def mask_array_for_one_d(self, az, elevs, cnos, mask, times):
#         """angles in radians"""
#         az = az[mask]
#         elevs = elevs[mask]
#         cnos = cnos[mask]
#         times = times[mask]
#         angles_rad = self.convert_angular_distance_from_center_beam(elevs, az)
#         angles_deg = angles_rad * 180 / np.pi
#         return angles_deg, cnos, times
#
#     def get_angles_for_one_d_prof(self, sat_name, shift=True):
#         """
#         Gives you a 1D profile for a satellite as a function of angular distance from beam and match times for
#         elevations and CNOs
#         :param sat_name - name of satellite
#         :param shift - whether or not to shift azimuths with shift_az
#
#         :return angles_deg, powers - theta and power of a one d beam profile for a given satellite
#         """
#         satellite = self.all_sat_dict[sat_name]
#         elevs = satellite["elevations"][0]
#         az = satellite["azimuths"][0]
#         cnos = satellite["cno"][0]
#         times_cno = satellite["times"][0]
#         times_elevs = satellite["elev_time"][0]
#         times_matched, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
#         is_sorted = lambda arr: np.all(arr[:-1] <= arr[1:])
#         print("Matched times is sorted" if is_sorted(times_matched) else "Array is not sorted")
#         if shift:  # shift > 180 to negative values
#             print("Shifted azimuths...")
#             az_deg = self.shift_az(az_deg)
#
#         elevs = elevs_deg * np.pi / 180
#         az = az_deg * np.pi / 180
#
#         ## REMOVE CNOS with values that are the same within threshold of 1000s
#
#
#         # TAKING ONLY ONE (ALT, AZ) & TIMESTAMP
#         unique_cnos, unique_indices_cnos = np.unique(cnos, return_index=True)
#         unique_times, unique_indices_times = np.unique(times_matched, return_index=True)
#         unique_az, unique_indices_az = np.unique(az, return_index=True)
#         unique_elevs, unique_indices_elevs = np.unique(elevs, return_index=True)
#         # arrays = [unique_indices_az, unique_indices_times]
#
#         # unique_location_indices = reduce(np.intersect1d, arrays) # applies intersect consecutively
#         # unique_location_indices = unique_indices_times
#         # print(np.sum(np.isin(unique_indices_cnos, unique_location_indices)))
#
#         # times_matched, az, elevs, cnos, elevs_deg, az_deg = times_matched[unique_location_indices], az[unique_location_indices], elevs[unique_location_indices], cnos[unique_location_indices], elevs_deg[unique_location_indices], az_deg[unique_location_indices]
#
#         # masking
#         right_mask = az_deg > D3A_AZ_deg  # right
#         left_mask = az_deg < D3A_AZ_deg  # left
#
#         # this takes in radians
#         right_angles_deg, right_powers, right_times = self.mask_array_for_one_d(az, elevs, cnos, right_mask, times_matched)
#         left_angles_deg, left_powers, left_times = self.mask_array_for_one_d(az, elevs, cnos, left_mask, times_matched)
#
#         print("right times is sorted" if is_sorted(right_times) else "Array is not sorted")
#         print("left times is sorted" if is_sorted(left_times) else "Array is not sorted")
#
#         all_angles_deg = np.concatenate([-left_angles_deg, right_angles_deg])
#         all_powers = np.concatenate([left_powers, right_powers])
#         all_times = np.concatenate([left_times, right_times])
#
#         print("all times is sorted" if is_sorted(all_times) else "Array is not sorted")
#
#         # Get the indices that would sort the 'times' array
#         sorted_indices = np.argsort(all_times)
#
#         # Sort all arrays based on 'times'
#         all_angles_deg = all_angles_deg[sorted_indices]
#         all_powers = all_powers[sorted_indices]
#         all_times = all_times[sorted_indices]
#
#         split_all_times, split_all_powers, split_all_angles_deg = self.split_sat(all_times, all_powers, all_angles_deg, sat_name=sat_name)
#
#
#         return split_all_times, split_all_powers, split_all_angles_deg, all_times, all_powers, all_angles_deg
#
#     def convert_C_n_to_SNR_linear(self, c_n, bw=2e6):
#         # assuming nominal L1 bandwidth of 2Mhz
#         c_n = np.asarray(c_n)
#         bw = np.log10(bw)
#         s_n_db = c_n - bw
#         s_n_linear = 10**(s_n_db/10)
#         return s_n_linear
#
#     def convert_SNR_linear_C_n(self, SNR_linear, bw=2e6):
#         # assuming nominal L1 bandwidth of 2Mhz
#         s_n_db = 10 * np.log10(SNR_linear)
#         bw = np.log10(bw)
#         c_n = s_n_db + bw
#         return c_n
#
#     def within_threshold(self, arr, threshold=2):
#         # Ensure arr is a NumPy array
#         arr = np.asarray(arr)
#
#         # If arr is a scalar or has only one element, return True
#         if arr.ndim == 0 or len(arr) <= 1:
#             return True
#
#         # Compute pairwise absolute differences
#         diff_matrix = np.abs(arr[:, None] - arr)
#         # Check if all differences are within the threshold
#         return np.all(diff_matrix >= threshold)
#
#     # def average_C_N_0(self, grouped_c_n, grouped_angles, decimals_round=0, min_passes_for_average=2):
#     #     grouped_angles_round = [np.round(arr, decimals=decimals_round) for arr in grouped_angles] # rounded decimals
#     #
#     #     unique_elements_per_array = [np.unique(arr) for arr in grouped_angles_round]
#     #     unique_elements_per_array = np.concatenate(unique_elements_per_array).flatten()
#     #     unique_angles_deg = np.unique(unique_elements_per_array)
#     #     c_n_angle = np.empty(len(unique_angles_deg))
#     #     c_n_angle_min =  np.empty(len(unique_angles_deg)) # min for degree
#     #     c_n_angle_max =  np.empty(len(unique_angles_deg)) # min for degree
#     #
#     #     c_n_angle_all = [] # all C_Ns in group
#     #     std_angle = np.empty(len(unique_angles_deg))
#     #     snr_angle = np.empty(len(unique_angles_deg))
#     #     std_from_linear = np.empty(len(unique_angles_deg))
#     #     counts = np.empty(len(unique_angles_deg))
#     #
#     #     # for i, angle in enumerate(unique_angles_deg):
#     #     #     c_n_curr_angle = []
#     #     #     for y in range(len(grouped_c_n)):
#     #     #         # group unique indices
#     #     #         mask = grouped_angles_round[y] == angle
#     #     #         c_n = grouped_c_n[y][mask]
#     #     #         if len(c_n) < 1:  # empty list
#     #     #             continue
#     #     #         c_n_curr_angle.append(c_n)
#     #     #
#     #     #     if len(c_n_curr_angle) < min_passes_for_average:
#     #     #         # c_n_curr_angle.append([0])
#     #     #         c_n_angle[i] = np.nan
#     #     #         std_angle[i] = np.nan
#     #     #         snr_angle[i] = np.nan
#     #     #         std_from_linear[i] = np.nan
#     #     #         c_n_angle_min[i] = np.nan
#     #     #         c_n_angle_max[i] = np.nan
#     #     #         c_n_angle_all.append([])
#     #     #         continue
#     #     #     ## average cn and std
#     #     #     counts[i] = len(c_n_curr_angle)
#     #     #     c_n_angle_all.append(c_n_curr_angle) # saving all C_Ns
#     #     #
#     #     #     c_n_curr_angle = [item for sublist in c_n_curr_angle for item in sublist] # flattening ragged list
#     #     #     c_n_angle_min[i] = np.min(c_n_curr_angle) # getting min
#     #     #     c_n_angle_max[i] = np.max(c_n_curr_angle) # getting max
#     #     #     c_n_angle[i] = np.mean(c_n_curr_angle) # mean of c_n_angle in logspace
#     #     #     std_angle[i] = np.std(c_n_curr_angle) # weird std of c_n_angle in logspace
#     #     #
#     #     #     ## average linear SNR
#     #     #     s_n_curr_angle = self.convert_C_n_to_SNR_linear(c_n_curr_angle) # convert from C_n to linear SNR
#     #     #     # mean_s_n_curr_angle = np./(s_n_curr_angle) # take mean of SNR
#     #     #     snr_angle[i] = mean_s_n_curr_angle # mean of SNR in linear space
#     #     #
#     #     #     ## linear SNR standard deviation
#     #     #     s_n_curr_angle_std = np.std(s_n_curr_angle)
#     #     #     if s_n_curr_angle_std == 0:
#     #     #         s_n_curr_angle_std = 1
#     #     #     s_n_curr_angle_std = 10 * np.log10(s_n_curr_angle_std) # convert back into db
#     #     #     std_from_linear[i] = s_n_curr_angle_std # convert back to log space
#     #     return unique_angles_deg, c_n_angle, c_n_angle_all, std_angle, snr_angle, std_from_linear, counts, c_n_angle_max, c_n_angle_min
#
#
#     def split_sat(self, times, powers, thetas, sat_name):
#         diff_times = np.diff(times) # getting difference between adjacent elements
#         is_sorted = lambda arr: np.all(arr[:-1] <= arr[1:])
#         assert is_sorted(times) # ensuring time is sorted so chunking works
#
#         gap_threshold = np.max(diff_times) * 0.7  # threshold for splitting chunks at 90% maximum
#
#         chunk_indices = np.where(diff_times > gap_threshold)[0] + 1
#         print("CHECK OUTPLOT PLOTS TO MAKE SURE CHUNKING IS HAPPENING CORRECTLY")
#         time_chunks = np.split(times, chunk_indices)
#         theta_chunks = np.split(thetas, chunk_indices)
#         power_chunks = np.split(powers, chunk_indices)
#
#         if sat_name == "E15": # just keeping chunks but saving by eye from plots
#             time_chunks = [time_chunks[1]]
#             theta_chunks = [theta_chunks[1]]
#             power_chunks = [power_chunks[1]]
#
#         if sat_name == "E36": # just keeping chunks but saving by eye from plots
#             time_chunks = [time_chunks[1]]
#             theta_chunks = [theta_chunks[1]]
#             power_chunks = [power_chunks[1]]
#
#         # Plot each chunk with a different color
#         plt.close()
#         plt.figure(figsize=(10, 6))
#         colors = plt.cm.viridis(np.linspace(0, 1, len(time_chunks)))  # Using a colormap for colors
#         print(sat_name)
#         print(len(time_chunks))
#         for i, (t_chunk, theta_chunk) in enumerate(zip(time_chunks, theta_chunks)):
#             print(f'Chunk {i + 1}')
#             plt.scatter(t_chunk, theta_chunk, color=colors[i], label=f'Chunk {i + 1}')
#
#         # Label the plot
#         plt.xlabel('Times')
#         plt.ylabel('Thetas')
#         plt.title('Color-Coded Chunks of Data')
#         plt.legend()
#         plt.savefig(f"chunks/{sat_name}_chunking.png")
#
#         return time_chunks, power_chunks, theta_chunks
#
#
#     def plot_panel_sats(self, rows=6, cols=3, start=0, end=19, chosen=True, figsize=(8, 12), offset=True,
#                         want_theta=False, want_time=False, rms=False, airy_disk=True, create_latex_table=True):
#         var = 0
#         print("num close_sats", len(self.close_sat_beams))
#         if chosen:
#             close_sat_beams = self.chosen_list
#         else:
#             close_sat_beams = self.close_sat_beams[start:end] # local to function version, NO SELF
#
#         sat_names = []
#         peak_sigmas = []
#         peak_powers = []
#         sigma_at_max = []
#         counts_min_sigmas = []
#         min_max_power_range = []
#         diff = []
#         diff_at_max = []
#         fig, axes = plt.subplots(rows, cols, figsize=figsize)
#         # Create a colormap to match the theta plot
#
#         for i in range(rows):
#             for j in range(cols):
#                 sat_name = close_sat_beams[var]
#                 times_grouped, powers_grouped, angles_grouped, _, _, _ = self.get_angles_for_one_d_prof(sat_name,
#                                                                                                         shift=True)
#                 if create_latex_table:
#                     decimals_round = 100 # NO BINNING ?
#                 else:
#                     decimals_round = 1
#                 unique_angles_deg, avg_powers, all_powers, avg_sigmas, snr_angle, std_from_linear, counts, power_angle_max, power_angle_min = (
#                     self.average_C_N_0(powers_grouped, angles_grouped, decimals_round=decimals_round, min_passes_for_average=2))
#                 cmap = cm.get_cmap("viridis", 6)  # Specify 6 distinct colors
#                 colors = [cmap(i) for i in range(6)]  # Extract colors for each pass
#                 if want_theta:
#                     if not np.all(np.isnan(avg_powers)):
#                         sat_names.append(sat_name)
#                         peak_sigmas.append(np.nanmax(avg_sigmas))
#                         peak_powers.append(np.nanmax(avg_powers))
#                         index = np.nanargmax(avg_powers)
#                         sigma_at_max.append(avg_sigmas[index])
#                         print()
#                         diff_at_max.append(power_angle_max[index] - power_angle_min[index])
#                         diff_all = []
#                         for z, power_chunk in enumerate(all_powers):
#                             if len(power_chunk) > 0:
#                                 power_chunk_flattened = np.concatenate(power_chunk).ravel()
#                                 diff_curr = np.nanmax(power_chunk_flattened)-np.nanmin(power_chunk_flattened)
#                                 diff_all.append(diff_curr)
#                         # print("all_powers_flattened_index")
#                         # print(all_powers_flattened_index)
#                         # min_max_power_range.append((min(all_powers_flattened_index), max(all_powers_flattened_index)))
#                         diff.append(min(diff_all))
#                         print()
#                     # Define colormap and normalization
#                     # Convert to asymmetric errors
#                     # lower_errors = np.where(avg_sigmas < 0, np.abs(avg_sigmas), 0)  # Positive values (lower bounds)
#                     # upper_errors = np.where(avg_sigmas > 0, avg_sigmas, 0)  # Negative values (upper bounds)
#
#                     # Combine into asymmetric yerr
#                     max_min_powers = (avg_powers-power_angle_min, power_angle_max-avg_powers)
#                     # Create the scatter plot with error bars
#
#                     sc_color = axes[i][j].scatter(
#                         unique_angles_deg, avg_powers,
#                         label=sat_name, s=5, cmap=cmap, norm=Normalize(vmin=1, vmax=6),
#                         c=counts, alpha=0.9)
#
#                     axes[i][j].errorbar(
#                         unique_angles_deg,  # X data
#                         avg_powers,  # Y data
#                         yerr=max_min_powers,  # Error values
#                         label=sat_name,  # Label for the plot
#                         fmt='none',  # No markers for data points
#                         capsize=0,  # No caps on the error bars,
#                         zorder=-100,
#                         ecolor='#D3D3D3'
#                     )
#                     # Add colorbar
#                     cbar = plt.colorbar(sc_color, ax=axes[i][j])
#                     cbar.set_label('Counts')
#                     # Set colorbar ticks to integers between 0 and 6
#                     cbar.set_ticks([1, 2, 3, 4, 5, 6])
#                     cbar.set_ticklabels([1, 2, 3, 4, 5, 6])
#
#                     sat_constellation = self.all_sat_dict[sat_name]['sig_type'][0][0]
#
#                     if "L1" in sat_constellation:
#                         obs_freq = 1575e6 #L1 GPS
#                     else:
#                         print("OBS FREQUENCY UNKNOWN.")
#
#                     if airy_disk:
#                         theta, intensity = self.airy_disk_pattern(3e8/obs_freq, 3, max_rad=5 * np.pi/180, max_I=np.max(snr_angle))
#                         theta *= 180.0 / np.pi
#                         axes[i][j].plot(theta, 10 * np.log10(intensity), c="k", label="Airy disk", alpha=0.5)
#                     axes[i][j].set_xlabel(r"${\theta \rm ~ [deg]}$")
#                     axes[i][j].set_ylabel(r'Binned ${C/N_0 \rm ~ [dB-Hz]}$')
#                     axes[i][j].set_title(rf"${sat_name}$")
#                     axes[i][j].set_xlim((-90, 90))
#                     axes[i][j].set_ylim((20, 65))
#                 elif want_time:
#                     pass_num = 0
#                     for times, powers in zip(times_grouped, powers_grouped):
#                         if pass_num > 5 or len(times) < 500:
#                             continue
#                         times -= np.min(times)
#                         times /= 3600
#                         if offset:
#                             print(pass_num)
#                             k = pass_num + 1 # avoiding 0
#                             axes[i][j].scatter(times, powers + 20 * k, label=f"Pass {k}", c=colors[int(pass_num)], s=0.5)
#                         else:
#                             axes[i][j].scatter(times, powers, c="k", s=0.05)
#                         pass_num += 1
#                     # Retrieve handles and labels from the scatter plot
#                     handles, labels = axes[i][j].get_legend_handles_labels()
#
#                     # Create a new legend with increased marker sizes
#                     axes[i][j].legend(handles, labels, markerscale=5, loc='upper right')
#                     axes[i][j].set_xlabel(r'Time [hr]')
#                     axes[i][j].set_ylabel(r'Offset ${C/N_0 \rm ~ [dB-Hz]}$')
#                     axes[i][j].set_title(rf"${sat_name}$")
#                     # get_max_times = np.concatenate(times_grouped) / 3600
#                     # axes[i][j].set_xlim((0, np.max(get_max_times)))
#                 var += 1
#         if rms:
#             fig.savefig(self.plot_dir + f"RMS_panel_{self.panel_plot_num}.png", dpi=300)
#         else:
#             if want_time:
#                 fig.tight_layout()
#                 fig.savefig(self.plot_dir + f"time_panel_{self.panel_plot_num}.png", dpi=300)
#                 plt.close(fig)
#             elif want_theta:
#                 fig.tight_layout()
#                 fig.savefig(self.plot_dir + f"theta_panel_{self.panel_plot_num}.png", dpi=300)
#                 plt.close(fig)
#             else:
#                 print("You didn't specify theta, time, or RMS.")
#         self.panel_plot_num += 1
#
#         if create_latex_table:
#             # Start constructing the LaTeX table string
#             latex_table = r"\begin{table}[ht]\n\centering\n\begin{tabular}{|c|c|c|c|c|}\n\hline\n"
#             latex_table += "Satellite & Max - Min $C/N_0$ [db-Hz] & Peak $C/N_0$  [db-Hz] & $\sigma$ at Peak $C/N_0$  [db-Hz] & Max $\sigma$  [db-Hz] \\\\ \\hline\n"
#
#             # Single loop to fill in table rows
#             for idx, (sat, diff, sigma_max, max_sigma, peak_power) in enumerate(zip(sat_names, diff_all, sigma_at_max, peak_sigmas, peak_powers)):
#                 latex_table += f"{sat} & {diff:.2f} & {peak_power:.2f} & {sigma_max:.2f}& {max_sigma:.2f} \\\\ \\hline\n"
#
#             # Finish the LaTeX table
#             latex_table += r"\end{tabular}\n\caption{Satellite Data Table}\n\end{table}"
#
#             # Output LaTeX table string
#             print(latex_table)
#         print(counts_min_sigmas)
#         print(min_max_power_range)
#         print(diff)
#
#         # plt.close()
#         # Create the plot
#         plt.figure(figsize=(8, 5))
#         plt.bar(sat_names, peak_sigmas, color='skyblue')
#
#         # Turn satellite names horizontally
#         plt.xticks(rotation=90)
#
#         # Labels and title
#         plt.xlabel("Satellite Name")
#         plt.ylabel(r"Maximum $\sigma$ Value")
#         # Show plot
#         plt.tight_layout()  # Adjusts plot to fit into figure area
#         plt.savefig(self.plot_dir + "max_sigmas.png")
#
#     ## AIRY DISK
#
#     def airy_disk_pattern(self, wavelength, aperture_radius, max_I, max_rad=np.pi/4):
#         # Calculate wave number
#         k = 2 * np.pi / wavelength
#
#         theta = np.linspace(-max_rad, max_rad, int(1e4))
#
#         # Calculate intensity pattern
#         intensity = max_I * (2 * special.j1(k * aperture_radius * np.sin(theta)) / (
#                 k * aperture_radius * np.sin(theta))) ** 2
#         intensity[theta == 0] = max_I  # Handle the singularity at theta = 0
#
#         return theta, intensity
#
#     def compare_sats(self, list_sats):
#         for sat in list_sats:
#             satellite = self.all_sat_dict[sat]
#             elevs = satellite["elevations"][0]
#             az = satellite["azimuths"][0]
#             cnos = satellite["cno"][0]
#             times_cno = satellite["times"][0]
#             times_elevs = satellite["elev_time"][0]
#             times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
#             plt.scatter(times_elevs, cnos, label=sat, s=0.1)
#             plt.title(f"L{self.mask_frequency}")
#         plt.legend()
#         plt.savefig(f"../plots/compare_sats_{self.mask_frequency}.png", dpi=400)
#
#     def overplot_days(self, sat_name):
#         satellite = self.all_sat_dict[sat_name]
#         times_cno = satellite["times"][0]
#         times_elevs = satellite["elev_time"][0]
#         cnos = satellite["cno"][0]
#         elevs = satellite["elevations"][0]
#         az = satellite["azimuths"][0]
#         times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
#
#         # below I'm splicing the array when the difference between times in subsequent elements is greater than 10k
#         # this splicing tolerance might have to be changed
#         diff_times = np.diff(times_elevs)
#         indices = np.argwhere(diff_times > 10000).flatten()
#         print(indices)
#         past_index = 0
#         i = 0
#         while i < len(indices):
#             next_index = indices[i]
#             print(next_index)
#             times = times_elevs[past_index:next_index]
#             times = times - times_elevs[past_index+1]
#             if len(times) < 250:
#                 past_index = indices[i]
#                 i += 1
#                 continue
#             plt.scatter(times, cnos[past_index:next_index] + 10*i, label=f"Day {i}", s=0.5)
#             past_index = indices[i]
#             i += 1
#
#         plt.ylabel(r'${C/N_0 \rm ~ [dB-Hz]}$')
#         plt.legend()
#         plt.title(sat_name)
#         plt.xlim(0, 25000)
#         plt.savefig(self.plot_dir + "overplotted_offset.png")
#
#     def plot_one_d_prof(self, sat_names, with_sim=False):
#         """
#         Plots 1D beam profile in time and angle
#         """
#         fig, axes = plt.subplots()
#
#         for sat_name in sat_names:
#             split_all_times, split_all_powers, split_all_angles_deg, all_times, all_powers, all_angles_deg = self.get_angles_for_one_d_prof(sat_name, shift=False)
#             axes.scatter(all_angles_deg, all_powers, s=0.5, label=sat_name, c="k")
#
#             plt.xlim(-90, 90)
#             # plt.close()
#             # cnos = satellite["cno"][0]
#             # satellite = self.all_sat_dict[sat_name]
#             # times_cno = satellite["times"][0]
#             # plt.scatter(times_cno, cnos, c="k", s=0.1)
#             # plt.xlabel(r"${\rm Time ~ [s]}$")
#             # plt.ylabel(r'${C/N_0 \rm ~ [dB-Hz]}$')
#             # plt.title(sat_name)
#             # plt.savefig(self.plot_dir + sat_name + "_time_indiv.png", bbox_inches='tight')
#         plt.xlabel(r"${\theta \rm ~ [deg]}$")
#         plt.ylabel(r'${C/N_0 \rm ~ [dB-Hz]}$')
#         plt.legend()
#         plt.title(sat_name[:3])
#         plt.savefig(self.plot_dir + sat_name + "_power_indiv_new.png", bbox_inches='tight', dpi=300)
#
#     @staticmethod
#     def get_bessel_0(ka_sintheta):
#         return special.j0(ka_sintheta)
#
#     @staticmethod
#     def get_bessel_1(ka_sintheta):
#         return special.j1(ka_sintheta)
#
#     @staticmethod
#     def convert_angular_distance_from_center_beam(alt, az):
#         """
#         Finds the angular distance from the beam center for a given alt and az using the dot product (lazy trig).
#         """
#         x_beam = np.cos(D3A_AZ) * np.cos(D3A_ALT)
#         y_beam = np.sin(D3A_AZ) * np.cos(D3A_ALT)
#         z_beam = np.sin(D3A_ALT)
#
#         x_sat = np.cos(az) * np.cos(alt)
#         y_sat = np.sin(az) * np.cos(alt)
#         z_sat = np.sin(alt)
#
#         cos_ang = x_beam*x_sat + y_beam*y_sat + z_beam*z_sat
#         ang = np.arccos(cos_ang)
#         return ang
#
#     @staticmethod
#     def convert_P(C_N, nothing=True):
#         """
#         This function converts the receiver's C/N_0 measurement into a classical power measurement in dBw
#         Convert C/N_0 to power in dBw.
#         """
#         if nothing:
#             return C_N
#         N_sys = 8.5  # dB
#         Tant = 30  # K
#         P = C_N + 10 * np.log10(Tant + 290 * (10 ** (N_sys / 10) - 1)) - 228.6  # last i power
#         return P
#
#
#     @staticmethod
#     def match_elevs(times_cno, times_elevs, cnos, elevs, az):
#         """
#         Septentrio returns cnos and postiions of different lengths. This function just finds their intersection.
#         :returns intersected times_elevs, cnos, elevs, az
#         """
#
#         is_sorted = lambda arr: np.all(arr[:-1] <= arr[1:])
#         # assert(is_sorted(times_elevs))
#         # assert(is_sorted(times_cno))
#
#         indices = np.intersect1d(times_elevs, times_cno, return_indices=True)
#         elev_indices = indices[1]
#         cnos_indices = indices[2]
#
#         times_elevs = times_elevs[elev_indices]
#         elevs = elevs[elev_indices]
#         az = az[elev_indices]
#         cnos = cnos[cnos_indices]
#         return times_elevs, cnos, elevs, az

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

        # beam_map.plot_panel_sats(start=0, end=8, rows=4, cols=2, chosen=True, want_time=True, want_theta=False, figsize=(8,12), rms=False, offset=True, airy_disk=False)
        beam_map.plot_panel_sats(start=0, end=8, rows=4, cols=2, chosen=True, want_time=False, want_theta=True, figsize=(8,12), rms=False, offset=False, airy_disk=False)

        # beam_map.make_plot(show_power=True, all_sat=True, plot_title=r"${\rm Approx. 72 ~ hours ~ of ~ GNSS ~ satellite ~ tracks ~ at ~ D3A}$", filename=f"{day}_all_sat.png")
        beam_map.make_healpix_plot(all_sat=True, plot_title=r"${\rm Approx. 72 ~ hours ~ of ~ GNSS ~ satellite ~ tracks ~ at ~ D3A}$", filename=f"{day}_all_sat.png")
