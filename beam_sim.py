from pyuvdata import UVBeam
import matplotlib.pyplot as plt
import numpy as np
beam = UVBeam()
# beam.peak_normalize()

beam.read_beamfits('/Users/sabrinaberger/mitiono/new_beam_data_2024/beam.fits', use_future_array_shapes=True)
power_beam = beam.efield_to_power(inplace=False)
print(beam.axis1_array)
# plt.plot(beam.data_array[0, 0, -1, :, 0])
# plt.plot(beam.data_array[0, 0, -1, :, 0])
plt.semilogy(beam.axis2_array * 180 / np.pi,  power_beam.data_array[0, 0, 20, :, 0])
# plt.semilogy(-power_beam.axis2_array * 180 / np.pi,  power_beam.data_array[0, 0, -1, :, 0])

# plt.plot(beam.axis1_array, power_beam[0, 0, -1, :, 0])
# plt.plot(-beam.axis1_array, power_beam[0, 0, -1, :, 0])

plt.show()
# the first dimension (of size 2) is the E-field vector components (phi and theta or theta and phi, I can never remember the order)
# the second dimension (of size 2) is the feed polarization, x and y, where x is aligned with East-West, and y will be identical, only rotated 90 degrees.
# the third dimension (of size 49) is the frequency dimension. You can see the frequency array (in Hz) with beam.freq_array
# the fourth dimension (of size 181) is the theta dimension. You can see the theta array (in rad) with beam.axis1_array
# the fifth dimension (of size 360) is the phi dimension. You can see the phi array (in rad) with beam.axis2_array


# def fix_hist_single_sat(self, sat_name, i_pol, i_freq, dB_offset_sim=0, dB_offset_data=0, zenith_rot=0, ns_rot=0, ew_rot=0,
#                     norm_power=True, shift=True):
#     """
#     1D profile where y = power, x = distance from pointing center
#     """
#
#     satellite = self.all_sat_dict[sat_name]
#     elevs = satellite["elevations"][0]
#     az = satellite["azimuths"][0]
#     cnos = satellite["cno"][0]
#     times_cno = satellite["times"][0]
#     times_elevs = satellite["elev_time"][0]
#     times_elevs, cnos, elevs_deg, az_deg = self.match_elevs(times_cno, times_elevs, cnos, elevs, az)
#
#     if shift: # shift > 180 to negative values
#         az_deg = np.asarray([-(360 - x) if x > 180 else x for x in az_deg])
#
#     elevs = elevs_deg * np.pi / 180
#     az = az_deg * np.pi / 180
#
#     # masking
#     mask_az = az_deg < D3A_AZ_deg  # left
#     mask = mask_az
#     left_az = az[mask]
#     left_elevs = elevs[mask]
#     left_cnos = cnos[mask]
#     left_angles_rad = self.convert_angular_distance_from_center_beam(left_elevs, left_az)
#     left_angles_deg = left_angles_rad * 180 / np.pi
#     left_angles_deg = -left_angles_deg
#     left_powers = self.convert_P(left_cnos)
#     mask_az = az_deg > D3A_AZ_deg  # right
#     mask = mask_az
#     right_az = az[mask]
#     right_elevs = elevs[mask]
#     right_cnos = cnos[mask]
#     right_angles_rad = self.convert_angular_distance_from_center_beam(right_elevs, right_az)
#     right_angles_deg = right_angles_rad * 180 / np.pi
#     right_powers = self.convert_P(right_cnos)
#
#     # hot check for min
#     print(right_angles_deg[np.argmin(right_powers)])
#     print(180 * right_elevs[np.argmin(right_powers)] / np.pi)
#     print(180 * right_az[np.argmin(right_powers)] / np.pi)
#
#     # The beam is currently defined in theta and phi,
#     # with the beam's boresight at beam.theta = 0 degrees,
#     # and its horizon is at beam.theta = 90 degrees
#
#     # There are three issues to fix:
#     # (1) The satellite positions are defined in azimuth and elevation (azel),
#     #     so we would like to either re-express beam's theta/phi coordinates in
#     #     terms of azel, or re-express the satellite positions in terms of
#     #     theta and phi.
#     # (2) The beam boresight is at an elevation of D3A_ALT_deg = 81.41,
#     #     which, in phi/theta coordinates, is 90-81.41 = 8.59 deg.
#     #     We would like to rotate the beam so that it is pointed at 81.41
#     #     degrees of elevation, or at theta = 8.59 deg. Note that since
#     #     the beam's resolution is 0.5 degrees, this will need to be rounded
#     #     to theta = 8.5 deg.
#     # (3) The beam may need to be rotated around its boresight (likely 45, 135,
#     #     225, or 315 degrees, to account for the feed's orientation).
#
#     # We can readily fix issues (2) and (3) when calling the beam, with
#     # the optional variables "zenith_rot," "ns_rot," and "ew_rot," e.g.
#     # beam = CSTBeam(beams_folder,zenith_rot = 45, ns_rot = 9.5)
#     # The variable zenith_rot is the rotation around boresight,
#     # which is performed first. In other words, before tilting the dish
#     # away from zenith, we look at it from above, and rotate it by
#     # the amount given by zenith_rot (in degrees), to account for feed rotation.
#     # The variables ns_rot and ew_rot are the rotation around a tilting axis.
#     # In other words, we tilt the dish away from zenith. A positive
#     # ns_rot tilts it towards the south, and negative, towards
#     # the north. Similarly for ew_rot, although D3A/6 dishes don't tilt
#     # in that direction.
#
#     # To load the beams, you only need to specify the folder.
#     # The zenith_rot, ns_rot, and ew_rot variables are optional.
#     # In our case, I added those variables to the hist_single_sat() function,
#     # so that we can easily change the orientation of the beam when
#     # calling that function.
#
#     beams_folder = '/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/Current Research/mitiono/notebook_and_beams/new_new_beams/'
#     new_beam = CSTBeam(beams_folder)
#     beam = new_beam.rotate(zenith_rot=D3A_ALT_deg, ns_rot=90-D3A_AZ_deg)
#
#     # Now, let's solve problem (1).
#     # We would have two options: converting the beam to azel,
#     # or converting the satellite to theta/phi.
#     # I chose to try the latter, because once we got the theta/phi
#     # coordinates, it's going to be easy to retrieve the corresponding
#     # beam indices.
#
#     # Satellite positions in theta/phi:
#     # Thankfully, we're using the "physics" theta/phi convention (not the
#     # "maths" one), so the conversion is (almost) trivial.
#     left_sat_theta = np.pi / 2 - left_elevs
#     left_sat_phi = left_az
#     right_sat_theta = np.pi / 2 - right_elevs
#     right_sat_phi = right_az
#
#     # # Convert to degrees, and don't forget to
#     # # round to the nearest beam.theta_step & beam.phi_step
#     left_sat_theta_deg = beam.theta_step * np.round((180 / np.pi) * left_sat_theta / beam.theta_step)
#     left_sat_phi_deg = beam.phi_step * np.round((180 / np.pi) * left_sat_phi / beam.phi_step)
#     right_sat_theta_deg = beam.theta_step * np.round((180 / np.pi) * right_sat_theta / beam.theta_step)
#     right_sat_phi_deg = beam.phi_step * np.round((180 / np.pi) * right_sat_phi / beam.phi_step)
#
#     # Now we want to find what points on the beam correspond to those
#     # coordinates we just found (sat_theta_deg and sat_phi_deg).
#     # The way I do it is to use beam.phi_step and beam.theta_step.
#     # Say a point is at theta = 9, the corresponding index in the
#     # theta axis will be theta / beam.theta_step = 9 / 0.5 = 18.
#     i_sat_theta_left = (left_sat_theta_deg // beam.theta_step).astype('int')
#     i_sat_phi_left = (left_sat_phi_deg // beam.phi_step).astype('int')
#     i_sat_left = (i_sat_phi_left, i_sat_theta_left)
#     i_sat_theta_right = (right_sat_theta_deg // beam.theta_step).astype('int')
#     i_sat_phi_right = (right_sat_phi_deg // beam.phi_step).astype('int')
#     i_sat_right = (i_sat_phi_right, i_sat_theta_right)
#     # Then, the corresponding beam cut will be beam.directivity[i_pol,i_freq,i_sat_phi,i_sat_theta,:]
#
#     # # Plotting time!
#     fig, ax = plt.subplots(1, 1, figsize=[10, 8])
#
#     # Plot the satellite
#     max_sat_power = max(left_powers.max(), right_powers.max())
#     # I added a variable called norm_power, that is at True by default.
#     # If True, the norm_power variable makes it so that the plotted power
#     # (both the data and simulated beam) is max subtracted in dB
#     # (i.e. normalized by max, in linear scale).
#     # If False, then the dB_offset_data and dB_offset_beam variables can be
#     # used to adjust the height of either curve.
#     if norm_power:
#         ax.scatter(left_angles_deg, left_powers - max_sat_power, c="k", s=0.1)
#         ax.scatter(right_angles_deg, right_powers - max_sat_power, c="k", s=0.1,)
#     else:
#         ax.scatter(left_angles_deg, left_powers - 59.875 + dB_offset_data, c="k", s=0.1)
#         ax.scatter(right_angles_deg, right_powers - 59.875 + dB_offset_data, c="k", s=0.1)
#     ax.scatter(left_angles_deg, left_powers, c="k", s=0.1)
#     ax.scatter(right_angles_deg, right_powers, c="k", s=0.1)
#
#     # Plot the simulated beam
#     line_props = {
#         'linestyle': '-',
#         'color': 'r',
#         'linewidth': 2,
#         'alpha': 1,
#         'color': 'b'
#     }
#     max_beam_power = 10 * np.log10(
#         max(beam.directivity[i_pol, i_freq][i_sat_left].max(), beam.directivity[i_pol, i_freq][i_sat_right].max()))
#     if norm_power:
#         ax.plot(left_angles_deg, 10 * np.log10(beam.directivity[i_pol, i_freq][i_sat_left]) - max_beam_power,
#                 **line_props, label='Simulated beam at {:.3g} GHz'.format(beam.freqs[i_freq]))
#         ax.plot(right_angles_deg, 10 * np.log10(beam.directivity[i_pol, i_freq][i_sat_right]) - max_beam_power,
#                 **line_props)
#
#     else:
#         ax.plot(left_angles_deg, 10 * np.log10(beam.directivity[i_pol, i_freq][i_sat_left]) + dB_offset_sim,
#                 **line_props, label='Simulated beam at {:.3g} GHz'.format(beam.freqs[i_freq]))
#         ax.plot(right_angles_deg, 10 * np.log10(beam.directivity[i_pol, i_freq][i_sat_right]) + dB_offset_sim,
#                 **line_props)
#
#     ax.set_xlim([-20, 20])
#     print(beam.freqs)
#     ax.set_xlabel(r"$\theta$, the angle from beam center [deg]")
#     ax.set_ylabel("Log Power [uncalibrated]")
#     ax.set_title(f"Galileo Satellte: {sat_name}")
#     ax.legend(loc='lower center')
#     # print("Saving single satellite plot at " +
#     #       f"/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/Current Research/mitiono/plots/" + f"sim_beam_{beam.freqs[i_freq]}_{sat_name}.png".format(beam.freqs[i_freq], sat_name))
#     fig.savefig(f"test_sim_beam_{sat_name}.png", dpi=300, bbox_inches='tight')