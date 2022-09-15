#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 17:02:11 2022.
@author: vincent
"""

import matplotlib.pyplot as plt
import numpy as np
import copy
from scipy.special import jn


class CSTBeam():
    """
    Class containing beams as simulated by CST.

    This class take a folder where a beam has been converted to the
    "all_*.npy" formats, and loads it.
    """

    def __init__(self, beams_folder, load_comps=False, load_axratios=False):
        self.directivity = np.load(beams_folder + 'all_directivity.npy')
        self.freqs = np.load(beams_folder + 'all_freqs.npy')
        self.phi = np.load(beams_folder + 'all_phi.npy')
        self.theta = np.load(beams_folder + 'all_theta.npy')
        self.phi_step = self.phi[1][0]  # = 0.5 (deg)
        self.theta_step = self.theta[0][1]
        self.freq_min = self.freqs[0]
        self.freq_step = self.freqs[1] - self.freqs[0]
        self.gains = np.max(self.directivity, (2, 3))
        self.wl = 2.99792458e8 / (self.freqs * 1e9)  # wavelength, in meters
        self.A_e = (self.gains * self.wl ** 2) / 4 / np.pi
        if load_comps:
            self.thetaphi_comps = np.load(beams_folder + 'all_thetaphi_comps.npy')
        if load_axratios:
            self.axratios = np.load(beams_folder + 'all_axratios.npy')

    def rotate(self, zenith_rot=0, ns_rot=0, ew_rot=0):
        """
        Rotate the beam.

        This function rotates the beam by the desired angles.
        """
        # Make a copy of the beam
        new_beam = copy.deepcopy(self)

        # Make sure the rotations are at the right resolution
        ns_rot = new_beam.theta_step * np.round(ns_rot / new_beam.theta_step)
        ew_rot = new_beam.theta_step * np.round(ew_rot / new_beam.theta_step)
        zenith_rot = new_beam.phi_step * np.round(zenith_rot / new_beam.phi_step)

        # I do the rotations using this method:
        # https://stla.github.io/stlapblog/posts/RotationSphericalCoordinates.html

        a_x = ew_rot * (np.pi / 180)
        a_y = ns_rot * (np.pi / 180)
        a_z = zenith_rot * (np.pi / 180)

        R_x = np.array([[np.cos(a_x / 2), -1j * np.sin(a_x / 2)], [-1j * np.sin(a_x / 2), np.cos(a_x / 2)]])
        R_y = np.array([[np.cos(a_y / 2), -np.sin(a_y / 2)], [np.sin(a_y / 2), np.cos(a_y / 2)]])
        R_z = np.array([[np.exp(-1j * a_z / 2), 0], [0, np.exp(1j * a_z / 2)]])

        for R in [R_z, R_x, R_y]:  # this determines the order of the rotations
            if not np.allclose(R, np.eye(2)):
                t = new_beam.theta * (np.pi / 180)
                p = new_beam.phi * (np.pi / 180)
                psi = np.array([np.cos(t / 2), np.exp(1j * p) * np.sin(t / 2)])
                psi[0] = R[0, 0] * psi[0] + R[0, 1] * psi[1]
                psi[1] = R[1, 0] * psi[0] + R[1, 1] * psi[1]
                new_theta = (180 / np.pi) * 2 * np.arctan2(np.abs(psi[1]), np.abs(psi[0]))
                new_phi = (180 / np.pi) * (np.angle(psi[1]) - np.angle(psi[0]))
                # Reindex
                i_new_theta = (new_theta / new_beam.theta_step).astype('int')
                i_new_phi = (new_phi / new_beam.phi_step).astype('int')
                new_beam.directivity = new_beam.directivity[:, :, i_new_phi, i_new_theta]

        return new_beam

    def plot_1d(self, freq=1.0, phi_cut=0, projection=None, i_pol=0, dB=True, norm_max=True, airy=False, r=3,
                show_ylabel=False, airy_alt=False, ax=None):
        """
        Plot a 1D beam cut .

        This function plots the beam in 1D, in given projection ('rectilinear' or 'polar'),
        along the cut phi_cut, at frequency index i_freq, for polarization pol.
        It returns the figure and axis.
        """
        line_props = {
            'linestyle': '-',
            'linewidth': 2,
            'alpha': 1
        }

        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        if type(phi_cut) == list:
            phi_cut = np.array(phi_cut)
        else:
            phi_cut = np.array([phi_cut])

        i_freq = np.argmin(np.abs(self.freqs - freq))
        if freq not in self.freqs:
            print('Frequency {:.3g} GHz was not simulated. Using the closest ({:.3g} GHz) instead.'.format(freq,
                                                                                                           self.freqs[
                                                                                                               i_freq]))
            print('A list of the simulated frequencies is stored as beam_name.freqs.')

        phi_cut_index = (phi_cut // self.phi_step).astype(int)
        phi_cut_index_2 = ((phi_cut + 180) // self.phi_step).astype(int)
        if ax != None:
            if projection != None:
                print('Both an \'ax\' and a \'projection\' arguments were passed.')
                print('The projection from the \'ax\' argument overrides the \'projection\' argument.')
            projection = ax.name
        if projection == 'rectilinear':
            figsize = [12, 8]
            radians = False
        if projection == 'polar':
            figsize = [16, 16]
            radians = True
        if ax == None:
            fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw={'projection': projection})
        rad_scale = 1 + (radians * (np.pi / 180 - 1))
        ylabel = r'$I(\theta)$ (arbitrary units)'

        power = self.directivity
        if norm_max:
            power = np.moveaxis(np.moveaxis(power, [0, 1], [2, 3]) / np.max(power, axis=(2, 3)), [0, 1], [2, 3])
        if dB:
            power = 10 * np.log10(power)
            ylabel = r'$I(\theta)$ (dB)'
        if projection == 'polar':
            ax.set_theta_zero_location("N")
        for i_phi_cut in range(phi_cut_index.shape[0]):
            ax.plot(rad_scale * self.theta[phi_cut_index[i_phi_cut]], power[i_pol, i_freq][phi_cut_index[i_phi_cut]][:],
                    colors[i_phi_cut], label=airy * 'Simulated beam ' + r'$\phi$ = {:.4g}$^\circ$'.format(
                    self.phi[phi_cut_index[i_phi_cut]][0]), zorder=1, **line_props)
            ax.plot(-rad_scale * self.theta[phi_cut_index_2[i_phi_cut]],
                    power[i_pol, i_freq][phi_cut_index_2[i_phi_cut]][:], colors[i_phi_cut], zorder=1, **line_props)

            if airy:
                airy_res = 0.1
                airy_theta = np.arange(-180, 180, airy_res)
                if dB:
                    I0 = 10 ** (np.max(power[i_pol, i_freq][i_phi_cut]) / 10)
                else:
                    I0 = np.max(power[i_pol, i_freq][i_phi_cut])

                if airy_alt:
                    k = 2 * r / self.wl[i_freq]
                    airy_func = lambda I0, k, theta: I0 * (np.sinc(theta * np.pi / 180 * k)) ** 2
                    airy = airy_func(I0, k, airy_theta)

                else:
                    k = 2 * np.pi / self.wl[i_freq]
                    airy_func = lambda I0, k, r, theta: I0 * (
                                2 * jn(1, k * r * np.sin(theta * np.pi / 180)) / k / r / np.sin(
                            theta * np.pi / 180)) ** 2
                    airy = airy_func(I0, k, r, airy_theta)

                if dB:
                    airy = 10 * np.log10(airy)
                ax.plot(rad_scale * airy_theta, airy, label='Airy pattern with same max', zorder=0, color='#eeaa00')

        ax.set_xlabel(r'$\theta$')
        if show_ylabel:
            ax.set_ylabel(ylabel)
        if projection == 'rectilinear':
            ax.set_title(
                '{:.3g} GHz, $\phi$ = {:.1f}$^\circ$'.format(self.freqs[i_freq], self.phi[phi_cut_index[i_phi_cut]][0]))
        if projection == 'polar':
            ax.set_title('{:.3g} GHz'.format(self.freqs[i_freq]))
        ax.legend(loc='lower center')

    def plot_2d(self, freq, mode='uv', i_pol=0, dB=True, front_cutoff=90, back_cutoff=-1, norm_max=False):
        """
        Plot a 2D beam (uv).

        This function plots the beam in a 2D uv plot.
        I need to check whether my definition of uv coordinates is good.
        """
        i_freq = np.argmin(np.abs(self.freqs - freq))
        if freq not in self.freqs:
            print('Frequency {:.3g} GHz was not simulated. Using the closest ({:.3g} GHz) instead.'.format(freq,
                                                                                                           self.freqs[
                                                                                                               i_freq]))
            print('A list of the simulated frequencies is stored as beam_name.freqs.')

        if mode == 'uv':
            u = np.sin(self.theta * np.pi / 180) * np.cos(self.phi * np.pi / 180)
            v = np.sin(self.theta * np.pi / 180) * np.sin(self.phi * np.pi / 180)
            if back_cutoff == -1:
                to_plot = np.unique(np.where(self.theta <= front_cutoff)[1])
                print('Plotting from theta = 0 to theta = {:.3g} degrees.'.format(front_cutoff))
                if front_cutoff > 90:
                    print('Careful, uv projection misbehaves if plotting more than 90 degrees of the beam.')
            else:
                to_plot = np.unique(np.where(self.theta >= back_cutoff)[1])
                print('Plotting from theta = {:.3g} to theta = 180 degrees.'.format(back_cutoff))
                if back_cutoff < 90:
                    print('Careful, uv projection misbehaves if plotting more than 90 degrees of the beam.')
            figsize = [12, 12]

        if mode == 'cart':
            figsize = [12, 8]
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        power = self.directivity
        if norm_max:
            power = np.moveaxis(np.moveaxis(power, [0, 1], [2, 3]) / np.max(power, axis=(2, 3)), [0, 1], [2, 3])
        if dB:
            power = 10 * np.log10(power)
            cb_label = 'dB'
        else:
            cb_label = 'Intensity (arbitrary units)'
            if norm_max:
                'Intensity (max normalized)'
        vmin = np.min(power[i_pol, i_freq, :, :])
        vmax = np.max(power[i_pol, i_freq, :, :])

        if mode == 'uv':
            image2d = ax.pcolormesh(u[:, to_plot], v[:, to_plot], power[i_pol, i_freq][:, to_plot], shading='gouraud',
                                    vmin=vmin, vmax=vmax)
            ax.set_aspect('equal')
            ax.set_xlim([-1, 1])
            ax.set_ylim([-1, 1])
            ax.set_xlabel('u')
            ax.set_ylabel('v')
        if mode == 'cart':
            image2d = ax.imshow(power[i_pol, i_freq], vmin=vmin, vmax=vmax, extent=[0, 180, 0, 360])

            ax.set_aspect('auto')
            ax.set_xlabel(r'$\theta$ (deg)')
            ax.set_ylabel(r'$\phi$ (deg)')
        fig.colorbar(image2d, label=cb_label, ax=ax)
        return fig, ax

    def get_e_A(self, r):
        """
        Get the aperture efficiency.

        This function gets the aperture efficiency of the beam from the max gain.
        """
        A_p = np.pi * r ** 2
        return self.A_e / A_p

    def get_fractional_power(self, theta_min=0, theta_max=180, phi_min=0, phi_max=360):
        """
        Get the fractional power.

        This function gets the power within a solid angle of the beam, normalized
        by the total output power.
        """
        scaling = np.abs(np.sin(self.theta * np.pi / 180))
        total_power = np.sum(scaling * self.directivity, axis=(2, 3))
        i_phi_min = int(phi_min // self.phi_step)
        i_phi_max = int(phi_max // self.phi_step)
        i_theta_min = int(theta_min // self.theta_step)
        i_theta_max = int(theta_max // self.theta_step)
        power_within_angles = np.sum((scaling * self.directivity)[:, :, i_phi_min:i_phi_max, i_theta_min:i_theta_max],
                                     axis=(2, 3))
        return power_within_angles / total_power

    def get_beamwidth(self, phi_cut, dB_threshold=3, deg=True):
        """
        Get the 3dB beamwidth.

        This function gets the 3dB beamwidth at a given phi cut.
        """
        i_phi_cut = int(phi_cut // self.phi_step)
        i_phi_cut_2 = int((phi_cut + 180) // self.phi_step)

        power = 10 * np.log10(self.directivity)
        beam_max = np.amax(power, axis=(2, 3))
        power_subtracted = np.moveaxis(
            np.subtract(np.moveaxis(power[:, :, i_phi_cut], 2, 0), (beam_max - dB_threshold)), 0, 2)
        i_beamwidth_right = np.argmax(power_subtracted <= 0, axis=2)
        power_subtracted = np.moveaxis(
            np.subtract(np.moveaxis(power[:, :, i_phi_cut_2], 2, 0), (beam_max - dB_threshold)), 0, 2)
        i_beamwidth_left = np.argmax(power_subtracted <= 0, axis=2)
        beamwidth = np.zeros([2, self.freqs.shape[0]])
        for i_pol in [0, 1]:
            beamwidth_right = self.theta[i_phi_cut, :][i_beamwidth_right[i_pol, :]]
            beamwidth_left = self.theta[i_phi_cut, :][i_beamwidth_left[i_pol, :]]
            beamwidth[i_pol, :] = beamwidth_right + beamwidth_left
        if not deg:
            beamwidth *= (np.pi / 180)
        return beamwidth