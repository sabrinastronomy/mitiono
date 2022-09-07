#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 17:02:11 2022

@author: vincent
"""

import matplotlib.pyplot as plt
import numpy as np

class CSTBeam():
    """
    This class take a folder where a beam has been converted to the
    "all_*.npy" formats, loads it, an 
    """
    def __init__(self,beams_folder, zenith_rot = 0,ns_rot = 0, ew_rot = 0, load_comps=False,load_axratios=False):
        self.beams_dir = np.load(beams_folder+'all_beams_dir.npy')
        self.freqs = np.load(beams_folder+'all_freqs.npy')
        self.phi = np.load(beams_folder + 'all_phi.npy')
        self.theta = np.load(beams_folder + 'all_theta.npy')
        self.phi_step = self.phi[1][0] # = 0.5 (deg)
        self.theta_step = self.theta[0][1]
        self.freq_min = self.freqs[0]
        self.freq_step = self.freqs[1] - self.freqs[0]
        if load_comps:
            self.thetaphi_comps = np.load(beams_folder + 'all_thetaphi_comps.npy')
        if load_axratios:
            self.axratios = np.load(beams_folder + 'all_axratios.npy')
        
        # Make sure the rotations are at the right resolution
        ns_rot = self.theta_step*np.round(ns_rot/self.theta_step)
        ew_rot = self.theta_step*np.round(ew_rot/self.theta_step)
        zenith_rot = self.phi_step*np.round(zenith_rot/self.phi_step)
        
        # I do the rotations using this method:
        # https://stla.github.io/stlapblog/posts/RotationSphericalCoordinates.html
        
        a_x = ew_rot * (np.pi/180)
        a_y = ns_rot * (np.pi/180)
        a_z = zenith_rot * (np.pi/180)
                
        R_x = np.array([ [np.cos(a_x/2) , -1j*np.sin(a_x/2)] , [-1j*np.sin(a_x/2) , np.cos(a_x/2)] ])
        R_y = np.array([ [ np.cos(a_y/2) , -np.sin(a_y/2) ] , [ np.sin(a_y/2) , np.cos(a_y/2) ] ])
        R_z = np.array([ [ np.exp(-1j * a_z/2) , 0 ] , [ 0 , np.exp(1j * a_z/2) ] ])


        
        for R in [R_z,R_x,R_y]: # this determines the order of the rotations
            t = self.theta * (np.pi/180)
            p = self.phi * (np.pi/180)
            psi = np.array([np.cos(t/2), np.exp(1j * p) * np.sin(t/2)])
            psi[0] = R[0,0] * psi[0] + R[0,1] * psi[1]
            psi[1] = R[1,0] * psi[0] + R[1,1] * psi[1]
            new_theta = (180/np.pi) * 2*np.arctan(np.abs(psi[1])/np.abs(psi[0]))
            new_phi = (180/np.pi) * (np.angle(psi[1]) - np.angle(psi[0]))
            # Reindex
            i_new_theta = (new_theta / self.theta_step).astype('int')
            i_new_phi = (new_phi / self.phi_step).astype('int')
            self.beams_dir = self.beams_dir[:,:,i_new_phi,i_new_theta]
        

         
    def plot_1d_cart(self,phi_cut,i_freq,pol=0, dB = True):
        """
        This function plots the beam in 1D in cartesian coordinates,
        along the cut phi_cut, at frequency index i_freq, for polarization pol.
        """
        line_props = {
        'linestyle': '-',
        'color': 'r',
        'linewidth': 2,
        'alpha': 1
        }
        
        i_phi_cut = int(phi_cut//self.phi_step)
        i_phi_cut_2 = int((phi_cut+180) // self.phi_step)
        fig,ax = plt.subplots(1,1,figsize=[20,10])
        ylabel = r'$I(\theta)$ (arbitrary units)'
        if pol==-1:
            #print('Plotting (pol_0 ^ 2 + pol_1 ^ 2) ^ 0.5')
            power = (self.beams_dir[0,:,:,:] ** 2. + self.beams_dir[1,:,:,:] ** 2.) ** 0.5
            if dB:
                power = 10 * np.log10(power)
                ylabel = r'$I(\theta)$ (dB)'
            ax.plot(self.theta[i_phi_cut,:] , power[i_freq][i_phi_cut][:], **line_props)
            # Need to plot the other side too, but the beam is perfectly symmetric in my simulation
            # so this is somewhat redundant  information.
            ax.plot(-self.theta[i_phi_cut_2] , power[i_freq][i_phi_cut][:], **line_props)
        else:
            #print('Plotting pol_{:d}'.format(pol))
            power = self.beams_dir
            if dB:
                power = 10 * np.log10(power)
                ylabel = r'$I(\theta)$ (dB)'
            ax.plot(self.theta[i_phi_cut,:] , power[pol,i_freq][i_phi_cut][:], **line_props)
            # Need to plot the other side too, but the beam is perfectly symmetric in my simulation
            # so this is somewhat redundant  information.
            ax.plot(-self.theta[i_phi_cut_2] , power[pol,i_freq][i_phi_cut][:], **line_props)
            
        ax.set_xlabel('Theta (deg)')
        ax.set_ylabel(ylabel)
        ax.set_title('Freq = {:.1f} GHz, $\phi$ = {:.1f}$^\circ$'.format(self.freqs[i_freq],self.phi[i_phi_cut][0]))
        ax.set_xlim([-20,20])

    def plot_1d_polar(self,phi_cut,i_freq,pol=0, dB = True):
        """
        This function plots the beam in 1D in polar coordinates,
        along the cut phi_cut, at frequency index i_freq, for polarization pol.
        """
        line_props = {
        'linestyle': '-',
        'color': 'r',
        'linewidth': 2,
        'alpha': 1
        }
        
        i_phi_cut = int(phi_cut//self.phi_step)
        i_phi_cut_2 = int((phi_cut+180) // self.phi_step)
        fig, ax = plt.subplots(1,1,figsize=[15,15],subplot_kw={'projection': 'polar'})
        ax.set_theta_zero_location("N")
        ylabel = r'$I(\theta)$ (arbitrary units)'
        if pol==-1:
            #print('Plotting (pol_0 ^ 2 + pol_1 ^ 2) ^ 0.5')
            power = (self.beams_dir[0,:,:,:] ** 2. + self.beams_dir[1,:,:,:] ** 2.) ** 0.5
            if dB:
                power = 10 * np.log10(power)
                ylabel = r'$I(\theta)$ (dB)'
            ax.plot((np.pi / 180)*self.theta[i_phi_cut,:] , power[i_freq][i_phi_cut][:], **line_props)
            # Need to plot the other side too, but the beam is perfectly symmetric in my simulation
            # so this is somewhat redundant  information.
            ax.plot(-(np.pi / 180)*self.theta[i_phi_cut_2] , power[i_freq][i_phi_cut][:], **line_props)
        else:
            #print('Plotting pol_{:d}'.format(pol))
            power = self.beams_dir
            if dB:
                power = 10 * np.log10(power)
                ylabel = r'$I(\theta)$ (dB)'
            ax.plot((np.pi / 180)*self.theta[i_phi_cut,:] , power[pol,i_freq][i_phi_cut][:], **line_props)
            # Need to plot the other side too, but the beam is perfectly symmetric in my simulation
            # so this is somewhat redundant  information.
            ax.plot(-(np.pi / 180)*self.theta[i_phi_cut_2] , power[pol,i_freq][i_phi_cut][:], **line_props)
            
        ax.set_xlabel('Theta (deg)')
        ax.set_ylabel(ylabel)
        ax.set_title('Freq = {:.1f} GHz, $\phi$ = {:.1f}$^\circ$'.format(self.freqs[i_freq],self.phi[i_phi_cut][0]))
        
    def plot_2d_uv(self, i_freq, i_pol, dB = True, front_cutoff = 90, back_cutoff = -1):
        """
        This function plots the beam in a 2D uv plot.
        I need to check whether my definition of uv coordinates is good.
        """
        u = np.sin(self.theta*np.pi/180) * np.cos(self.phi*np.pi/180)
        v = np.sin(self.theta*np.pi/180) * np.sin(self.phi*np.pi/180)
        if back_cutoff == -1:
            to_plot = np.unique(np.where(self.theta<=front_cutoff)[1])
        else:
            to_plot = np.unique(np.where(self.theta>=back_cutoff)[1])

        

        fig,ax = plt.subplots(1,1,figsize=[10,10])
        if dB:
            vmin = np.min(10*np.log10(self.beams_dir[i_pol,i_freq,:,:]))
            vmax = np.max(10*np.log10(self.beams_dir[i_pol,i_freq,:,:]))
            cb_label = 'dB'
            image2d = ax.pcolormesh(u[:,to_plot],v[:,to_plot],10*np.log10(self.beams_dir[i_pol,i_freq][:,to_plot]),shading='auto',vmin = vmin, vmax = vmax)
        else:
            vmin = np.min(self.beams_dir[i_pol,i_freq,:,:])
            vmax = np.max(self.beams_dir[i_pol,i_freq,:,:])
            cb_label = 'Intensity (arbitrary units)'
            image2d = ax.pcolormesh(u[:,to_plot],v[:,to_plot],self.beams_dir[i_pol,i_freq][:,to_plot],shading='auto',vmin = vmin, vmax = vmax)
        
        ax.set_aspect('equal')
        ax.set_xlim([-1,1])
        ax.set_ylim([-1,1])
        ax.set_xlabel('u')
        ax.set_ylabel('v')
        fig.colorbar(image2d,label = cb_label, ax=ax)
    
    def plot_2d_cart(self, i_freq, i_pol, dB = True):
        """
        This function plots the beam in a 2D cartesian plot.
        """
        fig,ax = plt.subplots(1,1,figsize=[20,10])
        
        if dB:
            vmin = np.min(10*np.log10(self.beams_dir[i_pol,i_freq,:,:]))
            vmax = np.max(10*np.log10(self.beams_dir[i_pol,i_freq,:,:]))
            cb_label = 'dB'
            image2d = ax.imshow(10*np.log10(self.beams_dir[i_pol,i_freq]),vmin = vmin, vmax = vmax, extent = [0 , 180, 0 , 360])
        else:
            vmin = np.min(self.beams_dir[i_pol,i_freq,:,:])
            vmax = np.max(self.beams_dir[i_pol,i_freq,:,:])
            cb_label = 'Intensity (arbitrary units)'
            image2d = ax.imshow(self.beams_dir[i_pol,i_freq],vmin = vmin, vmax = vmax, extent = [0 , 180, 0 , 360])
        
            
        ax.set_aspect('auto')
        ax.set_xlabel(r'$\theta$')
        ax.set_ylabel(r'$\phi$')
        fig.colorbar(image2d, label=cb_label, ax=ax)
