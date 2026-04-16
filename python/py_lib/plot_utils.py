"""
==========================================================================
PLOT_UTILS - Betaflight Controller Analysis Visualization Class
==========================================================================
Purpose: 
This class provides visualization methods for analyzing quadcopter 
flight controller performance and tuning data.

Author: [Janick Dort, Yuri Bianchi, Dario Jurietti]
Supervisor: [Michael Peter]
Date: [28.02.2026]
==========================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import control as ct
from typing import List, Dict, Optional
from matplotlib.gridspec import GridSpec


class PlotUtils:
    """Visualization utilities for flight analysis."""
    
    def __init__(self, legend: bool = True):
        """
        Initialize plot utilities.
        
        Parameters
        ----------
        legend : bool
            Whether to show legends on plots
        """
        self.do_insert_legends = legend
        self.linewidth = 1.0
        self.axis_names = ['Roll', 'Pitch', 'Yaw']
        
        # Set matplotlib style
        plt.style.use('seaborn-v0_8-darkgrid')
        plt.rcParams['figure.figsize'] = (12, 8)
        plt.rcParams['lines.linewidth'] = self.linewidth
    
    # ======================================================================
    # Figure Flight Data
    # ======================================================================
    def plot_flight_data(self, df, groups: List[Dict], fig_title: str = ''):
        """   
        Parameters
        ----------
        df : FlightData
            Flight data object
        groups : List[Dict]
            List of subplot specifications
        fig_title : str
            Figure title
        """
        n_subplots = len(groups)
        fig, axes = plt.subplots(n_subplots, 1, figsize=(8, 6))
        fig.canvas.manager.set_window_title(fig_title)
        
        if n_subplots == 1:
            axes = [axes]
        
        if fig_title:
            fig.suptitle(fig_title, fontsize=14, fontweight='bold')
        
        for i, (ax, group) in enumerate(zip(axes, groups)):
            idx_list = group['idx']
            
            # Get legend entries
            if 'legend' in group:
                legend_entries = group['legend']
            else:
                legend_entries = [f'Signal {k}' for k in range(len(idx_list))]
            
            # Plot data
            for k, idx in enumerate(idx_list):
                ax.plot(df.time, df.data[:, idx], label=legend_entries[k])
            
            ax.grid(True)
            ax.set_xlim([df.time[0], df.time[-1]])
            ax.set_ylabel(group['ylabel'])
            
            if self.do_insert_legends:
                ax.legend(loc='upper right')
            
            if i == n_subplots - 1:
                ax.set_xlabel('Time [s]')
        
        plt.tight_layout()
    
    # ======================================================================
    # Figure Evaluation Time
    # ======================================================================
    def plot_eval_time(self, time: np.ndarray):
        """
        Parameters
        ----------
        time : np.ndarray
            Time vector
        """
        delta_time_us = np.diff(time) * 1.0e6
        
        fig, ax = plt.subplots(figsize=(8, 6))
        fig.canvas.manager.set_window_title('Evaluation Time Analysis')
        
        ax.plot(time[:-1], delta_time_us)
        
        # Mean line
        mean_val = np.mean(delta_time_us)
        ax.axhline(mean_val, color='k', linestyle='--', linewidth=self.linewidth)
        ax.text(time[0], mean_val, 'Mean',
                ha='left', va='bottom',
                fontsize=9, color='k')
        
        ax.grid(True)
        ax.set_title(f'Mean: {np.mean(delta_time_us):.2f} μs, '
                    f'Median: {np.median(delta_time_us):.2f} μs, '
                    f'Std: {np.std(delta_time_us):.2f} μs')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Ts log (μs)')
        ax.set_xlim([0, time[-1]])
        
        plt.tight_layout()
    
    # ======================================================================
    # Figure Spectra
    # ======================================================================
    def plot_gyro_spectra(self, fa):
        """
        Parameters
        ----------
        fa : FlightAnalyzer
            Flight analyzer object
        """
        fig, axes = plt.subplots(3, 1, figsize=(8, 6))
        fig.canvas.manager.set_window_title('Gyro Spectra')
        
        # Unfiltered gyro spectra
        axes[0].plot(fa.freq_spectra, fa.spectra[:, 0:3])
        axes[0].grid(True, which='both', alpha=0.3)
        axes[0].set_ylabel('Gyro (deg/s)')
        axes[0].set_yscale('log')
        axes[0].set_title('Unfiltered gyro magnitude spectra')
        axes[0].legend(['Roll', 'Pitch', 'Yaw'], loc='upper right')
        
        # Filtered gyro spectra
        axes[1].plot(fa.freq_spectra, fa.spectra[:, 3:6])
        axes[1].grid(True, which='both', alpha=0.3)
        axes[1].set_ylabel('Gyro (deg/s)')
        axes[1].set_yscale('log')
        axes[1].set_title('Filtered (ADC) gyro magnitude spectra')
        
        if self.do_insert_legends:
            axes[1].legend(['Roll', 'Pitch', 'Yaw'], loc='upper right')
        
        # Axis sum spectra
        axes[2].plot(fa.freq_spectra, fa.spectra[:, 6:9])
        axes[2].grid(True, which='both', alpha=0.3)
        axes[2].set_ylabel('AxisSum')
        axes[2].set_xlabel('Frequency (Hz)')
        axes[2].set_yscale('log')
        axes[2].set_title('Axis sum spectra')
        axes[2].legend(['Roll', 'Pitch', 'Yaw'], loc='upper right')
        
        # Set limits
        nyq = 1 / (2 * fa.Ts_log)
        for ax in axes:
            ax.set_xlim([0, nyq])
            try:
                # ax.set_ylim([1e-3, 1e1])
                ax.autoscale(enable=True, axis='y', tight=False)
            except:
                pass
        
        plt.tight_layout()
    
    # ======================================================================
    # Figure Spectrogram
    # ======================================================================
    def plot_spectrogram(self, fa, num_spectrograms: int = 3):
        """
        Parameters
        ----------
        fa : FlightAnalyzer
            Flight analyzer object
        num_spectrograms : int
            Number of axes to plot
        """
        fig = plt.figure(figsize=(8, 6))
        fig.canvas.manager.set_window_title('Gyro Spectrograms')
        axes_labels = ['Roll', 'Pitch', 'Yaw']
        c_lim = [5e-2, 3e0]
        
        # Unfiltered spectrograms
        for i in range(num_spectrograms):
            ax = plt.subplot(2, 3, i + 1)
            
            im = ax.pcolormesh(
                fa.freq_spectrogram[i],
                fa.throttle_all[i],
                fa.spectrogram_unf[i],
                shading='auto',
                cmap='jet',
                norm=plt.matplotlib.colors.LogNorm(vmin=c_lim[0], vmax=c_lim[1])
            )
            
            if i == 0:
                ax.set_ylabel('Throttle (%)')
            ax.set_xlabel('Frequency (Hz)')
            ax.set_title(f'{axes_labels[i]} – without filter')
            ax.set_ylim([0, 100])
            plt.colorbar(im, ax=ax)
        
        # Filtered spectrograms
        for i in range(num_spectrograms):
            ax = plt.subplot(2, 3, i + 4)
            
            im = ax.pcolormesh(
                fa.freq_spectrogram[i],
                fa.throttle_all[i],
                fa.spectrogram_fil[i],
                shading='auto',
                cmap='jet',
                norm=plt.matplotlib.colors.LogNorm(vmin=c_lim[0], vmax=c_lim[1])
            )
            
            if i == 0:
                ax.set_ylabel('Throttle (%)')
            ax.set_xlabel('Frequency (Hz)')
            ax.set_title(f'{axes_labels[i]} – with filter')
            ax.set_ylim([0, 100])
            plt.colorbar(im, ax=ax)
        
        plt.suptitle('Gyro Spectrograms', fontsize=14, fontweight='bold')
        plt.tight_layout()
    
    # ======================================================================
    # Figure Bode Plot Plant Coherence
    # ======================================================================
    def plot_bode_plant(self, td, ind_ax: int, custom_label: str = 'Signal'):
        """
        Parameters
        ----------
        td : GyroCtrlTuning
            Tuning data object
        ind_ax : int
            Axis index
        custom_label : str
            Label for plot title
        """
        axis_name = self.axis_names[ind_ax]
        
        fig = plt.figure(figsize=(8, 6))
        fig.canvas.manager.set_window_title(f'Bode Plot {custom_label} - {axis_name}')
        gs = GridSpec(3, 1, height_ratios=[2, 2, 1])
        
        # Magnitude plot
        ax1 = fig.add_subplot(gs[0])
        mag, phase, omega = ct.bode(td.P[ind_ax], plot=False)
        freq_hz = omega / (2 * np.pi)
        ax1.semilogx(freq_hz, 20 * np.log10(mag), 'k', linewidth=self.linewidth)
        ax1.grid(True, which='both', alpha=0.3)
        ax1.set_ylabel('Magnitude (dB)')
        ax1.set_title(f'Plant P ({custom_label}) - {axis_name}')
        
        # Phase plot
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        phase_deg = phase * 180 / np.pi
        phase_wrapped = np.mod(phase_deg + 180, 360) - 180  # Wrap to [-180, 180]
        ax2.semilogx(freq_hz, phase_wrapped, 'k', linewidth=self.linewidth)
        ax2.grid(True, which='both', alpha=0.3)
        ax2.set_ylabel('Phase (deg)')
        
        # Coherence plot
        ax3 = fig.add_subplot(gs[2], sharex=ax1)
        coh_data = np.abs(np.squeeze(td.Coh[ind_ax].fresp))
        coh_freq = td.Coh[ind_ax].omega / (2 * np.pi)
        ax3.semilogx(coh_freq, coh_data, 'k', linewidth=self.linewidth)
        ax3.grid(True, which='both', alpha=0.3)
        ax3.set_ylabel('Coherence')
        ax3.set_ylim([0, 1])
        ax3.set_xlabel('Frequency (Hz)')
        
        plt.tight_layout()
    
    # ======================================================================
    # Figure Gang of Four
    # ======================================================================
    def plot_gang_of_four(self, td, label: str = 'Signal'):
        """
        Parameters
        ----------
        td : GyroCtrlTuning
            Tuning data object
        label : str
            Label for plots
        """
        ind_ax = td.ind_ax
        axis_name = self.axis_names[ind_ax]

        fig, axes = plt.subplots(2, 2, figsize=(8, 6))
        fig.canvas.manager.set_window_title(f'Gang of Four - {axis_name}')
        fig.suptitle(f'Gang of Four {label} - {axis_name}', fontsize=14, fontweight='bold')

        def plot_frd(ax, systems, title, xlabel=None):
            for sys, style, lbl in systems:
                mag, _, omega = ct.bode(sys, plot=False)
                freq_hz = omega / (2 * np.pi)
                mag_db = 20 * np.log10(mag.squeeze())
                ax.semilogx(freq_hz, mag_db, style, linewidth=self.linewidth, label=lbl)
            ax.grid(True, which='both', alpha=0.3)
            ax.set_ylabel('Magnitude (dB)')
            ax.set_title(title)
            if xlabel:
                ax.set_xlabel(xlabel)
            if self.do_insert_legends:
                ax.legend(loc='best')

        # Tracking T
        plot_frd(axes[0, 0], [
            (td.CloLoAan.T,    '--', 'actual'),
            (td.CloLoAanNew.T, '-',  'new'),
            (td.T[ind_ax],     ':',  'measured'),
        ], title='Tracking T')

        # Sensitivity S
        plot_frd(axes[0, 1], [
            (td.CloLoAan.S,    '--', 'actual'),
            (td.CloLoAanNew.S, '-',  'new'),
        ], title='Sensitivity S')

        # Controller Effort SC
        plot_frd(axes[1, 0], [
            (td.CloLoAan.SC,    '--', 'actual'),
            (td.CloLoAanNew.SC, '-',  'new'),
        ], title='Controller Effort SC', xlabel='Frequency (Hz)')

        # Compliance SP
        plot_frd(axes[1, 1], [
            (td.CloLoAan.SP,    '--', 'actual'),
            (td.CloLoAanNew.SP, '-',  'new'),
        ], title='Compliance SP', xlabel='Frequency (Hz)')

        plt.tight_layout()

    
    # ======================================================================
    # Figure Step Response
    # ======================================================================
    def plot_step_response(self, td, label: str = 'Signal', ylab: str = 'Response'):
        """
        Parameters
        ----------
        td : GyroCtrlTuning
            Tuning data object
        label : str
            Label for plots
        ylab : str
            Y-axis label
        """
        axis_name = self.axis_names[td.ind_ax]
        
        fig, axes = plt.subplots(2, 1, figsize=(8, 6))
        fig.canvas.manager.set_window_title(f'Step Response - {axis_name}')
        fig.suptitle(f'Step Response {label} - {axis_name}', fontsize=14, fontweight='bold')
        
        # Tracking
        ax = axes[0]
        ax.plot(td.step_time, td.step_resp_tra[:, 0], '--', label='actual calculated')
        ax.plot(td.step_time, td.step_resp_tra[:, 1], '-', label='new calculated')
        ax.plot(td.step_time, td.step_resp_tra[:, 2], ':', label='measured')
        ax.grid(True)
        ax.set_ylabel(ylab)
        ax.set_title(f'Tracking T - {axis_name}')
        if self.do_insert_legends:
            ax.legend(loc='best')
        ax.set_xlim([0, 0.5])
        
        # Compliance
        ax = axes[1]
        ax.plot(td.step_time, td.step_resp_com[:, 0], '--', label='actual calculated')
        ax.plot(td.step_time, td.step_resp_com[:, 1], '-', label='new calculated')
        ax.grid(True)
        ax.set_ylabel(ylab)
        ax.set_xlabel('Time (s)')
        ax.set_title('Compliance SP')
        if self.do_insert_legends:
            ax.legend(loc='best')
        ax.set_xlim([0, 0.5])
        
        plt.tight_layout()
        plt.show()