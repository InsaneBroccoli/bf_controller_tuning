"""
==========================================================================
ANGLE CTRL TUNING - Betaflight Controller Analysis (Angle tuning class)
==========================================================================
Betaflight Controller Tuning Analysis Script
Purpose: Calculations for angle tuning

Author: [Janick Dort, Yuri Bianchi, Dario Jurietti]
Supervisor: [Michael Peter]
Date: [28.02.2026]
==========================================================================
"""

import numpy as np
import control as ct
import copy

from scipy.signal.windows import hann
from typing import Optional, List

from .gyro_ctrl_tuning import get_ind_eval

from .pidtuninglib import (
    header_info,
    ClosedLoop,
    apply_rotfiltfilt,
    estimate_frequency_response,
    downsample_frd,
    calculate_step_response_from_frd,
    calculate_closed_loop_angle,
    get_filter,
)


class AngleCtrlTuning:
    """Angle controller tuning and analysis."""

    def __init__(
        self,
        data: np.ndarray,
        ind,
        Ts_log: float,
        para: header_info,
        Ts_cntr: float,
        gyro_tuning
    ):

        self.data = data
        self.ind = ind
        self.Ts_log = Ts_log
        self.para = para
        self.Ts_cntr = Ts_cntr
        self.gyro_tuning = gyro_tuning

        self.TargetAngle = None

        # Transfer functions
        self.T: List[ct.FRD] = []
        self.P: List[ct.FRD] = []
        self.C: List[ct.FRD] = []
        self.T_gy: List[ct.FRD] = []

        # Coherence
        self.Coh_T: List[ct.FRD] = []
        self.Coh_P: List[ct.FRD] = []
        self.Coh_C: List[ct.FRD] = []

        # Controller
        self.C_ana: Optional[ct.FRD] = None
        self.C_ana_new: Optional[ct.FRD] = None

        # Axis
        self.ind_ax: Optional[int] = None

        # Helpers
        self.omega_bode: Optional[np.ndarray] = None
        self.Nest: Optional[int] = None
        self.step_time: Optional[np.ndarray] = None

        # Closed-loop results
        self.CL_ana: Optional[ClosedLoop] = None
        self.CL_ana_new: Optional[ClosedLoop] = None

        # Step response
        self.step_resp_tra: Optional[np.ndarray] = None
        self.step_resp_com: Optional[np.ndarray] = None

    def calculate_angle_trans(
        self,
        Nestfatra: float,
        koverlaptra: float
    ) -> "AngleCtrlTuning":

        # =============================================================
        # Analysis window parameters
        # =============================================================

        self.Nest = int(np.round(Nestfatra / self.Ts_log))
        Noverlap = int(np.floor(koverlaptra * self.Nest))
        window = hann(self.Nest, sym=False)

        # =============================================================
        # Excitation filter
        # =============================================================

        Dlp = np.sqrt(3) / 2
        wlp = 2 * np.pi * 10

        Glp_cont = ct.tf(
            [wlp**2],
            [1, 2 * Dlp * wlp, wlp**2]
        )

        Glp = ct.c2d(
            Glp_cont,
            self.Ts_log,
            method="tustin"
        )

        # =============================================================
        # Preallocate, angle controller only roll and pitch
        # =============================================================

        n_axes = 2

        self.T = [None] * n_axes
        self.P = [None] * n_axes
        self.C = [None] * n_axes
        self.T_gy = [None] * n_axes

        self.Coh_T = [None] * n_axes
        self.Coh_P = [None] * n_axes
        self.Coh_C = [None] * n_axes

        sinarg_full = self.data[:, self.ind.sinarg].copy()

        # =============================================================
        # Process roll and pitch
        # =============================================================

        for ind_axis in range(n_axes):

            ind_eval = get_ind_eval(
                self.data[:, self.ind.sinarg],
                self.data[:, self.ind.gyroADC[ind_axis]]
            )

            sinarg_ax = sinarg_full.copy()
            sinarg_ax[~ind_eval] = 0

            # ---------------------------------------------------------
            # Input signal: Target Angle
            # ---------------------------------------------------------

            w = self.data[:, self.ind.angleTarget[ind_axis]]

            inp = apply_rotfiltfilt(
                Glp,
                sinarg_ax,
                w[:, np.newaxis]
            )

            # ---------------------------------------------------------
            # Output signal: Current Angle
            # ---------------------------------------------------------

            y = self.data[:, self.ind.currentAngle[ind_axis]]

            out_y = apply_rotfiltfilt(
                Glp,
                sinarg_ax,
                y[:, np.newaxis]
            )

            T_ax, C_T_ax, _, _ = estimate_frequency_response(
                inp[ind_eval].ravel(),
                out_y[ind_eval].ravel(),
                window,
                Noverlap,
                self.Nest,
                self.Ts_log
            )

            # ---------------------------------------------------------
            # Plant calculation: target angle to gyro
            # ---------------------------------------------------------

            v = self.data[:, self.ind.gyroADC[ind_axis]]

            out_v = apply_rotfiltfilt(
                Glp,
                sinarg_ax,
                v[:, np.newaxis]
            )

            G_wv_ax, C_G_wv_ax, _, _ = estimate_frequency_response(
                inp[ind_eval].ravel(),
                out_v[ind_eval].ravel(),
                window,
                Noverlap,
                self.Nest,
                self.Ts_log
            )

            P_angle_ax = T_ax / G_wv_ax

            # ---------------------------------------------------------
            # Controller calculation: target angle to setpoint
            # ---------------------------------------------------------

            c = self.data[:, self.ind.setpoint[ind_axis]]

            out_c = apply_rotfiltfilt(
                Glp,
                sinarg_ax,
                c[:, np.newaxis]
            )

            G_wc_ax, C_G_wc_ax, _, _ = estimate_frequency_response(
                inp[ind_eval].ravel(),
                out_c[ind_eval].ravel(),
                window,
                Noverlap,
                self.Nest,
                self.Ts_log
            )

            Cp_ax = G_wc_ax / (1 - T_ax)

            # ---------------------------------------------------------
            # Inner gyro loop transfer function
            # ---------------------------------------------------------

            T_gy_ax = G_wv_ax / G_wc_ax

            # ---------------------------------------------------------
            # Store results
            # ---------------------------------------------------------

            self.T[ind_axis] = T_ax
            self.P[ind_axis] = P_angle_ax
            self.C[ind_axis] = Cp_ax
            self.T_gy[ind_axis] = T_gy_ax

            self.Coh_T[ind_axis] = C_T_ax
            self.Coh_P[ind_axis] = C_T_ax * C_G_wv_ax
            self.Coh_C[ind_axis] = C_T_ax * C_G_wc_ax

            self.omega_bode = P_angle_ax.omega

        return self

    def calculate_new_controller(
        self,
        ind_ax: int,
        P_new: float,
        default_parameters: bool
    ) -> "AngleCtrlTuning":

        # =============================================================
        # Angle controller filter
        # =============================================================

        angle_lpf_hz = 50

        # =============================================================
        # Use default Betaflight angle P if requested
        # =============================================================

        if default_parameters:
            P_new = self.para.get_list("levelPID")[0]

        P_old = self.para.get_list("levelPID")[0]

        print("   used P parameter Angle Control are:")
        print(f"      P: {P_old}")

        print("   new P parameter Angle Control are:")
        print(f"      P: {P_new}")

        # =============================================================
        # New controller
        # =============================================================

        C_P_Angle = P_new * 0.1

        freq_hz = self.T[ind_ax].omega / (2 * np.pi)
        omega = 2 * np.pi * freq_hz

        response_new = C_P_Angle * np.ones_like(omega)

        C_P_Angle_frd = ct.FRD(
            response_new,
            omega,
            dt=self.Ts_cntr
        )

        Gf_ana = get_filter("pt3", angle_lpf_hz, self.Ts_cntr)[0]
        C_Angle_new = C_P_Angle_frd * Gf_ana

        self.C_ana_new = downsample_frd(
            C_Angle_new,
            self.Ts_log,
            freq_hz
        )

        # =============================================================
        # Old controller
        # =============================================================

        C_P_Angle_old = P_old * 0.1

        response_old = C_P_Angle_old * np.ones_like(omega)

        C_P_Angle_frd_old = ct.FRD(
            response_old,
            omega,
            dt=self.Ts_cntr
        )

        Gf_ana_old = get_filter("pt3", angle_lpf_hz, self.Ts_cntr)[0]
        C_Angle_ana = C_P_Angle_frd_old * Gf_ana_old

        self.C_ana = downsample_frd(
            C_Angle_ana,
            self.Ts_log,
            freq_hz
        )

        self.ind_ax = ind_ax

        return self

    def get_tuning_data(self) -> "AngleCtrlTuning":

        gyro = self.gyro_tuning

        # =============================================================
        # Closed loop calculation
        # =============================================================

        self.CL_ana = calculate_closed_loop_angle(
            self.C_ana,
            gyro.T[self.ind_ax],
            gyro.P[self.ind_ax],
            self.P[self.ind_ax],
            gyro.CL_ana.C
        )

        self.CL_ana_new = calculate_closed_loop_angle(
            self.C_ana_new,
            gyro.T[self.ind_ax],
            gyro.P[self.ind_ax],
            self.P[self.ind_ax],
            gyro.CL_ana.C
        )

        # =============================================================
        # Step response analysis
        # =============================================================

        f_max = 500

        T_mean = (
            0.1 * np.array([-1, 1])
            + (self.Nest * self.Ts_log) / 2
        )

        self.step_time = (
            np.arange(self.Nest)
            * self.Ts_log
        )

        step_resp = np.column_stack([
            calculate_step_response_from_frd(
                copy.deepcopy(self.CL_ana.T),
                f_max
            ),
            calculate_step_response_from_frd(
                copy.deepcopy(self.CL_ana_new.T),
                f_max
            ),
            calculate_step_response_from_frd(
                copy.deepcopy(self.T[self.ind_ax]),
                f_max
            )
        ])

        idx_mean = (
            (self.step_time > T_mean[0])
            & (self.step_time < T_mean[1])
        )

        step_resp_mean = np.mean(
            step_resp[idx_mean, :],
            axis=0
        )

        self.step_resp_tra = step_resp / step_resp_mean

        # =============================================================
        # Disturbance rejection
        # =============================================================

        self.step_resp_com = np.column_stack([
            calculate_step_response_from_frd(
                copy.deepcopy(self.CL_ana.SP),
                f_max
            ),
            calculate_step_response_from_frd(
                copy.deepcopy(self.CL_ana_new.SP),
                f_max
            )
        ])

        step_resp_mean = np.mean(
            self.step_resp_com[idx_mean, :],
            axis=0
        )

        self.step_resp_com = (
            self.step_resp_com
            - step_resp_mean
        )

        return self