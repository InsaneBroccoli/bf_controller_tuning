%==========================================================================
% Betaflight Controller Tuning Analysis Script
% Purpose: Analyzes flight logs and tunes PID controllers for a quadcopter
%
% Author: [Your Name]
% Date: [Current Date]
%==========================================================================

%% General

% Initialize workspace
clc, clear variables, close all
addpath(genpath('lib'));
addpath logs/
addpath class/

% Show Legends
do_insert_legends = true;

%% Load Data

% =========================================================================
%  Define path to *.bbl.csv file for the 
% =========================================================================

log_folder = '';
flight_folder = '20250907';
log_name = '20250907_apex5_00.bbl.csv';
file_path = fullfile(log_folder, flight_folder, log_name);

data_flight = flight_data(file_path);
% Load Data
data_flight = data_flight.get_data();

gyro_tuning = gyro_ctrl_tuning(data_flight.data, data_flight.ind, data_flight.Ts_log, ...
    data_flight.para, data_flight.Ts_cntr);

resolution_factor_tf = 2;    % Window length for spectral analysis (seconds)
overlap_tf = 0.9;              % Overlap factor for spectral analysis (0-1)
gyro_tuning = gyro_tuning.calculate_transfer_func(resolution_factor_tf, overlap_tf);

%% Flight Analyser
analysis_flight = flight_analyser(data_flight.data, data_flight.ind, data_flight.Ts_log);

% Data for Spectra
resolution_factor_spectra = 2;    % Window length for spectral analysis (seconds)
overlap_spectra = 0.9;              % Overlap factor for spectral analysis (0-1)
analysis_flight = analysis_flight.calculate_spectra(resolution_factor_spectra, ...
    overlap_spectra);

% Data for Spectogram
resolution_factor_spectogram = 0.2;    % Window length for spectral analysis (seconds)
overlap_spectogram = 0.9;              % Overlap factor for spectral analysis (0-1)
analysis_flight = analysis_flight.calculate_spectogram(resolution_factor_spectogram, ...
    overlap_spectogram);

%% Tuninig Data

% =========================================================================
%  Axis Selection: 1: roll, 2: pitch, 3: yaw
% =========================================================================

ind_ax = 1;     % keep it now until plot_utils is finished

% I-term Relax on/off
do_compensate_iterm = true;

% Transfer Function
do_show_transferfunction = false; % Plant, PI & D Controller, PID Controller 
Nestfatra = 2.0;                 % Window length for spectral analysis (seconds)
koverlaptra = 0.9;               % Overlap factor for spectral analysis (0-1) 

% New and old parameters are the same
default_parameters = false; 

% =========================================================================
%  First flight: Parameters
% =========================================================================

% Configure filters to match Betaflight settings
% All frequencies are in Hz
% Filter types:
%   0: PT1 (First order lowpass)
%   1: BIQUAD (Second order)
%   2: PT2 (Second order lowpass) 
%   3: PT3 (Third order lowpass)
para_new.gyro_lpf            = 0;       % dono what this is
para_new.gyro_lowpass_hz     = 0;       % frequency of gyro lpf 1
para_new.gyro_soft_type      = 0;       % type of gyro lpf 1
para_new.gyro_lowpass_dyn_hz = [0, 0];  % dyn gyro lpf overwrites gyro_lowpass_hz
para_new.gyro_lowpass2_hz    = 800;     % frequency of gyro lpf 2
para_new.gyro_soft2_type     = 0;       % type of gyro lpf 2
para_new.gyro_notch_hz       = [0, 520];  % frequency of gyro notch 1 and 2
para_new.gyro_notch_cutoff   = [0, 448] % % Cutoff frequency gyro notch 1 and 2
para_new.dterm_lpf_hz        = 0;       % frequency of dterm lpf 1
para_new.dterm_filter_type   = 0;       % type of dterm lpf 1
para_new.dterm_lpf_dyn_hz    = [0, 0];  % dyn dterm lpf overwrites dterm_lpf_hz
para_new.dterm_lpf2_hz       = 130;     % frequency of dterm lpf 2
para_new.dterm_filter2_type  = 3;       % type of dterm lpf 2
para_new.dterm_notch_hz      = 235;       % frequency of dterm notch
para_new.dterm_notch_cutoff  = 202;     % Cutoff frequency dterm notch
para_new.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)

%--------------------------------------------------------------------------
%  tune your PIDs here
%--------------------------------------------------------------------------
% Adjust multipliers (1.0 *) to tune the response
% Higher P: more immediate response but possible oscillations
% Higher I: better steady-state tracking but slower response
% Higher D: better damping but noise sensitive

switch ind_ax
    case 1 % roll: [45, 80, 30, 0] - Betaflight standard PIDs
        P_new       = 1.0 * 46;
        I_new       = 1.0 * 74;
        D_new       = 1.0 * 30;
    case 2 % pitch: [47, 84, 34, 0] - Betaflight standard PIDs
        P_new       = 1.0 * 50;
        I_new       = 1.0 * 88;
        D_new       = 1.0 * 43;
    case 3 % yaw: [45, 80, 0, 0] - Betaflight standard PIDs
        P_new       = 1.0 * 35;
        I_new       = 1.0 * 70;
        D_new       = 1.0 * 3;
end

gyro_tuning = gyro_tuning.calculate_new_controller(ind_ax, P_new, I_new, D_new, ...
    default_parameters, para_new);
gyro_tuning = gyro_tuning.get_tuning_data(do_compensate_iterm);

plotter = plot_utils(data_flight, gyro_tuning, analysis_flight);

%% Plot Gyro Data
plotter.plot_Gyro_Signal(do_insert_legends);
plotter.plot_Overview(do_insert_legends);
plotter.plot_Eval_Time();

%% Plot Flight Analyser

plotter.plot_Gyro_spectra(do_insert_legends);
plotter.plot_Spectogram(3);

%% Plot Tuning Data
plotter.plot_Bode_Plant(ind_ax);
plotter.plot_Bode_Contr(ind_ax, do_insert_legends);
plotter.plot_Gang_of_Four(do_insert_legends);
plotter.plot_Step_Response(do_insert_legends);