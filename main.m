%==========================================================================
% MAIN - Betaflight Controller Analysis Main File
%==========================================================================
% Betaflight Controller Tuning Analysis Script
% Purpose: Analyzes flight logs and tunes PID controllers for a quadcopter
%
% Author: [Janick Dort, Yuri Bianchi, Dario Jurietti]
% Supervisor: [Michael Peter]
% Date: [25.11.2025]
%==========================================================================

%% General

% Initialize workspace
clc, clear variables, close all
addpath(genpath('lib'));
addpath logs/
addpath class/

% Show Legends
do_insert_legends = true;

% Create Plotter Class
plotter = plot_utils(do_insert_legends);

%% Loading Data

% =========================================================================
%  Define path to *.bbl.csv file for the 
% =========================================================================

log_folder = 'logs';
flight_folder = '20260302';
log_name = '20260302_flipmini.csv';
file_path = fullfile(log_folder, flight_folder, log_name);

df = flight_data(file_path);
% Load Data
df = df.get_data();


roll = 1; pitch = 2; yaw = 3;
group1 = {
    struct( ...
        'idx', [df.ind.setpoint(roll), df.ind.gyroUnfilt(roll), df.ind.gyroADC(roll)], ...
        'ylabel', 'Roll Rate [deg/s]' )

    struct( ...
        'idx', [df.ind.setpoint(pitch), df.ind.gyroUnfilt(pitch), df.ind.gyroADC(pitch)], ...
        'ylabel', 'Pitch Rate [deg/s]' )

    struct( ...
        'idx', [df.ind.setpoint(yaw), df.ind.gyroUnfilt(yaw), df.ind.gyroADC(yaw)], ...
        'ylabel', 'Yaw Rate [deg/s]' )
};

plotter.plotFlightData(df, group1,'Flight Gyro Data');

 group2 = {
        struct('idx', df.ind.gyroUnfilt(1:3), ...
               'ylabel', 'Gyro (deg/s)', ...
               'legend', ["Roll","Pitch","Yaw"])

        struct('idx', df.ind.axisSum(1:3), ...
               'ylabel', 'AxisSum', ...
               'legend', ["Roll","Pitch","Yaw"])

        struct('idx', df.ind.motor, ...
               'ylabel', 'Motors', ...
               'legend', ["Motor 1","Motor 2","Motor 3","Motor 4"])

        struct('idx', df.ind.setpoint(4), ...
               'ylabel', 'Throttle', ...
               'legend', "setpoint")
    };

plotter.plotFlightData(df, group2,'Flight Overview');

 group3 = {
        struct('idx', [df.ind.currentAngle(1:2), df.ind.heading(1:2)], ...
               'ylabel', 'Current Angle [°]', ...
               'legend', ["Roll Current","Pitch Current","Roll Heading","Pitch Heading"])

        struct('idx', df.ind.angleTarget(1:2), ...
               'ylabel', 'Target Angle [°]', ...
               'legend', ["Roll","Pitch"])

    };

plotter.plotFlightData(df, group3,'Target Overview');
plotter.plot_Eval_Time(df.time);

%% Get Bodeplots

gyro_tuning = gyro_ctrl_tuning(df.data, df.ind, df.Ts_log, ...
    df.para, df.Ts_cntr);

angle_tuning = angle_ctrl_tuning(df,gyro_tuning);

resolution_factor_tf = 2;    % Window length for spectral analysis (seconds)
overlap_tf = 0.9;              % Overlap factor for spectral analysis (0-1)
gyro_tuning = gyro_tuning.calculate_transfer_func(resolution_factor_tf, overlap_tf);

resolution_factor_at = 15;    % Window length for spectral analysis (seconds)
overlap_at = 0.9;              % Overlap factor for spectral analysis (0-1)
angle_tuning = angle_tuning.calculate_Angle_trans(overlap_at, overlap_at);

% Options (gyro_tuning/angle_tuning, roll/pitch/yaw, 'Gyro'/'Angle', ...
% 'Plant'/'Complementary Sensitivity'/'Controller')

plotter.plot_Bode_Plant(gyro_tuning, roll, 'Gyro', 'Plant');
plotter.plot_Bode_Plant(angle_tuning, roll, 'Angle', 'Complementary Sensitivity');

%% Flight Analyser
analysis_flight = flight_analyzer(df.data, df.ind, df.Ts_log);

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

%% Plot Flight Analyser
plotter.plot_Spectogram(analysis_flight, 3);
plotter.plot_Gyro_spectra(analysis_flight);

%% Gyro Tuning Data

% =========================================================================
%  Axis Selection: 1: roll, 2: pitch, 3: yaw
% =========================================================================

ind_ax = 1;     % keep it now until plot_utils is finished

% I-term Relax on/off
do_compensate_iterm = true;

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
para_new.gyro_lpf            = 0;       % if lpf is static
para_new.gyro_lowpass_hz     = 0;       % frequency of gyro lpf 1
para_new.gyro_soft_type      = 0;       % type of gyro lpf 1
para_new.gyro_lowpass_dyn_hz = [0, 0];  % dyn gyro lpf overwrites gyro_lowpass_hz
para_new.gyro_lowpass2_hz    = 800;     % frequency of gyro lpf 2
para_new.gyro_soft2_type     = 0;       % type of gyro lpf 2
para_new.gyro_notch_hz       = [0, 0];  % frequency of gyro notch 1 and 2

% Gyro notch cutoff: either computed from D and center frequency or set manually
% para_new.gyro_notch_cutoff   = get_fcut_from_D_and_fcenter([0.0, 0.00], para_new.gyro_notch_hz); % damping of gyro notch 1 and 
para_new.gyro_notch_cutoff   = [0, 0]; % % Cutoff frequency gyro notch 1 and 2

para_new.dterm_lpf_hz        = 0;       % frequency of dterm lpf 1
para_new.dterm_filter_type   = 0;       % type of dterm lpf 1
para_new.dterm_lpf_dyn_hz    = [0, 0];  % dyn dterm lpf overwrites dterm_lpf_hz
para_new.dterm_lpf2_hz       = 102;     % frequency of dterm lpf 2
para_new.dterm_filter2_type  = 3;       % type of dterm lpf 2
para_new.dterm_notch_hz      = 0;       % frequency of dterm notch

% D Term notch cutoff:either computed from D and center frequency or set manually
% para_new.dterm_notch_cutoff  = get_fcut_from_D_and_fcenter(0.00, para_new.dterm_notch_hz); % damping of dterm notch
para_new.dterm_notch_cutoff  = 0;     % Cutoff frequency dterm notch
para_new.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)

%--------------------------------------------------------------------------
%  tune your PIDs here
%--------------------------------------------------------------------------
% Adjust multipliers PID values to tune the response
% Higher P: more immediate response but possible oscillations
% Higher I: better steady-state tracking but slower response
% Higher D: better damping but noise sensitive

switch ind_ax
    case 1 % Roll PID values [default: 45, 80, 30]
        P_new       = 38;
        I_new       = 80;
        D_new       = 27;
    case 2 % Pitch PID values [default: 47, 84, 34]
        P_new       = 49;
        I_new       = 96;
        D_new       = 35;
    case 3 % Yaw PID values [default: 45, 80, 0]
        P_new       = 45;
        I_new       = 86;
        D_new       = 1;
end

gyro_tuning = gyro_tuning.calculate_new_controller(ind_ax, P_new, I_new, D_new, ...
    default_parameters, para_new);
gyro_tuning = gyro_tuning.get_tuning_data(do_compensate_iterm);

% plotter.plot_Bode_Contr(ind_ax, do_insert_legends);
plotter.plot_Gang_of_Four(gyro_tuning,  'Gyro');
plotter.plot_Step_Response(gyro_tuning,  'Gyro', 'Gyro (deg/sec)');

%% Angle Tuning Data

% Old PT3
para.angle_lpf_hz = 50;     % frequency of Angle lpf (PT3)

% PT3 Angle Control
para_new.angle_lpf_hz = 50;     % frequency of Angle lpf (PT3)
P_Angle = 100;

angle_tuning = angle_tuning.calculate_new_controller(ind_ax, P_Angle, ...
    default_parameters, para_new, para);
