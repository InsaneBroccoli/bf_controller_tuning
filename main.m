clc, clear variables
addpath(genpath('lib_michi'))
addpath(genpath('lib'));
addpath logs/
addpath utils/
addpath class/
%%
% Choose an axis: 1: roll, 2: pitch, 3: yaw
ind_ax = 1;     % keep it now until plot_utils is finished


% -------------------------------------------------------------------------

% Define path to *.bbl.csv file for the first flight
log_folder1 = '';
flight_folder1 = '20250907';
log_name1 = '20250907_flipmini_00.bbl.csv';

file_path1 = fullfile(log_folder1, flight_folder1, log_name1);

% Define path to *.bbl.csv file for the second flight

second_flight = false;   % Set on true when you want to compare two flights

log_folder2 = '';
flight_folder2 = '20250908';
log_name2 = '20250908_flipmini_00.bbl.csv';

file_path2 = fullfile(log_folder2, flight_folder2, log_name2);

% Evaluation parameters
do_compensate_iterm  = false;
do_show_spec_figures = true;
do_insert_legends    = true;


pu = plot_utils;
pu.second_flight = second_flight;
pu.do_insert_legends = do_insert_legends;
pu.linewidth = 1.2;
pu.ind_ax = ind_ax;




% Defines
set(cstprefs.tbxprefs, 'MagnitudeUnits', 'abs');
set(cstprefs.tbxprefs, 'FrequencyUnits', 'Hz');
set(cstprefs.tbxprefs, 'UnwrapPhase', 'Off');
set(cstprefs.tbxprefs, 'Grid', 'On');

linewidth = 1.2;
set(0, 'defaultAxesColorOrder', get_my_colors);


% Bodeoptions
opt = bodeoptions('cstprefs');

pu.opt = opt;

% Parameter first flight
% type: 0: PT1, 1: BIQUAD, 2: PT2, 3: PT3
para_new.gyro_lpf            = 0;       % dono what this is
para_new.gyro_lowpass_hz     = 0;       % frequency of gyro lpf 1
para_new.gyro_soft_type      = 0;       % type of gyro lpf 1
para_new.gyro_lowpass_dyn_hz = [0, 0];  % dyn gyro lpf overwrites gyro_lowpass_hz
para_new.gyro_lowpass2_hz    = 800;     % frequency of gyro lpf 2
para_new.gyro_soft2_type     = 0;       % type of gyro lpf 2
para_new.gyro_notch_hz       = [0, 0]; % frequency of gyro notch 1 and 2
para_new.gyro_notch_cutoff   = get_fcut_from_D_and_fcenter([0.00, 0.00], para_new.gyro_notch_hz); % damping of gyro notch 1 and 2
para_new.dterm_lpf_hz        = 0;       % frequency of dterm lpf 1
para_new.dterm_filter_type   = 0;       % type of dterm lpf 1
para_new.dterm_lpf_dyn_hz    = [0, 0];  % dyn dterm lpf overwrites dterm_lpf_hz
para_new.dterm_lpf2_hz       = 120;     % frequency of dterm lpf 2
para_new.dterm_filter2_type  = 3;       % type of dterm lpf 2
para_new.dterm_notch_hz      = 0;     % frequency of dterm notch
para_new.dterm_notch_cutoff  = get_fcut_from_D_and_fcenter(0.00, para_new.dterm_notch_hz); % damping of dterm notch
para_new.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)
switch ind_ax
    case 1 % roll: [33, 52, 26, 0]
        P_new       = 0.4 * 33;
        I_ratio_new = 1.0 * 52/52;
        D_new       = 0.025 * 26;
    case 2 % pitch: [58, 98, 44, 0]
        P_new       = 1.0 * 58;
        I_ratio_new = 1.0 * 98/98;
        D_new       = 1.0 * 44;
    case 3 % yaw: [42, 65, 3, 0]
        P_new       = 1.0 * 42;
        I_ratio_new = 1.0 * 65/65;
        D_new       = 1.0 * 3;
end

flight1 = main_code.main_class(file_path1, para_new, ind_ax, do_compensate_iterm, ...
    do_show_spec_figures, do_insert_legends, opt, P_new, I_ratio_new, D_new);
flight1 = flight1.run();

if second_flight
    flight2 = main_code.main_class(file_path2, para_new, ind_ax, do_compensate_iterm, ...
        do_show_spec_figures, do_insert_legends, opt, P_new, I_ratio_new, D_new);
    flight2 = flight2.run();
    pu.plotGyroSignals(flight1, flight2);
    %pu.plotGyroSpectra(flight1, flight2); to be done!
    %pu.plotOverview (flight1, flight2); to be done!
else
    pu.plotGyroSignals(flight1);
    pu.plotGyroSpectra(flight1);
    pu.plotOverview (flight1);
    pu.plotBode(flight1);
end

