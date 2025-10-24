clc, clear variables
addpath(genpath('lib_michi'))
addpath(genpath('lib')); % fügt lib/ und auch alle unterordner hinzu
addpath logs/
addpath utils/
%%
% Choose an axis: 1: roll, 2: pitch, 3: yaw
ind_ax = 1;

% -------------------------------------------------------------------------

% Define quad and path to *.bbl.csv file
log_folder = '';
flight_folder = '20250907';
log_name = '20250907_aosmini_00.bbl.csv';

file_path = fullfile(log_folder, flight_folder, log_name);

% Evaluation parameters
do_compensate_iterm  = false;
do_show_spec_figures = true;
do_insert_legends    = true;

% Defines
set(cstprefs.tbxprefs, 'MagnitudeUnits', 'abs');
set(cstprefs.tbxprefs, 'FrequencyUnits', 'Hz');
set(cstprefs.tbxprefs, 'UnwrapPhase', 'Off');
set(cstprefs.tbxprefs, 'Grid', 'On');

linewidth = 1.2;
set(0, 'defaultAxesColorOrder', get_my_colors);
pos_bode = [0.1514, 0.5838-0.2, 0.7536, 0.3472+0.2; ... % this is a bit hacky
            0.1514, 0.1100    , 0.7536, 0.1917    ];

% Bodeoptions
opt = bodeoptions('cstprefs');


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

drone1 = main_code.main_class(file_path, para_new, ind_ax, do_compensate_iterm, ...
    do_show_spec_figures, do_insert_legends, opt, P_new, I_ratio_new, D_new);
drone1.run();