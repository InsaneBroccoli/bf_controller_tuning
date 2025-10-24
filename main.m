clc, clear variables
addpath(genpath('lib_michi'))
addpath(genpath('lib')); % fügt lib/ und auch alle unterordner hinzu
addpath logs/
addpath utils/
%%
% Choose an axis: 1: roll, 2: pitch, 3: yaw
ind_ax = 1;

% -------------------------------------------------------------------------

% Define path to *.bbl.csv file for the first flight
log_folder1 = '';
flight_folder1 = '20250907';
log_name1 = '20250907_flipmini_00.bbl.csv';

file_path1 = fullfile(log_folder1, flight_folder1, log_name1);

% Define path to *.bbl.csv file for the second flight

second_flight = true;   % Set on true when you want to compare two flights

log_folder2 = '';
flight_folder2 = '20250908';
log_name2 = '20250908_flipmini_00.bbl.csv';

file_path2 = fullfile(log_folder2, flight_folder2, log_name2);


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


% Bodeoptions
opt = bodeoptions('cstprefs');


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
end

%% Do plots

figure(1)

% --- Roll ---
ax(1) = subplot(311); hold(ax(1),'on')
plot(ax(1), flight1.gyroData.time, flight1.gyroData.roll)
if do_insert_legends, legend('setpoint F1', 'gyro F1', 'gyroADC F1', 'location', 'best'), end
if second_flight
    plot(ax(1), flight2.gyroData.time, flight2.gyroData.roll)
    if do_insert_legends, legend('setpoint F1', 'gyro F1', 'gyroADC F1', ...
            'setpoint F2', 'gyro F2', 'gyroADC F2', 'location', 'best'), end
end
grid(ax(1),'on'); ylabel(ax(1),'Roll (deg/sec)')
title(ax(1),'Gyro Signals')


% --- Pitch ---
ax(2) = subplot(312); hold(ax(2),'on')
plot(ax(2), flight1.gyroData.time, flight1.gyroData.pitch)
if second_flight
    plot(ax(2), flight2.gyroData.time, flight2.gyroData.pitch)
end
legend(ax(2),'off')
grid(ax(2),'on'); ylabel(ax(2),'Pitch (deg/sec)')


% --- Yaw ---
ax(3) = subplot(313); hold(ax(3),'on')
plot(ax(3), flight1.gyroData.time, flight1.gyroData.yaw)
if second_flight
    plot(ax(3), flight2.gyroData.time, flight2.gyroData.yaw)
end
legend(ax(3),'off')
grid(ax(3),'on'); ylabel(ax(3),'Yaw (deg/sec)'); xlabel(ax(3),'Time (sec)')

% Achsen verlinken & Limits setzen
linkaxes(ax,'x');

if second_flight
    xmax = max([flight1.gyroData.time; flight2.gyroData.time]);
else
    xmax = max(flight1.gyroData.time);
end
xlim(ax, [0 xmax]);

% Styling
set(findall(gcf,'type','line'), 'linewidth', linewidth)
