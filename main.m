clc, clear variables
addpath(genpath('lib'));
addpath logs/
addpath class/
%%
% Choose an axis: 1: roll, 2: pitch, 3: yaw
ind_ax = 1;     % keep it now until plot_utils is finished

% -------------------------------------------------------------------------

% Define path to *.bbl.csv file for the first flight
log_folder1 = 'logs';
flight_folder1 = '20250907';
log_name1 = '20250907_flipmini_00.bbl.csv';

file_path1 = fullfile(log_folder1, flight_folder1, log_name1);

% Define path to *.bbl.csv file for the second flight

second_flight = false;   % Set on true when you want to compare two flights

log_folder2 = '';
flight_folder2 = '20250908';
log_name2 = '20250908_flipmini_00.bbl.csv';

file_path2 = fullfile(log_folder2, flight_folder2, log_name2);

% Show Legend
do_insert_legends    = true;    % Set on true when you want to see the legende

% Define if Item Relax should be on
do_compensate_iterm  = false;   % Set on true if you want to compensate with item relax

% Parameters of Spectogramfigures
do_show_spec_figures = true;
Nestfaspec = 0.2;
koverlapspec = 0.9;

% Show Plots 
do_show_flight_track = true;
do_show_tuning_figures = true;

% Parameter of Transferfunction
do_show_transferfunction = true;
Nestfatra = 2.0;
koverlaptra = 0.9;

pu = plot_utils;
pu.second_flight = second_flight;
pu.do_insert_legends = do_insert_legends;
pu.linewidth = 1.2;
pu.ind_ax = ind_ax;


% Defines
set(cstprefs.tbxprefs,'MagnitudeUnits','abs');
set(cstprefs.tbxprefs,'FrequencyUnits','Hz');
set(cstprefs.tbxprefs,'UnwrapPhase','Off');
set(cstprefs.tbxprefs,'Grid','On');

% ---- Central Option for Bodeplots ----
opt = bodeoptions('cstprefs');    % start with cstprefs
opt.MagScale      = 'log';        % Magnitude logarithmic
opt.PhaseUnits    = 'deg';        % Phase in degrees
opt.PhaseWrapping = 'on';         % Phase in the area [-180,180]    
opt.YLimMode      = {'auto'; 'manual'};  % Autoscale for Mag and Manuel for Phase
opt.YLim          = {[-180 180]};   % Area of Phase
opt.Grid          = 'on';           % Grid ob

pu.opt = opt;

% New and old parameters are the same
default_parameters = false; 

% Parameter first flight
% type: 0: PT1, 1: BIQUAD, 2: PT2, 3: PT3
para_new1.gyro_lpf            = 0;       % dono what this is
para_new1.gyro_lowpass_hz     = 0;       % frequency of gyro lpf 1
para_new1.gyro_soft_type      = 0;       % type of gyro lpf 1
para_new1.gyro_lowpass_dyn_hz = [0, 0];  % dyn gyro lpf overwrites gyro_lowpass_hz
para_new1.gyro_lowpass2_hz    = 800;     % frequency of gyro lpf 2
para_new1.gyro_soft2_type     = 0;       % type of gyro lpf 2
para_new1.gyro_notch_hz       = [0, 0]; % frequency of gyro notch 1 and 2
para_new1.gyro_notch_cutoff   = get_fcut_from_D_and_fcenter([0.00, 0.00], para_new1.gyro_notch_hz); % damping of gyro notch 1 and 2
para_new1.dterm_lpf_hz        = 0;       % frequency of dterm lpf 1
para_new1.dterm_filter_type   = 0;       % type of dterm lpf 1
para_new1.dterm_lpf_dyn_hz    = [0, 0];  % dyn dterm lpf overwrites dterm_lpf_hz
para_new1.dterm_lpf2_hz       = 120;     % frequency of dterm lpf 2
para_new1.dterm_filter2_type  = 3;       % type of dterm lpf 2
para_new1.dterm_notch_hz      = 0;     % frequency of dterm notch
para_new1.dterm_notch_cutoff  = get_fcut_from_D_and_fcenter(0.00, para_new1.dterm_notch_hz); % damping of dterm notch
para_new1.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)
switch ind_ax
    case 1 % roll: [33, 52, 26, 0]
        P_new1       = 0.4 * 33;
        I_ratio_new1 = 1.0 * 52/52;
        D_new1       = 0.025 * 26;
    case 2 % pitch: [58, 98, 44, 0]
        P_new1       = 1.0 * 58;
        I_ratio_new1 = 1.0 * 98/98;
        D_new1       = 1.0 * 44;
    case 3 % yaw: [42, 65, 3, 0]
        P_new1       = 1.0 * 42;
        I_ratio_new1 = 1.0 * 65/65;
        D_new1       = 1.0 * 3;
end

% Parameter seconde flight
% type: 0: PT1, 1: BIQUAD, 2: PT2, 3: PT3
para_new2.gyro_lpf            = 0;       % dono what this is
para_new2.gyro_lowpass_hz     = 0;       % frequency of gyro lpf 1
para_new2.gyro_soft_type      = 0;       % type of gyro lpf 1
para_new2.gyro_lowpass_dyn_hz = [0, 0];  % dyn gyro lpf overwrites gyro_lowpass_hz
para_new2.gyro_lowpass2_hz    = 800;     % frequency of gyro lpf 2
para_new2.gyro_soft2_type     = 0;       % type of gyro lpf 2
para_new2.gyro_notch_hz       = [0, 0]; % frequency of gyro notch 1 and 2
para_new2.gyro_notch_cutoff   = get_fcut_from_D_and_fcenter([0.00, 0.00], para_new2.gyro_notch_hz); % damping of gyro notch 1 and 2
para_new2.dterm_lpf_hz        = 0;       % frequency of dterm lpf 1
para_new2.dterm_filter_type   = 0;       % type of dterm lpf 1
para_new2.dterm_lpf_dyn_hz    = [0, 0];  % dyn dterm lpf overwrites dterm_lpf_hz
para_new2.dterm_lpf2_hz       = 120;     % frequency of dterm lpf 2
para_new2.dterm_filter2_type  = 3;       % type of dterm lpf 2
para_new2.dterm_notch_hz      = 0;     % frequency of dterm notch
para_new2.dterm_notch_cutoff  = get_fcut_from_D_and_fcenter(0.00, para_new2.dterm_notch_hz); % damping of dterm notch
para_new2.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)
switch ind_ax
    case 1 % roll: [33, 52, 26, 0]
        P_new2       = 1.0 * 33;
        I_ratio_new2 = 1.0 * 52/52;
        D_new2       = 1.0 * 26;
    case 2 % pitch: [58, 98, 44, 0]
        P_new2       = 1.0 * 58;
        I_ratio_new2 = 1.0 * 98/98;
        D_new2       = 1.0 * 44;
    case 3 % yaw: [42, 65, 3, 0]
        P_new2       = 1.0 * 42;
        I_ratio_new2 = 1.0 * 65/65;
        D_new2       = 1.0 * 3;
end

flight1 = main_class(file_path1, para_new1, ind_ax, do_compensate_iterm, ...
    P_new1, I_ratio_new1, D_new1, Nestfaspec,koverlapspec, Nestfatra, koverlaptra ...
    ,default_parameters);
flight1 = flight1.run();

if second_flight
    flight2 = main_class(file_path2 ,para_new2, ind_ax, do_compensate_iterm, ...
        P_new2, I_ratio_new2, D_new2, Nestfaspec,koverlapspec,Nestfatra, koverlaptra, ...
        default_parameters);
    flight2 = flight2.run();
end
%% Plots

if second_flight
    pu.plotevaltime(flight1, flight2);
    if do_show_flight_track
        pu.plotGyroSignals(flight1, flight2);
        pu.plotOverview (flight1, flight2);
    end
    if do_show_tuning_figures
        pu.plotStepResp(flight1,flight2);
        pu.plotGangofFour(flight1, flight2);
    end

    if do_show_transferfunction
        pu.plotBode(flight1, flight2);
        pu.plotCPIDBode(flight1, flight2);
        pu.plotController(flight1, flight2);
    end
    if do_show_spec_figures
        pu.plotGyroSpectra(flight1, flight2);
        pu.plotspectogram(flight1, flight2);
    end

else
     pu.plotevaltime(flight1);
    if do_show_flight_track
        pu.plotGyroSignals(flight1);
        pu.plotOverview (flight1);
    end
    if do_show_tuning_figures
        pu.plotStepResp(flight1);
        pu.plotGangofFour(flight1);
    end
    if do_show_transferfunction
        pu.plotBode(flight1);
        pu.plotCPIDBode(flight1);
        pu.plotController(flight1);
    end
   
    if do_show_spec_figures
        pu.plotGyroSpectra(flight1);
        pu.plotspectogram(flight1);
    end
end