%% Simulation Frequenzgang Drohen

clc, clear variables
addpath lib/
s = tf('s');

%% Laden von Drohnenflüge als Referenz
flight_folder = '20250907';

quad = 'aosmini';
log_name = '20250907_aosmini_00.bbl.csv';

% quad = 'apex5';
% log_name = '20250907_apex5_00.bbl.csv';

% quad = 'flipmini';
% log_name = '20250907_flipmini_00.bbl.csv';

% Choose an axis: 1: roll, 2: pitch, 3: yaw
ind_ax = 1;

file_path = fullfile(flight_folder, log_name);
[para, Nheader, ind, ind_cntr] = extract_header_information(file_path);

% Read the data
%  - If its the first time from the .csv and save a mat, otherwise the
%    .mat. This increases load speed significantly.
tic
try
   load([file_path(1:end-8), '.mat'])
catch exception
   data = readmatrix(file_path, 'NumHeaderLines', Nheader);
   save([file_path(1:end-8), '.mat'], 'data');
end
[Ndata, Nsig] = size(data)
toc
% Expand index
ind.axisSumPI = ind_cntr + (1:3);
ind.sinarg = ind.debug(1);


%% Filterparameter
switch quad
    case 'aosmini'
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
                P_new       = 33;
                I_ratio_new = 52/52;
                D_new       = 26;
            case 2 % pitch: [58, 98, 44, 0]
                P_new       = 58;
                I_ratio_new = 98/98;
                D_new       = 44;
            case 3 % yaw: [42, 65, 3, 0]
                P_new       = 42;
                I_ratio_new = 65/65;
                D_new       = 3;
        end
    case 'apex5'
        % type: 0: PT1, 1: BIQUAD, 2: PT2, 3: PT3
        para_new.gyro_lpf            = 0;       % dono what this is
        para_new.gyro_lowpass_hz     = 0;       % frequency of gyro lpf 1
        para_new.gyro_soft_type      = 0;       % type of gyro lpf 1
        para_new.gyro_lowpass_dyn_hz = [0, 0];  % dyn gyro lpf overwrites gyro_lowpass_hz
        para_new.gyro_lowpass2_hz    = 800;     % frequency of gyro lpf 2
        para_new.gyro_soft2_type     = 0;       % type of gyro lpf 2
        para_new.gyro_notch_hz       = [0, 520]; % frequency of gyro notch 1 and 2
        para_new.gyro_notch_cutoff   = get_fcut_from_D_and_fcenter([0.00, 0.15], para_new.gyro_notch_hz); % damping of gyro notch 1 and 2
        para_new.dterm_lpf_hz        = 0;       % frequency of dterm lpf 1
        para_new.dterm_filter_type   = 0;       % type of dterm lpf 1
        para_new.dterm_lpf_dyn_hz    = [0, 0];  % dyn dterm lpf overwrites dterm_lpf_hz
        para_new.dterm_lpf2_hz       = 130;     % frequency of dterm lpf 2
        para_new.dterm_filter2_type  = 3;       % type of dterm lpf 2
        para_new.dterm_notch_hz      = 235;     % frequency of dterm notch
        para_new.dterm_notch_cutoff  = get_fcut_from_D_and_fcenter(0.15, para_new.dterm_notch_hz); % damping of dterm notch
        para_new.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)
        switch ind_ax
            case 1 % roll: [49, 83, 33, 0]
                P_new       = 49;
                I_ratio_new = 83/83;
                D_new       = 33;
            case 2 % pitch: [61, 103, 39, 0]
                P_new       = 61;
                I_ratio_new = 103/103;
                D_new       = 39;
            case 3 % yaw: [42, 104, 3, 0]
                P_new       = 42;
                I_ratio_new = 104/104;
                D_new       = 3;
        end
    case 'flipmini'
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
        para_new.dterm_lpf2_hz       = 140;     % frequency of dterm lpf 2
        para_new.dterm_filter2_type  = 3;       % type of dterm lpf 2
        para_new.dterm_notch_hz      = 0;     % frequency of dterm notch
        para_new.dterm_notch_cutoff  = get_fcut_from_D_and_fcenter(0.00, para_new.dterm_notch_hz); % damping of dterm notch
        para_new.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)
        switch ind_ax
            case 1 % roll: [46, 74, 30, 0]
                P_new       = 46;
                I_ratio_new = 74/74;
                D_new       = 30;
            case 2 % pitch: [71, 118, 47, 0]
                P_new       = 71;
                I_ratio_new = 118/118;
                D_new       = 47;
            case 3 % yaw: [35, 70, 3, 0]
                P_new       = 35;
                I_ratio_new = 70/70;
                D_new       = 3;
        end
    otherwise
        warning(' no valid quad selected');
end

%% PID Vektor
pid_scale = [get_pid_scale(ind_ax), 1];

PID_new(1) = P_new * pid_scale(1);                  % Proportionalanteil

%KI von KP abhängig, daher rückrechnung nötig
fI = PID(2) / (2 * pi * PID(1));        % alte Integralfrequenz aus altem PID
fI_new = fI * I_ratio_new;              % gewünschte skalierte Integralfrequenz

PID_new(2) = 2 * pi * PID_new(1) * fI_new;          % Integralanteil
%Differentialanteil unabhängig
PID_new(3) = D_new * pid_scale(3);                  % Differentialanteil
PID_new(4) = 0;

%% Generelle Gyro-Tiefpassfilter

%Berechnung LPF2
omega2 = 2*pi*para_new.gyro_lowpass_hz;
G_LPF2 = 1 / (1+ s/omega2);

%Berechnung D-spezifischer Tiefpassfilter PT3
omegad = 2*pi*para_new.dterm_lpf2_hz;
G_LPFD = (1 / (1+ s/omegad))^3;

%Berechnung Notchfilter

