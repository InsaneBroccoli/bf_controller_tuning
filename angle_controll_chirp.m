%% Start
clc, clear variables, close all
addpath(genpath('lib'));
addpath logs/

% Add file information
log_folder = 'logs';
flight_folder = '20251212';
log_name = '06_20251212_OvershootExpress.TXT.csv';
file_path = fullfile(log_folder, flight_folder, log_name);


% --- Load and Process Flight Log Data ---
[para, Nheader, ind, ind_cntr] = extract_header_information(file_path);

% Load data from CSV or cached MAT file for faster processing
% robust MAT path
[folder, base, ~] = fileparts(file_path);
mat_path = fullfile(folder, [base '.mat']);

try
    S = load(mat_path);
    data = S.data;
catch
    data = readmatrix(file_path, 'NumHeaderLines', Nheader);
    save(mat_path, 'data');
end

Ts= para.looptime * 1.0e-6;             % Gyro loop
Ts_cntr = para.pid_process_denom * Ts;        % Control loop
Ts_log  = para.frameIntervalPDenom * Ts_cntr; % Logging loop   

% Get a useable Timevector
time = (data(:,ind.time) - data(1,ind.time)) * 1.0e-6;

% Define numbers of axis
axis = 1;       %roll = 1; pitch = 2;


%% Creation of Chirpsignal
sinarg_idx = ind.debug(1);
sinarg_ax  = data(:, sinarg_idx) / 5e3;   % phase-like signal

switch axis
    case 1
        A = para.chirp_amplitude_roll;
    case 2
        A = para.chirp_amplitude_pitch;
    case 3
        A = para.chirp_amplitude_yaw;
    otherwise
        error('Unknown axis.');
end

% Detect active chirp interval
ind_eval = get_ind_eval(sinarg_ax, data(:, ind.gyroADC(axis)));
idx = ind_eval;

% Reconstruct chirp input
chirp = zeros(size(sinarg_ax));
chirp(ind_eval) = A * cos(sinarg_ax(ind_eval));

figure(11)
plot(time, chirp);


%%

% Check if only the chirpsignal is not zero
figure(2)
subplot(411)
plot(time(idx), tdata(idx, currentAngle(axis)), '-b');
hold on
plot(time(idx), tdata(idx, ind.heading(axis)),'-r');grid on;
legend('Current Angle','Heading','Location','best');
title('Compare Current Anlge and Heading');
xlabel('Time [s]'); ylabel('deg [°]');

subplot(412)
plot(time(idx), tdata(idx,angleRate(axis)),'-b');
hold on
plot(time(idx), tdata(idx, ind.setpoint(axis)),'-r');grid on;
legend('Angle Rate','Setpoint','Location','best');
title('Compare Angle Rate and Current Setpoint');
xlabel('Time [s]'); ylabel('deg/s [°/s]');

subplot(413)
plot(time(idx), tdata(idx,angleTarget(axis)),'-b');grid on;
title('Target Angle');
xlabel('Time [s]'); ylabel('deg [°]');

subplot(414)
plot(time(idx), tdata(idx,ind.gyroADC(axis)),'-b');grid on;
title('Gyro ADC');
xlabel('Time [s]'); ylabel('deg/s  [°/s]');
%% Lets Measure the Transferfunction and such kind of things

% Pararmeters
Nest = round(15 / Ts_log);               % Window length in samples
Noverlap = floor(0.9 * Nest);           % Overlap between windows
window   = hann(Nest, 'periodic');  % Hanning window for analysis

% Design linear filter for zero phase excitation filter (apply_rotfiltfilt)
Dlp = sqrt(3) / 2;    % Damping ratio
wlp = 2 * pi * 10;    % Cutoff frequency [rad/s]
Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), ...    % Discrete filter
        Ts_log, 'tustin');                   % Using Tustin transform

% Calculation of the Transferfunction of the System
w = tdata(:, angleTarget(axis));
inp = apply_rotfiltfilt(Glp, sinarg_ax, w);

y = tdata(:, currentAngle(axis));
out = apply_rotfiltfilt(Glp, sinarg_ax, y);

[T_ax, C_T_ax] = estimate_frequency_response(inp(idx), ...
    out(idx), window, Noverlap, Nest, Ts_log);

% Calculation of the plant
v = tdata(:, ind.gyroADC(axis));
out_v = apply_rotfiltfilt(Glp, sinarg_ax, v);

[G_wv, C_G_wv] = estimate_frequency_response(inp(idx), ...
    out_v(idx), window, Noverlap, Nest, Ts_log);

P_angle = T_ax / G_wv;

[G_vy, C_G_vy] = estimate_frequency_response(out_v(idx), ...
    out(idx), window, Noverlap, Nest, Ts_log);

% Calculation of the controller
c = tdata(:, ind.setpoint(axis));
out_c = apply_rotfiltfilt(Glp, sinarg_ax, c);

[G_wc,C_G_wc] = estimate_frequency_response(inp(idx), out_c(idx), ...
    window, Noverlap, Nest, Ts_log);

Cp_ax = G_wc / (1 - T_ax);

% Transferfunction of the Gyro loop
T_gy = G_wv / G_wc;

% Direct methode
[G_cv,C_G_cv] = estimate_frequency_response(out_c(idx), out_v(idx), ...
    window, Noverlap, Nest, Ts_log);

%% Lets Calculate the Transferfunctions

% Get the Controller
C_P_Angle = 50;   % 0.1 is scaling in Betaflight
f  = T_ax.Frequency*2*pi;            % Frequency

C_P_Angle_frd = frd( C_P_Angle * 0.1 *ones(size(f)), f, Ts_cntr );
fc = 50;  % Hz  (dein "Angle PT3")
angle_fcut_ana = fc;
Gf_ana = get_filter('pt3', fc, Ts_cntr);  % should be discrete

C_Angle_ana = C_P_Angle_frd * Gf_ana;
C_Angle_ana = downsample_frd(C_Angle_ana, Ts_log, T_ax.Frequency);

% Calculate the Transferfunction with this Controller
T_ana = (C_Angle_ana * T_gy * P_angle)/(1+ C_Angle_ana * T_gy * P_angle);

%% Lets start Tuning

% Add the Parameters
P_new = 100;            % Gain 
angle_lpf_hz = 50;      % Cut off Frequency

C_P_Angle_new_frd = frd( P_new * 0.1 *ones(size(f)), f, Ts_cntr );
Gf_new = get_filter('pt3', angle_lpf_hz, Ts_cntr);  % should be discrete

C_Angle_new = C_P_Angle_new_frd * Gf_ana;
C_Angle_new = downsample_frd(C_Angle_new, Ts_log, T_ax.Frequency);

% Calculate the new Transferfunction with the new Parameters
T_new = (C_Angle_new * T_gy * P_angle)/(1+ C_Angle_new * T_gy * P_angle);

%% Lets Plot the figures
% Bodeplot of the hole system
figure(3)
bode(T_ax, T_ana, T_new);grid on;
legend('Measured','Calculated','New Calculated','Location','best')
title('Measured Transferfunction System');

% Bodeplot of the plant angle
figure(4)
bode(P_angle, G_vy); grid on;
legend('Indirect','Direct','Location','best');
title('Measured Plant Angle');

% Bodeplot of the Controller
figure(5)
bode(Cp_ax, C_Angle_ana); grid on;
legend('Indirect','Simulation','Location','best');
title('Measured Controller Plant');

% Bodeplot of the Gyro Loop
figure(6)
bode(T_gy, G_cv); grid on;
title('Measured Transferfunction Gyro Loop')

%% Calculate Measured Step Response
fmax = 600;     % Not sure which value I should choose so YOLO

step_time = (0:Nest-1).*Ts_log;                 % Time vector [s]
T_mean = 0.1 * [-1, 1] + (Nest * Ts_log) / 2;    % Analysis window [s]

step_resp = [calculate_step_response_from_frd(T_ax , fmax), ...    % Used parameters
             calculate_step_response_from_frd(T_ana, fmax), ...    % New parameters
             calculate_step_response_from_frd(T_new, fmax)];  

step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2),:));
step_resp_meas = step_resp ./ step_resp_mean;



figure(7)
plot(step_time, step_resp_meas);grid on;
title('Step Response Measured')
xlabel('Time [s]'); ylabel('Angle [°]');
