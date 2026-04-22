%% Initialize workspace
clc;
clear variables;
close all;

addpath(genpath('lib'));
addpath('logs');

%% File paths
log_folder    = 'logs';
flight_folder1 = '20260408';
flight_folder2 = '20260420';

log_name1 = 'Flipmini_P100.TXT.csv';
log_name2 = '20260420_flipmini_P80.TXT.csv';

file_path1 = fullfile(log_folder, flight_folder1, log_name1);
file_path2 = fullfile(log_folder, flight_folder2, log_name2);

%% Load and process header information
[para1, Nheader1, ind1, ind_cntr1] = extract_header_information(file_path1);
[para2, Nheader2, ind2, ind_cntr2] = extract_header_information(file_path2);

%% Load data from CSV or cached MAT file
[folder1, base1, ~] = fileparts(file_path1);
mat_path1 = fullfile(folder1, [base1 '.mat']);

[folder2, base2, ~] = fileparts(file_path2);
mat_path2 = fullfile(folder2, [base2 '.mat']);

try
    S = load(mat_path1);
    if isfield(S, 'flight1')
        flight1 = S.flight1;
    elseif isfield(S, 'data')
        flight1 = S.data;
    else
        error('No valid variable found in MAT file for flight1.');
    end
catch
    flight1 = readmatrix(file_path1, 'NumHeaderLines', Nheader1);
    save(mat_path1, 'flight1');
end

try
    S = load(mat_path2);
    if isfield(S, 'flight2')
        flight2 = S.flight2;
    elseif isfield(S, 'data')
        flight2 = S.data;
    else
        error('No valid variable found in MAT file for flight2.');
    end
catch
    flight2 = readmatrix(file_path2, 'NumHeaderLines', Nheader2);
    save(mat_path2, 'flight2');
end

%% Sampling times
Ts      = para1.looptime * 1.0e-6;                 % Gyro loop
Ts_cntr = para1.pid_process_denom * Ts;            % Control loop
Ts_log  = para1.frameIntervalPDenom * Ts_cntr;     % Logging loop

%% Time vector
time1 = (flight1(:, ind1.time) - flight1(1, ind1.time)) * 1.0e-6;
time2 = (flight2(:, ind2.time) - flight2(1, ind2.time)) * 1.0e-6;

%% Define axis
axis = 1;   % roll = 1, pitch = 2

%% Debug indices
sinarg1       = ind1.debug(1);
currentAngle1 = [ind1.debug(2), ind1.debug(5)];
angleTarget1  = [ind1.debug(3), ind1.debug(6)];
angleRate1    = [ind1.debug(4), ind1.debug(7)];

sinarg2       = ind2.debug(1);
currentAngle2 = [ind2.debug(2), ind2.debug(5)];
angleTarget2  = [ind2.debug(3), ind2.debug(6)];
angleRate2    = [ind2.debug(4), ind2.debug(7)];

%% Scale flight 1 data
flight1(:, sinarg1)              = flight1(:, sinarg1) / 5e3;
flight1(:, currentAngle1(axis))  = flight1(:, currentAngle1(axis)) * 0.1;
flight1(:, ind1.setpoint(axis))  = flight1(:, ind1.setpoint(axis)) * 0.1;
flight1(:, angleRate1(axis))     = flight1(:, angleRate1(axis)) * 0.1;
flight1(:, angleTarget1(axis))   = flight1(:, angleTarget1(axis)) * 0.1;
flight1(:, ind1.heading(axis))   = flight1(:, ind1.heading(axis)) * 100;
flight1(:, ind1.gyroADC(axis))   = flight1(:, ind1.gyroADC(axis)) * 0.1;

%% Scale flight 2 data
flight2(:, sinarg2)              = flight2(:, sinarg2) / 5e3;
flight2(:, currentAngle2(axis))  = flight2(:, currentAngle2(axis)) * 0.1;
flight2(:, ind2.setpoint(axis))  = flight2(:, ind2.setpoint(axis)) * 0.1;
flight2(:, angleRate2(axis))     = flight2(:, angleRate2(axis)) * 0.1;
flight2(:, angleTarget2(axis))   = flight2(:, angleTarget2(axis)) * 0.1;
flight2(:, ind2.heading(axis))   = flight2(:, ind2.heading(axis)) * 100;
flight2(:, ind2.gyroADC(axis))   = flight2(:, ind2.gyroADC(axis)) * 0.1;

%% Evaluation masks
ind_eval1 = get_ind_eval(flight1(:, sinarg1), flight1(:, ind1.gyroADC(axis)));
ind_eval2 = get_ind_eval(flight2(:, sinarg2), flight2(:, ind2.gyroADC(axis)));

sinarg_ax1 = flight1(:, sinarg1);
sinarg_ax1(~ind_eval1) = 0;

sinarg_ax2 = flight2(:, sinarg2);
sinarg_ax2(~ind_eval2) = 0;

idx1 = find(ind_eval1);
idx2 = find(ind_eval2);

%% Frequency response estimation settings
Nest     = round(2 / Ts_log);          % Window length in samples
Noverlap = floor(0.9 * Nest);          % Overlap
window   = hann(Nest, 'periodic');     % Hann window

%% Excitation filter
Dlp = sqrt(3) / 2;                     % Damping ratio
wlp = 2 * pi * 10;                     % Cutoff frequency [rad/s]

Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');

%% =========================
%% Flight 1
%% =========================

% Closed loop transfer function
w1   = flight1(:, angleTarget1(axis));
inp1 = apply_rotfiltfilt(Glp, sinarg_ax1, w1);

y1   = flight1(:, currentAngle1(axis));
out1 = apply_rotfiltfilt(Glp, sinarg_ax1, y1);

[T_ax1, C_T_ax1] = estimate_frequency_response(inp1(idx1), out1(idx1), ...
    window, Noverlap, Nest, Ts_log);

% Plant
v1     = flight1(:, ind1.gyroADC(axis));
out_v1 = apply_rotfiltfilt(Glp, sinarg_ax1, v1);

[G_wv1, C_G_wv1] = estimate_frequency_response(inp1(idx1), out_v1(idx1), ...
    window, Noverlap, Nest, Ts_log);

P_angle1 = T_ax1 / G_wv1;

% Controller
c1     = flight1(:, ind1.setpoint(axis));
out_c1 = apply_rotfiltfilt(Glp, sinarg_ax1, c1);

[G_wc1, C_G_wc1] = estimate_frequency_response(inp1(idx1), out_c1(idx1), ...
    window, Noverlap, Nest, Ts_log);

Cp_ax1 = G_wc1 / (1 - T_ax1);

% Gyro loop transfer function
T_gy1 = G_wv1 / G_wc1;

[T_gy1, C_G_wc1] = estimate_frequency_response(out_c1(idx1), out_v1(idx1), ...
    window, Noverlap, Nest, Ts_log);

%% =========================
%% Flight 2
%% =========================

% Closed loop transfer function
w2   = flight2(:, angleTarget2(axis));
inp2 = apply_rotfiltfilt(Glp, sinarg_ax2, w2);

y2   = flight2(:, currentAngle2(axis));
out2 = apply_rotfiltfilt(Glp, sinarg_ax2, y2);

[T_ax2, C_T_ax2] = estimate_frequency_response(inp2(idx2), out2(idx2), ...
    window, Noverlap, Nest, Ts_log);

% Plant
v2     = flight2(:, ind2.gyroADC(axis));
out_v2 = apply_rotfiltfilt(Glp, sinarg_ax2, v2);

[G_wv2, C_G_wv2] = estimate_frequency_response(inp2(idx2), out_v2(idx2), ...
    window, Noverlap, Nest, Ts_log);

P_angle2 = T_ax2 / G_wv2;

% Controller
c2     = flight2(:, ind2.setpoint(axis));
out_c2 = apply_rotfiltfilt(Glp, sinarg_ax2, c2);

[G_wc2, C_G_wc2] = estimate_frequency_response(inp2(idx2), out_c2(idx2), ...
    window, Noverlap, Nest, Ts_log);

Cp_ax2 = G_wc2 / (1 - T_ax2);

% Gyro loop transfer function
T_gy2 = G_wv2 / G_wc2;

%% Analytical transfer function
C_P_Angle = 100;                                % Betaflight gain, scaled by 0.1
f = T_ax1.Frequency * 2 * pi;                  % Angular frequency [rad/s]

C_P_Angle_frd = frd(C_P_Angle * 0.1 * ones(size(f)), f, Ts_cntr);

fc = 50;                                       % Hz
Gf_ana = get_filter('pt3', fc, Ts_cntr);

C_Angle_ana = C_P_Angle_frd * Gf_ana;
C_Angle_ana = downsample_frd(C_Angle_ana, Ts_log, T_ax1.Frequency);

% Analytical closed-loop TF based on flight 1
T_ana1 = (C_Angle_ana * T_gy1 * P_angle1) / (1 + C_Angle_ana * T_gy1 * P_angle1);
T_ana2 = (C_Angle_ana * T_gy2 * P_angle2) / (1 + C_Angle_ana * T_gy2 * P_angle2);

%% Tuning
P_new = 100;               % New gain
angle_lpf_hz = 50;         % New cutoff frequency

C_P_Angle_new_frd = frd(P_new * 0.1 * ones(size(f)), f, Ts_cntr);
Gf_new = get_filter('pt3', angle_lpf_hz, Ts_cntr);

C_Angle_new = C_P_Angle_new_frd * Gf_new;
C_Angle_new = downsample_frd(C_Angle_new, Ts_log, T_ax1.Frequency);

T_new1 = (C_Angle_new * T_gy1 * P_angle1) / (1 + C_Angle_new * T_gy1 * P_angle1);
T_new2 = (C_Angle_new * T_gy2 * P_angle2) / (1 + C_Angle_new * T_gy2 * P_angle2);

%% Step responses
fmax = 400;

step_time = (0:Nest-1) .* Ts_log;
T_mean = 0.1 * [-1, 1] + (Nest * Ts_log) / 2;

step_resp = [ ...
    calculate_step_response_from_frd(T_ax2,  fmax), ...
    calculate_step_response_from_frd(T_new2, fmax), ...
    calculate_step_response_from_frd(T_ax1, fmax)];

step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2), :), 1);
step_resp_meas = step_resp ./ step_resp_mean;

%% Plot
figure(1);
plot(step_time, step_resp_meas, 'LineWidth', 1.2);
grid on;
title('Measured Step Response');
xlim([0 1])
xlabel('Time [s]');
ylabel('Angle [deg]');
legend('Measured P=80', 'Calculated P=100', 'Measured P=100', ...
    'Location', 'best');

figure(2)
bode(Cp_ax1,C_Angle_ana, C_Angle_new);