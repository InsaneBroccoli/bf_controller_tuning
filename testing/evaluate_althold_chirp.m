% Initialize workspace
clc, clear variables, close all
addpath("../lib/");
addpath(genpath("../logs/"));

%  Chrip Amplitude 100
% log_folder = '../logs';
% flight_folder = '20260323';
% log_name = 'A100.TXT.csv';
% file_path = fullfile(log_folder, flight_folder, log_name);

%  Chrip Amplitude 200
log_folder = '../logs';
flight_folder = '20260323';
log_name = 'A200.TXT.csv';
file_path = fullfile(log_folder, flight_folder, log_name);

%  Chrip Amplitude 250
% log_folder = '../logs';
% flight_folder = '20260323';
% log_name = 'A250.TXT.csv';
% file_path = fullfile(log_folder, flight_folder, log_name);

%  Chrip Amplitude 300
% log_folder = '../logs';
% flight_folder = '20260323';
% log_name = 'A300.TXT.csv';
% file_path = fullfile(log_folder, flight_folder, log_name);

%  Chrip Amplitude 400
% log_folder = '../logs';
% flight_folder = '20260323';
% log_name = 'A400.TXT.csv';
% file_path = fullfile(log_folder, flight_folder, log_name);

% Chrip Amplitude 800
% log_folder = '../logs';
% flight_folder = '20260317';
% log_name = 'LOG077.TXT.csv';
% file_path = fullfile(log_folder, flight_folder, log_name);

% log_folder = '../logs';
% flight_folder = '20260317';
% log_name = 'LOG070.TXT.csv';
% file_path = fullfile(log_folder, flight_folder, log_name);

% log_folder = '../logs';
% flight_folder = '20260316';
% log_name = 'LOG070.TXT.csv';
% file_path = fullfile(log_folder, flight_folder, log_name);

% --- Load and Process Flight Log Data ---
[para, Nheader, ind, ind_cntr] = extract_header_information(file_path);

% Load data from CSV or cached MAT file for faster processing
[folder, base, ~] = fileparts(file_path);
mat_path = fullfile(folder, base + ".mat");

try
  S = load(mat_path);
  data = S.data;
catch
  data = readmatrix(file_path, "NumHeaderLines", Nheader);
  save(mat_path, "data");
end

Ts      = para.looptime * 1.0e-6;             % Gyro loop period [s]
Ts_cntr = para.pid_process_denom * Ts;        % Control loop period [s]
Ts_log  = para.frameIntervalPDenom * Ts_cntr; % Logging period [s]

% Get a usable time vector
time = (data(:, ind.time) - data(1, ind.time)) * 1.0e-6;

% Extract debug signals
sinarg   = data(:, ind.debug(1)) / 5e3;
exec_c   = data(:, ind.debug(2));
fchirp   = data(:, ind.debug(3));
meas_alt = data(:, ind.debug(4));
set_alt  = data(:, ind.debug(5));

if (sinarg(1) ~= 0)
  disp("altered data")
  sinarg = fix_signal(sinarg);
  set_alt = fix_signal(set_alt);
end

idx = find(sinarg ~= 0, 1);
meas_alt = fix_startalt(meas_alt, idx);

% Extract motor signals and compute mean throttle
motor1 = data(:, ind.motor(1));
motor2 = data(:, ind.motor(2));
motor3 = data(:, ind.motor(3));
motor4 = data(:, ind.motor(4));
motor_mean = mean([motor1, motor2, motor3, motor4], 2);

% --- Figure 1: Raw signal inspection (full flight) ---
figure(1)
subplot(311)
plot(time, set_alt); grid on;
title('Set Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(312)
plot(time, meas_alt); grid on;
title('Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(313)
plot(time, motor1, time, motor2, time, motor3, time, motor4); grid on;
legend('M1', 'M2', 'M3', 'M4', 'Location', 'best');
title('Motor Commands');
xlabel('Time [s]'); ylabel('Throttle');
ylim([0 1000]);

% --- Select chirp evaluation window ---
idx = get_ind_eval(sinarg, meas_alt);

sinarg_ax = sinarg;
sinarg_ax(~idx) = 0;

% --- Figure 2: Chirp window signals ---
figure(2)
subplot(411)
plot(time(idx), set_alt(idx)); grid on;
title('Set Altitude (chirp window)');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(412)
plot(time(idx), meas_alt(idx)); grid on;
title('Measured Altitude (chirp window)');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(413)
plot(time(idx), motor_mean(idx)); grid on;
title('Mean Motor Throttle (chirp window)');
xlabel('Time [s]'); ylabel('Throttle');

subplot(414)
plot(time(idx), sinarg(idx)); grid on;
title('Chirp Signal');
xlabel('Time [s]'); ylabel('Amplitude');

%% Estimate Transfer Functions

% Welch parameters
Nest     = round(15 / Ts_log);       % Window length [samples]
Noverlap = floor(0.9 * Nest);        % Overlap between windows
window   = hann(Nest, 'periodic');   % Hanning window

% Low-pass filter for rotating-frame filtering
Dlp = sqrt(3) / 2;    % Damping ratio
wlp = 2 * pi * 10;    % Cutoff frequency [rad/s]
Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');

% Apply rotating-frame filter to all signals
inp   = apply_rotfiltfilt(Glp, sinarg_ax, set_alt);    % w: altitude setpoint
out   = apply_rotfiltfilt(Glp, sinarg_ax, meas_alt);   % y: altitude measurement
out_v = apply_rotfiltfilt(Glp, sinarg_ax, motor_mean); % v: mean motor throttle (plant input)

% Closed-loop altitude transfer function: set_alt → meas_alt
[T_ax, C_T_ax] = estimate_frequency_response(inp(idx), out(idx), window, Noverlap, Nest, Ts_log);

% Direct plant: motor_mean → meas_alt
[G_vy, C_G_vy] = estimate_frequency_response(out_v(idx), out(idx), window, Noverlap, Nest, Ts_log);

% Reference to plant input: set_alt → motor_mean
[G_wv, C_G_wv] = estimate_frequency_response(inp(idx), out_v(idx), window, Noverlap, Nest, Ts_log);

% Indirect estimates
P_alt = T_ax / G_wv;           % Indirect plant
C_alt = G_wv / (1 - T_ax);    % Indirect controller

%% Analytical Model

% Althold controller: P gain * low-pass filter
C_P_Alt = 15;       % Althold P gain
fc_ana  = 80;       % Filter cutoff [Hz]

f = T_ax.Frequency * 2*pi;
C_P_Alt_frd = frd(C_P_Alt * ones(size(f)), f, Ts_cntr);
Gf_ana = get_filter('pt1', fc_ana, Ts_cntr);

C_alt_ana = C_P_Alt_frd * Gf_ana;
C_alt_ana = downsample_frd(C_alt_ana, Ts_log, T_ax.Frequency);

T_ana = (C_alt_ana * G_vy) / (1 + C_alt_ana * G_vy);

%% Tuning

P_new  = 30;    % New P gain
fc_new = 80;    % New filter cutoff [Hz]

C_P_Alt_new_frd = frd(P_new * ones(size(f)), f, Ts_cntr);
Gf_new = get_filter('pt1', fc_new, Ts_cntr);

C_alt_new = C_P_Alt_new_frd * Gf_new;
C_alt_new = downsample_frd(C_alt_new, Ts_log, T_ax.Frequency);

T_new = (C_alt_new * G_vy) / (1 + C_alt_new * G_vy);

%% Plot Results

% Figure 3: Closed-loop system comparison
figure(3)
bode(T_ax, T_ana, T_new); grid on;
legend('Measured', 'Analytical', 'New', 'Location', 'best');
title('Altitude Hold Transfer Function');

% Figure 4: Plant (indirect vs direct)
figure(4)
bode(P_alt, G_vy); grid on;
legend('Indirect', 'Direct', 'Location', 'best');
title('Measured Plant (Throttle \rightarrow Altitude)');

% Figure 5: Controller (indirect vs analytical)
figure(5)
bode(C_alt, C_alt_ana); grid on;
legend('Indirect', 'Analytical', 'Location', 'best');
title('Altitude Controller');

% Figure 6: Coherence
figure(6)
bode(C_T_ax, C_G_vy, C_G_wv); grid on;
legend('T_{ax}', 'G_{vy}', 'G_{wv}', 'Location', 'best');
title('Coherence');

%% Step Response

fmax = 10;    % Hz — althold bandwidth is much lower than the rate loop

step_time = (0:Nest-1) .* Ts_log;
T_mean = 0.1 * [-1, 1] + (Nest * Ts_log) / 2;

step_resp = [calculate_step_response_from_frd(T_ax,  fmax), ...
  calculate_step_response_from_frd(T_ana, fmax), ...
  calculate_step_response_from_frd(T_new, fmax)];

step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2), :));
step_resp_norm = step_resp ./ step_resp_mean;

figure(7)
plot(step_time, step_resp_norm); grid on;
legend('Measured', 'Analytical', 'New', 'Location', 'best');
title('Step Response (Normalized)');
xlabel('Time [s]'); ylabel('Altitude [cm]');
