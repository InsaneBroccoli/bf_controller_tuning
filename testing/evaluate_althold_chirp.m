% Initialize workspace
clc, clear variables, close all
addpath("../lib/");
addpath(genpath("../logs/"));

% Add file information
log_folder = '../logs';
flight_folder = '20260331';
log_name = 'LOG090.TXT.csv';
file_path = fullfile(log_folder, flight_folder, log_name);

% --- Load and Process Flight Log Data ---
[para, Nheader, ind, ind_cntr] = extract_header_information(file_path);

% Load data from CSV or cached MAT file for faster processing
[folder, base, ~] = fileparts(file_path);
mat_path = fullfile(folder, [base '.mat']);

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

% Define debug signals (see althold-chirp-debug-signals.md)
%   [0] sinarg           - chirp window          (x 5e3)
%   [1] exc              - excitation             (x 1000)
%   [2] fchirp           - chirp frequency        (x 100)
%   [3] measuredAlt      - measured altitude       [cm]
%   [4] targetAlt        - target altitude         [cm]
%   [5] throttleOut      - controller output       (x 1000)
%   [6] verticalVelocity - vertical speed          (x 10)
%   [7] throttleOffset   - throttle offset
sinarg       = ind.debug(1);
meas_alt     = ind.debug(4);
set_alt      = ind.debug(5);
throttle_out = ind.debug(6);
vertical_v   = ind.debug(7);

% --- Figure 1: Raw signal inspection (full flight) ---
figure(1)
subplot(311)
plot(time, data(:, set_alt), '-b');
hold on
plot(time, data(:, meas_alt), '-r'); grid on;
legend('Target Altitude', 'Measured Altitude', 'Location', 'best');
title('Compare Target and Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(312)
plot(time, data(:, throttle_out) / 1000, '-b'); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle');

subplot(313)
plot(time, data(:, vertical_v) / 10, '-b'); grid on;
title('Vertical Speed');
xlabel('Time [s]'); ylabel('Speed [cm/s]');


% Scale data into working copy
tdata(:,:) = data(:,:);

tdata(:, sinarg)       = tdata(:, sinarg) / 5e3;
tdata(:, meas_alt)     = tdata(:, meas_alt);           % already in cm
tdata(:, set_alt)      = tdata(:, set_alt);             % already in cm
tdata(:, throttle_out) = tdata(:, throttle_out) / 1000;
tdata(:, vertical_v)   = tdata(:, vertical_v) / 10;

% Fix initial values if needed
if (tdata(1, sinarg) ~= 0)
  disp("altered data")
  tdata(:, sinarg)  = fix_signal(tdata(:, sinarg));
  tdata(:, set_alt) = fix_signal(tdata(:, set_alt));
end

% first_nonzero = find(tdata(:, sinarg) ~= 0, 1);
% tdata(:, meas_alt) = fix_startalt(tdata(:, meas_alt), first_nonzero);

% Create logical mask for desired chirp window
ind_eval = get_ind_eval(tdata(:, sinarg), tdata(:, meas_alt));
idx = ind_eval;

sinarg_ax = tdata(:, sinarg);
sinarg_ax(~ind_eval) = 0;

% --- Figure 2: Chirp window signals ---
figure(2)
subplot(311)
plot(time(idx), tdata(idx, set_alt), '-b');
hold on
plot(time(idx), tdata(idx, meas_alt), '-r'); grid on;
legend('Target Altitude', 'Measured Altitude', 'Location', 'best');
title('Compare Target and Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(312)
plot(time(idx), tdata(idx, throttle_out), '-b'); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle');

subplot(313)
plot(time(idx), tdata(idx, vertical_v), '-b'); grid on;
title('Vertical Speed');
xlabel('Time [s]'); ylabel('Speed [cm/s]');


%% Estimate Transfer Functions

% Welch parameters
Nest     = round(15 / Ts_log);          % Window length [samples]
Noverlap = floor(0.9 * Nest);           % Overlap between windows
window   = hann(Nest, 'periodic');      % Hanning window

% Low-pass filter for rotating-frame filtering (apply_rotfiltfilt)
Dlp = sqrt(3) / 2;    % Damping ratio
wlp = 2 * pi * 10;    % Cutoff frequency [rad/s]
Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');

% Apply rotating-frame filter to all signals
w = tdata(:, set_alt);
inp = apply_rotfiltfilt(Glp, sinarg_ax, w);

y = tdata(:, meas_alt);
out = apply_rotfiltfilt(Glp, sinarg_ax, y);

% Closed-loop: set_alt -> meas_alt
[T_ax, C_T_ax] = estimate_frequency_response(inp(idx), ...
  out(idx), window, Noverlap, Nest, Ts_log);

% Controller output path: set_alt -> throttle_out
c = tdata(:, throttle_out);
out_c = apply_rotfiltfilt(Glp, sinarg_ax, c);

[G_wc, C_G_wc] = estimate_frequency_response(inp(idx), out_c(idx), ...
  window, Noverlap, Nest, Ts_log);

% Indirect estimates
P_alt  = T_ax / G_wc;          % Plant: throttle -> altitude
C_alt  = G_wc / (1 - T_ax);   % Controller

% Velocity path: set_alt -> vertical_v
v = tdata(:, vertical_v);
out_v = apply_rotfiltfilt(Glp, sinarg_ax, v);

[G_wv, C_G_wv] = estimate_frequency_response(inp(idx), ...
  out_v(idx), window, Noverlap, Nest, Ts_log);

% Plant decomposition via velocity
P_vel = G_wv / G_wc;           % Throttle -> velocity (indirect)

% Direct plant: throttle_out -> meas_alt
[G_cy, C_G_cy] = estimate_frequency_response(out_c(idx), ...
  out(idx), window, Noverlap, Nest, Ts_log);

% Direct throttle -> velocity
[G_cv, C_G_cv] = estimate_frequency_response(out_c(idx), out_v(idx), ...
  window, Noverlap, Nest, Ts_log);

%% Analytical Model

% Althold controller: P gain * low-pass filter
C_P_Alt = 15;       % Althold P gain
fc_ana  = 15;       % Filter cutoff [Hz]

f = T_ax.Frequency * 2*pi;
C_P_Alt_frd = frd(C_P_Alt * ones(size(f)), f, Ts_cntr);
Gf_ana = get_filter('pt1', fc_ana, Ts_cntr);

C_alt_ana = C_P_Alt_frd * Gf_ana;
C_alt_ana = downsample_frd(C_alt_ana, Ts_log, T_ax.Frequency);

% Closed-loop using direct plant (single-loop structure)
T_ana = (C_alt_ana * G_cy) / (1 + C_alt_ana * G_cy);

%% Tuning

P_new  = 30;    % New P gain
fc_new = 30;    % New filter cutoff [Hz]

C_P_Alt_new_frd = frd(P_new * ones(size(f)), f, Ts_cntr);
Gf_new = get_filter('pt1', fc_new, Ts_cntr);

C_alt_new = C_P_Alt_new_frd * Gf_new;
C_alt_new = downsample_frd(C_alt_new, Ts_log, T_ax.Frequency);

T_new = (C_alt_new * G_cy) / (1 + C_alt_new * G_cy);

%% Plot Results

% Figure 3: Closed-loop system comparison
figure(3)
bode(T_ax, T_ana, T_new); grid on;
legend('Measured', 'Analytical', 'New', 'Location', 'best');
title('Altitude Hold Transfer Function');

% Figure 4: Plant (indirect vs direct)
figure(4)
bode(P_alt, G_cy); grid on;
legend('Indirect', 'Direct', 'Location', 'best');
title('Measured Plant (Throttle -> Altitude)');

% Figure 5: Controller (indirect vs analytical)
figure(5)
bode(C_alt, C_alt_ana); grid on;
legend('Indirect', 'Analytical', 'Location', 'best');
title('Altitude Controller');

% Figure 6: Velocity plant (indirect vs direct)
figure(6)
bode(P_vel, G_cv); grid on;
legend('Indirect', 'Direct', 'Location', 'best');
title('Measured Plant (Throttle -> Velocity)');

%% Step Response

fmax = 10;    % Hz -- althold bandwidth is much lower than the rate loop

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
