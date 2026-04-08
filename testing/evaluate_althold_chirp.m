%==========================================================================
% ALTHOLD CHIRP TUNING - Betaflight Altitude Hold Controller Analysis
%==========================================================================
% Purpose: Frequency response estimation and controller tuning for
%          altitude hold using chirp excitation signals.
%
% Signal definitions (matching gyro_ctrl_tuning convention):
%   w: reference input (targetAlt)
%   y: system output (measuredAlt)
%   u: control input (throttleOut, total PID output)
%   v: velocity feedback (verticalVelocity, used by D-term)
%
% See althold-chirp-debug-signals.md and althold_diagram.jpeg for details.
%==========================================================================

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

% Debug signal column indices (see althold-chirp-debug-signals.md)
%   [0] sinarg           (x 5e3)
%   [3] measuredAlt      [cm]       -> y
%   [4] targetAlt        [cm]       -> w
%   [5] throttleOut      (x 1000)   -> u
%   [6] verticalVelocity (x 10)     -> v
sinarg          = ind.debug(1);
meas_alt        = ind.debug(4);
set_alt         = ind.debug(5);
throttle_out    = ind.debug(6);
vertical_v      = ind.debug(7);
throttle_offset = ind.debug(8);

% --- Figure 1: Raw signal inspection (full flight) ---
figure(1)
subplot(311)
plot(time, data(:, set_alt), '-b'); hold on
plot(time, data(:, meas_alt), '-r'); grid on;
legend('Target Altitude', 'Measured Altitude', 'Location', 'best');
title('Compare Target and Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(312)
plot(time, data(:, throttle_out), '-b'); grid on;
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
tdata(:, throttle_out) = tdata(:, throttle_out);
tdata(:, vertical_v)   = tdata(:, vertical_v);

% Fix initial values if needed
if (tdata(1, sinarg) ~= 0)
  disp("altered data")
  tdata(:, sinarg)  = fix_signal(tdata(:, sinarg));
  tdata(:, set_alt) = fix_signal(tdata(:, set_alt));
end

% Create logical mask for desired chirp window
ind_eval = get_ind_eval(tdata(:, sinarg), tdata(:, meas_alt));
idx = ind_eval;

sinarg_ax = tdata(:, sinarg);
sinarg_ax(~ind_eval) = 0;

% --- Figure 2: Chirp window signals ---
figure(2)
subplot(311)
plot(time(idx), tdata(idx, set_alt), '-b'); hold on
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

% ----- Input signal: w (filtered setpoint) -----
w = tdata(:, set_alt);
inp = apply_rotfiltfilt(Glp, sinarg_ax, w);

% ----- Output y: filtered altitude measurement -----
y = tdata(:, meas_alt);
out_y = apply_rotfiltfilt(Glp, sinarg_ax, y);

% T: w -> y (complementary sensitivity)
[T, C_T] = estimate_frequency_response(inp(idx), ...
  out_y(idx), window, Noverlap, Nest, Ts_log);

% ----- Control signal u: filtered throttle output -----
u = tdata(:, throttle_out);
out_u = apply_rotfiltfilt(Glp, sinarg_ax, u);

% Guw: w -> u (reference to total controller output)
[Guw, C_Guw] = estimate_frequency_response(inp(idx), out_u(idx), ...
  window, Noverlap, Nest, Ts_log);

% ----- Velocity signal v: filtered vertical velocity -----
v = tdata(:, vertical_v);
out_v = apply_rotfiltfilt(Glp, sinarg_ax, v);

% Gvw: w -> v (reference to velocity)
[Gvw, C_Gvw] = estimate_frequency_response(inp(idx), ...
  out_v(idx), window, Noverlap, Nest, Ts_log);

% Indirect plant: P_gef = T / Guw (throttle -> altitude, includes filters)
P_gef = T / Guw;

% Total controller (indirect)
C_total = Guw / (1 - T);

% Direct estimates
% G_vy: v -> y (velocity -> altitude)
[G_vy, C_G_vy] = estimate_frequency_response(out_v(idx), ...
  out_y(idx), window, Noverlap, Nest, Ts_log);

% G_cy: u -> y (direct plant: throttle -> altitude)
[G_cy, C_G_cy] = estimate_frequency_response(out_u(idx), ...
  out_y(idx), window, Noverlap, Nest, Ts_log);

% G_cv: u -> v (direct: throttle -> velocity)
[G_cv, C_G_cv] = estimate_frequency_response(out_u(idx), out_v(idx), ...
  window, Noverlap, Nest, Ts_log);

% Coherence product
Coh = C_T * C_Guw;

% Frequency vector for Bode plots
omega_bode = 2*pi * P_gef.Frequency;


%% Analytical Model

% Althold PI+D controller (D on measurement, see althold_diagram.jpeg)
%   Forward: Kp * error  (+ Ki * integral(error))
%   D feedback: Kd * dBoost * PT2(dy/dt), dBoost = 1
Kp_alt = 15 * 0.01;        % P gain
Ki_alt = 15 * 0.003;         % I gain
Kd_alt = 15;         % D gain (TODO: set correct value)
fc_pt2 = 80;        % PT2 filter cutoff [Hz] (D-term)

% Build PI and D controllers (matching calculate_controllers structure)
PID_alt = [Kp_alt, Ki_alt, Kd_alt];
Gf_p = ss(tf(1, 1, Ts_cntr));                    % No P-term filter for althold
[Cpi_ana, Cd_ana] = calculate_controllers(PID_alt, Gf_p, Ts_cntr);

% D-term filter: PT2 (see althold_diagram.jpeg)
Gd_ana = get_filter('pt2', fc_pt2, Ts_cntr);
Cd_ana = Cd_ana * Gd_ana;

% No gyro filter in althold (measurement goes directly to controller)
Gf_ana = ss(tf(1, 1, Ts_cntr));

% Downsample to logging rate
Cpi_ana = downsample_frd(Cpi_ana, Ts_log, P_gef.Frequency);
Cd_ana  = downsample_frd(Cd_ana,  Ts_log, P_gef.Frequency);
Gf_ana  = downsample_frd(Gf_ana,  Ts_log, P_gef.Frequency);

% Plant (P = P_gef / Gf, but Gf = 1 for althold)
P = P_gef / Gf_ana;

% Closed-loop (2-DOF structure matching Betaflight)
CL_ana = calculate_closed_loop(Cpi_ana, tf(1,1,Ts_log), P, Gf_ana, Cd_ana);


%% Tuning

Kp_new = 30 * 0.01;        % New P gain
Ki_new = 30 * 0.003;         % New I gain
Kd_new = 30;         % New D gain (TODO)
fc_new = 80;        % New PT2 cutoff [Hz]

PID_new = [Kp_new, Ki_new, Kd_new];
[Cpi_new, Cd_new] = calculate_controllers(PID_new, Gf_p, Ts_cntr);

Gd_new = get_filter('pt2', fc_new, Ts_cntr);
Cd_new = Cd_new * Gd_new;

Gf_new = ss(tf(1, 1, Ts_cntr));

Cpi_new = downsample_frd(Cpi_new, Ts_log, P_gef.Frequency);
Cd_new  = downsample_frd(Cd_new,  Ts_log, P_gef.Frequency);
Gf_new  = downsample_frd(Gf_new,  Ts_log, P_gef.Frequency);

CL_new = calculate_closed_loop(Cpi_new, tf(1,1,Ts_log), P, Gf_new, Cd_new);


%% Plot Results

% Figure 3: Closed-loop comparison
figure(3)
bode(T, CL_ana.T, CL_new.T); grid on;
legend('Measured', 'Analytical', 'New', 'Location', 'best');
title('Altitude Hold Transfer Function');

% Figure 4: Plant (indirect vs direct)
figure(4)
bode(P_gef, G_cy); grid on;
legend('Indirect', 'Direct', 'Location', 'best');
title('Measured Plant (Throttle -> Altitude)');

% Figure 5: Controller (indirect vs analytical)
figure(5)
bode(C_total, Cpi_ana); grid on;
legend('Indirect', 'Analytical', 'Location', 'best');
title('Altitude Controller');

% Figure 6: Velocity plant (indirect vs direct)
figure(6)
bode(Gvw / Guw, G_cv); grid on;
legend('Indirect', 'Direct', 'Location', 'best');
title('Measured Plant (Throttle -> Velocity)');


%% Step Response

fmax = 10;    % Hz -- althold bandwidth is much lower than the rate loop

step_time = (0:Nest-1).' * Ts_log;
T_mean = 0.1 * [-1, 1] + (Nest * Ts_log) / 2;

% Tracking performance
step_resp = [calculate_step_response_from_frd(CL_ana.T, fmax), ...
  calculate_step_response_from_frd(CL_new.T, fmax), ...
  calculate_step_response_from_frd(T,        fmax)];

step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2), :));
step_resp_tra = step_resp ./ step_resp_mean;

figure(7)
plot(step_time, step_resp_tra); grid on;
legend('Analytical', 'New', 'Measured', 'Location', 'best');
title('Step Response - Tracking (Normalized)');
xlabel('Time [s]'); ylabel('Altitude [cm]');

% Disturbance rejection
step_resp_d = [calculate_step_response_from_frd(CL_ana.SP, fmax), ...
  calculate_step_response_from_frd(CL_new.SP, fmax)];

step_resp_d_mean = mean(step_resp_d(step_time > T_mean(1) & step_time < T_mean(2), :));
step_resp_d = step_resp_d - step_resp_d_mean;

figure(8)
plot(step_time, step_resp_d); grid on;
legend('Analytical', 'New', 'Location', 'best');
title('Step Response - Disturbance Rejection');
xlabel('Time [s]'); ylabel('Altitude [cm]');
