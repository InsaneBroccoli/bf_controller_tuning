% Initialize workspace
clc, clear variables, close all
addpath("../lib/");
addpath(genpath("../logs/"));

% Add file information
log_folder = '../logs';
flight_folder = '20260402';
log_name = 'LOG093.TXT.csv';
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

time = (data(:, ind.time) - data(1, ind.time)) * 1e6;

sinarg          = data(:,ind.debug(1)) / 5e3;
meas_alt        = data(:,ind.debug(4));
set_alt         = data(:,ind.debug(5));
throttle_out    = data(:,ind.debug(6));
vertical_v      = data(:,ind.debug(7));
throttle_offset = data(:,ind.debug(8));

time            = time(1:20:end);
sinarg          = sinarg(1:20:end);
meas_alt        = meas_alt(1:20:end);
set_alt         = set_alt(1:20:end);
throttle_out    = throttle_out(1:20:end);
vertical_v      = vertical_v(1:20:end);
throttle_offset = throttle_offset(1:20:end);

figure(1)
subplot(411)
plot(time, meas_alt, '-r'); hold on
plot(time, set_alt, '-b'); grid on;
legend('Target Altitude', 'Measured Altitude', 'Location', 'best');
title('Compare Target and Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(412)
plot(time, throttle_out); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle');

subplot(413)
plot(time, throttle_out); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle');

subplot(414)
plot(time, vertical_v); grid on;
title('Vertical Speed');
xlabel('Time [s]'); ylabel('Speed [cm/s]');


idx = get_ind_eval(sinarg, meas_alt);
meas_alt = fix_offset(meas_alt, idx);

time_c            = time(idx);
meas_alt_c        = meas_alt(idx);
set_alt_c         = set_alt(idx);
throttle_out_c    = throttle_out(idx);
vertical_v_c      = vertical_v(idx);
throttle_offset_c = throttle_offset(idx);


figure(2)
subplot(411)
plot(time_c, meas_alt_c, '-r'); hold on
plot(time_c, set_alt_c, '-b'); grid on;
legend('Target Altitude', 'Measured Altitude', 'Location', 'best');
title('Compare Target and Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(412)
plot(time_c, throttle_out_c); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle');

subplot(413)
plot(time_c, throttle_out_c); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle');

subplot(414)
plot(time_c, vertical_v_c); grid on;
title('Vertical Speed');
xlabel('Time [s]'); ylabel('Speed [cm/s]');

Ts_log  = 1 / 100;
Ts_cntr = Ts_log;   % althold control loop runs at logging rate


%% Estimate Transfer Functions

% Welch parameters
Nest     = round(15 / Ts_log);
Noverlap = floor(0.9 * Nest);
window   = hann(Nest, 'periodic');

% Low-pass filter for rotating-frame filtering
Dlp = sqrt(3) / 2;
wlp = 2 * pi * 10;
Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');

% Build sinarg mask (zero outside chirp window)
sinarg_ax = sinarg;
sinarg_ax(~idx) = 0;

% Rotate-filter input and output signals (full length, then crop to idx)
inp   = apply_rotfiltfilt(Glp, sinarg_ax, set_alt);
out_y = apply_rotfiltfilt(Glp, sinarg_ax, meas_alt);
out_u = apply_rotfiltfilt(Glp, sinarg_ax, throttle_out);

% T: w -> y (complementary sensitivity)
[T, C_T] = estimate_frequency_response(inp(idx), out_y(idx), ...
  window, Noverlap, Nest, Ts_log);

% Guw: w -> u (reference to throttle output)
[Guw, C_Guw] = estimate_frequency_response(inp(idx), out_u(idx), ...
  window, Noverlap, Nest, Ts_log);

% Indirect plant: P_gef = T / Guw
P_gef = T / Guw;

% Coherence product
Coh = C_T * C_Guw;

% Frequency vector for Bode plots
omega_bode = 2 * pi * P_gef.Frequency;


%% Analytical Model

% Althold PI+D controller (D on measurement)
Kp_alt = 15 * 0.01;
Ki_alt = 15 * 0.003;
Kd_alt = 15 * 0.01;
fc_pt2 = 1;

PID_alt = [Kp_alt, Ki_alt, Kd_alt];
Gf_p = ss(tf(1, 1, Ts_cntr));
[Cpi_ana, Cd_ana] = calculate_controllers(PID_alt, Gf_p, Ts_cntr);

% D-term filter: PT2
Gd_ana = get_filter('pt2', fc_pt2, Ts_cntr);
Cd_ana = Cd_ana * Gd_ana;

% No gyro filter in althold
Gf_ana = ss(tf(1, 1, Ts_cntr));

% Downsample to logging rate
Cpi_ana = downsample_frd(Cpi_ana, Ts_log, P_gef.Frequency);
Cd_ana  = downsample_frd(Cd_ana,  Ts_log, P_gef.Frequency);
Gf_ana  = downsample_frd(Gf_ana,  Ts_log, P_gef.Frequency);

% Plant (Gf = 1 for althold, so P = P_gef)
P = P_gef / Gf_ana;

% Closed-loop
CL_ana = calculate_closed_loop(Cpi_ana, tf(1, 1, Ts_log), P, Gf_ana, Cd_ana);


%% Tuning

Kp_new = 12 * 0.01;
Ki_new = 12 * 0.003;
Kd_new = 12 * 0.01;
fc_new = 3;

PID_new = [Kp_new, Ki_new, Kd_new];
[Cpi_new, Cd_new] = calculate_controllers(PID_new, Gf_p, Ts_cntr);

Gd_new = get_filter('pt2', fc_new, Ts_cntr);
Cd_new = Cd_new * Gd_new;

Gf_new = ss(tf(1, 1, Ts_cntr));

Cpi_new = downsample_frd(Cpi_new, Ts_log, P_gef.Frequency);
Cd_new  = downsample_frd(Cd_new,  Ts_log, P_gef.Frequency);
Gf_new  = downsample_frd(Gf_new,  Ts_log, P_gef.Frequency);

CL_new = calculate_closed_loop(Cpi_new, tf(1, 1, Ts_log), P, Gf_new, Cd_new);


%% Plot Results

% Figure 3: Closed-loop Bode comparison
figure(3)
bode(T, CL_ana.T, CL_new.T); grid on;
legend('Measured', 'Analytical', 'New', 'Location', 'best');
title('Altitude Hold Transfer Function');

% Figure 4: Plant with coherence
figure(4)
subplot(211)
bodemag(P_gef); grid on;
title('Indirect Plant (Throttle -> Altitude)');
subplot(212)
semilogx(P_gef.Frequency, squeeze(abs(Coh.ResponseData))); grid on;
title('Coherence'); xlabel('Frequency [Hz]'); ylabel('Coherence');
ylim([0 1]);


%% Step Response

fmax = 10;

step_time = (0:Nest-1).' * Ts_log;
T_mean = 0.1 * [-1, 1] + (Nest * Ts_log) / 2;

% Tracking
step_resp = [calculate_step_response_from_frd(CL_ana.T, fmax), ...
  calculate_step_response_from_frd(CL_new.T, fmax), ...
  calculate_step_response_from_frd(T, fmax)];

step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2), :));
step_resp_tra = step_resp ./ step_resp_mean;

figure(5)
plot(step_time, step_resp_tra); grid on;
legend('Analytical', 'New', 'Measured', 'Location', 'best');
title('Step Response - Tracking (Normalized)');
xlabel('Time [s]'); ylabel('Altitude [cm]');

% Disturbance rejection
step_resp_d = [calculate_step_response_from_frd(CL_ana.SP, fmax), ...
  calculate_step_response_from_frd(CL_new.SP, fmax)];

step_resp_d_mean = mean(step_resp_d(step_time > T_mean(1) & step_time < T_mean(2), :));
step_resp_d = step_resp_d - step_resp_d_mean;

figure(6)
plot(step_time, step_resp_d); grid on;
legend('Analytical', 'New', 'Location', 'best');
title('Step Response - Disturbance Rejection');
xlabel('Time [s]'); ylabel('Altitude [cm]');
