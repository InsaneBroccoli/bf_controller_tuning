
% Initialize workspace
clc, clear variables, close all
addpath("../lib/");
addpath(genpath("../logs/"));

% Add file information
log_folder = '../logs';
flight_folder = '20260423';
log_name = 'LOG012.TXT.csv';
file_path = fullfile(log_folder, flight_folder, log_name);

% --- Load and Process Flight Log Data ---
[para, Nheader, ind] = extract_header_information(file_path);

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

% Decimate the full matrix once (blackbox runs at the FC loop rate; althold
% debug values update at 100 Hz, so 20:1 brings the effective rate to ~100 Hz
% while keeping a single decimation factor in one place).
data = data(1:20:end, :);

% ind.time column is in microseconds (Betaflight blackbox convention).
time = (data(:, ind.time) - data(1, ind.time)) * 1e-6;

% Undo firmware-side debug scaling so all signals are in physical units.
% Firmware logs (alt_hold_multirotor.c, autopilot_multirotor.c)
sinarg          = data(:,ind.debug(1)) / 5e3;
meas_alt        = data(:,ind.debug(2));
set_alt         = data(:,ind.debug(3));
meas_pi         = data(:,ind.debug(4));
meas_d          = data(:,ind.debug(5));
throttle_out    = data(:,ind.debug(6)) / 1e3;   % normalised throttle 0..1
vertical_v      = data(:,ind.debug(7)) / 10;    % cm/s
throttle_offset = data(:,ind.debug(8));         % PWM units

figure(1)
subplot(511)
plot(time, meas_alt, '-r'); hold on
plot(time, set_alt, '-b'); grid on;
legend('Measured Altitude', 'Target Altitude', 'Location', 'best');
title('Compare Target and Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(512)
plot(time, throttle_out); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle [0..1]');

subplot(513)
plot(time, vertical_v); grid on;
title('Vertical Speed');
xlabel('Time [s]'); ylabel('Speed [cm/s]');

subplot(514)
plot(time, throttle_offset); grid on;
title('Throttle Offset');
xlabel('Time [s]'); ylabel('Throttle [PWM units]');

subplot(515)
plot(time, sinarg); grid on;
title('Sinarg');

subplot(616)
plot(time, sinarg); grid on;
title('Sinarg');

sinarg = fix_signal(sinarg);
idx = get_ind_eval(sinarg, meas_alt);

figure(2)
plot(time, idx)

time_c            = time(idx);
meas_alt_c        = meas_alt(idx);
set_alt_c         = set_alt(idx);
throttle_out_c    = throttle_out(idx);
vertical_v_c      = vertical_v(idx);
throttle_offset_c = throttle_offset(idx);

figure(3)
subplot(411)
plot(time_c, meas_alt_c, '-r'); hold on
plot(time_c, set_alt_c, '-b'); grid on;
legend('Measured Altitude', 'Target Altitude', 'Location', 'best');
title('Compare Target and Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(412)
plot(time_c, throttle_out_c); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle [0..1]');

subplot(413)
plot(time_c, vertical_v_c); grid on;
title('Vertical Speed');
xlabel('Time [s]'); ylabel('Speed [cm/s]');

subplot(414)
plot(time_c, throttle_offset_c); grid on;
title('Throttle Offset');
xlabel('Time [s]'); ylabel('Throttle [PWM units]');

Ts_log  = 1 / 100;
Ts_cntr = Ts_log;   % althold control loop runs at logging rate


% %% Estimate Transfer Functions
%
% % Welch parameters
% Nest     = round(20 / (Ts_log));
% Noverlap = floor(0.9 * Nest);
% window   = hann(Nest, 'periodic');
%
% % Low-pass filter for rotating-frame filtering
% Dlp = sqrt(3) / 2;
% wlp = 2 * pi * 10;
% Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');
%
% % Build sinarg mask (zero outside chirp window)
% sinarg_ax = sinarg;
% sinarg_ax(~idx) = 0;
%
% % Rotate-filter input and output signals (full length, then crop to idx).
% % Use throttle_offset (debug[7], PWM units) instead of throttle_out so the
% % plant input units match the analytical controller output (PWM units),
% % bypassing the firmware's mincheck/PWM_RANGE rescaling and tilt boost.
% inp   = apply_rotfiltfilt(Glp, sinarg_ax, set_alt);
% out_y = apply_rotfiltfilt(Glp, sinarg_ax, meas_alt);
% out_u = apply_rotfiltfilt(Glp, sinarg_ax, throttle_offset);
%
% % T: w -> y (complementary sensitivity)
% [T, C_T] = estimate_frequency_response(inp(idx), out_y(idx), ...
%   window, Noverlap, Nest, Ts_log);
%
% % Guw: w -> u (reference to throttle offset, PWM units)
% [Guw, C_Guw] = estimate_frequency_response(inp(idx), out_u(idx), ...
%   window, Noverlap, Nest, Ts_log);
%
% % Indirect plant: P_gef = T / Guw  (throttle_offset [PWM] -> altitude [cm])
% P_gef = T / Guw;
%
% % Coherence product
% Coh = C_T * C_Guw;
%
% % Frequency vector for Bode plots
% omega_bode = 2 * pi * P_gef.Frequency;
%
%
% %% Analytical Model
%
% % Althold PID controller (autopilot_multirotor.c::altitudeControl):
% %   altitudeP = (target - measured) * Kp           -> P on error
% %   altitudeI += (target - measured) * Ki * dt     -> I on error
% %   altitudeD = verticalVelocity     * Kd * filter -> D on measurement
% % F-term is omitted: targetVelocity is not updated by altChirpUpdate, so
% % Kf * targetVelocity = 0 throughout the chirp.
% % iTermRelax (gain x0.1 above 200 cm error) and the +-200 PWM I clamp are
% % non-linear and not modelled; keep altChirpAmpl <= ~150 cm to stay linear.
% Kp_alt = 15 * 0.01;
% Ki_alt = 15 * 0.003;
% Kd_alt = 15 * 0.01;
% fc_pt2 = 1;
%
% % Discrete integrator from firmware: I(k) = I(k-1) + Ki*dt*e(k)
% Ki_int_ana = Ki_alt * Ts_cntr * tf([1 0], [1 -1], Ts_cntr);
%
% % Backward-difference D on measurement, then PT2 (matches altitude_d_lpf):
% %   Kd/Ts * (1 - z^-1) = Kd/Ts * (z - 1)/z
% Cd_ana_raw = ss(Kd_alt / Ts_cntr * tf([1 -1], [1 0], Ts_cntr));
% Gd_pt2 = get_filter('pt2', fc_pt2, Ts_cntr);
%
% % 2-DOF split for "P+I on error, D on measurement":
% %   u = (Kp + Ki/s) * r  -  (Kp + Ki/s + Kd*s*PT2) * y
% Co_ana = ss(Kp_alt + Ki_int_ana);
% Gd_ana = Kp_alt + Ki_int_ana + Cd_ana_raw * Gd_pt2;
%
% % Alias for calculate_closed_loop: Co -> Cpi_ana, Gd -> Cd_ana
% Cpi_ana = Co_ana;
%
% % No gyro filter in althold
% Gf_ana = ss(tf(1, 1, Ts_cntr));
%
% % Downsample to logging rate
% Cpi_ana = downsample_frd(Cpi_ana, Ts_log, P_gef.Frequency);
% Cd_ana  = downsample_frd(Gd_ana,  Ts_log, P_gef.Frequency);
% Gf_ana  = downsample_frd(Gf_ana,  Ts_log, P_gef.Frequency);
%
% % Plant (Gf = 1 for althold, so P = P_gef)
% P = P_gef / Gf_ana;
%
% % Closed-loop
% CL_ana = calculate_closed_loop(Cpi_ana, tf(1, 1, Ts_log), P, Gf_ana, Cd_ana);
%
%
% %% Tuning
%
% Kp_new = 12 * 0.01;
% Ki_new = 12 * 0.003;
% Kd_new = 12 * 0.01;
% fc_new = 3;
%
% % Same 2-DOF split as the analytical model: P+I on error, D on measurement
% Ki_int_new = Ki_new * Ts_cntr * tf([1 0], [1 -1], Ts_cntr);
% Cd_new_raw = ss(Kd_new / Ts_cntr * tf([1 -1], [1 0], Ts_cntr));
% Gd_pt2_new = get_filter('pt2', fc_new, Ts_cntr);
%
% Co_new = ss(Kp_new + Ki_int_new);
% Gd_new = Kp_new + Ki_int_new + Cd_new_raw * Gd_pt2_new;
%
% Cpi_new = Co_new;
%
% Gf_new = ss(tf(1, 1, Ts_cntr));
%
% Cpi_new = downsample_frd(Cpi_new, Ts_log, P_gef.Frequency);
% Cd_new  = downsample_frd(Gd_new,  Ts_log, P_gef.Frequency);
% Gf_new  = downsample_frd(Gf_new,  Ts_log, P_gef.Frequency);
%
% CL_new = calculate_closed_loop(Cpi_new, tf(1, 1, Ts_log), P, Gf_new, Cd_new);
%
%
% %% Plot Results
%
% % Figure 3: Closed-loop Bode comparison
% figure(3)
% bode(T, CL_ana.T, CL_new.T); grid on;
% legend('Measured', 'Analytical', 'New', 'Location', 'best');
% title('Altitude Hold Transfer Function');
%
% % Figure 4: Plant with coherence
% figure(4)
% subplot(211)
% bodemag(P_gef); grid on;
% title('Indirect Plant (Throttle -> Altitude)');
% subplot(212)
% semilogx(P_gef.Frequency, squeeze(abs(Coh.ResponseData))); grid on;
% title('Coherence'); xlabel('Frequency [Hz]'); ylabel('Coherence');
% ylim([0 1]);
%
%
% %% Gang of Four
%
% % S/SC/SP have no measured counterpart in this 2-DOF setup (T + S != 1,
% % see lib/calculate_closed_loop.m), so those panels are analytical vs new only.
% figure(7)
%
% ax_gof(1) = subplot(2, 2, 1);
% bodemag(CL_ana.T, CL_new.T, T); grid on;
% title('Tracking T  (w \rightarrow y)');
% legend('Analytical', 'New', 'Measured', 'Location', 'best');
%
% ax_gof(2) = subplot(2, 2, 2);
% bodemag(CL_ana.S, CL_new.S); grid on;
% title('Sensitivity S  (1 / (1+L))');
% legend('Analytical', 'New', 'Location', 'northwest');
%
% ax_gof(3) = subplot(2, 2, 3);
% bodemag(CL_ana.SC, CL_new.SC); grid on;
% title('Controller Effort SC  (n \rightarrow u)');
% legend('Analytical', 'New', 'Location', 'northwest');
%
% ax_gof(4) = subplot(2, 2, 4);
% bodemag(CL_ana.SP, CL_new.SP); grid on;
% title('Compliance SP  (d \rightarrow y)');
% legend('Analytical', 'New', 'Location', 'southwest');
%
% linkaxes(ax_gof, 'x');
% sgtitle('Gang of Four - Altitude Hold');
%
%
% %% Step Response
%
% fmax = 10;
%
% step_time = (0:Nest-1).' * Ts_log;
% T_mean = 0.1 * [-1, 1] + (Nest * Ts_log) / 2;
%
% % Tracking
% step_resp = [calculate_step_response_from_frd(CL_ana.T, fmax), ...
%   calculate_step_response_from_frd(CL_new.T, fmax), ...
%   calculate_step_response_from_frd(T, fmax)];
%
% step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2), :));
% step_resp_tra = step_resp ./ step_resp_mean;
%
% figure(5)
% plot(step_time, step_resp_tra); grid on;
% legend('Analytical', 'New', 'Measured', 'Location', 'best');
% title('Step Response - Tracking (Normalized)');
% xlabel('Time [s]'); ylabel('Altitude [cm]');
% xlim([0 10]);
%
% % Disturbance rejection
% step_resp_d = [calculate_step_response_from_frd(CL_ana.SP, fmax), ...
%   calculate_step_response_from_frd(CL_new.SP, fmax)];
%
% step_resp_d_mean = mean(step_resp_d(step_time > T_mean(1) & step_time < T_mean(2), :));
% step_resp_d = step_resp_d - step_resp_d_mean;
%
% figure(6)
% plot(step_time, step_resp_d); grid on;
% legend('Analytical', 'New', 'Location', 'best');
% title('Step Response - Disturbance Rejection');
% xlabel('Time [s]'); ylabel('Altitude [cm]');
% xlim([0 10]);
