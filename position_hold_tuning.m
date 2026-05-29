
%% Initialize workspace
clc, clear, close all
addpath("./lib/");
addpath(genpath("./logs/"));

pos_bode = [0.1514, 0.5838-0.2, 0.7536, 0.3472+0.2; ...
            0.1514, 0.1100    , 0.7536, 0.1917    ];

set(cstprefs.tbxprefs, 'MagnitudeUnits', 'dB');
set(cstprefs.tbxprefs, 'FrequencyUnits', 'Hz');
set(cstprefs.tbxprefs, 'PhaseUnits',     'deg');
set(cstprefs.tbxprefs, 'UnwrapPhase',    'Off');
set(cstprefs.tbxprefs, 'Grid',           'On');

opt = bodeoptions('cstprefs');
opt.MagScale      = 'linear';
opt.PhaseWrapping = 'on';

linewidth = 1.5;

%% Load File

% Add file information
log_folder = './logs';
flight_folder = '20260529';
log_name = '20260528_6_inch_1.TXT.csv';
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


%% Get Log / Control Loop Time

% Decimate the full matrix once (blackbox runs at the FC loop rate; althold
% debug values update at 100 Hz, so 20:1 brings the effective rate to ~100 Hz
% while keeping a single decimation factor in one place).
Ts_log_raw = para.looptime * 1e-6 * para.pid_process_denom * para.frameIntervalPDenom;
dec        = round((1/Ts_log_raw) / 10);
data       = data(1:dec:end, :);
Ts_log     = Ts_log_raw * dec;
Ts_cntr    = Ts_log;   % poshold control loop runs at logging rate (GPS-limited)

% ind.time column is in microseconds (Betaflight blackbox convention).
time = (data(:, ind.time) - data(1, ind.time)) * 1e-6;

%% Scaling Debugs

% Undo firmware-side debug scaling so all signals are in physical units
gps_error       = data(:,ind.debug(1)) / 10;    % gps error
angle_target    = data(:,ind.debug(2)) / 10;    % target angle
chirp           = data(:,ind.debug(3)) / 10;    % cm (injected position setpoint)
current_angle   = data(:,ind.debug(4)) / 10;    % current angle in BF [deg]
pid_sum_EF      = data(:,ind.debug(5)) / 10;    % PID_sum Earth Frame
sinarg          = data(:,ind.debug(6)) / 5e3;   % Injected Chirp Signal
active_axis     = data(:,ind.debug(7)) * 2;     % Low = LON/ROLL, High = LAT/PITCH
pidDA_limit     = data(:,ind.debug(8)) / 10;    % deg

% calculate current and target position from gps error and chirp
target_position  = chirp;
current_position = chirp - gps_error;

%% DA Limit Plot

% Plot pidDA Limit, to see if clipping occurs
figure(10)
plot(time, pidDA_limit); grid on;
yline(35, '--', sprintf('D+A vector saturation limit'), ...
      'LineWidth', 2, 'LabelHorizontalAlignment', 'right', ...
      'LabelVerticalAlignment', 'bottom');
title('Combined D+A Vector length with 35° tilt Limit');
legend('combined D+A vector length', Location='northwest');
xlabel('Time [s]'); ylabel('Angle [deg]');
xlim([time(1) time(end)]);

%% Overview Plots

figure(1)
subplot(211)
plot(time, gps_error, '-b'); hold on;
grid on;
title('GPS error');
xlabel('Time [s]');
ylabel('Error [cm]');
xlim([time(1) time(end)]);

% Angle Plots
subplot(212)
plot(time, pid_sum_EF, 'r'); hold on;
title('Angles');

plot(time, angle_target, '-b'); 
plot(time, current_angle, 'g');
grid on;
legend('PID Sum EF', 'Target Angle BF', 'Current Angle BF')
xlabel('Time [s]'); ylabel('Angle [deg]');
xlim([time(1) time(end)]);

% Chirp Plots
figure(2)
subplot(211)
plot(time, sinarg, '-m'); hold on;
plot(time, active_axis, '-k');grid on
legend('sinarg', 'axis (low = roll, high = pitch)');
title('Chirp Excitation');
ylabel('Excitation Signal');
xlim([time(1) time(end)]);
ylim([-1 7]);

subplot(212)
plot(time, chirp, '-m'); %hold on;
% plot(time, chirp_inst_freq * 50, '-k');
grid on
legend('chirp exc');
xlabel('Time [s]'); ylabel('Amplitude [cm]');
xlim([time(1) time(end)]);
ylim([-100 100]);

% Find the active chirp window
idx = get_ind_eval(sinarg, chirp);

% -------------------------------------------------------------------------
% Axis index masks
%   active_axis < 1  -> Roll  (LON)
%   active_axis >= 1 -> Pitch (LAT)
% -------------------------------------------------------------------------
idx_roll  = idx & (active_axis < 1);
idx_pitch = idx & (active_axis >= 1);



%% Estimate Transfer Functions

% Welch parameters
frame = 15;
Nest     = round(frame / (Ts_log));
Noverlap = floor(0.9 * Nest);
window   = hann(Nest, 'periodic');

% Low-pass filter for rotating-frame filtering
Dlp = sqrt(3) / 2;
wlp = 2 * pi * 1;
Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');

% Build sinarg mask (zero outside chirp window)
sinarg_ax = sinarg;
sinarg_ax(~idx) = 0;

% =========================================================================
% Rotating-frame filtering -- separate per axis
% =========================================================================

% --- Roll ---
sinarg_roll = sinarg_ax .* double(active_axis < 1);
inp_r   = apply_rotfiltfilt(Glp, sinarg_roll, target_position);
out_y_r = apply_rotfiltfilt(Glp, sinarg_roll, current_position);
out_u_r = apply_rotfiltfilt(Glp, sinarg_roll, pid_sum_EF);
out_v_r = apply_rotfiltfilt(Glp, sinarg_roll, angle_target);
out_q_r = apply_rotfiltfilt(Glp, sinarg_roll, current_angle);

% --- Pitch ---
sinarg_pitch = sinarg_ax .* double(active_axis >= 1);
inp_p   = apply_rotfiltfilt(Glp, sinarg_pitch, target_position);
out_y_p = apply_rotfiltfilt(Glp, sinarg_pitch, current_position);
out_u_p = apply_rotfiltfilt(Glp, sinarg_pitch, pid_sum_EF);
out_v_p = apply_rotfiltfilt(Glp, sinarg_pitch, angle_target);
out_q_p = apply_rotfiltfilt(Glp, sinarg_roll, current_angle);


%% Estimation Transferfunction T
% current position - target position

[T_roll, C_T_roll] = estimate_frequency_response(inp_r(idx_roll), ...
    out_y_r(idx_roll), window, Noverlap, Nest, Ts_log);

[T_pitch, C_T_pitch] = estimate_frequency_response(inp_p(idx_roll), ...
    out_y_p(idx_roll), window, Noverlap, Nest, Ts_log);

f_bode = squeeze(C_T_roll.Frequency);

C_T_roll = squeeze(C_T_roll.ResponseData);
C_T_pitch = squeeze(C_T_pitch.ResponseData);
omega_bode = 2*pi*f_bode;

figure(3)
ax(1) = subplot('Position', pos_bode(1,:));
bode(ax(1), T_roll, 'k',T_pitch, '-r', omega_bode, opt);
legend('Roll', 'Pitch', 'Location','best')
title('Bode Plot Transfer Function');

ax(2) = subplot('Position', pos_bode(2,:));
semilogx(ax(2), f_bode, C_T_roll, 'k', f_bode,C_T_pitch,'-r', 'LineWidth', linewidth);
grid(ax(2), 'on');
ylabel(ax(2), 'Coherence');
xlabel(ax(2), 'Frequency [Hz]');
ylim(ax(2), [0 1]);
xlim(ax(2), [f_bode(2,:) max(f_bode)]);

linkaxes(ax, 'x');

%% Get the Controller PID
% PID sum - target position

[Guw_roll, C_Guw_roll] = estimate_frequency_response(inp_r(idx_roll), ...
    out_u_r(idx_roll),  window, Noverlap, Nest, Ts_log);

[Guw_pitch, C_Guw_pitch] = estimate_frequency_response(inp_p(idx_pitch), ...
    out_u_p(idx_pitch),  window, Noverlap, Nest, Ts_log);

Cpid_roll = Guw_roll / (1 - T_roll);
Cpid_pitch = Guw_pitch / (1 - T_pitch);

C_Guw_roll = squeeze(C_Guw_roll.ResponseData);
C_Guw_pitch = squeeze(C_Guw_pitch.ResponseData);

C_P_roll = C_Guw_roll.* C_T_roll;
C_P_pitch = C_Guw_pitch.* C_T_pitch;

freq_vector = Guw_roll.Frequency;

% Analytical PID (standard parameters)
P_cur = 30;
I_cur = 30;
D_cur = 30;
A_cur = 30;
fc_pt1_cur = 0.8;

[Cpid_ana] = calculate_poshold_controller( ...
  P_cur,I_cur, D_cur, A_cur, fc_pt1_cur, Ts_cntr, Ts_log, freq_vector);

figure(4)
bode(Cpid_roll, 'k', omega_bode, opt); hold on;
bode(Cpid_pitch, 'r', omega_bode, opt); hold on;
bode(Cpid_ana, 'b', omega_bode, opt);
xlim([f_bode(2,:) max(f_bode)]);
legend('measured roll', 'measured pitch' , 'analytical');
title('Bode Plot PID Controller');

%% Estimation Plant

P_total_roll = T_roll / Guw_roll;    % Estimation Plant with T_Angle included
P_total_pitch = T_pitch / Guw_pitch;    % Estimation Plant with T_Angle included

figure(5)
ax(1) = subplot('Position', pos_bode(1,:));
bode(ax(1), P_total_roll, 'k',P_total_pitch, 'r', omega_bode, opt);
legend('Roll', 'Pitch', 'Location','best')
title('Bode Plot Plant with T-Angle');

ax(2) = subplot('Position', pos_bode(2,:));
semilogx(ax(2), f_bode, C_P_roll, 'k', f_bode, C_P_pitch, 'r' ,'LineWidth', linewidth);
grid(ax(2), 'on');
ylabel(ax(2), 'Coherence');
xlabel(ax(2), 'Frequency [Hz]');
ylim(ax(2), [0 1]);
xlim(ax(2), [f_bode(2,:) max(f_bode)]);

linkaxes(ax, 'x');

% Estimate Plant (current_angle -> current_position)
[Gqw_roll, C_Gqw_roll] = estimate_frequency_response(inp_r(idx_roll), ...
    out_q_r(idx_roll), window, Noverlap, Nest, Ts_log);
[Gqw_pitch, C_Gqw_pitch] = estimate_frequency_response(inp_p(idx_pitch), ...
    out_q_p(idx_pitch), window, Noverlap, Nest, Ts_log);

C_Gqw_roll = squeeze(C_Gqw_roll.ResponseData);
C_Gqw_pitch = squeeze(C_Gqw_pitch.ResponseData);

P_roll = T_roll / Gqw_roll;      % Just Plant (currentAngle-->currentPosition)
P_pitch = T_pitch / Gqw_pitch;      % Just Plant (currentAngle-->currentPosition)

% Plot
figure(55)
ax(1) = subplot('Position', pos_bode(1,:));
bode(ax(1), P_roll, 'k', P_pitch, 'r', omega_bode, opt);
title('Plant Decomposition');
legend('Roll', 'Pitch');

ax(2) = subplot('Position', pos_bode(2,:));
semilogx(ax(2), f_bode, C_Gqw_roll, 'k', 'LineWidth', linewidth); hold on;
semilogx(ax(2), f_bode, C_Gqw_pitch, 'r', 'LineWidth', linewidth);
grid(ax(2), 'on');
ylabel(ax(2), 'Coherence');
xlabel(ax(2), 'Frequency [Hz]');
ylim(ax(2), [0 1]);
xlim(ax(2), [f_bode(2,:) max(f_bode)]);

linkaxes(ax, 'x');

%% Get closed Loop Data

CL_ana_roll = calculate_closed_loop(Cpid_ana, tf(1,1,Ts_log), P_total_roll, tf(1,1,Ts_log), tf(0,1));
CL_ana_pitch = calculate_closed_loop(Cpid_ana, tf(1,1,Ts_log), P_total_pitch, tf(1,1,Ts_log), tf(0,1));

%% New Controller
default_parameters = false;

if default_parameters
  P_new = P_cur;
  I_new = I_cur;
  D_new = D_cur;
  A_new = A_cur;
  fc_pt1_new = fc_pt1_cur;
else
  P_new = 50;
  I_new = 36;
  D_new = 33;
  A_new = 25;
  fc_pt1_new = 0.8;
end

[Cpid_ana_new] = calculate_poshold_controller( ...
  P_new, I_new, D_new, A_new, fc_pt1_new, Ts_cntr, Ts_log, P_total_pitch.Frequency);

CL_new_roll = calculate_closed_loop(Cpid_ana_new, tf(1,1,Ts_log), P_total_roll, tf(1,1,Ts_log), tf(0,1));
CL_new_pitch = calculate_closed_loop(Cpid_ana_new, tf(1,1,Ts_log), P_total_pitch, tf(1,1,Ts_log), tf(0,1));

%% Gang of Four Plot Roll

figure(6)
bode(T_roll, CL_ana_roll.T, omega_bode, opt);
xlim(ax(1), [0.0667 5]);
title('Comparison Measured to Analytical Transfer Function - Roll');
grid on;
legend('Measured','Calculated');

figure(7)
ax(1) = subplot(2,2,1);
bodemag(ax(1), CL_ana_roll.T, CL_new_roll.T, T_roll, omega_bode, opt);
title('Tracking T');
legend('Actual','New','Measured','Location','best');
grid on;

ax(2) = subplot(2,2,2);
bodemag(ax(2), CL_ana_roll.S, CL_new_roll.S, omega_bode, opt);
title('Sensitivity S')
legend('Actual','New','Location','northwest')
grid on

ax(3) = subplot(2,2,3);
bodemag(ax(3), CL_ana_roll.SC, CL_new_roll.SC, omega_bode, opt);
title('Controller Effort SC');
legend('Actual','New','Location','northwest');
grid on;

ax(4) = subplot(2,2,4);
bodemag(ax(4), CL_ana_roll.SP, CL_new_roll.SP, omega_bode, opt);
title('Compliance SP');
legend('Actual','New','Location','southwest');
grid on;

linkaxes(ax,'x');
xlim(ax(1), [0.0667 5]);
sgtitle('Gang of Four - Position Hold - Roll');

%% Gang of Four Plot Pitch

figure(8)
bode(T_pitch, CL_ana_pitch.T, omega_bode, opt);
xlim(ax(1), [0.0667 5]);
title('Comparison Measured to Analytical Transfer Function - Pitch');
grid on;
legend('Measured','Calculated');

figure(9)
ax(1) = subplot(2,2,1);
bodemag(ax(1), CL_ana_pitch.T, CL_new_pitch.T, T_pitch, omega_bode, opt);
title('Tracking T');
legend('Actual','New','Measured','Location','best');
grid on;

ax(2) = subplot(2,2,2);
bodemag(ax(2), CL_ana_pitch.S, CL_new_pitch.S, omega_bode, opt);
title('Sensitivity S')
legend('Actual','New','Location','northwest')
grid on

ax(3) = subplot(2,2,3);
bodemag(ax(3), CL_ana_pitch.SC, CL_new_pitch.SC, omega_bode, opt);
title('Controller Effort SC');
legend('Actual','New','Location','northwest');
grid on;

ax(4) = subplot(2,2,4);
bodemag(ax(4), CL_ana_pitch.SP, CL_new_pitch.SP, omega_bode, opt);
title('Compliance SP');
legend('Actual','New','Location','southwest');
grid on;

linkaxes(ax,'x');
xlim(ax(1), [0.0667 5]);
sgtitle('Gang of Four - Position Hold - Pitch');

%% Get Step Response

f_max = 2.5;
T_mean = 0.1 * [-1, 1] + (Nest * Ts_log) / 2;

step_time = (0:Nest-1).' * Ts_log;

step_resp = [calculate_step_response_from_frd(CL_new_pitch.T,     f_max), ...
             calculate_step_response_from_frd(CL_new_roll.T, f_max), ...
             calculate_step_response_from_frd(T_pitch, f_max), ...
             calculate_step_response_from_frd(T_roll,            f_max)];

% Normalization
step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2), :));
step_resp = step_resp ./ step_resp_mean;

figure(10)
plot(step_time, step_resp)
grid on;
xlabel('Time [s]');
ylabel('Position [cm]');
title('Step Response Position Hold');
legend('New Pitch', 'New Roll', 'Measured Pitch','Measured Roll' , 'location', 'best');
xlim([0, frame/2]);
ylim([-0.3, 1.8]);
