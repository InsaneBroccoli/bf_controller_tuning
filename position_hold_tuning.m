
% Initialize workspace
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

% Add file information
log_folder = './logs';
flight_folder = '20260527';
log_name = '20260527_6_inch_1.csv';
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
Ts_log_raw = para.looptime * 1e-6 * para.pid_process_denom * para.frameIntervalPDenom;
dec        = round((1/Ts_log_raw) / 100);
data       = data(1:dec:end, :);
Ts_log     = Ts_log_raw * dec;
Ts_cntr    = Ts_log;   % althold control loop runs at logging rate

% ind.time column is in microseconds (Betaflight blackbox convention).
time = (data(:, ind.time) - data(1, ind.time)) * 1e-6;

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

inp   = apply_rotfiltfilt(Glp, sinarg_ax, target_position);
out_q = apply_rotfiltfilt(Glp, sinarg_ax, current_angle);
out_y = apply_rotfiltfilt(Glp, sinarg_ax, current_position);
out_u = apply_rotfiltfilt(Glp, sinarg_ax, pid_sum_EF);

out_v = apply_rotfiltfilt(Glp, sinarg_ax, angle_target);

%% Estimation Transferfunction T
% current position - target position

[T, C_T] = estimate_frequency_response(inp(idx), out_y(idx), ...
  window, Noverlap, Nest, Ts_log);
C_T_data = squeeze(C_T.ResponseData);
f_bode = squeeze(C_T.Frequency);
omega_bode = 2*pi*f_bode;

figure(3)
ax(1) = subplot('Position', pos_bode(1,:));
bode(ax(1), T, 'k', omega_bode, opt);
title('Bode Plot Transfer Function');

ax(2) = subplot('Position', pos_bode(2,:));
semilogx(ax(2), f_bode, C_T_data, 'k', 'LineWidth', linewidth);
grid(ax(2), 'on');
ylabel(ax(2), 'Coherence');
xlabel(ax(2), 'Frequency [Hz]');
ylim(ax(2), [0 1]);
xlim(ax(2), [f_bode(2,:) max(f_bode)]);

linkaxes(ax, 'x');

%% Get the Controller PID
% PID sum - target position

[Guw, C_Guw] = estimate_frequency_response(inp(idx), out_u(idx), ...
  window, Noverlap, Nest, Ts_log);

Cpid = Guw / (1 - T);

C_uw_data = squeeze(C_Guw.ResponseData);
C_P = C_uw_data.* C_T_data;

freq_vector = Guw.Frequency;

% Analytical PID (standard parameters)
P_cur = 30;
I_cur = 30;
D_cur = 30;
A_cur = 30;
fc_pt1_cur = 0.8;

[Cpid_ana] = calculate_poshold_controller( ...
  P_cur,...
  I_cur,...
  D_cur,...
  A_cur,...
  fc_pt1_cur,...
  Ts_cntr,...
  Ts_log,...
  freq_vector);

figure(4)
bode(Cpid, 'k', omega_bode, opt); hold on;
bode(Cpid_ana, 'r', omega_bode, opt);
xlim([f_bode(2,:) max(f_bode)]);
legend('measured', 'analytical');
title('Bode Plot PID Controller');

%% Estimation Plant

P = T / Guw;

figure(5)
ax(1) = subplot('Position', pos_bode(1,:));
bode(ax(1), P, 'k', omega_bode, opt);
title('Bode Plot Plant');

ax(2) = subplot('Position', pos_bode(2,:));
semilogx(ax(2), f_bode, C_P, 'k', 'LineWidth', linewidth);
grid(ax(2), 'on');
ylabel(ax(2), 'Coherence');
xlabel(ax(2), 'Frequency [Hz]');
ylim(ax(2), [0 1]);
xlim(ax(2), [f_bode(2,:) max(f_bode)]);

linkaxes(ax, 'x');

% Estimate Outer Plant (current_angle -> current_position)
[Gqw, C_Gqw] = estimate_frequency_response(inp(idx), out_q(idx), ...
  window, Noverlap, Nest, Ts_log);

C_qw_data = squeeze(C_Gqw.ResponseData);

P_outer = T / Gqw;

% Plot
figure(55)
ax(1) = subplot('Position', pos_bode(1,:));
bode(ax(1), P, 'k', P_outer, 'b', omega_bode, opt);
title('Plant Decomposition');
legend('P total (T/Guw)', 'P outer (kinematics)');

ax(2) = subplot('Position', pos_bode(2,:));
semilogx(ax(2), f_bode, C_P,          'k', 'LineWidth', linewidth); hold on;
semilogx(ax(2), f_bode, C_qw_data,    'r', 'LineWidth', linewidth);
grid(ax(2), 'on');
ylabel(ax(2), 'Coherence');
xlabel(ax(2), 'Frequency [Hz]');
ylim(ax(2), [0 1]);
xlim(ax(2), [f_bode(2,:) max(f_bode)]);
legend(ax(2), 'C_P total', 'C_P outer');

linkaxes(ax, 'x');

%% 
% angle_target - current position

% [T3, C_T3] = estimate_frequency_response(out_v(idx), out_v(idx), ...
%   window, Noverlap, Nest, Ts_log);
% C_T3_data = squeeze(C_T3.ResponseData);
% 
% figure(6)
% opt = bodeoptions('cstprefs');
% opt.MagUnits = 'dB';
% opt.MagScale = 'linear';
% opt.PhaseUnits = 'deg';
% opt.FreqUnits = 'Hz';
% opt.PhaseWrapping = 'on';
% bode(T3, 'k', omega_bode, opt);
% title('Angle Target - gpsSol.llh (current position)');


%% Get closed Loop Data

CL_ana = calculate_closed_loop(Cpid_ana, tf(1,1,Ts_log), P, tf(1,1,Ts_log), tf(0,1));

%% New Controller
default_parameters = false;

if default_parameters
  P_new = P_cur;
  I_new = I_cur;
  D_new = D_cur;
  A_new = A_cur;
  fc_pt1_new = fc_pt1_cur;
else
  P_new = 42;
  I_new = 80;
  D_new = 55;
  A_new = 10;
  fc_pt1_new = 0.8;
end

[Cpid_ana_new] = calculate_poshold_controller( ...
  P_new,...
  I_new,...
  D_new,...
  A_new,...
  fc_pt1_new,...
  Ts_cntr,...
  Ts_log,...
  P.Frequency);

CL_ana_new = calculate_closed_loop(Cpid_ana_new, tf(1,1,Ts_log), P, tf(1,1,Ts_log), tf(0,1));

figure(6)
bode(T, CL_ana.T, omega_bode, opt);
xlim([3e-2 100]);
title('Comparison Measured to Analytical Transfer Function');
grid on;
legend('Measured','Calculated');

figure(7)
ax(1) = subplot(2,2,1);
bodemag(ax(1), CL_ana.T, CL_ana_new.T, T, omega_bode, opt);
title('Tracking T');
legend('Actual','New','Measured','Location','best');
grid on;

ax(2) = subplot(2,2,2);
bodemag(ax(2), CL_ana.S, CL_ana_new.S, omega_bode, opt);
title('Sensitivity S')
legend('Actual','New','Location','northwest')
grid on

ax(3) = subplot(2,2,3);
bodemag(ax(3), CL_ana.SC, CL_ana_new.SC, omega_bode, opt);
title('Controller Effort SC');
legend('Actual','New','Location','northwest');
grid on;

ax(4) = subplot(2,2,4);
bodemag(ax(4), CL_ana.SP, CL_ana_new.SP, omega_bode, opt);
title('Compliance SP');
legend('Actual','New','Location','southwest');
grid on;

linkaxes(ax,'x');
xlim(ax(1), [3e-2 100]);
sgtitle('Gang of Four - Position Hold');

%% Get Step Response

f_max = 2.5;

% Nur sauberen Mittelteil anzeigen (wrap-around am Ende abschneiden)
t_plot_max = 10;  % s

% Normierung auf eingeschwungenen Bereich verschieben (weg von t_mid)
T_mean = [6, 9];  % s – dort sind analytical curves eingeschwungen

step_time = (0:Nest-1).' * Ts_log;

step_resp = [calculate_step_response_from_frd(CL_ana.T,     f_max), ...
             calculate_step_response_from_frd(CL_ana_new.T, f_max), ...
             calculate_step_response_from_frd(T,            f_max)];

% Normierung
step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2), :));
step_resp = step_resp ./ step_resp_mean;

% Plot nur bis t_plot_max
idx_plot = step_time <= t_plot_max;

figure(8)
plot(step_time(idx_plot), step_resp(idx_plot,:))
grid on;
xlabel('Time [s]');
ylabel('Position [norm.]');
title('Step Response Position Hold');
legend('actual', 'new', 'measured', 'location', 'best');
xlim([0, t_plot_max]);
ylim([-0.3, 1.8]);
