
% Initialize workspace
clc, clear variables, close all
addpath("../lib/");
addpath(genpath("../logs/"));

pos_bode = [0.1514+0.05, 0.5838-0.2, 0.7536, 0.3472+0.2; ...
                           0.1, 0.08,       0.7536+0.05, 0.2017]

% Add file information
log_folder = '../logs';
flight_folder = '20260423';
log_name = 'LOG005.TXT.csv';
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

% Derive sample times from the log header. Blackbox frame period is
% looptime * pid_process_denom * frameIntervalPDenom; althold updates at
% 100 Hz, so we decimate down to that rate and run the controller analysis
% at the same Ts.
Ts_log_raw = para.looptime * 1e-6 * para.pid_process_denom * para.frameIntervalPDenom;
dec        = round((1/Ts_log_raw) / 100);
data       = data(1:dec:end, :);
Ts_log     = Ts_log_raw * dec;
Ts_cntr    = Ts_log;   % althold control loop runs at logging rate

% ind.time column is in microseconds (Betaflight blackbox convention).
time = (data(:, ind.time) - data(1, ind.time)) * 1e-6;

% Undo firmware-side debug scaling so all signals are in physical units.
% Firmware logs (alt_hold_multirotor.c, autopilot_multirotor.c)
sinarg          = data(:,ind.debug(1)) / 5e3;
meas_alt        = data(:,ind.debug(2));
set_alt         = data(:,ind.debug(3));
meas_pi         = data(:,ind.debug(4));
throttle_out    = data(:,ind.debug(6)) / 1e3;   % normalised throttle 0..1
vertical_v      = data(:,ind.debug(7)) / 10;    % cm/s
throttle_offset = data(:,ind.debug(8));         % PWM units


%% Overview Plot

figure(1)
subplot(411)
plot(time, meas_alt, '-r'); hold on
plot(time, set_alt, '-b'); grid on;
legend('Measured Altitude', 'Target Altitude', 'Location', 'best');
title('Compare Target and Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(412)
plot(time, throttle_out); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle [0..1]');

subplot(413)
plot(time, vertical_v); grid on;
title('Vertical Speed');
xlabel('Time [s]'); ylabel('Speed [cm/s]');

subplot(414)
plot(time, throttle_offset); grid on;
title('Throttle Offset');
xlabel('Time [s]'); ylabel('Throttle [PWM units]');

sinarg = fix_signal(sinarg);
idx = get_ind_eval(sinarg, meas_alt);

%% Estimate Transfer Functions

% Welch parameters
frame    = 12.5;
Nest     = round(frame / (Ts_log));
Noverlap = floor(0.9 * Nest);
window   = hann(Nest, 'periodic');

% Low-pass filter for rotating-frame filtering
Dlp = sqrt(3) / 2;
wlp = 2 * pi * 10;
Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');

% Build sinarg mask (zero outside chirp window)
sinarg_ax = sinarg;
sinarg_ax(~idx) = 0;

% Rotate-filter input and output signals (full length, then crop to idx).
% Use throttle_offset (debug[7], PWM units) instead of throttle_out so the
% plant input units match the analytical controller output (PWM units),
% bypassing the firmware's mincheck/PWM_RANGE rescaling and tilt boost.
inp   = apply_rotfiltfilt(Glp, sinarg_ax, set_alt);
out_y = apply_rotfiltfilt(Glp, sinarg_ax, meas_alt);
out_u = apply_rotfiltfilt(Glp, sinarg_ax, throttle_offset);
out_v = apply_rotfiltfilt(Glp, sinarg_ax, meas_pi);

%% Estimation Transferfunction T
[T, C_T] = estimate_frequency_response(inp(idx), out_y(idx), ...
  window, Noverlap, Nest, Ts_log);
C_T_data = squeeze(C_T.ResponseData);
f_bode = squeeze(C_T.Frequency);
omega_bode = 2*pi*f_bode;

linewidth = 1.5;

figure(2)
ax(1) = subplot('Position', pos_bode(1,:));
opt = bodeoptions('cstprefs');
opt.MagUnits = 'dB';
opt.PhaseUnits = 'deg';
opt.FreqUnits = 'Hz';
opt.PhaseWrapping = 'on';
opt.PhaseWrappingBranch = -180;
opt.YLim = {[-60 20], [-180 180]};

bode(ax(1), T, 'k', omega_bode, opt)
title('Bode Plot Transfer Function')
grid(ax(1), 'on')

ax(2) = subplot('Position', pos_bode(2,:));
semilogx(ax(2), f_bode, C_T_data, 'k', 'LineWidth', linewidth)
grid(ax(2), 'on')
ylabel(ax(2), 'Coherence')
xlabel(ax(2), 'Frequency [Hz]')
ylim(ax(2), [0 1])
xlim(ax(2), [min(f_bode) max(f_bode)])

linkaxes(ax, 'x')

%% Estimation Plant

[Guw, C_uw] = estimate_frequency_response(inp(idx), out_u(idx), ...
  window, Noverlap, Nest, Ts_log);

P = T / Guw;

C_uw_data = squeeze(C_uw.ResponseData);
C_P = C_uw_data.* C_T_data;

figure(3)
ax(1) = subplot('Position', pos_bode(1,:));
opt = bodeoptions('cstprefs');
opt.MagUnits = 'dB';
opt.MagScale = 'linear';
opt.PhaseUnits = 'deg';
opt.FreqUnits = 'Hz';
opt.PhaseWrapping = 'on';
opt.YLim = {[-60 50], [-180 180]};

bode(ax(1), P, 'k', omega_bode, opt)   % bode bekommt weiter rad/s
title('Bode Plot Plant')
grid(ax(1), 'on')

ax(2) = subplot('Position', pos_bode(2,:));
semilogx(ax(2), f_bode, C_P, 'k', 'LineWidth', linewidth)
grid(ax(2), 'on')
ylabel(ax(2), 'Coherence')
xlabel(ax(2), 'Frequency [Hz]')
ylim(ax(2), [0 1])
xlim(ax(2), [min(f_bode) max(f_bode)])

linkaxes(ax, 'x')

%% Get the Controller PI

%Measured Transfer function PI
[Gvw, C_Gvw] = estimate_frequency_response(inp(idx), out_v(idx), ...
  window, Noverlap, Nest, Ts_log);

Cpi = Gvw / (1 - T);

% Analytical PI + D (current tune)
P_cur = 20;
I_cur = 20;
D_cur = 18;
fc_pt2_cur = 1;
[Cpi_ana, Cd_ana] = calculate_althold_controllers( ...
  P_cur,...
  I_cur,...
  D_cur,...
  fc_pt2_cur,...
  Ts_cntr,...
  Ts_log,...
  P.Frequency);

% Iterm Relax
% Cpi_com = Cpi / Cpi_ana;
% Cpi_ana = Cpi_ana * Cpi_com;

figure(4)
opt = bodeoptions('cstprefs');
opt.MagUnits = 'dB';
opt.MagScale = 'linear';
opt.PhaseUnits = 'deg';
opt.FreqUnits = 'Hz';
opt.PhaseWrapping = 'on';
bode(Cpi, Cpi_ana, opt);grid on;
title('Bode Plot PI Controller')
legend('Measured','Calculated')

%% Calculation of the D-Controller

%Measured Transfer function PI
Cd  = (Gvw - Guw) / T;
figure(5)
opt = bodeoptions('cstprefs');
opt.MagUnits = 'dB';
opt.MagScale = 'linear';
opt.PhaseUnits = 'deg';
opt.FreqUnits = 'Hz';
opt.PhaseWrapping = 'on';
bode(Cd, Cd_ana, opt);grid on;
title('Bode Plot DT2 Controller')
legend('Measured','Calculated')

%% Get Closed Loop Data

CL_ana= calculate_closed_loop(Cpi_ana, tf(1,1,Ts_log), P, tf(1,1,Ts_log), Cd_ana);

%% New Controller (proposed tune)

% When true, the "new" controller mirrors the current analytical one — useful
% as a sanity check on the first pass. Flip to false and edit P_new/I_new/
% D_new/fc_pt2_new below to propose a different tune.
default_parameters = true;

if default_parameters
  P_new = P_cur;
  I_new = I_cur;
  D_new = D_cur;
  fc_pt2_new = fc_pt2_cur;
else
  P_new = 22;
  I_new = 22;
  D_new = 20;
  fc_pt2_new = 1;
end

[Cpi_ana_new, Cd_ana_new] = calculate_althold_controllers( ...
  P_new,...
  I_new,...
  D_new,...
  fc_pt2_new,...
  Ts_cntr,...
  Ts_log,...
  P.Frequency);

CL_ana_new = calculate_closed_loop(Cpi_ana_new, tf(1,1,Ts_log), P, tf(1,1,Ts_log), Cd_ana_new);

figure(6)
opt = bodeoptions('cstprefs');
opt.MagUnits = 'dB';
opt.MagScale = 'linear';
opt.PhaseUnits = 'deg';
opt.FreqUnits = 'Hz';
opt.PhaseWrapping = 'on';
bode(T, CL_ana.T, omega_bode, opt)
xlim([min(f_bode) 10])
title('Comparison Measured to Analytical Transfer Function')
grid on
legend('Measured','Calculated')


figure(7)
opt = bodeoptions('cstprefs');
opt.MagUnits = 'dB';
opt.MagScale = 'linear';
opt.PhaseUnits = 'deg';
opt.FreqUnits = 'Hz';
opt.PhaseWrapping = 'on';

ax(1) = subplot(2,2,1);
bodemag(ax(1), CL_ana.T, CL_ana_new.T, T, omega_bode, opt);
title('Tracking T')
legend('Actual','New','Measured','Location','best')
grid on

ax(2) = subplot(2,2,2);
bodemag(ax(2), CL_ana.S, CL_ana_new.S, omega_bode, opt);
title('Sensitivity S')
legend('Actual','New','Location','northwest')
grid on

ax(3) = subplot(2,2,3);
bodemag(ax(3), CL_ana.SC, CL_ana_new.SC, omega_bode, opt);
title('Controller Effort SC')
legend('Actual','New','Location','northwest')
grid on

ax(4) = subplot(2,2,4);
bodemag(ax(4), CL_ana.SP, CL_ana_new.SP, omega_bode, opt);
title('Compliance SP')
legend('Actual','New','Location','southwest')
grid on

linkaxes(ax,'x');
xlim(ax(1), [min(f_bode) 10])
sgtitle('Gang of Four - Altitude Hold')

%% Get Step Response
% Step responses
f_max = 3;
T_mean = 0.1 * [-1, 1] + (Nest * Ts_log) / 2;
step_time = (0:Nest-1).'*Ts_log;

% Actual controller parameters
step_resp = [calculate_step_response_from_frd(CL_ana.T    , f_max), ...
  calculate_step_response_from_frd(CL_ana_new.T, f_max), ...
  calculate_step_response_from_frd(T           , f_max)];
step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2),:));
step_resp = step_resp ./ step_resp_mean;

figure(8)
plot(step_time, step_resp), grid on, ylabel('Altitude (cm)')
title('Step Response Altitude Hold')
legend('actual', 'new', 'measured', 'location', 'best')
ylim([-0.1 1.4]); xlim([0 frame/2]);

