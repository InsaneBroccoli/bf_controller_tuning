
% Initialize workspace
clc, clear variables, close all
addpath("../lib/");
addpath(genpath("../logs/"));

s = tf('s');

pos_bode = [0.1514, 0.5838-0.2, 0.7536, 0.3472+0.2; ... % this is a bit hacky
            0.1514, 0.1100    , 0.7536, 0.1917    ];
opt = bodeoptions('cstprefs');


% Add file information
log_folder = '../logs';
flight_folder = '20260420';
log_name = 'LOG000.TXT.csv';
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

idx = get_ind_eval(sinarg, meas_alt);


Ts_log  = 1 / 100;  % DO WE HAVE THIS IN THE LOG FILE?
Ts_cntr = Ts_log;   % althold control loop runs at logging rate SURE?

%% Estimate Transfer Functions

% Welch parameters
Nest     = round(15 / (Ts_log));
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
out_d = apply_rotfiltfilt(Glp, sinarg_ax, meas_d);

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
opt.MagScale = 'linear';
opt.PhaseUnits = 'deg';
opt.FreqUnits = 'Hz';
opt.PhaseWrapping = 'on';
opt.YLim = {[-60 20], [-180 180]};

bode(ax(1), T, 'k', omega_bode, opt)   % bode bekommt weiter rad/s
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
semilogx(ax(2), omega_bode, C_P, 'k', 'LineWidth', linewidth)
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

% Calculated Transfer function PI

Kp_alt = 15 * 0.01;
Ki_alt = 15 * 0.003;

Cpi_ana = ss(Kp_alt + Ki_alt*Ts_cntr*tf([1 0],[1 -1], Ts_cntr));

Cpi_ana  = downsample_frd(Cpi_ana , Ts_log, P.Frequency);

% Iterm Relax
Cpi_com = Cpi / Cpi_ana;
Cpi_ana = Cpi_ana * Cpi_com;

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
Cd  = Guw * Gvw / T * (1 / Guw - 1 / Gvw);

Kd_alt = 15 * 0.01;
fc_pt2 = 1;

Cd_only = ss( Kd_alt/Ts_cntr*tf([1 -1], [1 0], Ts_cntr) );
Gf_D_Term = get_filter('pt2', fc_pt2, Ts_cntr);

Cd_ana = Cd_only * Gf_D_Term;

Cd_ana  = downsample_frd(Cd_ana , Ts_log, P.Frequency);

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

figure(6)
opt = bodeoptions('cstprefs');
opt.MagUnits = 'dB';
opt.MagScale = 'linear';
opt.PhaseUnits = 'deg';
opt.FreqUnits = 'Hz';
opt.PhaseWrapping = 'on';
bode(T, CL_ana.T, omega_bode, opt)
xlim([1e-2 10])
title('Comparison Measured to Analytical Transfer Function')
grid on
legend('Measured','Calculated')

%% Get Step Response
% Step responses
f_max = 3;
T_mean = 0.1 * [-1, 1] + (Nest * Ts_log) / 2;
step_time = (0:Nest-1).'*Ts_log;

% Actual controller parameters
step_resp = [calculate_step_response_from_frd(CL_ana.T    , f_max), ...
             calculate_step_response_from_frd(T           , f_max)];
step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2),:));
step_resp = step_resp ./ step_resp_mean;

figure(7)
plot(step_time, step_resp), grid on, ylabel('Altitude (cm)')
title('Step Response Altitude Hold')
legend('actual', 'measured', 'location', 'best')
ylim([-0.5 1.8]); xlim([0 7.5]);

