%% Start
clc, clear variables, close all
addpath(genpath('lib'));
addpath logs/

% Add file information
log_folder = 'logs';
flight_folder = '20260302';
log_name = '20260302_flipmini.csv';
file_path = fullfile(log_folder, flight_folder, log_name);


% --- Load and Process Flight Log Data ---
[para, Nheader, ind, ind_cntr] = extract_header_information(file_path);

% Load data from CSV or cached MAT file for faster processing
% robust MAT path
[folder, base, ~] = fileparts(file_path);
mat_path = fullfile(folder, [base '.mat']);

try
    S = load(mat_path);
    data = S.data;
catch
    data = readmatrix(file_path, 'NumHeaderLines', Nheader);
    save(mat_path, 'data');
end

Ts= para.looptime * 1.0e-6;             % Gyro loop
Ts_cntr = para.pid_process_denom * Ts;        % Control loop
Ts_log  = para.frameIntervalPDenom * Ts_cntr; % Logging loop   

% Get a useable Timevector
time = (data(:,ind.time) - data(1,ind.time)) * 1.0e-6;

% Define numbers of axis
axis = 1;       %roll = 1; pitch = 2;

% Define debugs
sinarg = ind.debug(1);

% Scaling of data
data(:, sinarg) = data(:, sinarg) / 5e3;
data(:, ind.heading(1:3)) = data(:, ind.heading) * 100;
data(:, ind.setpoint(1:3)) = data(:, ind.setpoint(1:3)) /10;
data(:, ind.gyroADC(1:3)) = data(:, ind.gyroADC(1:3)) / 10;

idx = get_ind_eval(data(:,sinarg), data(:,ind.gyroADC(axis)));
sinarg_ax = data(:, sinarg);
sinarg_ax(~idx) = 0;

%% Check Signal during Chirp

figure(1)
subplot(311)
plot(time, sinarg_ax)

subplot(312)
plot(time(idx), data(idx, ind.gyroADC(axis)))

subplot(313)
plot(time(idx), data(idx, ind.setpoint(axis)))
xlabel('Time [s]')

%% Lets Measure the Transferfunction and such kind of things

% Pararmeters
Nest = round(15 / Ts_log);               % Window length in samples
Noverlap = floor(0.9 * Nest);           % Overlap between windows
window   = hann(Nest, 'periodic');  % Hanning window for analysis

% Design linear filter for zero phase excitation filter (apply_rotfiltfilt)
Dlp = sqrt(3) / 2;    % Damping ratio
wlp = 2 * pi * 10;    % Cutoff frequency [rad/s]
Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), ...    % Discrete filter
        Ts_log, 'tustin');                   % Using Tustin transform

% Calculation of the Transferfunction of the Gyro Loop
w = data(:, ind.setpoint(axis));
inp = apply_rotfiltfilt(Glp, sinarg_ax, w);

y = data(:, ind.gyroADC(axis));
out = apply_rotfiltfilt(Glp, sinarg_ax, y);

[T_gyro, C_T_gyro] = estimate_frequency_response(inp(idx), ...
    out(idx), window, Noverlap, Nest, Ts_log);

% Calculation of the Transferfunction of the Angle Plant

z = data(:, ind.heading(axis));
out_angle = apply_rotfiltfilt(Glp, sinarg_ax, z);

[P_angle, C_T_angle] = estimate_frequency_response(out(idx), ...
    out_angle(idx), window, Noverlap, Nest, Ts_log);

%% Plot Transfer Function + Coherence Plant Angle

figure(2)
opt = bode_plot_options('dB', 'linear', 'deg', 'Hz');

ax(1) = subplot(2,1,1);
bode(ax(1), P_angle, 'k', opt);
title(ax(1), 'Plant P (Angle)');

ax(2) = subplot(2,1,2);
f_Hz = P_angle.Frequency / (2*pi);

ax(2) = subplot(2,1,2);

opt_coh = bode_plot_options('abs', 'linear', 'deg', 'Hz');
bodemag(ax(2), C_T_angle, '-k', opt_coh);

ylabel(ax(2), 'Coherence (abs)');
ylim(ax(2), [0 1]);
title(ax(2), '');

linkaxes(ax, 'x');
set(findall(gcf, 'type', 'line'), 'LineWidth', 1.5);

%% Plot Transfer Function + Coherence Gyro Loop

figure(3)
opt = bode_plot_options('dB', 'linear', 'deg', 'Hz');

ax(1) = subplot(2,1,1);
bode(ax(1), T_gyro, 'k', opt);
title(ax(1), 'Plant P (Angle)');

ax(2) = subplot(2,1,2);
f_Hz = P_angle.Frequency / (2*pi);

ax(2) = subplot(2,1,2);

opt_coh = bode_plot_options('abs', 'linear', 'deg', 'Hz');
bodemag(ax(2), C_T_gyro, '-k', opt_coh);

ylabel(ax(2), 'Coherence (abs)');
ylim(ax(2), [0 1]);
title(ax(2), '');

linkaxes(ax, 'x');
set(findall(gcf, 'type', 'line'), 'LineWidth', 1.5);

%% Lets Calculate the Transferfunctions

% Get the Controller
C_P_Angle = 50;   % 0.1 is scaling in Betaflight
f  = T_gyro.Frequency*2*pi;            % Frequency

C_P_Angle_frd = frd( C_P_Angle * 0.1 *ones(size(f)), f, Ts_cntr );
fc = 50;  % Hz  (dein "Angle PT3")
angle_fcut_ana = fc;
Gf_ana = get_filter('pt3', fc, Ts_cntr);  % should be discrete

C_Angle_ana = C_P_Angle_frd * Gf_ana;
C_Angle_ana = downsample_frd(C_Angle_ana, Ts_log, T_gyro.Frequency);

% Calculate the Transferfunction with this Controller
T_ana = (C_Angle_ana * T_gyro * P_angle)/(1+ C_Angle_ana * T_gyro * P_angle);

%% Lets start Tuning

% Add the Parameters
P_new = 100;            % Gain 
angle_lpf_hz = 50;      % Cut off Frequency

C_P_Angle_new_frd = frd( P_new * 0.1 *ones(size(f)), f, Ts_cntr );
Gf_new = get_filter('pt3', angle_lpf_hz, Ts_cntr);  % should be discrete

C_Angle_new = C_P_Angle_new_frd * Gf_ana;
C_Angle_new = downsample_frd(C_Angle_new, Ts_log, T_gyro.Frequency);

% Calculate the new Transferfunction with the new Parameters
T_new = (C_Angle_new * T_gyro * P_angle)/(1+ C_Angle_new * T_gyro * P_angle);

%% Lets Plot the figures
% Bodeplot of the hole system
figure(4)
bode(T_ana, T_new);grid on;
legend('Calculated','New Calculated','Location','best')
title('Measured Transferfunction System');


