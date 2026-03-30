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

Ts      = 1 / 100;     % Althold loop period [s]
Ts_cntr = Ts;          % Control loop period [s]
Ts_log  = 1 / 2000;   % Logging period [s]

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

% Extract motor signals
motor1 = data(:, ind.motor(1));
motor2 = data(:, ind.motor(2));
motor3 = data(:, ind.motor(3));
motor4 = data(:, ind.motor(4));

% --- Figure 1: Raw signal inspection (full flight) ---
figure(1)
% subplot(411)
% plot(time, sinarg); grid on;
% title('Chirp Signal (sinarg)');
% xlabel('Time [s]'); ylabel('Amplitude');

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
subplot(311)
plot(time(idx), set_alt(idx)); grid on;
title('Set Altitude (chirp window)');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(312)
plot(time(idx), meas_alt(idx)); grid on;
title('Measured Altitude (chirp window)');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

% subplot(313)
% plot(time(idx), sinarg(idx)); grid on;
% title('Chirp Signal (chirp window)');
% xlabel('Time [s]'); ylabel('Amplitude');

subplot(313)
plot(time(idx), motor1(idx), time(idx), motor2(idx), time(idx), motor3(idx), time(idx), motor4(idx)); grid on;
legend('M1', 'M2', 'M3', 'M4', 'Location', 'best');
title('Motor Commands');
xlabel('Time [s]'); ylabel('Throttle');

%% Estimate Transfer Function

% Welch parameters
Nest     = round(10 / Ts_log);       % Window length [samples]
Noverlap = floor(0.9 * Nest);        % Overlap between windows
window   = hann(Nest, 'periodic');   % Hanning window

% Low-pass filter for rotating-frame filtering
Dlp = sqrt(3) / 2;    % Damping ratio
wlp = 2 * pi * 10;    % Cutoff frequency [rad/s]
Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');

% Apply rotating-frame filter to input and output
inp = apply_rotfiltfilt(Glp, sinarg_ax, set_alt);
out = apply_rotfiltfilt(Glp, sinarg_ax, meas_alt);

% Estimate frequency response
[T_ax, C_T_ax] = estimate_frequency_response(inp(idx), out(idx), window, Noverlap, Nest, Ts_log);


%% Plot Results
figure(3)
bode(T_ax); grid on;
title('Measured Altitude Hold Transfer Function');

