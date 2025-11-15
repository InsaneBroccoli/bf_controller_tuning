% =========================================================================
% Simple CSV Reader - Extract Setpoint and GyroUnfilt Data
% No external functions required
% =========================================================================

clc, clear variables
addpath('./lib');
%% 1. DEFINE FILE PATH
file_path = './20250908/20250908_flipmini_00.bbl.csv'

%% 2. READ THE ENTIRE FILE:
fprintf('Reading CSV file...\n');
tic

% Open file
fid = fopen(file_path, 'r');

% Skip header lines - find where data starts
% Header ends where we see a line starting with a number (loopIteration)
line_count = 0;
while ~feof(fid)
    line = fgetl(fid);
    line_count = line_count + 1;
    
    % Check if line starts with "loopIteration" - this is the column header
    if startsWith(line, '"loopIteration"')
        % Parse column names
        column_names = strsplit(line, ',');
        break;
    end
end

% Find the column indices we need
col_time = find(strcmp(column_names, '"time"'));
col_setpoint_roll = find(strcmp(column_names, '"setpoint[0]"'));
col_gyro_unfilt_roll = find(strcmp(column_names, '"gyroUnfilt[0]"'));
col_debug_chirp = find(strcmp(column_names, '"debug[0]"'));

fclose(fid);

% Read data matrix starting after header
try
    load([file_path(1:end-8), '.mat'])
catch exception
    data = readmatrix(file_path, 'NumHeaderLines', line_count);
    save([file_path(1:end-8), '.mat'], 'data')
end

toc
fprintf('Data loaded: %d rows\n', size(data, 1));

%% 3. EXTRACT THE DATA

% Time vector (convert from microseconds to seconds)
time = (data(:, col_time) - data(1, col_time)) * 1.0e-6;

% Extract chirp signal
chirp = data(:, col_debug_chirp);

% Detect chirp start and end (transition from 0 to non-zero)
chirp_binary = chirp > 0;
chirp_diff = diff([0; chirp_binary]);  % prepend 0 to align indices
chirp_start_idx = find(chirp_diff == 1);  % rising edge
chirp_end_idx = find(chirp_diff == -1);   % falling edge

% Get first chirp signal
chirp_start = chirp_start_idx(1);
chirp_end = chirp_end_idx(1);

setpoint_roll = data(:, col_setpoint_roll);
gyro_unfilt_roll = data(:, col_gyro_unfilt_roll);

% Get data over chirp
time_over_chirp = time(chirp_start:chirp_end);
setpoint_over_chirp = setpoint_roll(chirp_start:chirp_end);
gyro_over_chirp = gyro_unfilt_roll(chirp_start:chirp_end);

fprintf('\nExtracted %d samples (%.2f seconds of data)\n', length(time), time(end));

%% 4. PROCESS DATA
% Calculate sampling time
Ts = 125 * 1.e-6;
Ts_cntr = 2 * Ts;
Ts_log = 2 * Ts_cntr;

% Parameters
Nest     = round(2 / Ts_log);
koverlap = 0.9;
Noverlap = floor(koverlap * Nest);
window   = hann(Nest, 'periodic');

[G, C, ~, ~] = estimate_frequency_response(setpoint_over_chirp, gyro_over_chirp, window, Noverlap, Nest, Ts_log);

f_max_hz = 200
step_resp = calculate_step_response_from_frd(G, f_max_hz);

%% 5. PLOT DATA 
figure(1)
subplot(211)
plot(time, setpoint_roll)
title("setpoint")
subplot(212)
plot(time, gyro_unfilt_roll)
title("gyro")

figure(2)
subplot(211)
plot(time_over_chirp, setpoint_over_chirp)
title("setpoint")
subplot(212)
plot(time_over_chirp, gyro_over_chirp)
title("gyro")

figure(3)
plot(G, C);
xlabel("Frequency (Hz)");
title("Frequency Response");
xlim([0 , f_max_hz])
grid on;

figure(4)
plot(step_resp);
title("Step Response")
xlabel("time (ms)");
xlim([0 , f_max_hz])
grid on;
fprintf("\n=== Processing Complete ===\n");
