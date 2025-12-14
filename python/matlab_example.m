% =========================================================================
% Flight Data Analysis - Extract and Process Chirp Signal Data
% =========================================================================

clc, clear variables
addpath('../lib');

%% 1. DEFINE FILE PATH
file_path = '../20250908/20250908_flipmini_00.bbl.csv';

%% 2. READ THE ENTIRE FILE
fprintf('Reading CSV file...\n');
tic

% Find line where column headers start (contains "loopIteration")
fid = fopen(file_path, 'r');

line_count = 0;
while ~feof(fid)
    line = fgetl(fid);
    line_count = line_count + 1;

    if startsWith(line, '"loopIteration"')
        column_names = strsplit(line, ',');
        break;
    end
end

% Find column indices for analysis
col_time            = find(strcmp(column_names, '"time"'));
col_setpoint_roll   = find(strcmp(column_names, '"setpoint[0]"'));
col_gyro_unfilt_roll = find(strcmp(column_names, '"gyroUnfilt[0]"'));
col_debug_chirp     = find(strcmp(column_names, '"debug[0]"'));

fclose(fid);

% Read data matrix, using cache if available
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

% Detect chirp start and end (rising and falling edges)
chirp = data(:, col_debug_chirp);
chirp_binary = chirp > 0;
chirp_diff = diff([0; chirp_binary]);
chirp_start_idx = find(chirp_diff == 1);  % Rising edge
chirp_end_idx   = find(chirp_diff == -1); % Falling edge

chirp_start = chirp_start_idx(1);
chirp_end   = chirp_end_idx(1);

setpoint_roll    = data(:, col_setpoint_roll);
gyro_unfilt_roll = data(:, col_gyro_unfilt_roll);

% Python-style slicing: end index is EXCLUSIVE (use chirp_end-1)
time_over_chirp     = time(chirp_start:chirp_end-1);
setpoint_over_chirp = setpoint_roll(chirp_start:chirp_end-1);
gyro_over_chirp     = gyro_unfilt_roll(chirp_start:chirp_end-1);

fprintf('\nExtracted %d samples (%.2f seconds of data)\n', length(time), time(end));

%% 4. PROCESS DATA
% Define sampling periods for control loop and logging
Ts      = 125 * 1.e-6;   % Base sampling period (125 μs)
Ts_cntr = 2 * Ts;        % Control loop period
Ts_log  = 2 * Ts_cntr;   % Logging period

% Welch method parameters for frequency response estimation
Nest     = round(2 / Ts_log);       % FFT window size
koverlap = 0.9;
Noverlap = round(koverlap * Nest);  % Window overlap (round, not floor)
window   = hann(Nest, 'periodic');  % Matches Python's sym=False

% Estimate frequency response from input/output data
[G, C, ~, ~] = estimate_frequency_response(setpoint_over_chirp, gyro_over_chirp, window, Noverlap, Nest, Ts_log);

% Calculate step response from frequency response data
f_max_hz = 200;
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

figure(3); clf

% Configure Bode plot
opts = bodeoptions;
opts.FreqUnits      = 'Hz';
opts.MagUnits       = 'dB';
opts.PhaseUnits     = 'deg';
opts.PhaseWrapping  = 'off';  
opts.Grid           = 'on';

f_min_hz = max(0.5, time_over_chirp(2) - time_over_chirp(1));
f_max_hz = 200;

wmin = 2*pi*f_min_hz;
wmax = 2*pi*f_max_hz;

h = bodeplot(G, {wmin, wmax}, opts);
title("FREQUENCY RESPONSE");

% Plot step response
t_step = (0:numel(step_resp)-1) * Ts_log;

figure(4)
plot(t_step, step_resp);
title("Step Response")
xlabel("time (s)");
xlim([0 , 0.2])
grid on;

fprintf("\n=== Processing Complete ===\n");

% Save results for validation against Python implementation
save('./testing/output.mat', 't_step', 'step_resp')
