% =========================================================================
% Simple CSV Reader - Extract Setpoint and GyroUnfilt Data
% No external functions required
% =========================================================================

clc, clear variables
addpath ../

%% 1. DEFINE FILE PATH
file_path = '20250908/20250908_flipmini_00.bbl.csv';

%% 2. READ THE ENTIRE FILE
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
col_setpoint_pitch = find(strcmp(column_names, '"setpoint[1]"'));
col_setpoint_yaw = find(strcmp(column_names, '"setpoint[2]"'));
col_gyro_unfilt_roll = find(strcmp(column_names, '"gyroUnfilt[0]"'));
col_gyro_unfilt_pitch = find(strcmp(column_names, '"gyroUnfilt[1]"'));
col_gyro_unfilt_yaw = find(strcmp(column_names, '"gyroUnfilt[2]"'));
col_debug_chirp = find(strcmp(column_names, '"debug[0]"'));

fprintf('Found column headers at line %d\n', line_count);
fprintf('Column indices:\n');
fprintf('  time: %d\n', col_time);
fprintf('  debug[1] (chirp): %d\n', col_debug_chirp);
fprintf('  setpoint: [%d, %d, %d]\n', col_setpoint_roll, col_setpoint_pitch, col_setpoint_yaw);
fprintf('  gyroUnfilt: [%d, %d, %d]\n', col_gyro_unfilt_roll, col_gyro_unfilt_pitch, col_gyro_unfilt_yaw);

% Now read the actual data (skip the header lines)
fclose(fid);

% Read data matrix starting after header
try
  load([file_path(1, end - 8), '.mat'])
catch exeption
  data = readmatrix(file_path, 'NumHeaderLines', line_count);
  save([file_path(1, end - 8), '.mat'], 'data')
end

toc
fprintf('Data loaded: %d rows\n', size(data, 1));

%% 3. EXTRACT THE DATA

% Time vector (convert from microseconds to seconds)
time = (data(:, col_time) - data(1, col_time)) * 1.0e-6;

% extract chirp signal
chirp = data(:, col_debug_chirp);
chirp_active = chirp > 0;

% Detect chirp start (transition from 0 to non-zero)
chirp_diff = diff([0; chirp > 0]);  % prepend 0 to align indices
chirp_start_idx = find(chirp_diff == 1);  % rising edge
chirp_end_idx = find(chirp_diff == -1);   % falling edge

time_active = time(chirp_active);

% Setpoint data (deg/sec)
setpoint_roll = data(:, col_setpoint_roll);
setpoint_pitch = data(:, col_setpoint_pitch);
setpoint_yaw = data(:, col_setpoint_yaw);

% Gyro unfiltered data (deg/sec)
gyro_unfilt_roll = data(:, col_gyro_unfilt_roll);
gyro_unfilt_pitch = data(:, col_gyro_unfilt_pitch);
gyro_unfilt_yaw = data(:, col_gyro_unfilt_yaw);

u_roll = setpoint_roll(chirp_start_idx:chirp_end_idx);

y_roll = gyro_unfilt_roll(chirp_start_idx:chirp_end_idx);

fprintf('\nExtracted %d samples (%.2f seconds of data)\n', length(time), time(end));
fprintf('Chirp active for %d samples (%.2f seconds)\n', ...
        sum(chirp_active), sum(chirp_active) * (time(2) - time(1)));

%% 4. FFT

U_roll  = fft(u_roll);
U_pitch = fft(setpoint_pitch);
U_yaw   = fft(setpoint_yaw);

Y_roll  = fft(y_roll);
Y_pitch = fft(gyro_unfilt_pitch);
Y_yaw   = fft(gyro_unfilt_yaw);

%% 5. Y/U

G_roll  = Y_roll ./ U_roll;
G_pitch = Y_pitch ./ U_pitch;
G_yaw   = Y_yaw ./ U_yaw;

dt = time(2) - time(1);
fs = 1 / dt;

h_roll = real(ifft(G_roll));
h_pitch = real(ifft(G_pitch));
h_yaw = real(ifft(G_yaw));

Ts = 1 / fs;
sys_roll = tf(h_roll', 1, Ts);
sys_pitch = tf(h_pitch', 1, Ts);
sys_yaw = tf(h_yaw', 1, Ts);
%%% 7. BODE PLOTS
%
figure('Name', 'Bode Plot - Roll', 'NumberTitle', 'off');
bode(sys_roll);
title('Roll Transfer Function');
grid on;
%
%figure('Name', 'Bode Plot - Pitch', 'NumberTitle', 'off');
%bode(sys_pitch);
%title('Pitch Transfer Function');
%grid on;
%
%figure('Name', 'Bode Plot - Yaw', 'NumberTitle', 'off');
%bode(sys_yaw);
%title('Yaw Transfer Function');
%grid on;
%
%%% 8. STEP RESPONSE PLOTS
%
figure('Name', 'Step Response - Roll', 'NumberTitle', 'off');
step(sys_roll);
title('Roll Step Response');
grid on;
%
%figure('Name', 'Step Response - Pitch', 'NumberTitle', 'off');
%step(sys_pitch);
%title('Pitch Step Response');
%grid on;
%
%figure('Name', 'Step Response - Yaw', 'NumberTitle', 'off');
%step(sys_yaw);
%title('Yaw Step Response');
%grid on;
%
