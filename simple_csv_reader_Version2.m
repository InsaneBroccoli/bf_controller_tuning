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
fprintf('  debug[0] (chirp): %d\n', col_debug_chirp);
fprintf('  setpoint: [%d, %d, %d]\n', col_setpoint_roll, col_setpoint_pitch, col_setpoint_yaw);
fprintf('  gyroUnfilt: [%d, %d, %d]\n', col_gyro_unfilt_roll, col_gyro_unfilt_pitch, col_gyro_unfilt_yaw);

% Now read the actual data (skip the header lines)
fclose(fid);

% Read data matrix starting after header
try
    load([file_path(1:end-4), '.mat'])
catch exception
    data = readmatrix(file_path, 'NumHeaderLines', line_count);
    save([file_path(1:end-4), '.mat'], 'data')
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

% Handle edge case where chirp is still active at end of recording
if ~isempty(chirp_start_idx) && (isempty(chirp_end_idx) || chirp_end_idx(end) < chirp_start_idx(end))
    chirp_end_idx = [chirp_end_idx; length(time)];
end

fprintf('\nFound %d chirp events:\n', length(chirp_start_idx));
for i = 1:length(chirp_start_idx)
    fprintf('  Chirp %d: Start at %.3f s (idx %d), End at %.3f s (idx %d), Duration: %.3f s\n', ...
            i, time(chirp_start_idx(i)), chirp_start_idx(i), ...
            time(chirp_end_idx(i)), chirp_end_idx(i), ...
            time(chirp_end_idx(i)) - time(chirp_start_idx(i)));
end

% Setpoint data (deg/sec)
setpoint_roll = data(:, col_setpoint_roll);
setpoint_pitch = data(:, col_setpoint_pitch);
setpoint_yaw = data(:, col_setpoint_yaw);

% Gyro unfiltered data (deg/sec)
gyro_unfilt_roll = data(:, col_gyro_unfilt_roll);
gyro_unfilt_pitch = data(:, col_gyro_unfilt_pitch);
gyro_unfilt_yaw = data(:, col_gyro_unfilt_yaw);

fprintf('\nExtracted %d samples (%.2f seconds of data)\n', length(time), time(end));

%% 4. PROCESS EACH CHIRP EVENT

% Calculate sampling time
dt = median(diff(time));  % Use median for robustness
fs = 1 / dt;
fprintf('\nSampling rate: %.2f Hz (dt = %.6f s)\n', fs, dt);

% Process each chirp event
for chirp_num = 1:length(chirp_start_idx)
    fprintf('\n--- Processing Chirp Event %d ---\n', chirp_num);
    
    % Extract data for this chirp period
    idx_start = chirp_start_idx(chirp_num);
    idx_end = chirp_end_idx(chirp_num);
    
    % Extract time-aligned data
    time_chirp = time(idx_start:idx_end) - time(idx_start);  % Zero-referenced time
    
    u_roll = setpoint_roll(idx_start:idx_end);
    u_pitch = setpoint_pitch(idx_start:idx_end);
    u_yaw = setpoint_yaw(idx_start:idx_end);
    
    y_roll = gyro_unfilt_roll(idx_start:idx_end);
    y_pitch = gyro_unfilt_pitch(idx_start:idx_end);
    y_yaw = gyro_unfilt_yaw(idx_start:idx_end);
    
    % Remove DC offset
    u_roll = u_roll - mean(u_roll);
    u_pitch = u_pitch - mean(u_pitch);
    u_yaw = u_yaw - mean(u_yaw);
    
    y_roll = y_roll - mean(y_roll);
    y_pitch = y_pitch - mean(y_pitch);
    y_yaw = y_yaw - mean(y_yaw);
    
    fprintf('Chirp duration: %.3f s (%d samples)\n', time_chirp(end), length(time_chirp));
    
    %% 5. PLOT RAW INPUT/OUTPUT SIGNALS
    
    figure('Name', sprintf('Chirp %d - Raw Signals', chirp_num), 'NumberTitle', 'off');
    
    subplot(3,1,1)
    plot(time_chirp, u_roll, 'b', 'LineWidth', 1.5); hold on;
    plot(time_chirp, y_roll, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Rate (deg/s)');
    title(sprintf('Roll - Chirp Event %d', chirp_num));
    legend('Setpoint (Input)', 'Gyro Unfilt (Output)');
    grid on;
    
    subplot(3,1,2)
    plot(time_chirp, u_pitch, 'b', 'LineWidth', 1.5); hold on;
    plot(time_chirp, y_pitch, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Rate (deg/s)');
    title(sprintf('Pitch - Chirp Event %d', chirp_num));
    legend('Setpoint (Input)', 'Gyro Unfilt (Output)');
    grid on;
    
    subplot(3,1,3)
    plot(time_chirp, u_yaw, 'b', 'LineWidth', 1.5); hold on;
    plot(time_chirp, y_yaw, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Rate (deg/s)');
    title(sprintf('Yaw - Chirp Event %d', chirp_num));
    legend('Setpoint (Input)', 'Gyro Unfilt (Output)');
    grid on;
    
    %% 6. ESTIMATE TRANSFER FUNCTION USING TFESTIMATE
    
    % Window settings for spectral estimation
    window = hann(min(512, floor(length(u_roll)/4)));
    noverlap = floor(length(window)/2);
    nfft = max(1024, 2^nextpow2(length(u_roll)));
    
    % Estimate transfer functions
    [H_roll, f] = tfestimate(u_roll, y_roll, window, noverlap, nfft, fs);
    [H_pitch, ~] = tfestimate(u_pitch, y_pitch, window, noverlap, nfft, fs);
    [H_yaw, ~] = tfestimate(u_yaw, y_yaw, window, noverlap, nfft, fs);
    
    %% 7. PLOT FREQUENCY RESPONSE (BODE-STYLE)
    
    figure('Name', sprintf('Chirp %d - Frequency Response', chirp_num), 'NumberTitle', 'off');
    
    % Magnitude
    subplot(2,1,1)
    semilogx(f, 20*log10(abs(H_roll)), 'b', 'LineWidth', 1.5); hold on;
    semilogx(f, 20*log10(abs(H_pitch)), 'r', 'LineWidth', 1.5);
    semilogx(f, 20*log10(abs(H_yaw)), 'g', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title(sprintf('Transfer Function Magnitude - Chirp Event %d', chirp_num));
    legend('Roll', 'Pitch', 'Yaw');
    grid on;
    
    % Phase
    subplot(2,1,2)
    semilogx(f, unwrap(angle(H_roll))*180/pi, 'b', 'LineWidth', 1.5); hold on;
    semilogx(f, unwrap(angle(H_pitch))*180/pi, 'r', 'LineWidth', 1.5);
    semilogx(f, unwrap(angle(H_yaw))*180/pi, 'g', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
    title('Transfer Function Phase');
    legend('Roll', 'Pitch', 'Yaw');
    grid on;
    
    %% 8. COMPUTE IMPULSE RESPONSE AND STEP RESPONSE
    
    % Convert frequency response to impulse response
    h_roll = real(ifft([H_roll; conj(flipud(H_roll(2:end-1)))]));
    h_pitch = real(ifft([H_pitch; conj(flipud(H_pitch(2:end-1)))]));
    h_yaw = real(ifft([H_yaw; conj(flipud(H_yaw(2:end-1)))]));
    
    % Truncate to reasonable length (e.g., first half)
    n_impulse = min(500, floor(length(h_roll)/2));
    h_roll = h_roll(1:n_impulse);
    h_pitch = h_pitch(1:n_impulse);
    h_yaw = h_yaw(1:n_impulse);
    
    % Compute step response (cumulative sum of impulse response)
    step_roll = cumsum(h_roll) * dt;
    step_pitch = cumsum(h_pitch) * dt;
    step_yaw = cumsum(h_yaw) * dt;
    
    time_step = (0:n_impulse-1)' * dt;
    
    %% 9. PLOT STEP RESPONSES
    
    figure('Name', sprintf('Chirp %d - Step Response', chirp_num), 'NumberTitle', 'off');
    
    subplot(3,1,1)
    plot(time_step, step_roll, 'b', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Roll Step Response - Chirp Event %d', chirp_num));
    grid on;
    
    subplot(3,1,2)
    plot(time_step, step_pitch, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Pitch Step Response - Chirp Event %d', chirp_num));
    grid on;
    
    subplot(3,1,3)
    plot(time_step, step_yaw, 'g', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Yaw Step Response - Chirp Event %d', chirp_num));
    grid on;
    
    %% 10. PLOT IMPULSE RESPONSES
    
    figure('Name', sprintf('Chirp %d - Impulse Response', chirp_num), 'NumberTitle', 'off');
    
    subplot(3,1,1)
    plot(time_step, h_roll, 'b', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Roll Impulse Response - Chirp Event %d', chirp_num));
    grid on;
    
    subplot(3,1,2)
    plot(time_step, h_pitch, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Pitch Impulse Response - Chirp Event %d', chirp_num));
    grid on;
    
    subplot(3,1,3)
    plot(time_step, h_yaw, 'g', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Yaw Impulse Response - Chirp Event %d', chirp_num));
    grid on;
end

fprintf('\n=== Processing Complete ===\n');