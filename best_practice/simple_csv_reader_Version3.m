% =========================================================================
% Chirp Analysis for Flight Controller - Estimate FRF and Step Response
% Similar structure to serial_stream_pi_controller.py
% =========================================================================

clc, clear variables
addpath ../

%% 1. DEFINE FILE PATH AND LOAD DATA
file_path = '../20250908/20250908_flipmini_00.bbl.csv';

fprintf('Reading CSV file...\n');
tic

% Open file and find header
fid = fopen(file_path, 'r');

% Skip header lines - find where data starts
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

% Close and read data matrix
fclose(fid);

% Read data matrix starting after header (load from .mat if available)
mat_file_path = [file_path(1:end-4), '.mat'];
try
    load(mat_file_path, 'data')
    fprintf('Loaded data from .mat file\n');
catch
    data = readmatrix(file_path, 'NumHeaderLines', line_count);
    save(mat_file_path, 'data')
    fprintf('Saved data to .mat file\n');
end

toc
fprintf('Data loaded: %d rows\n', size(data, 1));

%% 2. EXTRACT AND ORGANIZE DATA

% Time vector (convert from microseconds to seconds)
time = (data(:, col_time) - data(1, col_time)) * 1.0e-6;

% Extract chirp signal
chirp = data(:, col_debug_chirp);

% Detect chirp start and end (transition from 0 to non-zero)
chirp_binary = chirp > 0;
chirp_diff = diff([0; chirp_binary]);
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

%% 3. EVALUATE TIME - VERIFY SAMPLING RATE

Ts = median(diff(time));  % Use median for robustness
fs = 1 / Ts;

figure('Name', 'Sampling Time Analysis', 'NumberTitle', 'off');
plot(time(1:end-1), diff(time) * 1e6, 'b');
grid on;
title(sprintf('Mean %.0f μs, Std. %.0f μs, Med. dT = %.0f μs', ...
              mean(diff(time) * 1e6), std(diff(time) * 1e6), median(diff(time) * 1e6)));
xlabel('Time (sec)');
ylabel('dTime (μs)');
xlim([0, time(end-1)]);
ylim([0, 1.2 * max(diff(time) * 1e6)]);

fprintf('\nSampling rate: %.2f Hz (Ts = %.6f s)\n', fs, Ts);

%% 4. FUNCTION: ESTIMATE FRF AND COHERENCE

function [freq, g, c] = estimate_frf_and_coherence(x, y, fs, window, nperseg, noverlap)
    % Estimate cross spectral density and power spectral densities
    [Pxy, freq] = cpsd(x, y, window, noverlap, nperseg, fs);
    Pxx = pwelch(x, window, noverlap, nperseg, fs);
    Pyy = pwelch(y, window, noverlap, nperseg, fs);
    
    % Calculate frequency response function
    g = Pxy ./ Pxx;
    
    % Calculate coherence
    c = abs(Pxy).^2 ./ (Pxx .* Pyy);
    
    % Truncate DC (freq=0) to avoid divide-by-zero issues
    freq = freq(2:end);
    g = g(2:end);
    c = c(2:end);
end

%% 5. FUNCTION: GET STEP RESPONSE FROM FRD

function step_resp = get_step_resp_from_frd(g, freq, f_max_hz, Ts)
    % Zero out above f_max_hz
    g(freq > f_max_hz) = 0;
    
    % Reconstruct DC (simulate symmetry at zero freq)
    g_dc = g(1);
    g = [g_dc; g];
    freq = [0; freq];
    
    % Construct full symmetric spectrum
    g_full = [g; conj(flipud(g(2:end-1)))];
    
    % Step response is cumulative sum of real part of IFFT
    step_resp = cumsum(real(ifft(g_full))) * Ts;
end

%% 6. PROCESS EACH CHIRP EVENT

% Select which chirp to analyze (change this to process different chirps)
chirp_num = 1;

if chirp_num > length(chirp_start_idx)
    error('Chirp number %d not found. Only %d chirps available.', chirp_num, length(chirp_start_idx));
end

fprintf('\n=== Processing Chirp Event %d ===\n', chirp_num);

% Extract data for this chirp period
idx_start = chirp_start_idx(chirp_num);
idx_end = chirp_end_idx(chirp_num);

% Extract time-aligned data
time_chirp = time(idx_start:idx_end) - time(idx_start);

% Input and output signals
inp_roll = setpoint_roll(idx_start:idx_end);
out_roll = gyro_unfilt_roll(idx_start:idx_end);

inp_pitch = setpoint_pitch(idx_start:idx_end);
out_pitch = gyro_unfilt_pitch(idx_start:idx_end);

inp_yaw = setpoint_yaw(idx_start:idx_end);
out_yaw = gyro_unfilt_yaw(idx_start:idx_end);

fprintf('Chirp duration: %.3f s (%d samples)\n', time_chirp(end), length(time_chirp));

%% 7. PLOT RAW SIGNALS (Similar to Python Figure 2 & 3)

figure('Name', sprintf('Chirp %d - Input and Output Signals', chirp_num), 'NumberTitle', 'off');

subplot(3, 1, 1)
plot(time_chirp, inp_roll, 'b', 'LineWidth', 1.5); hold on;
plot(time_chirp, out_roll, 'r', 'LineWidth', 1.5);
grid on;
ylabel('Rate (deg/s)');
title(sprintf('Roll - Chirp Event %d', chirp_num));
legend('Setpoint (Input)', 'Gyro Unfilt (Output)', 'Location', 'best');

subplot(3, 1, 2)
plot(time_chirp, inp_pitch, 'b', 'LineWidth', 1.5); hold on;
plot(time_chirp, out_pitch, 'r', 'LineWidth', 1.5);
grid on;
ylabel('Rate (deg/s)');
title('Pitch');
legend('Setpoint (Input)', 'Gyro Unfilt (Output)', 'Location', 'best');

subplot(3, 1, 3)
plot(time_chirp, inp_yaw, 'b', 'LineWidth', 1.5); hold on;
plot(time_chirp, out_yaw, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time (sec)');
ylabel('Rate (deg/s)');
title('Yaw');
legend('Setpoint (Input)', 'Gyro Unfilt (Output)', 'Location', 'best');

%% 8. FREQUENCY RESPONSE ESTIMATION

% Window settings for spectral estimation
Nest = round(2.0 / Ts);
koverlap = 0.9;
Noverlap = round(koverlap * Nest);
window = hann(Nest);

% Estimate FRF and coherence for each axis
[freq_roll, g_roll, c_roll] = estimate_frf_and_coherence(...
    diff(inp_roll), diff(out_roll), fs, window, Nest, Noverlap);

[freq_pitch, g_pitch, c_pitch] = estimate_frf_and_coherence(...
    diff(inp_pitch), diff(out_pitch), fs, window, Nest, Noverlap);

[freq_yaw, g_yaw, c_yaw] = estimate_frf_and_coherence(...
    diff(inp_yaw), diff(out_yaw), fs, window, Nest, Noverlap);

% Print DC-Gain
fprintf('\nMeasured DC-Gains:\n');
fprintf('  Roll:  %.2f dB\n', 20*log10(abs(g_roll(1))));
fprintf('  Pitch: %.2f dB\n', 20*log10(abs(g_pitch(1))));
fprintf('  Yaw:   %.2f dB\n', 20*log10(abs(g_yaw(1))));

%% 9. BODE PLOT (Similar to Python Figure 4)

figure('Name', sprintf('Chirp %d - Frequency Response (Bode)', chirp_num), 'NumberTitle', 'off');

% Magnitude
subplot(2, 1, 1)
semilogx(freq_roll, 20*log10(abs(g_roll)), 'b', 'LineWidth', 1.5); hold on;
semilogx(freq_pitch, 20*log10(abs(g_pitch)), 'r', 'LineWidth', 1.5);
semilogx(freq_yaw, 20*log10(abs(g_yaw)), 'g', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title(sprintf('Transfer Function Magnitude - Chirp Event %d', chirp_num));
legend('Roll', 'Pitch', 'Yaw', 'Location', 'best');

% Phase
subplot(2, 1, 2)
semilogx(freq_roll, unwrap(angle(g_roll))*180/pi, 'b', 'LineWidth', 1.5); hold on;
semilogx(freq_pitch, unwrap(angle(g_pitch))*180/pi, 'r', 'LineWidth', 1.5);
semilogx(freq_yaw, unwrap(angle(g_yaw))*180/pi, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
title('Transfer Function Phase');
legend('Roll', 'Pitch', 'Yaw', 'Location', 'best');

%% 10. COHERENCE PLOT (Similar to Python Figure 5)

figure('Name', sprintf('Chirp %d - Coherence', chirp_num), 'NumberTitle', 'off');
semilogx(freq_roll, c_roll, 'b', 'LineWidth', 1.5); hold on;
semilogx(freq_pitch, c_pitch, 'r', 'LineWidth', 1.5);
semilogx(freq_yaw, c_yaw, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Coherence (abs)');
title(sprintf('Coherence - Chirp Event %d', chirp_num));
legend('Roll', 'Pitch', 'Yaw', 'Location', 'best');
ylim([0, 1]);

%% 11. STEP RESPONSE (Similar to Python Figure 6)

% Maximum frequency for step response calculation
f_max = 80;  % Hz

% Calculate step responses
step_resp_roll = get_step_resp_from_frd(g_roll, freq_roll, f_max, Ts);
step_resp_pitch = get_step_resp_from_frd(g_pitch, freq_pitch, f_max, Ts);
step_resp_yaw = get_step_resp_from_frd(g_yaw, freq_yaw, f_max, Ts);

% Time vector for step response
step_time = (0:length(step_resp_roll)-1)' * Ts;

% Plot step responses
figure('Name', sprintf('Chirp %d - Step Response', chirp_num), 'NumberTitle', 'off');

subplot(3, 1, 1)
plot(step_time, step_resp_roll, 'b', 'LineWidth', 1.5);
grid on;
ylabel('Amplitude');
title(sprintf('Roll Step Response - Chirp Event %d', chirp_num));
xlim([0, 0.2]);

subplot(3, 1, 2)
plot(step_time, step_resp_pitch, 'r', 'LineWidth', 1.5);
grid on;
ylabel('Amplitude');
title('Pitch Step Response');
xlim([0, 0.2]);

subplot(3, 1, 3)
plot(step_time, step_resp_yaw, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Time (sec)');
ylabel('Amplitude');
title('Yaw Step Response');
xlim([0, 0.2]);

fprintf('\n=== Processing Complete ===\n');