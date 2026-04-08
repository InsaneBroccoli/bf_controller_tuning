% Initialize workspace
clc, clear variables, close all
addpath("../lib/");
addpath(genpath("../logs/"));

% Add file information
log_folder = '../logs';
flight_folder = '20260402';
log_name = 'LOG093.TXT.csv';
file_path = fullfile(log_folder, flight_folder, log_name);

% --- Load and Process Flight Log Data ---
[para, Nheader, ind, ind_cntr] = extract_header_information(file_path);

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

time = (data(:, ind.time) - data(1, ind.time)) * 1e6;

sinarg          = data(:,ind.debug(1)) / 5e3;
meas_alt        = data(:,ind.debug(4));
set_alt         = data(:,ind.debug(5));
throttle_out    = data(:,ind.debug(6));
vertical_v      = data(:,ind.debug(7));
throttle_offset = data(:,ind.debug(8));

time            = time(1:20:end);
sinarg          = sinarg(1:20:end);
meas_alt        = meas_alt(1:20:end);
set_alt         = set_alt(1:20:end);
throttle_out    = throttle_out(1:20:end);
vertical_v      = vertical_v(1:20:end);
throttle_offset = throttle_offset(1:20:end);

figure(1)
subplot(411)
plot(time, meas_alt, '-r'); hold on
plot(time, set_alt, '-b'); grid on;
legend('Target Altitude', 'Measured Altitude', 'Location', 'best');
title('Compare Target and Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(412)
plot(time, throttle_out); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle');

subplot(413)
plot(time, throttle_out); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle');

subplot(414)
plot(time, vertical_v); grid on;
title('Vertical Speed');
xlabel('Time [s]'); ylabel('Speed [cm/s]');


idx = get_ind_eval(sinarg, meas_alt);
meas_alt = fix_offset(meas_alt, idx);

time_c            = time(idx);
meas_alt_c        = meas_alt(idx);
set_alt_c         = set_alt(idx);
throttle_out_c    = throttle_out(idx);
vertical_v_c      = vertical_v(idx);
throttle_offset_c = throttle_offset(idx);


figure(2)
subplot(411)
plot(time_c, meas_alt_c, '-r'); hold on
plot(time_c, set_alt_c, '-b'); grid on;
legend('Target Altitude', 'Measured Altitude', 'Location', 'best');
title('Compare Target and Measured Altitude');
xlabel('Time [s]'); ylabel('Altitude offset [cm]');

subplot(412)
plot(time_c, throttle_out_c); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle');

subplot(413)
plot(time_c, throttle_out_c); grid on;
title('Throttle Out');
xlabel('Time [s]'); ylabel('Throttle');

subplot(414)
plot(time_c, vertical_v_c); grid on;
title('Vertical Speed');
xlabel('Time [s]'); ylabel('Speed [cm/s]');

Ts_log = 1 / 100;
