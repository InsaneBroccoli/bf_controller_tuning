% =========================================================================
% Position Hold Frequency Response Analysis
% Axis-aware version: separates Roll (LON) and Pitch (LAT) excitation
%
% axis_select controls which axes are shown in all plots:
%   'roll'  -> only Roll (LON)
%   'pitch' -> only Pitch (LAT)
%   'both'  -> Roll and Pitch overlaid
% =========================================================================

clc, clear, close all
addpath("../lib/");
addpath(genpath("../logs/"));

pos_bode = [0.1514, 0.5838-0.2, 0.7536, 0.3472+0.2; ...
            0.1514, 0.1100    , 0.7536, 0.1917    ];

linewidth = 1.5;

% -------------------------------------------------------------------------
% USER SETTINGS
% -------------------------------------------------------------------------
axis_select = 'both';   % 'roll' | 'pitch' | 'both'

% Colors: index 1 = Roll, index 2 = Pitch
colors  = {'k', 'b'};
markers = {'-', '--'};

% -------------------------------------------------------------------------
% File information
% -------------------------------------------------------------------------
log_folder  = '../logs';
flight_folder = '20260512';
log_name    = '20260512_6_inch_2.csv';
file_path   = fullfile(log_folder, flight_folder, log_name);

% =========================================================================
% Load and Process Flight Log Data
% =========================================================================
[para, Nheader, ind] = extract_header_information(file_path);

[folder, base, ~] = fileparts(file_path);
mat_path = fullfile(folder, [base '.mat']);

try
    S    = load(mat_path);
    data = S.data;
catch
    data = readmatrix(file_path, "NumHeaderLines", Nheader);
    save(mat_path, "data");
end

% Decimate: blackbox runs at FC loop rate; althold debug updates at ~100 Hz
Ts_log_raw = para.looptime * 1e-6 * para.pid_process_denom * para.frameIntervalPDenom;
dec        = round((1/Ts_log_raw) / 100);
data       = data(1:dec:end, :);
Ts_log     = Ts_log_raw * dec;
Ts_cntr    = Ts_log;

% Time vector (microseconds -> seconds)
time = (data(:, ind.time) - data(1, ind.time)) * 1e-6;

% -------------------------------------------------------------------------
% Undo firmware-side debug scaling
% -------------------------------------------------------------------------
gps_error       = data(:,ind.debug(1)) / 10;    % cm
angle_target    = data(:,ind.debug(2)) / 10;    % deg
chirp           = data(:,ind.debug(3)) / 10;    % cm (injected setpoint)
chirp_inst_freq = data(:,ind.debug(4)) / 100;   % Hz
pid_sum_EF      = data(:,ind.debug(5)) / 10;    % PID sum Earth Frame
sinarg          = data(:,ind.debug(6)) / 5e3;   % injected chirp signal
active_axis     = data(:,ind.debug(7)) * 2;     % 0 = Roll/LON, >0 = Pitch/LAT
pidDA_limit     = data(:,ind.debug(8)) / 10;    % deg

target_position  = chirp;
current_position = chirp - gps_error;

% -------------------------------------------------------------------------
% DA-limit / chirp frequency diagnostic plot
% -------------------------------------------------------------------------
figure(10)
plot(time, pidDA_limit); hold on
plot(time, chirp_inst_freq*10); grid on
legend('Combined D+A vector length', 'Instantaneous chirp frequency (scaled)', ...
    Location='northwest');
title('DA Limit Diagnostic');
xlabel('Time [s]');

% =========================================================================
% Overview Plots
% =========================================================================
figure(1)
subplot(211)
plot(time, gps_error, '-b'); grid on;
title('GPS Error');
xlabel('Time [s]'); ylabel('Error [cm]');

subplot(212)
plot(time, pid_sum_EF, 'r'); hold on;
plot(time, angle_target, '-b'); grid on;
title('PID Sum EF & Angle Target');
xlabel('Time [s]'); ylabel('Angle [deg]');

figure(2)
subplot(211)
plot(time, sinarg, '-m'); hold on;
plot(time, active_axis, '-k'); grid on;
legend('sinarg', 'axis (low=roll, high=pitch)');
title('Chirp Excitation');
ylabel('Excitation Signal');
ylim([-1 7]);

subplot(212)
plot(time, chirp, '-m'); hold on;
plot(time, chirp_inst_freq * 50, '-k'); grid on;
legend('chirp exc', 'inst. freq (scaled)');
xlabel('Time [s]'); ylabel('cm / scaled Hz');
ylim([-210 210]);

% =========================================================================
% Chirp window detection
% =========================================================================
idx = get_ind_eval(sinarg, chirp);

% -------------------------------------------------------------------------
% Axis index masks
%   active_axis < 1  -> Roll  (LON)
%   active_axis >= 1 -> Pitch (LAT)
% -------------------------------------------------------------------------
idx_roll  = idx & (active_axis < 1);
idx_pitch = idx & (active_axis >= 1);

fprintf('Chirp samples  total: %d | roll: %d | pitch: %d\n', ...
    sum(idx), sum(idx_roll), sum(idx_pitch));

% =========================================================================
% Welch / Rotfilter parameters
% =========================================================================
frame    = 20;
Nest     = round(frame / Ts_log);
Noverlap = floor(0.9 * Nest);
window   = hann(Nest, 'periodic');

Dlp = sqrt(3) / 2;
wlp = 2 * pi * 1;
Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');

% sinarg zeroed outside chirp window
sinarg_ax = sinarg;
sinarg_ax(~idx) = 0;

% =========================================================================
% Rotating-frame filtering -- separate per axis
% =========================================================================

% --- Roll ---
sinarg_roll = sinarg_ax .* double(active_axis < 1);
inp_r   = apply_rotfiltfilt(Glp, sinarg_roll, target_position);
out_y_r = apply_rotfiltfilt(Glp, sinarg_roll, current_position);
out_u_r = apply_rotfiltfilt(Glp, sinarg_roll, pid_sum_EF);
out_v_r = apply_rotfiltfilt(Glp, sinarg_roll, angle_target);

% --- Pitch ---
sinarg_pitch = sinarg_ax .* double(active_axis >= 1);
inp_p   = apply_rotfiltfilt(Glp, sinarg_pitch, target_position);
out_y_p = apply_rotfiltfilt(Glp, sinarg_pitch, current_position);
out_u_p = apply_rotfiltfilt(Glp, sinarg_pitch, pid_sum_EF);
out_v_p = apply_rotfiltfilt(Glp, sinarg_pitch, angle_target);

% =========================================================================
% Frequency response estimation -- per axis
% =========================================================================

% --- Roll ---
[T_r,   C_T_r]   = estimate_frequency_response(inp_r(idx_roll),  out_y_r(idx_roll), ...
    window, Noverlap, Nest, Ts_log);
[Guw_r, C_Guw_r] = estimate_frequency_response(inp_r(idx_roll),  out_u_r(idx_roll), ...
    window, Noverlap, Nest, Ts_log);

% --- Pitch ---
[T_p,   C_T_p]   = estimate_frequency_response(inp_p(idx_pitch), out_y_p(idx_pitch), ...
    window, Noverlap, Nest, Ts_log);
[Guw_p, C_Guw_p] = estimate_frequency_response(inp_p(idx_pitch), out_u_p(idx_pitch), ...
    window, Noverlap, Nest, Ts_log);

% Derived quantities
P_r = T_r / Guw_r;
P_p = T_p / Guw_p;

Cpid_r = Guw_r / (1 - T_r);
Cpid_p = Guw_p / (1 - T_p);

% Coherence data
C_T_r_data   = squeeze(C_T_r.ResponseData);
C_T_p_data   = squeeze(C_T_p.ResponseData);
C_Guw_r_data = squeeze(C_Guw_r.ResponseData);
C_Guw_p_data = squeeze(C_Guw_p.ResponseData);
C_P_r        = C_Guw_r_data .* C_T_r_data;
C_P_p        = C_Guw_p_data .* C_T_p_data;

% Shared frequency axis (from Roll; both axes use same Welch settings)
f_bode     = squeeze(C_T_r.Frequency);
omega_bode = 2*pi*f_bode;

% =========================================================================
% Axis selection helper
% =========================================================================
[Tsel, Psel, Cpid_sel, Coh_T, Coh_P, ax_labels] = ...
    select_axes(axis_select, ...
        T_r,   T_p,   P_r,    P_p,   Cpid_r, Cpid_p, ...
        C_T_r_data, C_T_p_data, C_P_r, C_P_p);

% =========================================================================
% Figure 3 -- Transfer Function T
% =========================================================================
opt_bode = make_bodeopts();

figure(3)
ax_f3(1) = subplot('Position', pos_bode(1,:));
for i = 1:numel(Tsel)
    bode(ax_f3(1), Tsel{i}, [colors{i} markers{i}], omega_bode, opt_bode);
    hold(ax_f3(1), 'on');
end
legend(ax_labels); grid(ax_f3(1), 'on');
title(ax_f3(1), 'Bode Plot Transfer Function T');

ax_f3(2) = subplot('Position', pos_bode(2,:));
for i = 1:numel(Tsel)
    semilogx(ax_f3(2), f_bode, Coh_T{i}, [colors{i} markers{i}], 'LineWidth', linewidth);
    hold(ax_f3(2), 'on');
end
legend(ax_labels); grid(ax_f3(2), 'on');
ylabel(ax_f3(2), 'Coherence'); xlabel(ax_f3(2), 'Frequency [Hz]');
ylim(ax_f3(2), [0 1]); xlim(ax_f3(2), [min(f_bode) max(f_bode)]);
linkaxes(ax_f3, 'x');

% =========================================================================
% Figure 4 -- PID Controller (measured vs analytical)
% =========================================================================

% Analytical PID (current tune)
P_cur      = 30;
I_cur      = 30;
D_cur      = 30;
A_cur      = 30;
fc_pt1_cur = 0.8;

freq_vector = T_r.Frequency;

Cpid_ana = calculate_poshold_controller( ...
    P_cur, I_cur, D_cur, A_cur, fc_pt1_cur, Ts_cntr, Ts_log, freq_vector);

figure(4)
for i = 1:numel(Cpid_sel)
    bode(Cpid_sel{i}, [colors{i} markers{i}], omega_bode, opt_bode); hold on;
end
bode(Cpid_ana, 'r-', omega_bode, opt_bode);
legend([ax_labels, {'Analytical'}]);
title('Bode Plot PID Controller');

% =========================================================================
% Figure 5 -- Plant P
% =========================================================================
figure(5)
ax_f5(1) = subplot('Position', pos_bode(1,:));
for i = 1:numel(Psel)
    bode(ax_f5(1), Psel{i}, [colors{i} markers{i}], omega_bode, opt_bode);
    hold(ax_f5(1), 'on');
end
legend(ax_labels); grid(ax_f5(1), 'on');
title(ax_f5(1), 'Bode Plot Plant P');

ax_f5(2) = subplot('Position', pos_bode(2,:));
for i = 1:numel(Psel)
    semilogx(ax_f5(2), f_bode, Coh_P{i}, [colors{i} markers{i}], 'LineWidth', linewidth);
    hold(ax_f5(2), 'on');
end
legend(ax_labels); grid(ax_f5(2), 'on');
ylabel(ax_f5(2), 'Coherence'); xlabel(ax_f5(2), 'Frequency [Hz]');
ylim(ax_f5(2), [0 1]); xlim(ax_f5(2), [min(f_bode) max(f_bode)]);
linkaxes(ax_f5, 'x');

% =========================================================================
% Closed-Loop (uses Roll plant as representative; adapt if needed)
% =========================================================================
CL_ana = calculate_closed_loop(Cpid_ana, tf(1,1,Ts_log), P_r, tf(1,1,Ts_log), tf(0,1));

% =========================================================================
% New Controller (proposed tune)
% =========================================================================
default_parameters = true;

if default_parameters
    P_new = P_cur;  I_new = I_cur;  D_new = D_cur;
    A_new = A_cur;  fc_pt1_new = fc_pt1_cur;
else
    P_new = 20;  
    I_new = 20;  
    D_new = 20;  
    A_new = 20;  
    fc_pt1_new = 1;
end

Cpid_ana_new = calculate_poshold_controller( ...
    P_new, I_new, D_new, A_new, fc_pt1_new, Ts_cntr, Ts_log, P_r.Frequency);

CL_ana_new = calculate_closed_loop(Cpid_ana_new, tf(1,1,Ts_log), P_r, tf(1,1,Ts_log), tf(0,1));

% =========================================================================
% Figure 6 -- Measured T vs Analytical CL
% =========================================================================
figure(6)
for i = 1:numel(Tsel)
    bode(Tsel{i}, [colors{i} markers{i}], omega_bode, opt_bode); hold on;
end
bode(CL_ana.T, 'r-', omega_bode, opt_bode);
xlim([1e-2 10]); grid on;
title('Comparison Measured T to Analytical Closed Loop');
legend([ax_labels, {'CL Analytical'}]);

% =========================================================================
% Figure 7 -- Gang of Four
% =========================================================================
figure(7)

% --- Tracking T ---
ax_gof(1) = subplot(2,2,1);
bodemag(ax_gof(1), CL_ana.T,     omega_bode, opt_bode); hold(ax_gof(1), 'on');
bodemag(ax_gof(1), CL_ana_new.T, omega_bode, opt_bode);
for i = 1:numel(Tsel)
    bodemag(ax_gof(1), Tsel{i}, [colors{i} markers{i}], omega_bode, opt_bode);
end
title(ax_gof(1), 'Tracking T');
legend(['Actual', 'New', ax_labels], 'Location', 'best');
grid(ax_gof(1), 'on');

% --- Sensitivity S ---
ax_gof(2) = subplot(2,2,2);
bodemag(ax_gof(2), CL_ana.S,     omega_bode, opt_bode); hold(ax_gof(2), 'on');
bodemag(ax_gof(2), CL_ana_new.S, omega_bode, opt_bode);
title(ax_gof(2), 'Sensitivity S');
legend({'Actual', 'New'}, 'Location', 'northwest');
grid(ax_gof(2), 'on');

% --- Controller Effort SC ---
ax_gof(3) = subplot(2,2,3);
bodemag(ax_gof(3), CL_ana.SC,     omega_bode, opt_bode); hold(ax_gof(3), 'on');
bodemag(ax_gof(3), CL_ana_new.SC, omega_bode, opt_bode);
title(ax_gof(3), 'Controller Effort SC');
legend({'Actual', 'New'}, 'Location', 'northwest');
grid(ax_gof(3), 'on');

% --- Compliance SP ---
ax_gof(4) = subplot(2,2,4);
bodemag(ax_gof(4), CL_ana.SP,     omega_bode, opt_bode); hold(ax_gof(4), 'on');
bodemag(ax_gof(4), CL_ana_new.SP, omega_bode, opt_bode);
title(ax_gof(4), 'Compliance SP');
legend({'Actual', 'New'}, 'Location', 'southwest');
grid(ax_gof(4), 'on');

linkaxes(ax_gof, 'x');
xlim(ax_gof(1), [1e-2 10]);
sgtitle('Gang of Four - Position Hold');

% =========================================================================
% Step Response placeholder
% =========================================================================


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [Tsel, Psel, Cpid_sel, Coh_T, Coh_P, labels] = select_axes( ...
        axis_select, ...
        T_r, T_p, P_r, P_p, Cpid_r, Cpid_p, ...
        C_T_r, C_T_p, C_P_r, C_P_p)
% SELECT_AXES  Return FRD cell arrays for the requested axis selection.
%
%   axis_select : 'roll' | 'pitch' | 'both'
%
%   Returns cell arrays with one entry per selected axis, plus matching
%   string labels for legends.

    switch lower(axis_select)
        case 'roll'
            Tsel     = {T_r};
            Psel     = {P_r};
            Cpid_sel = {Cpid_r};
            Coh_T    = {C_T_r};
            Coh_P    = {C_P_r};
            labels   = {'Roll'};

        case 'pitch'
            Tsel     = {T_p};
            Psel     = {P_p};
            Cpid_sel = {Cpid_p};
            Coh_T    = {C_T_p};
            Coh_P    = {C_P_p};
            labels   = {'Pitch'};

        case 'both'
            Tsel     = {T_r,   T_p};
            Psel     = {P_r,   P_p};
            Cpid_sel = {Cpid_r, Cpid_p};
            Coh_T    = {C_T_r, C_T_p};
            Coh_P    = {C_P_r, C_P_p};
            labels   = {'Roll', 'Pitch'};

        otherwise
            error('select_axes: axis_select must be ''roll'', ''pitch'', or ''both''.');
    end
end


function opt = make_bodeopts()
% MAKE_BODEOPTS  Return a standard bodeoptions struct used throughout.
    opt              = bodeoptions('cstprefs');
    opt.MagUnits     = 'dB';
    opt.MagScale     = 'linear';
    opt.PhaseUnits   = 'deg';
    opt.FreqUnits    = 'Hz';
    opt.PhaseWrapping = 'on';
    opt.Grid         = 'on';
end