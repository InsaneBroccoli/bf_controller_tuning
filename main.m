clc, clear variables
addpath lib_michi/
addpath(genpath('lib')); % fügt lib/ und auch alle unterordner hinzu
addpath logs/
addpath utils/
%%

config; % load axis_r, colors, lineWidth etc.

% -------------------------------------------------------------------------

% Extract header information
[para, Nheader, ind, ind_cntr] = extract_header_information(file_path);

%=data_io==================================================================

% Read the data
%  - If its the first time from the .csv and save a mat, otherwise the
%    .mat. This increases load speed significantly.
tic
try
   load([file_path(1:end-8), '.mat'])
catch exception
   data = readmatrix(file_path, 'NumHeaderLines', Nheader);
   save([file_path(1:end-8), '.mat'], 'data');
end
[Ndata, Nsig] = size(data)
toc

% Expand index
ind.axisSumPI = ind_cntr + (1:3);
ind.sinarg = ind.debug(1);

% Convert and evaluate time
time = (data(:,ind.time) - data(1,ind.time)) * 1.0e-6;
delta_time_mus = diff(time) * 1.0e6;

figure(99)
plot(time(1:end-1), delta_time_mus), grid on
title(sprintf('Mean: %0.2f mus, Median: %0.2f mus, Std: %0.2f mus\n', ...
      mean(delta_time_mus), ...
      median(delta_time_mus), ...
      std(delta_time_mus)))
xlabel('Time (sec)'), ylabel('Ts log (mus)')
xlim([0, time(end)])
set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)

% Unscale highResolutionGain
if para.blackbox_high_resolution
    blackbox_high_resolution_scale = 10.0;
    ind_bb_high_res = [ind.gyroADC, ind.gyroUnfilt, ind.rcCommand, ind.setpoint(1:3)];
    data(:, ind_bb_high_res) = 1.0 / blackbox_high_resolution_scale * data(:, ind_bb_high_res);
end

% Unscale and remap sinarg
sinargScale = 5.0e3;
data(:,ind.sinarg) = 1.0 / sinargScale * data(:,ind.sinarg);

% Assign negative sign for pid error
data(:,ind.axisError) = -data(:,ind.axisError);

% Create an additional entry for the pi sum
data = [data, data(:,ind.axisP) + data(:,ind.axisI)];

% Create different sampling times
Ts      = para.looptime * 1.0e-6;             % Gyro loop
Ts_cntr = para.pid_process_denom * Ts;        % Control loop
Ts_log  = para.frameIntervalPDenom * Ts_cntr; % Logging loop

% Get evaluation index where Chirp was active
ind_eval = get_ind_eval(data(:,ind.sinarg), data(:,ind.gyroADC(ind_ax)));
data(~ind_eval,ind.sinarg) = 0.0;
T_eval_tot = size(data(ind_eval,ind.sinarg), 1) * Ts_log

% Calculate average throttle
throttle_avg = median(data(ind_eval,ind.setpoint(4))) / 1.0e3;


%% show Gyro to select Teval and spectra (gyro and pid sum)

figure(1)
ax(1) = subplot(311);
plot(ax(1), time, data(:,[ind.setpoint(1), ind.gyroUnfilt(1), ind.gyroADC(1)])), grid on, ylabel('Roll (deg/sec)')
title('Gyro Signals')
if do_insert_legends, legend('setpoint', 'gyro', 'gyroADC', 'location', 'best'), end
ax(2) = subplot(312);
plot(ax(2), time, data(:,[ind.setpoint(2), ind.gyroUnfilt(2), ind.gyroADC(2)])), grid on, ylabel('Pitch (deg/sec)')
ax(3) = subplot(313);
plot(ax(3), time, data(:,[ind.setpoint(3), ind.gyroUnfilt(3), ind.gyroADC(3)])), grid on, ylabel('Yaw (deg/sec)'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, time(end)])
set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)


% Select data for spectra
data_for_spectra = data(:,[ind.gyroUnfilt, ...
                           ind.gyroADC, ...
                           ind.axisSum, ...
                           ind.setpoint(1:3)]);

% Parameters
Nest     = round(2.0 / Ts_log);
koverlap = 0.9;
Noverlap = floor(koverlap * Nest);
window   = hann(Nest, 'periodic');

%=Spectra==================================================================

[pxx, freq] = estimate_spectra(data_for_spectra, window, Noverlap, Nest, Ts_log);
spectra = sqrt(pxx); % power -> amplitude (dc needs to be scaled differently)

figure(2)
ax(1) = subplot(211);
plot(ax(1), freq, spectra(:, 1:6)), grid on, ylabel('Gyro (deg/sec)'), set(gca, 'YScale', 'log')
title('Magnitude Spectra')
if do_insert_legends, legend('gyro Roll', 'gyro Pitch', 'gyro Yaw', 'gyroADC Roll', 'gyroADC Pitch', 'gyroADC Yaw', 'location', 'best'), end
ax(2) = subplot(212);
plot(ax(2), freq, spectra(:, 7:9)), grid on, ylabel('AxisSum'), xlabel('Frequency (Hz)'), set(gca, 'YScale', 'log')
if do_insert_legends, legend('axisSum Roll', 'axisSum Pitch', 'axisSum Yaw', 'location', 'best'), end
linkaxes(ax), clear ax, axis([0 1/2/Ts_log 1e-3 1e1])
set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)


%%
%=Spectrogram==============================================================

if (do_show_spec_figures)

    % Parameters
    Nest     = round(0.2 / Ts_log);
    koverlap = 0.9;
    Noverlap = floor(koverlap * Nest);
    window   = hann(Nest, 'periodic');
    Nres     = floor(max(data(:,ind.setpoint(4))) / 1e1 / 2) % should give 40 at 80% throttle constrain

    c_lim = [5e-2 3e0];

    for spectrogram_nr = 1:3
        [pxx, freq, throttle] = estimate_spectrogram(data(:,ind.gyroUnfilt(spectrogram_nr)), ...
                                                     data(:,ind.setpoint(4)) / 10.0, ...
                                                     window, Noverlap, Nest, Nres, Ts_log);
        spectrograms = sqrt(pxx); % power -> amplitude (dc needs to be scaled differently)
        
        figure(22)
        sgtitle('Gyro Spectrograms')
        axes_labels = {'Roll', 'Pitch', 'Yaw'};
        subplot(230 + spectrogram_nr)
        qmesh = pcolor(freq, throttle, spectrograms);
        set(qmesh, 'EdgeColor', 'None');
        % xlabel('Frequency (Hz)')
        if spectrogram_nr == 1
            ylabel('Throttle (%)')
        end
        title([axes_labels{spectrogram_nr}, ' – ohne Filter'])
        % colorbar()
        colormap('jet')
        set(gca, 'ColorScale', 'log')
        clim(c_lim);
        ylim([0 100])
    end

   

    for spectrogram_nr = 1:3
        [pxx, freq, throttle] = estimate_spectrogram(data(:,ind.gyroADC(spectrogram_nr)), ...
                                                     data(:,ind.setpoint(4)) / 10.0, ...
                                                     window, Noverlap, Nest, Nres, Ts_log);
        spectrograms = sqrt(pxx); % power -> amplitude (dc needs to be scaled differently)
        
        figure(22)
        subplot(230 + spectrogram_nr + 3)
        qmesh = pcolor(freq, throttle, spectrograms);
        set(qmesh, 'EdgeColor', 'None');
        xlabel('Frequency (Hz)')
        if spectrogram_nr == 1
            ylabel('Throttle (%)')
        end
        title([axes_labels{spectrogram_nr}, ' – mit Filter'])
        % colorbar()
        colormap('jet')
        set(gca, 'ColorScale', 'log')
        clim(c_lim);
        ylim([0 100])
    end
end


%% Some relevant fligth data

figure(3)
ax(1) = subplot(411);
plot(ax(1), time, data(:,ind.gyroUnfilt)), grid on, ylabel('Gyro (deg/sec)')
ax(2) = subplot(412);
plot(ax(2), time, data(:,ind.axisSum)), grid on, ylabel('AxisSum')
ax(3) = subplot(413);
plot(ax(3), time, data(:,ind.motor)), grid on, ylabel('Motor')
ax(4) = subplot(414);
plot(ax(4), time, data(:,ind.setpoint(4))), grid on, ylabel('Throttle'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, time(end)])
set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)


%% 
%=frequency_response=======================================================

est = frequency_response.frequency_response_estimator(data, ind, ind_ax, ind_eval, Ts_log, para, throttle_avg, Ts_cntr);

freq_resp_result = est.estimate();

Nest = est.Nest;

% Index and frequency for bode plots
omega_bode = 2*pi*freq_resp_result.P.Frequency;

% figure(5) Compare analytical to estimated controllers
frequency_response.show_frequency_response(freq_resp_result.Cpi, freq_resp_result.Cd, ...
                                           freq_resp_result.Cpi_ana, freq_resp_result.Cd_ana, ...
                                           omega_bode, opt, do_insert_legends, linewidth, ...
                                           @expand_multiple_figure_nr, multp_fig_nr);



%% Plant and used controllers

figure(expand_multiple_figure_nr(4, multp_fig_nr))
ax(1) = subplot('Position', pos_bode(1,:));
opt.YLim = {[1e-4 1e2], [-180 180]}; opt.MagScale = 'log';
bode(ax(1), freq_resp_result.P / freq_resp_result.Gf_ana, 'k', omega_bode, opt), title('Plant P')
hold off, grid on
ax(2) = subplot('Position', pos_bode(2,:));
opt.YLimMode = {'auto'}; opt.MagScale = 'linear';
bodemag(ax(2), freq_resp_result.C_T * freq_resp_result.C_Guw, 'k', omega_bode, opt), title(''), ylabel('Coherence')
linkaxes(ax, 'x'), clear ax
set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)


%% 
%=flight_parameters========================================================

% New controller and filter parameters

tic

pid_axis = {'rollPID', 'pitchPID', 'yawPID'};

% PID parameters
fprintf('   used PID parameters are:\n');
fprintf(['      ', pid_axis{ind_ax}, ':  %d, %d, %d\n'], ...
    para.(pid_axis{ind_ax})(1:3));

% Inform user about parameters
para_used_fieldnames = fieldnames(freq_resp_result.para_used);
Npara_used = size(para_used_fieldnames, 1);
fprintf('   used parameters are:\n');
for i = 1:Npara_used
    fprintf(['      ', para_used_fieldnames{i},': %d\n'], eval(['round(', 'freq_resp_result.para_used.', para_used_fieldnames{i}, ');']));
end

% First create new parameters the same as the actual ones
para_new = para;

% You can use the following command to generate the text below for the 
% actual parameters
% get_switch_case_text_from_para(para)


% type: 0: PT1, 1: BIQUAD, 2: PT2, 3: PT3
para_new.gyro_lpf            = 0;       % dono what this is
para_new.gyro_lowpass_hz     = 0;       % frequency of gyro lpf 1
para_new.gyro_soft_type      = 0;       % type of gyro lpf 1

para_new.gyro_lowpass_dyn_hz = [0, 0];  % dyn gyro lpf overwrites gyro_lowpass_hz
para_new.gyro_lowpass2_hz    = 800;     % frequency of gyro lpf 2
para_new.gyro_soft2_type     = 0;       % type of gyro lpf 2

para_new.gyro_notch_hz       = [0, 0]; % frequency of gyro notch 1 and 2
para_new.gyro_notch_cutoff   = get_fcut_from_D_and_fcenter([0.00, 0.00], para_new.gyro_notch_hz); % damping of gyro notch 1 and 2

para_new.dterm_lpf_hz        = 0;       % frequency of dterm lpf 1
para_new.dterm_filter_type   = 0;       % type of dterm lpf 1
para_new.dterm_lpf_dyn_hz    = [0, 0];  % dyn dterm lpf overwrites dterm_lpf_hz
para_new.dterm_lpf2_hz       = 120;     % frequency of dterm lpf 2
para_new.dterm_filter2_type  = 3;       % type of dterm lpf 2

para_new.dterm_notch_hz      = 0;     % frequency of dterm notch
para_new.dterm_notch_cutoff  = get_fcut_from_D_and_fcenter(0.00, para_new.dterm_notch_hz); % damping of dterm notch

para_new.yaw_lpf_hz          = 200;     % frequency of yaw lpf (pt1)

switch ind_ax
    case 1 % roll: [33, 52, 26, 0]
        P_new       = 1.0 * 33;
        I_ratio_new = 1.0 * 52/52;
        D_new       = 1.0 * 26;
    case 2 % pitch: [58, 98, 44, 0]
        P_new       = 1.0 * 58;
        I_ratio_new = 1.0 * 98/98;
        D_new       = 1.0 * 44;
    case 3 % yaw: [42, 65, 3, 0]
        P_new       = 1.0 * 42;
        I_ratio_new = 1.0 * 65/65;
        D_new       = 1.0 * 3;
end

% Scale to new PID parameters
pid_scale = [get_pid_scale(ind_ax), 1];
PID_new(1) = P_new * pid_scale(1);
fI         = freq_resp_result.PID(2) / (2 * pi * freq_resp_result.PID(1)); % extract fn from initial parametrization
fI_new     = fI * I_ratio_new;
PID_new(2) = 2 * pi * PID_new(1) * fI_new;
PID_new(3) = D_new * pid_scale(3);
PID_new(4) = 0;

fprintf('   used fI is: %0.2f Hz\n\n', fI);

% New PID parameters
fprintf('   new PID parameters are:\n');
para_new.(pid_axis{ind_ax}) = round( PID_new ./ pid_scale);
para_new.(pid_axis{ind_ax}) = [para_new.(pid_axis{ind_ax})(1:3), ...
                               para_new.(pid_axis{ind_ax})(3), ...
                               para_new.(pid_axis{ind_ax})(4)];
fprintf(['      ', pid_axis{ind_ax}, ':  %d, %d, %d\n'], ...
    para_new.(pid_axis{ind_ax})(1:3));

[Cpi_ana_new, Cd_ana_new, Gf_ana_new, PID_new, para_used_new] = ...
    calculate_transfer_functions(para_new, ind_ax, throttle_avg, Ts_cntr);

% Inform user about new parameters
para_used_fieldnames_new = fieldnames(para_used_new);
Npara_used_new = size(para_used_fieldnames_new, 1);
fprintf('   new parameters are:\n');
for i = 1:Npara_used_new
    fprintf(['      ', para_used_fieldnames_new{i},': %d\n'], ...
        eval(['round(', 'para_used_new.', para_used_fieldnames_new{i}, ');']));
end

fprintf('   new used fI is: %0.2f Hz\n\n', fI_new);

% Downsample analytical controller transferfunction and convert to frd objects
if Gf_ana_new.Ts < Ts_log % by using Gf_ana.Ts we secure that we do this only once
    Gf_ana_new  = downsample_frd(Gf_ana_new , Ts_log, freq_resp_result.P.Frequency);
    Cpi_ana_new = downsample_frd(Cpi_ana_new, Ts_log, freq_resp_result.P.Frequency);
    Cd_ana_new  = downsample_frd(Cd_ana_new , Ts_log, freq_resp_result.P.Frequency);
end

%% =Controller Analysis====================================================

CL_ana     = controller_analysis.calculate_closed_loop(freq_resp_result.Cpi_ana    , tf(1,1,Ts_log), freq_resp_result.P / freq_resp_result.Gf_ana, freq_resp_result.Gf_ana    , freq_resp_result.Cd_ana    );
CL_ana_new = controller_analysis.calculate_closed_loop(Cpi_ana_new, tf(1,1,Ts_log), freq_resp_result.P / freq_resp_result.Gf_ana, Gf_ana_new, Cd_ana_new);
if do_compensate_iterm
    % Compensate only PI part
    Cpi_com = freq_resp_result.Cpi / freq_resp_result.Cpi_ana;
    CL_ana_      = calculate_closed_loop(freq_resp_result.Cpi_ana     * Cpi_com, tf(1,1,freq_resp_result.Ts_log), freq_resp_result.P / freq_resp_result.Gf_ana, freq_resp_result.Gf_ana    , freq_resp_result.Cd_ana    );
    CL_ana_new_  = calculate_closed_loop(Cpi_ana_new * Cpi_com, tf(1,1,freq_resp_result.Ts_log), P / Gf_ana, Gf_ana_new, Cd_ana_new);
    CL_ana.T     = CL_ana_.T;
    CL_ana_new.T = CL_ana_new_.T;
end

% 
% uncomment für packages und den ganzen teil fürs plotten unten löschen
% controller_analysis.show_controller_analysis(CL_ana, CL_ana_new, T, omega_bode, do_insert_legends);

% Closed-loop bode plots (gang of four)
figure(expand_multiple_figure_nr(6, multp_fig_nr))
ax(1) = subplot(221);
opt.YLim = {[1e-3 1e1], [-180 180]}; opt.MagScale = 'log';
bodemag(ax(1), CL_ana.T , CL_ana_new.T , freq_resp_result.T, omega_bode, opt), title('Tracking T')
if do_insert_legends, legend('actual', 'new', 'location', 'best'), end
ax(2) = subplot(222);
bodemag(ax(2), CL_ana.S , CL_ana_new.S , omega_bode, opt), title('Sensitivity S')
ax(3) = subplot(223);
opt.YLim = {[1e-2 1e2], [-180 180]};
bodemag(ax(3), CL_ana.SC, CL_ana_new.SC, omega_bode, opt), title('Controller Effort SC')
ax(4) = subplot(224);
opt.YLim = {[1e-3 1e1], [-180 180]};
bodemag(ax(4), CL_ana.SP, CL_ana_new.SP, omega_bode, opt), title('Compliance SP')
linkaxes(ax, 'x'), clear ax
set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)

% Step responses
f_max = min([para.dyn_notch_min_hz, para.gyro_rpm_notch_min]);
T_mean = 0.1 * [-1, 1] + (Nest * Ts_log) / 2;
step_time = (0:Nest-1).'*Ts_log;

% Actual controller parameters
step_resp = [calculate_step_response_from_frd(CL_ana.T    , f_max), ...
             calculate_step_response_from_frd(CL_ana_new.T, f_max), ...
             calculate_step_response_from_frd(freq_resp_result.T, f_max)];
step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2),:));
step_resp = step_resp ./ step_resp_mean;

figure(expand_multiple_figure_nr(7, multp_fig_nr))
ax(1) = subplot(211);
plot(ax(1), step_time, step_resp), grid on, ylabel('Gyro (deg/sec)')
title('Tracking T')
if do_insert_legends, legend('actual', 'new', 'location', 'best'), end
ylim([0 1.3])

% New controller parameters
step_resp = [calculate_step_response_from_frd(CL_ana.SP    , f_max), ...
             calculate_step_response_from_frd(CL_ana_new.SP, f_max)];
step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2),:));
step_resp = step_resp - step_resp_mean;

ax(2) = subplot(212);
plot(ax(2), step_time, step_resp), grid on
title('Compliance SP'), xlabel('Time (sec)'), ylabel('Gyro (deg/sec)')
ylim([-0.2 1.1])
linkaxes(ax, 'x'), clear ax, xlim([0 0.5])
set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)

% Controllers
figure(expand_multiple_figure_nr(8, multp_fig_nr))
opt.YLim = {[1e-1 1e2], [-180 180]};
bode(CL_ana.C, CL_ana_new.C, omega_bode, opt)
title('Controller C')
if do_insert_legends, legend('actual', 'new', 'location', 'best'), end
set(findall(gcf, 'type', 'line'), 'linewidth', linewidth)

toc
