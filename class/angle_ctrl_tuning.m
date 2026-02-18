%==========================================================================
% ANGLE CTRL TUNING - Betaflight Controller Analysis ANGLE TUNING CLASS
%==========================================================================
% Betaflight Controller Tuning Analysis Script
% Purpose: Calculation for Angle Tuning
%
% Author: [Janick Dort, Yuri Bianchi, Dario Jurietti]
% Supervisor: [Michael Peter]
% Date: [25.11.2025]

%==========================================================================
%  ADDITIONAL INFORMATION
% =============================================================
% Estimates and calculates frequency responses for system identification
% and controller analysis.


classdef angle_ctrl_tuning
    properties
        data_flight     % Get all data and parameters
        gyro_tuning     % Get values of the inner cascade
        TargetAngle
        % Neede in different function
        T
        Guw
        Gvw
        P
        Cpi
        Cd
        C_T
        C_Guw
        P_gef
        Coh
        omega_bode
    end
    methods

        function obj = angle_ctrl_tuning(data_flight, gyro_tuning)

            obj.data_flight = data_flight;
            obj.gyro_tuning = gyro_tuning;
        end

       function obj = get_angle_data(obj)
            dataf = obj.data_flight;
            angleLimit = 60;
        
            angleTarget = zeros(size(dataf.data,1), 2);   % prealloc
        
            for ind_axis = 1:2
                InvmaxRcRate = 1 / getMaxRcRate(ind_axis, dataf.para.rates_type, ...
                    dataf.para.rc_expo, dataf.para.rc_rates, dataf.para.rates);
        
                currentSetpoint = dataf.data(:, dataf.ind.setpoint(ind_axis));
                angleTarget(:, ind_axis) = angleLimit .* currentSetpoint .* InvmaxRcRate;
            end
        
            obj.TargetAngle = angleTarget;
       end
       function obj = calculate_Angle_trans(obj, Nestfatra, koverlaptra)
           dataf = obj.data_flight;

            % Analysis window parameters
            Nest     = round(Nestfatra / dataf.Ts_log);    % Window length in samples
            Noverlap = floor(koverlaptra * Nest);        % Overlap between windows
            window   = hann(Nest, 'periodic');               % Hanning window for analysis

            % Design linear filter for zero phase excitation filter (apply_rotfiltfilt)
            Dlp = sqrt(3) / 2;    % Damping ratio
            wlp = 2 * pi * 10;    % Cutoff frequency [rad/s]
            Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), ...    % Discrete filter
                    dataf.Ts_log, 'tustin');                   % Using Tustin transform

            % Preallocate cells for 2 axes
            n_axes = 2;
            obj.T   = cell(1, n_axes);
            obj.Guw = cell(1, n_axes);
            obj.Gvw = cell(1, n_axes);
            obj.P   = cell(1, n_axes);
            obj.Cpi = cell(1, n_axes);
            obj.Cd  = cell(1, n_axes);
            obj.C_T   = cell(1, n_axes);
            obj.C_Guw = cell(1, n_axes);
            obj.P_gef = cell(1, n_axes);
            obj.Coh = cell(1, n_axes);

            sinarg_full = dataf.data(:, dataf.ind.sinarg);  % Copy Data to adjust it

            for ind_axis = 1:n_axes       % Calculate Transferfunction for Roll, Pitch and Yaw
            
                % Check where the Chirp Signal was activated
                ind_eval = get_ind_eval(dataf.data(:,dataf.ind.sinarg), ...
                    dataf.data(:,dataf.ind.gyroADC(ind_axis)));

                sinarg_ax = sinarg_full;
                sinarg_ax(~ind_eval) = 0;

                % ----- Input signal: w (filtered setpoint for this axis) -----
                w = obj.TargetAngle(:, ind_axis);
                inp = apply_rotfiltfilt(Glp, sinarg_ax, w);

                % ----- Output y: Angle for this axis -----
                y = dataf.data(:, dataf.ind.heading(ind_axis));
                out_y = apply_rotfiltfilt(Glp, sinarg_ax, y);

                % Calculate complementary sensitivity (T) and input-output responses
                % T  , Gyw: Target Angle -> Current Angle
                [T_ax, C_T_ax] = estimate_frequency_response( ...
                    inp(ind_eval), out_y(ind_eval), window, Noverlap, ...
                    Nest, dataf.Ts_log);

                % Calculate control sensitivity (Represents total controller output response)
                % SCw, Guw: Target Angle -> Current Angular Rate
                u = dataf.data(:, dataf.ind.gyroADC);
                out_u = apply_rotfiltfilt(Glp, sinarg_ax, u);
                [Guw_ax, C_Guw_ax] = estimate_frequency_response( ...
                    inp(ind_eval), out_u(ind_eval), window, Noverlap, ...
                    Nest, dataf.Ts_log);

                % Calculate plant response (indirect method: P = T/Guw for better noise immunity)
                % P  , Gyu: Current Angular Rate -> Current Angle
                P = T_ax / Guw_ax;


                % Prepare frequency vector for Bode plotspara
                omega_bode_ax = 2*pi*T_ax.Frequency;    % Convert Hz to rad/s

                % ----- Store for this axis -----
                obj.T{ind_axis}   = T_ax;
                obj.C_T{ind_axis} = C_T_ax;
                obj.P{ind_axis} = P;
                obj.Coh{ind_axis} = C_T_ax * C_Guw_ax;
                obj.omega_bode = omega_bode_ax;
                
            end
       end

    end
end