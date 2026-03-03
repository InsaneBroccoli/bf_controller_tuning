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


classdef angle_ctrl_tuning < handle
    properties
        data_flight     % Get all data and parameters
        gyro_tuning     % Get values of the inner cascade
        TargetAngle
    
        % Needed in different functions
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
    
        % Filters / controller
        Gf_ana
        Gf_new
        Cp_new
        Cp_ana
    
        % store used cutoffs for step-response f_max
        angle_fcut_ana
        angle_fcut_new
    
        % analysis helpers
        omega_bode
        Nest
        step_time
    
        % closed-loop results
        T_ana
        T_ana_new
        Se_ana
        Se_ana_new
    
        % step response
        step_resp_tra
    end
    methods

        function obj = angle_ctrl_tuning(data_flight, gyro_tuning)

            obj.data_flight = data_flight;
            obj.gyro_tuning = gyro_tuning;
        end

       function obj = calculate_Angle_trans(obj, Nestfatra, koverlaptra)
            dataf = obj.data_flight;
        
            % Analysis window parameters
            obj.Nest   = round(Nestfatra / dataf.Ts_log);
            Noverlap   = floor(koverlaptra * obj.Nest);
            window     = hann(obj.Nest, 'periodic');
        
            % Excitation filter
            Dlp = sqrt(3) / 2;
            wlp = 2 * pi * 10;
            Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), dataf.Ts_log, 'tustin');
        
            % Preallocate cells for 2 axes
            n_axes  = 2;
            obj.P   = cell(1, n_axes);
            obj.C_T = cell(1, n_axes);
            obj.Coh = cell(1, n_axes);
        
            % --- Define analysis filter once ---
            fc = 50;  % Hz  (dein "Angle PT3")
            obj.angle_fcut_ana = fc;
            obj.Gf_ana = get_filter('pt3', fc, dataf.Ts_cntr);  % should be discrete
            Cp_ana_ax = ss(dataf.para.levelPID(1)/10);
        
            sinarg_full = dataf.data(:, dataf.ind.sinarg);
        
            for ind_axis = 1:n_axes
        
                ind_eval = get_ind_eval( ...
                    dataf.data(:,dataf.ind.sinarg), ...
                    dataf.data(:,dataf.ind.gyroADC(ind_axis)));
        
                sinarg_ax = sinarg_full;
                sinarg_ax(~ind_eval) = 0;
        
                % Input signal
                w = dataf.data(:, dataf.ind.gyroADC(ind_axis));
                inp = apply_rotfiltfilt(Glp, sinarg_ax, w);
        
                % Output signal (Angle)
                y = dataf.data(:, dataf.ind.heading(ind_axis))*100;
                out_y = apply_rotfiltfilt(Glp, sinarg_ax, y);
        
                % FRD estimate (this returns FRD on FFT-friendly grid)
                [P_ax, C_T_ax] = estimate_frequency_response( ...
                    inp(ind_eval), out_y(ind_eval), window, Noverlap, obj.Nest, dataf.Ts_log);
        
                obj.C_T{ind_axis} = C_T_ax;
                obj.Coh{ind_axis} = C_T_ax;
                obj.omega_bode    = 2*pi*P_ax.Frequency;

                  % Downsample to logging grid if necessary
                if obj.Gf_ana.Ts < dataf.Ts_log
                    obj.Gf_ana  = downsample_frd(obj.Gf_ana , dataf.Ts_log, P_ax.Frequency);
                    Cp_ana_ax = downsample_frd(Cp_ana_ax, dataf.Ts_log, P_ax.Frequency);
                end
                
                obj.Cp_ana{ind_axis} = Cp_ana_ax;
                % Remove filter from identified response (keep FRD)
                % obj.P{ind_axis} = P_ax / obj.Gf_ana;
                obj.P{ind_axis} = P_ax;
            end
       end

       function obj = calculate_new_controller(obj, P_new, default_parameters, para_new)

            dataf = obj.data_flight;
        
            if default_parameters
                P_new = dataf.para.levelPID(1);
            end
        
            fprintf('   used P parameter Angle Control are:\n');
            fprintf('      P: %d\n', dataf.para.levelPID(1));
          
            % Store new controller as discrete gain on Ts_log
            Cp = tf(P_new/10, 1, dataf.Ts_log);   % diskrete TF
            Cp = ss(Cp);                         % Umwandlung in Zustandsraum
            obj.Cp_new = Cp;
        
            % New angle filter (PT3) on Ts_log
            obj.Gf_new = get_filter('pt3', para_new.angle_lpf_hz, dataf.Ts_log);
        end


        function obj = get_tuning_data(obj, ind_ax)

            gt    = obj.gyro_tuning;
            dataf = obj.data_flight;
        
            % --- Axis pick ---
            Tg = gt.T{ind_ax};     % FRD (inner loop complementary sensitivity)
            Pg = obj.P{ind_ax};    % FRD (measured rate -> angle plant, filter removed)
            
            % --- Make FRDs consistent without resampling Tg/Pg (avoid Hz/rad/s traps) ---
            Tg_frd = Tg;  Tg_frd.Ts = dataf.Ts_log;
            Pg_frd = Pg;  Pg_frd.Ts = dataf.Ts_log;
            
            % Use the measured frequency grid (in Hz) and convert to rad/s ONLY for frd(sys,omega)
            freqHz = Pg_frd.Frequency(:);     % Hz
            omega  = 2*pi*freqHz;             % rad/s
            
            % Convert discrete models to FRD on the SAME omega grid
            Cp_new_frd = frd(obj.Cp_new, omega);  Cp_new_frd.Ts = dataf.Ts_log;
            Gf_new_frd = frd(obj.Gf_new, omega);  Gf_new_frd.Ts = dataf.Ts_log;
            
            % Old models also on the same grid (for fair comparison)
            Cp_ana_frd = frd(obj.Cp_ana{ind_ax}, omega);  Cp_ana_frd.Ts = dataf.Ts_log;
            Gf_ana_frd = frd(obj.Gf_ana, omega);         Gf_ana_frd.Ts = dataf.Ts_log;
            
            % --- OPTIONAL: if response looks "flipped", try a sign correction on Pg ---
            % (use phase check instead of guessing if you want)
            % If your step goes negative or phase looks wrong, uncomment:
            % Pg_frd = -Pg_frd;
            
            % --- Closed loop (old) ---
            L_ana  = Cp_ana_frd .* Tg_frd .* Gf_ana_frd .* Pg_frd;
            Cl_ana = L_ana / (1 + L_ana);
            S_ana  = 1 / (1 + L_ana);
            
            % --- Closed loop (new) using measured Pg (NOT Pg_test) ---
            L_new      = Cp_new_frd .* Tg_frd .* Gf_new_frd .* Pg_frd;
            Cl_ana_new = L_new / (1 + L_new);
            S_ana_new  = 1 / (1 + L_new);
            
            % Store
            obj.T_ana      = Cl_ana;
            obj.T_ana_new  = Cl_ana_new;
            obj.Se_ana     = S_ana;
            obj.Se_ana_new = S_ana_new;

        
            % --- f_max choice for PT3 + cap to FRD max ---
            f_nyq    = 1/(2*dataf.Ts_log);
            f_cut    = obj.angle_fcut_ana;           % use analysis filter cutoff
            f_frdmax = min([max(Pg.Frequency), max(Tg.Frequency)]);
        
            f_max = min([50, 0.25*f_nyq, f_frdmax]);
        
            % Time vector (consistent with Nest)
            T_mean = 0.1 * [-1, 1] + (obj.Nest * dataf.Ts_log) / 2;
            obj.step_time = (0:obj.Nest-1).' * dataf.Ts_log;
        
            % Step responses (FRD in -> your IFFT method works)
            step_resp = [ ...
                calculate_step_response_from_frd(Cl_ana    , f_max), ...
                calculate_step_response_from_frd(Cl_ana_new, f_max) ];
        
            % Normalize around mean window
            idx_mean = (obj.step_time > T_mean(1)) & (obj.step_time < T_mean(2));
            step_resp_mean = mean(step_resp(idx_mean,:), 1);
            obj.step_resp_tra = step_resp ./ step_resp_mean;
        end
    end
    
end