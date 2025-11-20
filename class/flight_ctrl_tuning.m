classdef flight_ctrl_tuning
    properties       
        % Raw flight data inputs

        time             % Time vector [s]
        data             % Flgiht data
        ind              % Columns
        para             % para
        Ts_log           % Logging Time [s]
        Ts_cntr

        % Neede in different function
         T
        Guw
        Gvw
        P
        Cpi
        Cd
        Gf_ana
        Cpi_ana
        Cd_ana
        PID
        para_used
        omega_bode
        C_T
        C_Guw

        % Default values
        linewidth            (1,1) double  = 1.2
            
    end

     methods
         function obj = flight_ctrl_tuning(data, ind, Ts_log, para, Ts_cntr)
            obj.data = data;
            obj.ind = ind;
            obj.Ts_log = Ts_log;
            obj.para = para;
            obj.Ts_cntr = Ts_cntr;


        end
        function obj = calculate_transfer_func(obj, Nestfatra, koverlaptra)
             % =============================================================
            %  Frequency response estimation and calculation
            % =============================================================
            % Estimates and calculates frequency responses for system identification
            % and controller analysis.
            %
            % Key components analyzed:
            % - T  (Complementary Sensitivity): Response from reference to output
            % - P  (Plant): Response from control input to output
            % - Cpi (PI Controller): Response from error to control input
            % - Cd (D Controller): Response from derivative to control input
            %
            % Signal definitions:
            % w: reference input (setpoint)
            % y: system output (gyro measurements)
            % u: control input (total PID output)
            % v: PI controller output
            
            % Get evaluation index where Chirp was active
            
            % Analysis window parameters
            Nest     = round(Nestfatra / obj.Ts_log);    % Window length in samples
            Noverlap = floor(koverlaptra * Nest);        % Overlap between windows
            window   = hann(Nest, 'periodic');               % Hanning window for analysis
            
            % Design linear filter for zero phase excitation filter (apply_rotfiltfilt)
            Dlp = sqrt(3) / 2;    % Damping ratio
            wlp = 2 * pi * 10;    % Cutoff frequency [rad/s]
            Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), ...    % Discrete filter
                    obj.Ts_log, 'tustin');                   % Using Tustin transform

            % Preallocate cells for 3 axes
            n_axes = 3;
            obj.T   = cell(1, n_axes);
            obj.Guw = cell(1, n_axes);
            obj.Gvw = cell(1, n_axes);
            obj.P   = cell(1, n_axes);
            obj.Cpi = cell(1, n_axes);
            obj.Cd  = cell(1, n_axes);
            obj.C_T   = cell(1, n_axes);
            obj.C_Guw = cell(1, n_axes);

            sinarg_full = obj.data(:, obj.ind.sinarg);  % Copy Data to adjust it

            for ind_ax = 1:n_axes       % Calculate Transferfunction for Roll, Pitch and Yaw
            
                ind_eval = get_ind_eval(obj.data(:,obj.ind.sinarg), obj.data(:,obj.ind.gyroADC(ind_ax)));
                
                sinarg_ax = sinarg_full;
                sinarg_ax(~ind_eval) = 0;

                % Calculate average throttle
                throttle_avg = median(obj.data(ind_eval,obj.ind.setpoint(4))) / 1.0e3;

                
                 % ----- Input signal: w (filtered setpoint for this axis) -----
                w = obj.data(:, obj.ind.setpoint(ind_ax));
                inp = apply_rotfiltfilt(Glp, sinarg_ax, w);
                
                % ----- Output y: filtered gyro for this axis -----
                y = obj.data(:, obj.ind.gyroADC(ind_ax));
                out_y = apply_rotfiltfilt(Glp, sinarg_ax, y);

                % Calculate complementary sensitivity (T) and input-output responses
                % T  , Gyw: w -> y
                [T_ax, C_T_ax] = estimate_frequency_response( ...
                    inp(ind_eval), out_y(ind_eval), window, Noverlap, Nest, obj.Ts_log);
                            
                % Calculate control sensitivity (Represents total controller output response)
                % SCw, Guw: w -> u
                u = obj.data(:, obj.ind.axisSum(ind_ax));
                out_u = apply_rotfiltfilt(Glp, sinarg_ax, u);
                [Guw_ax, C_Guw_ax] = estimate_frequency_response( ...
                    inp(ind_eval), out_u(ind_eval), window, Noverlap, Nest, obj.Ts_log);
                
                % Calculate PI controller response (v is the output of just the PI portion of controller)
                % Gvw: w -> v
                v = obj.data(:, obj.ind.axisSumPI(ind_ax));
                out_v = apply_rotfiltfilt(Glp, sinarg_ax, v);
                [Gvw_ax, ~] = estimate_frequency_response( ...
                    inp(ind_eval), out_v(ind_eval), window, Noverlap, Nest, obj.Ts_log);
                
                % Calculate plant response (indirect method: P = T/Guw for better noise immunity)
                % P  , Gyu: u -> y
                P_ax = T_ax / Guw_ax;
                               
                % Calculate controller frequency responses
                % Split into PI and D components for analysis
                Cpi_ax = Gvw_ax / (1 - T_ax);
                Cd_ax  = Guw_ax * Gvw_ax / T_ax * (1 / Guw_ax - 1 / Gvw_ax);

                % Prepare frequency vector for Bode plotspara
                omega_bode_ax = 2*pi*P_ax.Frequency;    % Convert Hz to rad/s
                

                % ----- Store for this axis -----
                obj.T{ind_ax}   = T_ax;
                obj.Guw{ind_ax} = Guw_ax;
                obj.Gvw{ind_ax} = Gvw_ax;
                obj.P{ind_ax}   = P_ax;
                obj.Cpi{ind_ax} = Cpi_ax;
                obj.Cd{ind_ax}  = Cd_ax;
                obj.omega_bode = omega_bode_ax;
                obj.C_T{ind_ax} = C_T_ax;
                obj.C_Guw{ind_ax} = C_Guw_ax;
            
                
                % Analytical controller transfer functions for this axis
                [Cpi_ana_ax, Cd_ana_ax, Gf_ana_ax, PID_ax, para_used_ax] = ...
                    calculate_transfer_functions(obj.para, ind_ax, throttle_avg, obj.Ts_cntr);
    
                 % Downsample to logging grid if necessary
                if Gf_ana_ax.Ts < obj.Ts_log
                    Gf_ana_ax  = downsample_frd(Gf_ana_ax , obj.Ts_log, P_ax.Frequency);
                    Cpi_ana_ax = downsample_frd(Cpi_ana_ax, obj.Ts_log, P_ax.Frequency);
                    Cd_ana_ax  = downsample_frd(Cd_ana_ax , obj.Ts_log, P_ax.Frequency);
                end
                 % Store results per axis
                obj.Gf_ana{ind_ax}  = Gf_ana_ax;
                obj.Cpi_ana{ind_ax} = Cpi_ana_ax;
                obj.Cd_ana{ind_ax}  = Cd_ana_ax;
                obj.PID{ind_ax} = PID_ax;
                obj.para_used{ind_ax} = para_used_ax;                
               
            end
        end

        function calculate_tuning_data()
            tic
            
            % Initialize axis names for reporting
            pid_axis = {'rollPID', 'pitchPID', 'yawPID'};

            % Handle default parameter case
            if obj.default_parameters
                obj.para_new = para;
                 pid_vec = para.(pid_axis{obj.ind_ax});

                 % Use existing parameters
                 obj.P_new = pid_vec(1);    % Current P gain
                 obj.I_new = pid_vec(2);       % Keep same I gain
                 obj.D_new = pid_vec(3);    % Current D gain
            end
            
            % Display used PID configuration
            fprintf('   used PID parameters are:\n');
            fprintf(['      ', pid_axis{obj.ind_ax}, ':  %d, %d, %d\n'], ...
                para.(pid_axis{obj.ind_ax})(1:3));
            
            % Display all used parameters
            para_used_fieldnames = fieldnames(para_used);
            Npara_used = size(para_used_fieldnames, 1);
            fprintf('   used parameters are:\n');
            for i = 1:Npara_used
                fprintf(['      ', para_used_fieldnames{i},': %d\n'], ...
                    eval(['round(', 'para_used.', para_used_fieldnames{i}, ');']));
            end
            
            % --- Calculate New PID Parameters ---
            % Get axis-specific scaling factors (roll/pitch/yaw have different gains)
            pid_scale = [get_pid_scale(obj.ind_ax), 1];    % Get axis scaling factors

            % Calculate new gains with scaling
            PID_new(1) = obj.P_new * pid_scale(1);         % Scale P gain
            fI         = PID(2) / (2 * pi * PID(1));       % Extract current I frequency
            PID_new(2) = obj.I_new *pid_scale(2);
            fI_new     = PID_new(2) / (2 * pi * PID_new(1));    % Scale I frequency
            PID_new(3) = obj.D_new * pid_scale(3);         % Scale D gain
            PID_new(4) = 0;                                % No feedforward
            
            % Display integral frequencies
            fprintf('   used fI is: %0.2f Hz\n\n', fI);
            
            % --- Update and Display New Parameters ---
            fprintf('   new PID parameters are:\n');

            % Round scaled PID values and store in parameter structure
            obj.para_new.(pid_axis{obj.ind_ax}) = round( PID_new ./ pid_scale);
            obj.para_new.(pid_axis{obj.ind_ax}) = [obj.para_new.(pid_axis{obj.ind_ax})(1:3), ...
                                           obj.para_new.(pid_axis{obj.ind_ax})(3), ...
                                           obj.para_new.(pid_axis{obj.ind_ax})(4)];

            % Display new PID configuration
            fprintf(['      ', pid_axis{obj.ind_ax}, ':  %d, %d, %d\n'], ...
                obj.para_new.(pid_axis{obj.ind_ax})(1:3));

            % --- Generate New Transfer Functions ---
            % Calculate analytical transfer functions with new parameters
            [Cpi_ana_new, Cd_ana_new, Gf_ana_new, PID_new, para_used_new] = ...
                calculate_transfer_functions(obj.para_new, obj.ind_ax, throttle_avg, Ts_cntr);
            
            % Display new parameter values
            para_used_fieldnames_new = fieldnames(para_used_new);
            Npara_used_new = size(para_used_fieldnames_new, 1);
            fprintf('   new parameters are:\n');
            for i = 1:Npara_used_new
                fprintf(['      ', para_used_fieldnames_new{i},': %d\n'], ...
                    eval(['round(', 'para_used_new.', para_used_fieldnames_new{i}, ');']));
            end
            
            fprintf('   new used fI is: %0.2f Hz\n\n', fI_new);
            
            % --- Process Transfer Functions ---
            % Downsample analytical controller transferfunction and convert to frd objects
            if Gf_ana_new.Ts < obj.Ts_log % by using Gf_ana.Ts we secure that we do this only once
                Gf_ana_new  = downsample_frd(Gf_ana_new , obj.Ts_log, P.Frequency);
                Cpi_ana_new = downsample_frd(Cpi_ana_new, obj.Ts_log, P.Frequency);
                Cd_ana_new  = downsample_frd(Cd_ana_new , obj.Ts_log, P.Frequency);
            end
            
            % --- Calculate Closed Loop Responses ---
            % Generate closed loop responses for both current and new parameters
            CL_ana     = calculate_closed_loop(Cpi_ana    , tf(1,1,obj.Ts_log), P / Gf_ana, Gf_ana    , Cd_ana    );
            CL_ana_new = calculate_closed_loop(Cpi_ana_new, tf(1,1,obj.Ts_log), P / Gf_ana, Gf_ana_new, Cd_ana_new);

            % Apply I-term compensation if enabled
            if obj.do_compensate_iterm
                % Compensate only PI part
                Cpi_com = Cpi / Cpi_ana;

                % Recalculate closed loop responses with compensation
                CL_ana_      = calculate_closed_loop(Cpi_ana     * Cpi_com, tf(1,1,obj.Ts_log), P / Gf_ana, Gf_ana    , Cd_ana    );
                CL_ana_new_  = calculate_closed_loop(Cpi_ana_new * Cpi_com, tf(1,1,obj.Ts_log), P / Gf_ana, Gf_ana_new, Cd_ana_new);

                % Update complementary sensitivity
                CL_ana.T     = CL_ana_.T;
                CL_ana_new.T = CL_ana_new_.T;
            end

            % Store final closed loop responses
            obj.CloLoAan = CL_ana;
            obj.CloLoAanNew = CL_ana_new;
            
            % =============================================================
            %  Step Response Analysis
            % =============================================================
            % Calculate and compare step responses for current and new parameters
            % Includes both tracking performance and disturbance rejection analysis
            %
            % Parameters:
            % - f_max: Maximum frequency from notch filter settings
            % - T_mean: Time window for response normalization
            % - Responses calculated for current and new parameters

            % --- Configure Analysis Parameters ---
            % Set frequency limit based on notch filter settings
            f_max = min([para.dyn_notch_min_hz, para.gyro_rpm_notch_min]);

            % Define time window for response analysis
            T_mean = 0.1 * [-1, 1] + (Nest * obj.Ts_log) / 2;    % Analysis window [s]
            step_time = (0:Nest-1).'*obj.Ts_log;                 % Time vector [s]
            
            % --- Tracking Performance Analysis ---
            % Calculate step responses for setpoint tracking
            step_resp = [calculate_step_response_from_frd(CL_ana.T    , f_max), ...    % Used parameters
                         calculate_step_response_from_frd(CL_ana_new.T, f_max), ...    % New parameters
                         calculate_step_response_from_frd(T           , f_max)];       % Measured response

            % Normalize responses around their mean values
            step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2),:));
            step_resp = step_resp ./ step_resp_mean;

            % Store tracking response results
            obj.step_resp_tra = step_resp;
            obj.step_resp_tim = step_time;
            
            % --- Disturbance Rejection Analysis ---
            % Calculate responses to disturbance inputs
            step_resp = [calculate_step_response_from_frd(CL_ana.SP    , f_max), ...    % Used parameters
                         calculate_step_response_from_frd(CL_ana_new.SP, f_max)];       % New parameters
            
            % Center responses around their mean
            step_resp_mean = mean(step_resp(step_time > T_mean(1) & step_time < T_mean(2),:));
            step_resp = step_resp - step_resp_mean;

            % Store disturbance response results
            obj.step_resp_com = step_resp;
                       
            toc
        end

        end


     end
end