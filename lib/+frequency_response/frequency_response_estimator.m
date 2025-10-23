
% Frequency response estimation and calculation

classdef frequency_response_estimator
    properties
        data
        ind
        ind_ax
        ind_eval
        Ts_log
        para
        throttle_avg
        Ts_cntr
        Nest
        Noverlap
        koverlap
        window
        Glp
    end

    methods

        function obj = frequency_response_estimator(data, ind, ind_ax, ind_eval, Ts_log, para, throttle_avg, Ts_cntr)
            obj.data = data;
            obj.ind = ind;
            obj.ind_ax = ind_ax;
            obj.ind_eval = ind_eval;
            obj.Ts_log = Ts_log;
            obj.para = para;
            obj.throttle_avg = throttle_avg;
            obj.Ts_cntr = Ts_cntr;
            
            % calculate base parameters
            obj.Nest = round(2 / obj.Ts_log);
            obj.koverlap = 0.9;
            obj.Noverlap = floor(obj.koverlap * obj.Nest);
            obj.window = hann(obj.Nest, 'periodic');

            % generate filters
            Dlp = sqrt(3)/2;
            wlp = 2 * pi * 10;
            obj.Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');
        end

        function results = estimate(obj)

            % T  , Gyw: w -> y
            inp = apply_rotfiltfilt(obj.Glp, obj.data(:,obj.ind.sinarg), obj.data(:,obj.ind.setpoint(obj.ind_ax)));
            out = apply_rotfiltfilt(obj.Glp, obj.data(:,obj.ind.sinarg), obj.data(:,obj.ind.gyroADC(obj.ind_ax)) );
            [T, C_T] = estimate_frequency_response(inp(obj.ind_eval), out(obj.ind_eval), obj.window, obj.Noverlap, obj.Nest, obj.Ts_log);

            % SCw, Guw: w -> u
            out = apply_rotfiltfilt(obj.Glp, obj.data(:,obj.ind.sinarg), obj.data(:,obj.ind.axisSum(obj.ind_ax)));
            [Guw, C_Guw] = estimate_frequency_response(inp(obj.ind_eval), out(obj.ind_eval), obj.window, obj.Noverlap, obj.Nest, obj.Ts_log);

            % Gvw: w -> v (v := u only from PI cntrl)
            out = apply_rotfiltfilt(obj.Glp, obj.data(:,obj.ind.sinarg), obj.data(:,obj.ind.axisSumPI(obj.ind_ax)));
            [Gvw, C_Gvw] = estimate_frequency_response(inp(obj.ind_eval), out(obj.ind_eval), obj.window, obj.Noverlap, obj.Nest, obj.Ts_log);

            % calculations
            % P  , Gyu: u -> y
            P = T / Guw;

            % Calculated controller frequency response estimates
            Cpi = Gvw / (1 - T);
            Cd  = Guw * Gvw / T * (1 / Guw - 1 / Gvw);

            % Downsample analytical controller transferfunction and convert to frd objects
            [Cpi_ana, Cd_ana, Gf_ana, PID, para_used] = calculate_transfer_functions(obj.para, obj.ind_ax, obj.throttle_avg, obj.Ts_cntr);

            if Gf_ana.Ts < obj.Ts_log % by using Gf_ana.Ts we secure that we do this only once
                Gf_ana  = downsample_frd(Gf_ana , obj.Ts_log, P.Frequency);
                Cpi_ana = downsample_frd(Cpi_ana, obj.Ts_log, P.Frequency);
                Cd_ana  = downsample_frd(Cd_ana , obj.Ts_log, P.Frequency);
            end

            % save results
            results = struct( ...
                'T', T, ...
                'C_T', C_T, ...
                'Guw', Guw, ...
                'C_Guw', C_Guw, ...
                'Gvw', Gvw, ...
                'C_Gvw', C_Gvw, ...
                'P', P, ...
                'Cpi', Cpi, ...
                'Cd', Cd, ...
                'Cpi_ana', Cpi_ana, ...
                'Cd_ana', Cd_ana, ...
                'Gf_ana', Gf_ana, ...
                'PID', PID, ...
                'para_used', para_used ...
            );
        end
    end
end
