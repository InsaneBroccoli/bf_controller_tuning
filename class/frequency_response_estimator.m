% FREQUENCY_RESPONSE_ESTIMATOR
% This class estimates the frequency response of PID controllers
% It estimates how your drone reacts to control inputs using flight logs.
% Helps you understand and tune your PID settings (especially PI and D)

classdef frequency_response_estimator
    properties
        % ====== Input & Data Properties ======
        data            % Logged data
        ind             % signal indices
        ind_ax          % Axis index (roll = 1, pitch = 2, yaw = 3)
        ind_eval        % Indices defining the evaluation time range
        Ts_log          % Sampling time of the logged data
    
        % ====== Controller and System Parameters ======
        para            % Controller parameter structure
        throttle_avg    % Average throttle value during logging
        Ts_cntr         % Controller sampling time
        data_sinarg     % 
        setpoint        % Setpoint for rotfiltfilt

        % ====== Generate filters ======
        Dlp = sqrt(3)/2;
        wlp = 2 * pi * 10;  % Cutoff frequency: 10 Hz
        Glp = c2d(tf(wlp^2, [1 2*Dlp*wlp wlp^2]), Ts_log, 'tustin');

        % ====== Calculate base parameters for spectral estimation ======
        Nest = round(2 / Ts_log);                   % Window size for 2 seconds
        koverlap = 0.9;                             % 90% overlap
        Noverlap = floor(obj.koverlap * obj.Nest);  % Overlap in samples
        window = hann(obj.Nest, 'periodic');        % Hann window
    
        % ====== Signal Indexing for System Identification ======
        axisSum         % Index or signal for summed control outputs (all axes)
        axisSumPI       % Index or signal for summed PI controller outputs
    
        % ====== Estimated Frequency-Domain Results ======
        T               % Tracking transfer function (w -> y)
        C_T             % Coherence of tracking transfer function T
        Guw             % Transfer function from input to control signal (w -> u)
        C_Guw           % Coherence of Guw
        Gvw             % Transfer function from input to PI controller output (w -> v)
        C_Gvw           % Coherence of Gvw
        P               % Identified plant transfer function
        Cpi             % Identified PI controller transfer function
        Cd              % Identified D controller transfer function
    
        % ====== Analytical (Model-Based) Transfer Functions ======
        Gf_ana          % Analytical filter response
        Cpi_ana         % Analytical PI controller transfer function
        Cd_ana          % Analytical D controller transfer function
    
        % ====== Additional Information ======
        omega_bode      % Frequency vector for Bode plot representation
        PID             % PID parameters obtained from analytical model
        para_used       % Parameter set used in analytical transfer function calculation
    end

    methods
        % =================================================================
        %  CONSTRUCTOR
        % =================================================================
        function obj = frequency_response_estimator(data, ind, ind_ax, ind_eval, Ts_log, para, throttle_avg, Ts_cntr, data_sinarg, setpoint)
            obj.data = data;
            obj.ind = ind;
            obj.ind_ax = ind_ax;
            obj.ind_eval = ind_eval;
            obj.Ts_log = Ts_log;
            obj.para = para;
            obj.throttle_avg = throttle_avg;
            obj.Ts_cntr = Ts_cntr;
            obj.data_sinarg = data_sinarg;
            obj.setpoint = setpoint;
        end

        % =================================================================
        %  APPLY ROTFILTFILT FUNCTION - Quadrature rotation + zero-phase filtering + back-rotation
        % =================================================================
        function xf = apply_rotfiltfilt(sinarg, x)

            % PURPOSE
            %   - Demodulate a narrowband component by complex rotation
            %   - Apply zero-phase IIR filtering to the in-phase and quadrature parts
            %   - Remodulate back to the original carrier to obtain a cleaned real signal
            %
            % INPUTS
            %   - Glp      IIR filter struct with fields G.num{1}, G.den{1} for filtfilt
            %   - sinarg   [N x 1] phase argument per sample in radians, e.g. 2*pi*f0*t
            %   - setpoint [N x M] real signal matrix, each column a channel
            %
            % OUTPUTS
            %   - xf       [N x M] real filtered signal after rotate–filter–inverse-rotate
            %
            % METHOD
            %   1) Remove column mean from x
            %   2) Form complex phasor p = exp(1j*sinarg)
            %   3) Rotate yR = y .* p and yQ = y .* conj(p) to isolate sidebands
            %   4) Apply zero-phase filtering yR = filtfilt(G.num{1}, G.den{1}, yR) and same for yQ
            %   5) Back-rotate and recombine xf = real((yR.*conj(p) + yQ.*p) * 0.5)
            %
            % NOTES
            %   - Uses filtfilt which requires a stable, causal IIR and typically Signal Processing Toolbox
            %   - sinarg must be scalar or match the first dimension of x
            %   - The 0.5 scaling is omitted here; relative scaling cancels in ratio-based analyses
            %   - Ensure G.den{1}(1) == 1 or normalize the filter coefficients before calling

            % Get signal dimensions
            [Nx, nx] = size(obj.setpoint);
            xf = zeros(Nx, nx);
            
            % Create complex phasor
            p = exp(1i * obj.sinarg);
            
            for i = 1:nx
                % Remove mean from signal
                y = x(:,i) - mean(x(:,i));
                
                % Rotate to isolate sidebands
                yR = y .* p;
                yQ = y .* conj(p);
                
                % Apply zero-phase filtering
                yR = filtfilt(obj.Glp.num{1}, obj.Glp.den{1}, yR);
                yQ = filtfilt(obj.Glp.num{1}, obj.Glp.den{1}, yQ);
                
                % Back-rotate and recombine
                xf(:,i) = real((yR.*conj(p) + yQ.*p) * 0.5);
            end
        end

        % =================================================================
        %  ESTIMATION FUNCTION
        % =================================================================
        function obj = estimate(obj)
            % T  , Gyw: w -> y
            inp = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.setpoint(ind_ax)));
            out = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.gyroADC(ind_ax)) );
            [T, C_T] = estimate_frequency_response(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, Ts_log);
            
            % SCw, Guw: w -> u
            out = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.axisSum(ind_ax)));
            [Guw, C_Guw] = estimate_frequency_response(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, Ts_log);
            
            %      Gvw: w -> v (v := u only from PI cntrl)
            out = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.axisSumPI(ind_ax)));
            [Gvw, C_Gvw] = estimate_frequency_response(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, Ts_log);
            
            % P  , Gyu: u -> y
            P = T / Guw;
            
            % % P  , Gyu: u -> y (direct measurement, results are slightly worse)
            % inp = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.axisSum(ind_ax)));
            % out = apply_rotfiltfilt(Glp, data(:,ind.sinarg), data(:,ind.gyroADC(ind_ax)));
            % [Pd, C_Pd] = estimate_frequency_response(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, Ts_log);
            
            % Calculated controller frequency response estimates
            Cpi = Gvw / (1 - T);
            Cd  = Guw * Gvw / T * (1 / Guw - 1 / Gvw);
            
            % Index and frequency for bode plots
            omega_bode = 2*pi*P.Frequency;
            
            
            %% Downsample analytical controller transferfunction and convert to frd objects
            
            [Cpi_ana, Cd_ana, Gf_ana, PID, para_used] = ...
                calculate_transfer_functions(para, ind_ax, throttle_avg, Ts_cntr);
            
            if Gf_ana.Ts < Ts_log % by using Gf_ana.Ts we secure that we do this only once
                Gf_ana  = downsample_frd(Gf_ana , Ts_log, P.Frequency);
                Cpi_ana = downsample_frd(Cpi_ana, Ts_log, P.Frequency);
                Cd_ana  = downsample_frd(Cd_ana , Ts_log, P.Frequency);
            end
        end
    end
end
