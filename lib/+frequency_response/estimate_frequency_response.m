
% =================================================================
%  FREQUENCY RESPONSE ESTIMATION FUNCTION
% =================================================================
function results = estimate_frequency_response(data, ind, ind_ax, ind_eval, Ts_log, para, throttle_avg, Ts_cntr, data_sinarg)
    % ESTIMATE_FREQUENCY_RESPONSE - Estimates frequency response of PID controllers
    %
    % INPUTS:
    %   data         - Logged flight data
    %   ind          - Signal indices structure
    %   ind_ax       - Axis index (1=roll, 2=pitch, 3=yaw)
    %   ind_eval     - Indices for evaluation time range
    %   Ts_log       - Sampling time of logged data
    %   para         - Controller parameters structure
    %   throttle_avg - Average throttle during logging
    %   Ts_cntr      - Controller sampling time
    %   data_sinarg  - Sine argument data
    %   setpoint     - Setpoint for filtering
    %
    % OUTPUTS:
    %   results - Structure containing all frequency response results with fields:
    %       .T            - Tracking transfer function (w -> y)
    %       .C_T          - Coherence of T
    %       .Guw          - Transfer function (w -> u)
    %       .C_Guw        - Coherence of Guw
    %       .Gvw          - Transfer function (w -> v)
    %       .C_Gvw        - Coherence of Gvw
    %       .P            - Plant transfer function
    %       .Cpi          - PI controller transfer function
    %       .Cd           - D controller transfer function
    %       .Cpi_ana      - Analytical PI controller
    %       .Cd_ana       - Analytical D controller
    %       .Gf_ana       - Analytical filter response
    %       .omega_bode   - Frequency vector for Bode plots
    %       .PID          - PID parameters
    %       .para_used    - Parameters used in calculation

    % ====== Create lowpass filter ======
    Glp = create_lowpass_filter(Ts_log);
    
    % ====== Calculate spectral estimation parameters ======
    [Nest, Noverlap, window] = calculate_spectral_params(Ts_log);
    
    % ====== T, Gyw: w -> y ======
    inp = apply_rotfiltfilt(Glp, data_sinarg, data(:, ind.setpoint(ind_ax)));
    out = apply_rotfiltfilt(Glp, data_sinarg, data(:, ind.gyroADC(ind_ax)));
    [T, C_T] = tfestimate(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, 1/Ts_log);
    
    % ====== SCw, Guw: w -> u ======
    out = apply_rotfiltfilt(Glp, data_sinarg, data(:, ind.axisSum(ind_ax)));
    [Guw, C_Guw] = tfestimate(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, 1/Ts_log);
    
    % ====== Gvw: w -> v (v := u only from PI controller) ======
    out = apply_rotfiltfilt(Glp, data_sinarg, data(:, ind.axisSumPI(ind_ax)));
    [Gvw, C_Gvw] = tfestimate(inp(ind_eval), out(ind_eval), window, Noverlap, Nest, 1/Ts_log);
    
    % ====== P, Gyu: u -> y ======
    P = T ./ Guw;
    
    % ====== Calculated controller frequency response ======
    Cpi = Gvw ./ (1 - T);
    Cd  = Guw .* Gvw ./ T .* (1 ./ Guw - 1 ./ Gvw);
    
    % ====== Frequency vector for bode plots ======
    omega_bode = 2*pi*P.Frequency;
    
    % ====== Analytical controller transfer functions ======
    [Cpi_ana, Cd_ana, Gf_ana, PID, para_used] = ...
        calculate_transfer_functions(para, ind_ax, throttle_avg, Ts_cntr);
    
    % Downsample analytical functions if needed
    if Gf_ana.Ts < Ts_log
        Gf_ana  = downsample_frd(Gf_ana,  Ts_log, P.Frequency);
        Cpi_ana = downsample_frd(Cpi_ana, Ts_log, P.Frequency);
        Cd_ana  = downsample_frd(Cd_ana,  Ts_log, P.Frequency);
    end
    
    % ====== Package all results into output structure ======
    results.T           = T;
    results.C_T         = C_T;
    results.Guw         = Guw;
    results.C_Guw       = C_Guw;
    results.Gvw         = Gvw;
    results.C_Gvw       = C_Gvw;
    results.P           = P;
    results.Cpi         = Cpi;
    results.Cd          = Cd;
    results.Cpi_ana     = Cpi_ana;
    results.Cd_ana      = Cd_ana;
    results.Gf_ana      = Gf_ana;
    results.omega_bode  = omega_bode;
    results.PID         = PID;
    results.para_used   = para_used;
end