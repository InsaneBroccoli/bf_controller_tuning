function Glp = create_lowpass_filter(Ts_log)
    % CREATE_LOWPASS_FILTER - Create discrete lowpass filter
    %
    % INPUTS:
    %   Ts_log - Sampling time
    %
    % OUTPUTS:
    %   Glp - Discrete transfer function
    
    Dlp = sqrt(3)/2;
    wlp = 2 * pi * 10;  % Cutoff frequency: 10 Hz
    
    % Create continuous-time filter
    Glp_cont = tf(wlp^2, [1, 2*Dlp*wlp, wlp^2]);
    
    % Convert to discrete time using Tustin
    Glp = c2d(Glp_cont, Ts_log, 'tustin');
end