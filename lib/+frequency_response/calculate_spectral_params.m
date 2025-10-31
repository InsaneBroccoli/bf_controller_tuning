function [Nest, Noverlap, window] = calculate_spectral_params(Ts_log, koverlap)
    % CALCULATE_SPECTRAL_PARAMS - Calculate parameters for spectral estimation
    %
    % INPUTS:
    %   Ts_log    - Sampling time
    %   koverlap  - Overlap factor (default: 0.9)
    %
    % OUTPUTS:
    %   Nest      - Window size in samples
    %   Noverlap  - Overlap in samples
    %   window    - Window function
    
    if nargin < 2
        koverlap = 0.9;  % Default 90% overlap
    end
    
    Nest = round(2 / Ts_log);                   % Window size for 2 seconds
    Noverlap = floor(koverlap * Nest);          % Overlap in samples
    window = hann(Nest, 'periodic');            % Hann window
end