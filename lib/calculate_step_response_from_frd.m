function step_resp = calculate_step_response_from_frd(G, f_max)
%CALCULATE_STEP_RESPONSE_FROM_FRD  Step response from frequency response data
%   step_resp = calculate_step_response_from_frd(G, f_max)
%
% PURPOSE
%   - Compute a time-domain step response approximation from an FRD object
%   - Uses inverse FFT of truncated, symmetrized frequency response data
%
% INPUTS
%   - G      frd object containing complex frequency response samples
%   - f_max  scalar cutoff frequency [Hz]; frequencies above are zeroed
%
% OUTPUTS
%   - step_resp [N x 1] real-valued step response approximation in time domain
%
% METHOD
%   1) Extract frequency response data g = G.ResponseData
%   2) Replace missing DC (NaN) with the second frequency point if needed
%   3) Zero out response above f_max
%   4) Construct a full symmetric spectrum [g; conj(flip(g))]
%   5) Inverse FFT to impulse response, then cumulative sum to step response
%
% NOTES
%   - Assumes G.Frequency is equally spaced from DC up to Nyquist
%   - f_max truncation removes high-frequency contributions (low-pass effect)
%   - Step response is approximate, resolution depends on frequency grid

    g = squeeze(G.ResponseData);
    if isnan(abs(g(1))) % TODO: interpolate based on point 2 and 3
        g(1) = g(2);
    end

    % Zero out above f_max_hz
    freq = G.Frequency;
    g(freq > f_max) = 0;

    % Construct full symmetric spectrum
    g_full = [g; conj(g(end-1:-1:2))];

    % Step response is cumulative sum of real part of IFFT
    step_resp = cumsum(real(ifft(g_full)));

end
