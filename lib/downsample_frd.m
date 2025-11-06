function Gds = downsample_frd(G, Ts, freq)
%DOWNSAMPLE_FRD  Re-sample an FRD onto a target frequency grid (Hz)
%   Gds = downsample_frd(G, Ts, freq)
%
% PURPOSE
%   - Evaluate an FRD object on a new frequency vector in Hz
%   - Preserve DC by guarding the 0-Hz point when the original response is singular or missing
%
% INPUTS
%   - G     frd object with existing frequency response
%   - Ts    sample time in seconds for the returned frd (use [] for continuous-time semantics)
%   - freq  [N x 1] target frequency vector in Hz, monotonic increasing and including 0
%
% OUTPUTS
%   - Gds   frd object evaluated on freq (Hz), with size matching G’s I/O
%
% METHOD
%   1) Evaluate G at 2*pi*freq(2:end) using freqresp (skip DC)
%   2) Prepend the first evaluated point to stand in for DC if needed
%   3) Construct a new frd with Units 'Hz' at sampling time Ts
%
% NOTES
%   - If G has a well-defined DC value you may wish to insert it explicitly instead of copying the first nonzero bin
%   - freq must start at 0 for the DC placeholder logic to make sense

    % The frequency response might be infinite as zero frequency (integrating)
    g = squeeze(freqresp(G, 2*pi*freq(2:end)));
    g = [g(1); g];

    Gds = frd(g, freq, Ts, 'Units', 'Hz');

end
