function [G, C, freq, Pavg] = estimate_frequency_response(inp, out, window, Noverlap, Nest, Ts, delta)
%ESTIMATE_FREQUENCY_RESPONSE  Welch single-sided X/Y spectra with amplitude-calibrated scaling
%   [G, C, freq, Pavg] = estimate_frequency_response(inp, out, window, Noverlap, Nest, Ts, delta)
%
% PURPOSE
%   - Estimate single-sided Welch-averaged spectra Suu, Syu, Syy for real input/output
%   - Return frequency response G = Syu ./ Suu and coherence C consistent with single-sided scaling
%
% INPUTS
%   - inp      [Ndata x 1] real input signal
%   - out      [Ndata x 1] real output signal
%   - window   [Nest x 1] analysis window, e.g. hann(Nest,'periodic')
%   - Noverlap integer overlap in samples, 0..Nest-1
%   - Nest     even integer segment/FFT length, required even in this implementation
%   - Ts       sample time in seconds, Fs = 1/Ts
%   - delta    nonnegative regularization added to Suu to stabilize division, default 0
%
% OUTPUTS
%   - G     frd single-input single-output frequency response estimate over freq [Hz]
%   - C     frd magnitude-squared coherence over freq [Hz]
%   - freq  [Npos x 1] frequency in Hz from 0 to Fs/2 inclusive, Npos = Nest/2 + 1
%   - Pavg  [Npos x 3] single-sided power spectra columns [Suu, Syu, Syy] averaged over segments
%
% SCALING
%   - Internal normalization divides by sum(window)/2 which pre-doubles all bins
%   - DC and Nyquist must not be doubled, so those two bins are divided by 4 in power for all spectra
%   - Interior bins remain doubled, giving correct single-sided scaling; common factors cancel in G and C
%
% METHOD
%   1) Remove global mean from inp and out
%   2) Segment with length Nest and overlap Noverlap, apply window
%   3) FFT each segment of inp and out, normalize by sum(window)/2
%   4) Form two-sided Suu, Syu, Syy, convert to one-sided and fix DC & Nyquist by dividing by 4
%   5) Average over segments to get Pavg, compute G and C, return frd with Units 'Hz'
%
% NOTES
%   - Assumes even Nest and real signals
%   - 50–90% overlap and a tapered window give smoother estimates
%   - Per-segment mean removal further suppresses drift but attenuates true DC
%   - delta prevents division by zero in low-excitation bands and does not affect coherence scaling

    if (~exist('delta','var') || isempty(delta))
        delta = 0.0;
    end

    % Assumptions for this minimal fix
    assert(mod(Nest,2)==0, 'This implementation assumes even Nest.');
    assert(numel(window)==Nest, 'window length must equal Nest.');

    % Global mean removal
    inp = inp - mean(inp);
    out = out - mean(out);

    Ndata = size(inp, 1);

    % Frequency axis (0 .. Fs/2), length Nfreq = Nest/2+1
    % df   = 1 / (Nest * Ts);
    % freq = (0:df:1/Ts-df).';
    % ind = freq <= 1 / (2 * Ts);
    % freq = freq(ind);
    fs = 1/Ts;
    freq = (0:Nest/2).' * (fs / Nest);
    Nfreq = length(freq);

    % Normalization (bakes single-sided doubling into all bins)
    % W = sum(window)/Nest/2  => dividing by Nest*W == dividing by (sum(window)/2)
    W = sum(window) / Nest / 2;

    % Columns: [Suu, Syu, Syy]
    Pavg = zeros(Nfreq, 3);

    Navg = 0;

    ind_start = 1;
    ind_end   = Nest;
    Ndelta    = Nest - Noverlap;

    while ind_end <= Ndata
        ind = ind_start:ind_end;

        inp_act = inp(ind);
        out_act = out(ind);

        % Optional per-segment mean removal (comment out to preserve DC)
        inp_act = inp_act - mean(inp_act);
        out_act = out_act - mean(out_act);

        inp_act = window .* inp_act;
        out_act = window .* out_act;

        U = fft(inp_act) / (Nest * W);
        Y = fft(out_act) / (Nest * W);

        % Two-sided spectra
        % Suu = U*conj(U), Syu = Y*conj(U), Syy = Y*conj(Y)
        Pact = [U.*conj(U), Y.*conj(U), Y.*conj(Y)];

        % One-sided conversion with DC & Nyquist fix (power /4)
        Pseg = Pact(1:Nfreq, :);
        Pseg(1,  :) = Pseg(1,  :) / 4; % DC
        Pseg(end,:) = Pseg(end,:) / 4; % Nyquist (exists since Nest is even)

        Pavg = Pavg + Pseg;
        Navg = Navg + 1;

        % Next segment
        ind_start = ind_start + Ndelta;
        ind_end   = ind_end   + Ndelta;
    end

    if Navg > 0
        Pavg = Pavg / Navg;
    end

    [g, c] = calc_freqresp_and_cohere(Pavg, delta);

    % Return frd with frequency in Hz
    G = frd(g, freq, Ts, 'Units', 'Hz');
    C = frd(c, freq, Ts, 'Units', 'Hz');

end

function [g, c] = calc_freqresp_and_cohere(P, delta)
%CALC_FREQRESP_AND_COHERE  Helper for G and coherence
%   Inputs
%     - P    [Npos x 3] columns [Suu, Syu, Syy], single-sided power
%     - delta scalar nonnegative regularization added to Suu
%   Outputs
%     - g    complex frequency response Syu ./ (Suu + delta)
%     - c    magnitude-squared coherence |Syu|^2 ./ (Suu .* Syy)

    P(:,1) = P(:,1) + delta;      % regularize Suu
    g = P(:,2) ./ P(:,1);         % Syu / Suu
    c = abs(P(:,2)).^2 ./ (P(:,1) .* P(:,3));

end
