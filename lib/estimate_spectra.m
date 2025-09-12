function [Pavg, freq] = estimate_spectra(inp, window, Noverlap, Nest, Ts)
%ESTIMATE_SPECTRA  Welch single-sided spectrum (power) with amplitude-calibrated scaling
%   [Pavg, freq] = estimate_spectra(inp, window, Noverlap, Nest, Ts)
%
% PURPOSE
%   - Compute a single-sided Welch-averaged spectrum for real signals
%   - Calibrated so that a sine of amplitude A at a bin appears as ~A in sqrt(Pavg)
%
% INPUTS
%   - inp      [Ndata x Nsignals] real data, each column one signal
%   - window   [Nest x 1] analysis window, e.g. hann(Nest,'periodic')
%   - Noverlap integer overlap in samples, 0..Nest-1
%   - Nest     even integer segment/FFT length, required even in this implementation
%   - Ts       sample time in seconds, Fs = 1/Ts
%
% OUTPUTS
%   - Pavg   [Npos x Nsignals] single-sided power spectrum averaged over segments
%   - freq   [Npos x 1] frequency in Hz from 0 to Fs/2 inclusive, Npos = Nest/2 + 1
%   - Use spectra = sqrt(Pavg) to obtain single-sided amplitude
%
% SCALING
%   - Internal normalization divides by sum(window)/2 which pre-doubles all bins
%   - DC and Nyquist must not be doubled, so those two bins are divided by 4 in power
%   - Interior bins remain doubled, yielding correct single-sided amplitude after sqrt
%
% METHOD
%   1) Remove global mean per column
%   2) Segment with length Nest and overlap Noverlap, apply window
%   3) FFT each segment, normalize by sum(window)/2
%   4) Form one-sided power, fix DC and Nyquist by dividing those bins by 4
%   5) Average over segments to get Pavg
%
% NOTES
%   - Assumes even Nest and real inputs
%   - 50–90% overlap and a tapered window give smoother estimates
%   - Per-segment mean removal suppresses drift but attenuates true DC
%   - Units carry through, e.g. if inp is deg/s then sqrt(Pavg) is deg/s
%
% PSD ALTERNATIVE
%   - For PSD, use RMS window normalization sum(window.^2), average power, account for bin width

    % Assumptions for this minimal fix
    assert(mod(Nest,2)==0, 'This implementation assumes even Nest.');
    assert(numel(window)==Nest, 'window length must equal Nest.');

    % Global mean removal
    inp = inp - mean(inp);

    [Ndata, Nsignals] = size(inp);

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

    Pavg = zeros(Nfreq, Nsignals);

    for i = 1:Nsignals
        Navg = 0;

        ind_start = 1;
        ind_end   = Nest;
        Ndelta    = Nest - Noverlap;

        while ind_end <= Ndata
            inp_act = inp(ind_start:ind_end, i);

            % Optional per-segment mean removal (comment out to preserve DC)
            inp_act = inp_act - mean(inp_act);

            inp_act = window .* inp_act;

            U    = fft(inp_act) / (Nest * W);
            Pact = U .* conj(U); % two-sided power

            % Take one-sided and fix DC & Nyquist (power /4)
            Pseg = Pact(1:Nfreq);
            Pseg(1)   = Pseg(1)   / 4; % DC
            Pseg(end) = Pseg(end) / 4; % Nyquist (exists since Nest is even)

            Pavg(:, i) = Pavg(:, i) + Pseg;
            Navg = Navg + 1;

            % Next segment
            ind_start = ind_start + Ndelta;
            ind_end   = ind_end   + Ndelta;
        end

        if Navg > 0
            Pavg(:, i) = Pavg(:, i) / Navg;
        end
    end

end
