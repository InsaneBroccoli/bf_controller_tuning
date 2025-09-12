function [Pavg, freq, y_axis] = estimate_spectrogram(inp, y, window, Noverlap, Nest, Nres, Ts)
%ESTIMATE_SPECTROGRAM  Welch single-sided spectrogram (power) with amplitude-calibrated scaling
%   [Pavg, freq, y_axis] = estimate_spectrogram(inp, y, window, Noverlap, Nest, Nres, Ts)
%
% PURPOSE
%   - Compute a single-sided Welch-averaged spectrogram for a real signal over a secondary coordinate y
%   - Calibrated so that interior bins are doubled and DC/Nyquist are corrected like in estimate_spectra
%
% INPUTS
%   - inp      [Ndata x 1] real input signal
%   - y        [Ndata x 1] secondary coordinate per sample, binned into Nres rows
%   - window   [Nest x 1] analysis window, e.g. hann(Nest,'periodic')
%   - Noverlap integer overlap in samples, 0..Nest-1
%   - Nest     even integer segment/FFT length, required even in this implementation
%   - Nres     integer number of y bins
%   - Ts       sample time in seconds, Fs = 1/Ts
%
% OUTPUTS
%   - Pavg   [Nres x Npos] single-sided power spectrogram averaged within each y bin
%   - freq   [Npos x 1] frequency in Hz from 0 to Fs/2 inclusive, Npos = Nest/2 + 1
%   - y_axis [Nres x 1] linearly spaced bin centers from min(y) to max(y)
%   - Use sqrt(Pavg) to obtain single-sided amplitude per y bin
%
% SCALING
%   - Internal normalization divides by sum(window)/2 which pre-doubles all bins
%   - DC and Nyquist must not be doubled, so those two bins are divided by 4 in power
%   - Interior bins remain doubled, yielding correct single-sided amplitude after sqrt
%
% METHOD
%   1) Define Nres linear bins across y
%   2) Segment inp with length Nest and overlap Noverlap, apply window
%   3) FFT each segment, normalize by sum(window)/2
%   4) Form one-sided power, fix DC and Nyquist by dividing those bins by 4
%   5) Accumulate segment spectra into y bins and average to get Pavg
%
% NOTES
%   - Assumes even Nest and real inp
%   - 50–90% overlap and a tapered window give smoother estimates
%   - Per-segment mean removal suppresses drift but attenuates true DC
%   - Units carry through, e.g. sqrt(Pavg) has the same amplitude units as inp

    % Assumptions for this minimal fix
    assert(mod(Nest,2)==0, 'This implementation assumes even Nest.');
    assert(numel(window)==Nest, 'window length must equal Nest.');

    % % Global mean removal
    % inp = inp - mean(inp);
    % y   = y   - mean(y);

    Ndata = size(inp, 1);

    % y-axis bins (linear)
    y_min = min(y);
    y_max = max(y);
    dy    = (y_max - y_min) / Nres;
    y_axis = (y_min:dy:y_max-dy).';

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

    Pavg = zeros(Nres, Nfreq);

    Navg = zeros(Nres, 1);

    ind_start = 1;
    ind_end   = Nest;
    Ndelta    = Nest - Noverlap;

    while ind_end <= Ndata

        % Extract segment
        seg_idx = ind_start:ind_end;
        inp_act = inp(seg_idx);

        % % Optional per-segment mean removal (comment out to preserve DC)
        % inp_act = inp_act - mean(inp_act);

        inp_act = window .* inp_act;

        U  = fft(inp_act) / (Nest * W);
        Pact = U .* conj(U); % two-sided power

        % Take one-sided and fix DC & Nyquist (power /4)
        Pseg = Pact(1:Nfreq).';
        Pseg(1)   = Pseg(1)   / 4; % DC
        Pseg(end) = Pseg(end) / 4; % Nyquist (exists since Nest is even)

        % Map y values in this segment to spectrogram rows
        y_seg = y(seg_idx);
        % linear bin index in [1..Nres], with safety clamp
        ind_y = round((y_seg - y_min) / max(y_max - y_min, eps) * (Nres - 1)) + 1;
        ind_y = max(1, min(Nres, ind_y));
        ind_y = sort(ind_y);

        % Count occurrences per unique row (kept close to your original approach)
        ind_y_count = zeros(size(ind_y));
        j = 1;
        for k = 1:length(ind_y)
            if ind_y(k) ~= ind_y(j)
                j = j + 1;
                ind_y(j) = ind_y(k);
            end
            ind_y_count(j) = ind_y_count(j) + 1;
        end

        % Accumulate weighted spectra into rows
        for r = 1:j
            Pavg(ind_y(r), :) = Pavg(ind_y(r), :) + ind_y_count(r) * Pseg;
            Navg(ind_y(r))    = Navg(ind_y(r))    + ind_y_count(r);
        end

        % Next segment
        ind_start = ind_start + Ndelta;
        ind_end   = ind_end   + Ndelta;
    end

    % Row-wise average
    ind = Navg ~= 0;
    Pavg(ind, :) = Pavg(ind, :) ./ Navg(ind);

    Pavg = smooth2d(Pavg);

end

function PavgSmoothed = smooth2d(Pavg)
%SMOOTH2D  Internal 2D smoothing helper
%   Inputs
%     - Pavg [Ny x Nf] matrix to smooth
%   Outputs
%     - PavgSmoothed [Ny x Nf] smoothed matrix using a weighted 3x3 kernel

    kernel = [1 3 1; ...
              3 5 3; ...
              1 3 1];
    kernel = kernel / sum(kernel(:));

    num = conv2(Pavg, kernel, 'same');
    den = conv2(ones(size(Pavg)), kernel, 'same'); % normalize at edges
    PavgSmoothed = num ./ den;

end
