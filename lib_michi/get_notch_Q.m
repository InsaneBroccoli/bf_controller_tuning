function Q = get_notch_Q(centerFreq, cutoffFreq)
%GET_NOTCH_Q  Compute quality factor Q for a digital notch filter
%   Q = get_notch_Q(centerFreq, cutoffFreq)
%
% PURPOSE
%   - Calculate the quality factor Q of a notch filter from its center
%     frequency and lower cutoff frequency
%
% INPUTS
%   - centerFreq  scalar center frequency f0 in Hz
%   - cutoffFreq  scalar lower cutoff frequency f1 in Hz
%
% OUTPUTS
%   - Q           scalar quality factor, Q = f0 / (f2 - f1) with
%                 f2 = f0^2 / f1
%
% NOTES
%   - cutoffFreq must be less than centerFreq
%   - Useful for designing digital notch filters where bandwidth is
%     defined by cutoff frequencies

    Q = centerFreq * cutoffFreq / (centerFreq * centerFreq - cutoffFreq * cutoffFreq);

end
