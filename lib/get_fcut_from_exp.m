function fcut = get_fcut_from_exp(dynLpfMin, dynLpfMax, expo, throttle)
%GET_FCUT_FROM_EXP  Dynamic low-pass cutoff frequency from exponential curve
%   fcut = get_fcut_from_exp(dynLpfMin, dynLpfMax, expo, throttle)
%
% PURPOSE
%   - Compute a throttle-dependent cutoff frequency for a dynamic low-pass filter
%   - Uses an exponential curve shaping to transition smoothly between min and max cutoff
%
% INPUTS
%   - dynLpfMin scalar minimum cutoff frequency [Hz] at throttle = 0
%   - dynLpfMax scalar maximum cutoff frequency [Hz] at throttle = 1
%   - expo      scalar curve shaping parameter, larger values push transition toward mid-throttle
%   - throttle  scalar or vector in [0..1], normalized throttle input
%
% OUTPUTS
%   - fcut scalar or vector cutoff frequency [Hz] corresponding to throttle
%
% METHOD
%   1) Scale expo by 0.1 (expo/10) to form expof
%   2) Apply shaping curve = throttle*(1-throttle)*expof + throttle
%   3) Interpolate linearly between dynLpfMin and dynLpfMax using curve
%
% NOTES
%   - throttle must be normalized in [0..1]
%   - expo = 0 results in a purely linear mapping between min and max
%   - Larger expo values create a stronger dip around mid-throttle

    expof = expo / 10.0;
    curve = throttle .* (1 - throttle) * expof + throttle;
    fcut = (dynLpfMax - dynLpfMin) * curve + dynLpfMin;

end
