function fcut = get_fcut_from_D_and_fcenter(D, fcent)
%GET_FCUT_FROM_D_AND_FCENTER  Convert damping ratio and center frequency to cutoff frequency
%   fcut = get_fcut_from_D_and_fcenter(D, fcent)
%
% PURPOSE
%   - Compute the effective cutoff frequency fcut for a second-order system
%     given damping ratio D and center frequency fcent
%
% INPUTS
%   - D     scalar or vector damping ratio(s), dimensionless
%   - fcent scalar or vector center frequency/frequencies in Hz
%
% OUTPUTS
%   - fcut  scalar or vector cutoff frequency/frequencies in Hz
%
% METHOD
%   1) Convert damping ratio D to quality factor Q = 1 / (2D)
%   2) Compute cutoff using analytical formula:
%        fcut = fcent / (2Q) * ( -1 + sqrt(1 + 4Q^2) )
%
% NOTES
%   - Supports scalar or element-wise vector computation
%   - Formula derived from standard second-order low-pass response

    Q = 1/2./D;
    fcut = fcent/2./Q.*(-1 + sqrt(1 + 4*Q.^2));

end
