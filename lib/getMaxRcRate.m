function angleRate = getMaxRcRate(axis, rates_type, rcExpoAxis, rcRatesAxis, ratesAxis)

%GETMAXRCRATE  Return the maximum RC rate (deg/s) for a given axis and rates type
%   maxRcRate = getMaxRcRate(axis, rates_type)
%
% PURPOSE
%   - Provide the peak achievable angular rate commanded by the RC mapping
%     for the selected axis (roll/pitch/yaw), depending on the Betaflight
%     "rates_type" (e.g., Betaflight / Raceflight / KISS / Actual / Quick).
%   - Used to scale/normalize setpoints or to convert normalized stick inputs
%     into physical rate limits for controller analysis.
%
% INPUTS
%   - axis        Axis selector (e.g., 0=Roll, 1=Pitch, 2=Yaw) or string
%                ('roll','pitch','yaw') depending on your project convention.
%   - rates_type  Betaflight rates type as numeric code (from log/CLI):
%                0=BETAFLIGHT, 1=RACEFLIGHT, 2=KISS, 3=ACTUAL, 4=QUICK
%
% OUTPUTS
%   - maxRcRate   Scalar maximum commanded angular rate for the given axis
%                in degrees per second [deg/s].
%
% METHOD
%   1) Decode/validate the axis selection (roll/pitch/yaw)
%   2) Decode the rates_type and select the corresponding rate model
%   3) Compute the maximum RC rate using the model-specific mapping parameters
%      (e.g., rcRate / superRate / expo or Actual: centerSensitivity/maxRate/expo)
%   4) Return the resulting maxRcRate (deg/s) for normalization/scaling
%
% NOTES
%   - The function assumes Betaflight’s standard rates_type enumeration:
%     0=Betaflight, 1=Raceflight, 2=KISS, 3=Actual, 4=Quick.
%   - If your analysis is performed in SI units, convert with:
%       maxRcRate_rad = deg2rad(maxRcRate);
%   - The required rate parameters must be available in your data/config
%     (e.g., from log header/CLI dump or stored parameter struct).


rcCommandf = 1;
rcCommandAbs = 1;

switch rates_type
    case 0  % Type BETAFLIGHT
  
    case 1  %Type RACEFLIGHT
        name = "RACEFLIGHT";
        
    case 2  %Type KISS
        name = "KISS";
        
    case 3  %Type ACTUAL
        expoParam     = rcExpoAxis(axis)/100;
        expComponent  = rcCommandAbs * (rcCommandf^5 * expoParam + rcCommandf * (1 - expoParam));

        centerSensitivity = rcRatesAxis(axis)*10;
        stickMovement     = max(0, ratesAxis(axis) * 10 - centerSensitivity);
        
        angleRate = rcCommandf * centerSensitivity + stickMovement * expComponent;  % deg/s
        
    case 4  %Type QUICK
        name = "QUICK";
        
    otherwise

end
