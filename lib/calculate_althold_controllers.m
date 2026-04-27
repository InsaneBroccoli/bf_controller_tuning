function [Cpi, Cd] = calculate_althold_controllers(P, I, D, fc_pt2, Ts_cntr, Ts_log, freq)
%CALCULATE_ALTHOLD_CONTROLLERS  Discrete-time analytical althold PI + D controllers.
%
% INPUTS
%   P, I, D    Betaflight CLI integer gains (e.g. 20, 20, 18)
%   fc_pt2     D-term PT2 cutoff [Hz]
%   Ts_cntr    althold control loop sample time [s] (typically 0.01)
%   Ts_log     sample time of the returned FRDs [s]
%   freq       target frequency vector [Hz] (e.g. P.Frequency)
%
% OUTPUTS
%   Cpi        FRD: PI controller, Kp + Ki*Ts*z/(z-1)
%   Cd         FRD: (Kd/Ts)*(z-1)/z cascaded with a PT2 filter at fc_pt2
%
% Scales AP_*_SCALE are firmware constants from
% src/main/flight/alt_hold_multirotor.c / autopilot_multirotor.c.

    AP_P_SCALE = 0.01;
    AP_I_SCALE = 0.003;
    AP_D_SCALE = 0.01;

    Kp = P * AP_P_SCALE;
    Ki = I * AP_I_SCALE;
    Kd = D * AP_D_SCALE;

    Cpi = ss(Kp + Ki*Ts_cntr*tf([1 0],[1 -1], Ts_cntr));
    Cpi = downsample_frd(Cpi, Ts_log, freq);

    Cd_only = ss( Kd/Ts_cntr*tf([1 -1], [1 0], Ts_cntr) );
    Gf_pt2  = get_filter('pt2', fc_pt2, Ts_cntr);
    Cd      = downsample_frd(Cd_only * Gf_pt2, Ts_log, freq);

end
