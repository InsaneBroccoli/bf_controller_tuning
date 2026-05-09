function [Cpid] = calculate_poshold_controller(P, I, D, A, fc_pt1, Ts_cntr, Ts_log, freq)
%CALCULATE_POSHOLD_CONTROLLER  Discrete-time analytical position hold PIDA controller.
%
% INPUTS
%   P, I, D, A  Betaflight CLI integer gains (e.g. 30, 30, 30, 30)
%   fc_pt1      velocity/acceleration PT1 filter cutoff [Hz] (positionCutoff * 0.01)
%   Ts_cntr     poshold control loop sample time [s]
%   Ts_log      sample time of the returned FRDs [s]
%   freq        target frequency vector [Hz]
%
%
% Scales AP_*_SCALE are firmware constants from
% src/main/flight/pos_hold_multirotor.c / autopilot_multirotor.c.

    AP_P_SCALE = 0.0012;
    AP_I_SCALE = 0.0001;
    AP_D_SCALE = 0.0015;
    AP_A_SCALE = 0.0008;

    Kp = P * AP_P_SCALE;
    Ki = I * AP_I_SCALE;
    Kd = D * AP_D_SCALE;
    Ka = A * AP_A_SCALE;

    z = tf('z');

    % PT1 filter for velocity and acceleration (vaLpfCutoff in firmware)
    Gpt1 = get_filter('pt1', fc_pt1, Ts_cntr);

    % Discrete differentiator: (z-1)/(Ts*z)
    diff_z = (z - 1) / (Ts_cntr * z);

    % PIDA controller:
    % P: proportional on error
    % I: forward Euler integrator
    % D: 1x differentiator + 1x PT1
    % A: 2x differentiator + 2x PT1 (acceleration)
    Cpid_tf = Kp ...
            + Ki * Ts_cntr * z / (z - 1) ...
            + Kd * diff_z * Gpt1 ...
            + Ka * diff_z^2 * Gpt1^2;

    % Complete controller with filter
    Cpid_full = Cpid_tf;
    Cpid = downsample_frd(Cpid_full, Ts_log, freq);
end