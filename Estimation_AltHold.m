clc; close all;
s = tf('s');

%% Drone parameters
m   = 0.255;        % kg
g   = 9.81;         % m/s^2
rho = 1.225;        % kg/m^3
CT  = 0.10;         % thrust coefficient
D   = 0.0889;       % m   (3.5 inch propeller)
% D = 0.0762;       % m   (3 inch propeller)

%% Motor parameters
R   = 0.08;         % Ohm
L   = 20e-6;        % H
Jm  = 3.0e-6;       % kg*m^2
k   = 0.0032;       % torque / back-EMF constant

%% Offsets

hoverOffset = 1275;
PWM_min = 1000;
PWM_max = 2000;

V_Bat = 16;

%% Hover values
T0_total = m*g;         % total thrust at hover
T0_motor = T0_total/4;  % thrust per motor at hover

% Hover speed per motor
n0 = sqrt(T0_motor/(CT*rho*D^4));   % rev/s
rpm0 = n0*60;                       % rpm

%% Plant
% Motor: U -> omega [rad/s]
P_Motor = k / (Jm*L*s^2 + Jm*R*s + k^2);

% Thrust/speed linearization for all 4 motors together
P_Alt_vel = (8*CT*rho*D^4*n0)/(m*s);

% Velocity -> position
P_Alt_pos = 1/s;

Td = 0.01;   % 100 Hz -> 10 ms

P_delay = exp(-Td*s);

%% PI Controllers

KI = 30;

KP_inner = 15*0.005;
KI_inner = KI*0.00012;
KP_outer = 15*0.04;
KI_outer = KI*0.0008;

C_inner = KP_inner + KI_inner/s;
C_outer = KP_outer + KI_outer/s;

%% Inner loop Estimation

L_inner = C_inner * P_Motor * P_Alt_vel * P_delay;

T_inner = L_inner / (1+ L_inner);

S_inner = 1 / (1 + L_inner);

figure(1)
margin(L_inner);grid on;
title('Open Loop Inner')

figure(2)
bode(S_inner);
title('Sensitivity Inner Loop');


%% Outer loop Estimation

L_outer = C_outer * T_inner * P_Alt_pos;

T_outer = L_outer / (1+L_outer);

S_outer = 1 / (1 + L_outer);


figure(3)
margin(L_outer);grid on;
title('Open Loop Outer')

figure(4)
bode(S_outer);
title('Sensitivity Outer Loop');


figure(5)
step(T_outer);