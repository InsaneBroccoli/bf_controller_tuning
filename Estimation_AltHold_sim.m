clc; close all; clear all;
s = tf('s');

%% Drone parameters

CT  = 0.15;         % thrust coefficient
m   = 0.300;          % kg
g   = 981;            % cm/s^2
rho = 1.225e-6;       % kg/cm^3
D   = 8.89;           % cm

b = 1e-5;   % Beispielwert, muss abgestimmt werden

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

%% PI Controllers

KI = 15;

KP_inner = 15*0.005;
KI_inner = KI*0.00012;
KP_outer = 15*0.04;
KI_outer = KI*0.0008;

C_inner = KP_inner + KI_inner/s;
C_outer = KP_outer + KI_outer/s;

%% Plant
% Motor: U -> omega [rad/s]
P_Motor = k / (Jm*L*s^2 + Jm*R*s + k^2);

Td = 0.005;   % 5 ms als effektives Delay
P_delay = exp(-Td*s);