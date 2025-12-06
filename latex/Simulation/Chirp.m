%% Chirp-Signal und Phase arg(t)

clear; close all; clc;

% Zeitbasis
Fs = 20000;            % Abtastrate
T  = 10;               % Dauer
t  = 0:1/Fs:T;

% Instantaneous frequency f(t) (linear chirp 50 Hz → 300 Hz)
f0 = 0.01;
f1 = 10;
f  = f0 + (f1 - f0) * t / T;

% Phase argument arg(t) = 2*pi * integral(f(t)*dt)
arg = 2*pi * cumtrapz(t, f);

% Exponentieller Frequenzverlauf f(t)
f_exp = f0 * (f1/f0).^(t/T);

% Phase arg(t) = 2*pi * integral f(t)
arg_exp = 2*pi * cumtrapz(t, f_exp);

% Wrapped phase modulo 2*pi
arg_wrapped = mod(arg, 2*pi);
arg_exp_wrapped = mod(arg_exp, 2*pi);

% Chirp-Signal x(t)
x = sin(arg);

%% Plotten
figure(1);
plot(t, x, 'LineWidth', 1.2);
title('Chirp Signal  x(t) = sin(arg(t))');
xlabel('Zeit [s]'); ylabel('Amplitude');

figure(2)
subplot(1,2,1)
plot(t, arg, 'LineWidth', 1.2);
title('Unwrapped Phase  arg(t)');
xlabel('Zeit [s]'); ylabel('Phase [rad]');

subplot(1,2,2)
plot(t, arg_wrapped, 'LineWidth', 1.2);
title('Wrapped Phase  mod(arg(t), 2\pi)');
xlabel('Zeit [s]'); ylabel('Phase [rad]');
ylim([0 2*pi]);

figure(3)
subplot(1,2,1)
plot(t, arg_exp, 'LineWidth', 1.2);
title('Unwrapped Phase  arg(t)');
xlabel('Zeit [s]'); ylabel('Phase [rad]');

subplot(1,2,2)
plot(t, arg_wrapped, 'LineWidth', 1.2);
title('Wrapped Phase  mod(arg(t), 2\pi)');
xlabel('Zeit [s]'); ylabel('Phase [rad]');
ylim([0 2*pi]);