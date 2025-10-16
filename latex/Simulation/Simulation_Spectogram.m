clc; clear variables

%% Sampling Time 
Ts      = 125 * 1.0e-6;        % Gyro loop
Ts_cntr = 2 * Ts;              % Control loop
Ts_log  = 2 * Ts_cntr;         % Logging loop (not used here, just info)
fs      = 1/Ts;

t = (0:Ts:10).';               % Time vector from 0s to 10s (column vector)
N = length(t);                 % number of samples

%% Thrust
w_thr = 5 * 2*pi;
thrust = 2 + 2*sin(w_thr*t);

%% Input Signal
omega = 10.1 * 2*pi();
u = 0.6*sin(2*pi*(10 + 1*t).*t) + 0.25*randn(size(t));    %Input signal
w = hann(N, 'periodic');          % window over the entry singal
u_win = u .* w;                   % signal with window

subplot(311)
plot(t,u);grid on;
xlim([0 2]);
xlabel('Time [s]');
title('Input Signal');

subplot(312)
plot(t,u_win);grid on;
xlim([0 2]);
xlabel('Time [s]');
title('Input Signal with Window');

subplot(313)
plot(t, thrust);grid on;
xlim([0 2]);
xlabel('Time [s]');
title('Thrust');

%% Spectrogram without Segments

U = fft(u);
