clc, clear variables

%% Sampling Time 
Ts      = 125 * 1.0e-6;        % Gyro loop
Ts_cntr = 2 * Ts;              % Control loop
Ts_log  = 2 * Ts_cntr;         % Logging loop (not used here, just info)
fs      = 1/Ts;

t = (0:Ts:10).';               % Time vector from 0s to 10s (column vector)
N = length(t);                 % number of samples

%% Input Signal
omega = 10.1 * 2*pi();
u = 0.6*sin(2*pi*(10 + 0.5*t).*t) + 0.25*randn(size(t));            %Input signal
w = hann(N, 'periodic');          % window over the entry singal
u_win = u .* w;                   % signal with window

figure(1)
subplot(211)
plot(t,u);grid on;
title('Input Signal without Window');

subplot(212)
plot(t,u_win);grid on;
title('Input Sigal with Window');

%% Spectrum of Input Signal

U = fft(u);
U_WIN = fft(u_win);
w_fac = sum(w) / 2;  

% Frequenzachse (Hz)
f = (0:N-1)*(fs/N);

% Calculation of spectra
A_noWin = abs(U(1:N/2)) / (N/2);    % New vector without Niquist
A_win   = abs(U_WIN(1:N/2)) / w_fac;   %New vector with upscaling because of win
f_half  = f(1:N/2);                 % New vector without Niquist

% Plot
figure(2)
subplot(2,1,1)
plot(f_half, A_noWin)
xlabel('Frequenz [Hz]')
ylabel('Amplitude')
title('Spektrum ohne Fenster')
grid on
xlim([0 100])

subplot(2,1,2)
plot(f_half, A_win)
xlabel('Frequenz [Hz]')
ylabel('Amplitude')
title('Spektrum mit Hann-Fenster')
grid on
xlim([0 100])

%% Split up in Segments

seg = 5;            %Number of segments
overlap = 0.5;      %Overlap in next segement x*100%
Nest = floor(N / (1+(seg-1)*(1-overlap)));  %Number of points in segment
window = hann(Nest, 'periodic');        % Window over signal
Nshift = round(Nest * (1 - overlap));   % Number of points to start next segment
freq = (0:Nest/2) * (fs / Nest);        

u_seg  = zeros(Nest, seg);      %Preparation array
U_SEG  = zeros(Nest, seg);
U_SEG_WIN  = zeros(Nest, seg);

for i = 1:seg
    start = (i-1)*Nshift + 1;  
    stop  = start + Nest - 1;  
    if stop > N
        break;   % Stop when segments are full
    end
    u_seg(:, i) = u(start:stop);    % Add Segment to array
    u_seg_win(:,i) = u(start:stop) .* window;   % Add segement to array with window

    U_SEG(:,i) = abs(fft(u_seg(:,i)));      % Absolut value of frequent component
    U_SEG_WIN(:,i) = abs(fft(u_seg_win(:,i)));
end

U_Com_seg = mean(U_SEG, seg);           % Meanvalue of frequent component over all segments
U_Com_win = mean(U_SEG_WIN, seg);

K = floor(Nest/2);                  % one sided spectrum without niquist
freq = (0:K-1).' * (fs/Nest);       % Kx1, frequency axis as column vector

% Scaling like full signal
A_seg = U_Com_seg(1:K)     / (Nest/2);
w_fac = sum(window) / 2;    %upscaling factor of hann window
A_win = U_Com_win(1:K) / w_fac;

figure(3)
subplot(211)
plot(t(1:Nest), u_seg(:,1));grid on;
title('Segment 1');
subplot(212)
plot(t(Nshift:(Nshift+Nest-1)),u_seg(:,2));grid on;
title('Segment 2');

figure(4)
subplot(2,1,1)
plot(freq, A_seg); grid on
xlabel('Frequenz [Hz]'); ylabel('Amplitude');
title('Segmentiertes Spektrum ohne Fenster'); xlim([0 100])

subplot(2,1,2)
plot(freq, A_win); grid on
xlabel('Frequenz [Hz]'); ylabel('Amplitude');
title('Segmentiertes Spektrum mit Hann-Fenster'); xlim([0 100])
