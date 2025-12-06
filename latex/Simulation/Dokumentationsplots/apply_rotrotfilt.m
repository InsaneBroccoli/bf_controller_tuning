%% --- Daten ---
t  = (0:1e-3:0.4).';           % Ts = 0.01 s  -> Fs = 100 Hz
Fs = 1/(t(2)-t(1));            % = 100 Hz
f  = 10;                       % Signal-Frequenz [Hz]

u     = sin(2*pi*f*t);         % Originales Sinussignal
A_l   = 0.4;                   % Rauschamplitude
l     = A_l*randn(size(t));    % Rauschsignal
utild = u + l;                 % verrauschtes Signal

%% --- Zeitbereich-Plot ---
figure(1); clf
subplot(1,2,1)
plot(t, u, 'LineWidth', 1.2); grid on
xlabel('Time [s]'); ylabel('Amplitude'); title('Original Signal')

subplot(1,2,2)
plot(t, utild, 'LineWidth', 1.2); grid on
xlabel('Time [s]'); ylabel('Amplitude'); title('Noisy Signal')

%% --- Rotation ins Basisband ---
p  = exp(1i*2*pi*f*t);         % zeitabhängiger Phasor e^{j2πft}
y  = utild - mean(utild);      % DC-Anteil entfernen
yR = y .* p;                   % Forward rotation  (e^{+j2πft})
yQ = y .* conj(p);             % Backward rotation (e^{-j2πft})

%% --- FFTs für Figure 2 (zweiseitig) ---
N = length(utild);
f_axis_full = (-N/2:N/2-1)*(Fs/N);

U_full  = fftshift(fft(utild));
YQ_full = fftshift(fft(yQ));
YR_full = fftshift(fft(yR));

U_mag_2s  = abs(U_full / N);
YQ_mag_2s = abs(YQ_full / N);
YR_mag_2s = abs(YR_full / N);

step  = 3;                                       % nur jeden 3. Punkt zeigen
ymax2 = 1.1*max([U_mag_2s; YQ_mag_2s; YR_mag_2s]);

%% --- FIGURE 2: Spektren (ungefiltert, zweiseitig) ---
figure(2); clf
subplot(1,3,1)
stem(f_axis_full(1:step:end), U_mag_2s(1:step:end), 'filled'); grid on
xlabel('Frequency [Hz]'); ylabel('Amplitude')
title('Spectrum (two-sided)')
xlim([-Fs/2, Fs/2]); ylim([0 ymax2])

subplot(1,3,2)
stem(f_axis_full(1:step:end), YQ_mag_2s(1:step:end), 'filled'); grid on
xlabel('Frequency [Hz]'); ylabel('Amplitude')
title('Backward Rotation')
xlim([-Fs/2, Fs/2]); ylim([0 ymax2])

subplot(1,3,3)
stem(f_axis_full(1:step:end), YR_mag_2s(1:step:end), 'filled'); grid on
xlabel('Frequency [Hz]'); ylabel('Amplitude')
title('Forward Rotation')
xlim([-Fs/2, Fs/2]); ylim([0 ymax2])

%% === Low-Pass im Basisband (Butterworth + filtfilt) für Figure 3 ===
fc = 2;                          % Grenzfrequenz [Hz] (anpassbar)
Wn = fc / (Fs/2);
[b,a] = butter(4, Wn);           % 4. Ordnung
yR_f = filtfilt(b, a, yR);       % forward-rot. gefiltert
yQ_f = filtfilt(b, a, yQ);       % backward-rot. gefiltert

% Zweiseitige Spektren (gleiche Darstellung wie Figure 2)
YQf_full    = fftshift(fft(yQ_f));
YRf_full    = fftshift(fft(yR_f));
YQf_mag_2s  = abs(YQf_full / N);
YRf_mag_2s  = abs(YRf_full / N);

ymax3 = 1.1*max([U_mag_2s; YQf_mag_2s; YRf_mag_2s]);

%% --- FIGURE 3: gleicher Aufbau wie Figure 2, aber gefiltert ---
figure(3); clf
subplot(1,3,1)
stem(f_axis_full(1:step:end), U_mag_2s(1:step:end), 'filled'); grid on
xlabel('Frequency [Hz]'); ylabel('Amplitude')
title('Spectrum (two-sided)')
xlim([-Fs/2, Fs/2]); ylim([0 ymax3])

subplot(1,3,2)
stem(f_axis_full(1:step:end), YQf_mag_2s(1:step:end), 'filled'); grid on
xlabel('Frequency [Hz]'); ylabel('Amplitude')
title(sprintf('Backward Rotation + LP)', fc))
xlim([-Fs/2, Fs/2]); ylim([0 ymax3])

subplot(1,3,3)
stem(f_axis_full(1:step:end), YRf_mag_2s(1:step:end), 'filled'); grid on
xlabel('Frequency [Hz]'); ylabel('Amplitude')
title(sprintf('Forward Rotation + LP', fc))
xlim([-Fs/2, Fs/2]); ylim([0 ymax3])

%% --- FIGURE 4: Rekonstruiertes Zeit-Signal nach LP ---
u_rec = real(0.5*(yR_f .* conj(p) + yQ_f .* p));

figure(4); clf
plot(t, u, 'k', 'LineWidth', 1.2); hold on
plot(t, utild, 'Color', [0.7 0.7 0.7]);
plot(t, u_rec, 'r', 'LineWidth', 1.2);
grid on
legend('Original','Noisy','Reconstructed (LP)','Location','best')
xlabel('Time [s]'); ylabel('Amplitude')
title(sprintf('Time Domain Reconstruction after Baseband LP (f_c = %.1f Hz)', fc))
