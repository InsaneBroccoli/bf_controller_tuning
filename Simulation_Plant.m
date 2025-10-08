clc, clear variables

%% Bodeplot Settings
% Preferences
set(cstprefs.tbxprefs, 'MagnitudeUnits', 'dB');
set(cstprefs.tbxprefs, 'FrequencyUnits', 'Hz');
set(cstprefs.tbxprefs, 'UnwrapPhase', 'Off');
set(cstprefs.tbxprefs, 'Grid', 'On');
set(groot,'defaultLineLineWidth',1.2);   % global for all new lines

% Bodeoptions
opt = bodeoptions('cstprefs');
opt.Xlim = { [0.1 1e3] };      % x-axis: 0.1 Hz to 1000 Hz

%% Sampling Time 
Ts      = 125 * 1.0e-6;        % Gyro loop
Ts_cntr = 2 * Ts;              % Control loop
Ts_log  = 2 * Ts_cntr;         % Logging loop (not used here, just info)
fs      = 1/Ts;

t = (0:Ts:10).';               % Time vector from 0s to 10s (column vector)
N = length(t);                 % number of samples

% frequency vector in Hz
f = (0:N-1)*(fs/N);

s = tf('s');
z = tf('z');

%% Chirp Signal (logarithmic)
f0 = 0.1; f1 = 500;            % Hz
u = chirp(t, f0, 10, f1, 'logarithmic', 0);   % numerical vector
U_In = timeseries(u, t);       % optional: for Simulink/Scopes

%% Interference Signals

% Input interference
A_l = 0.5;
omega_l = 2*pi*10;
l = A_l * sawtooth(omega_l*t);  % interference signal (numeric)
Int_In = timeseries(l, t);      % optional: for Simulink/Scopes

% Output interference
A_n = 0.5;
omega_n = 2*pi*20;
n = A_n * square(omega_n*t);    % interference signal (numeric)
Int_Out = timeseries(n, t);     % optional: for Simulink/Scopes

%% Combined Input Signals
% (For MATLAB calculations use numeric vectors; timeseries above only for Simulink)
u_tilde = u + l;

%% Force everything into column form
[t, u, l, n, u_tilde] = deal(t(:), u(:), l(:), n(:), u_tilde(:));

%% Transfer function Plant
P = (s+10) / ((s + 1)*(s+2));

%% Output Signal
y_tilde  = lsim(P, u_tilde, t);    % system response
y        = y_tilde + n;            % add output interference

%% Plot of Input
figure(1)
subplot(311)
plot(t,u);grid on;
title('Input u(t)');
subplot(312)
plot(t,l);grid on;
title('Input Interference Signal l(t)');
subplot(313)
plot(t,u_tilde);grid on;
title('Input u-tilde(t)');
sgtitle('Combination Input Signals')

%% Plot of Input
figure(2)
subplot(311)
plot(t,y_tilde);grid on;
title('Output y-tilde(t)');
subplot(312)
plot(t,n);grid on;
title('Input Interference Signal n(t)');
subplot(313)
plot(t,y);grid on;
title('Output y-tilde(t)');
sgtitle('Combination Output Signals')


%% Plot of Signals

figure(3)
subplot(411)
plot(t,u);grid on;
title('Input u(t)');
subplot(412)
plot(t,u_tilde);grid on;
title('Input u-tilde(t)');
subplot(413)
plot(t,y_tilde);grid on;
title('Output y-tilde(t)');
subplot(414)
plot(t,y);grid on;
title('Output y(t)');
sgtitle('Input and Output Signals')

%% Transferfunction without adjustments

U = fft(u);
Y = fft(y);
UTIL = fft(u_tilde);
YTIL = fft(y_tilde);
H1 = Y ./U;

% use only positive frequencies up to Nyquist
H1 = H1(1:floor(N/2));
f = f(1:floor(N/2));

% create FRD object for bode plot
Hfrd = frd(H1, 2*pi*f);   % bode expects rad/s
figure(4)
h = bodeplot(P, Hfrd, opt);
grid on
legend('P(s) (true)','Estimated (FFT)')
op = getoptions(h);
op.Title.String = '';
setoptions(h, op);
sgtitle('Transfer Function Y(s)/U(s)');

%% Transferfunctions with Signalenergy
Syy = Y .* conj(Y);
Suu = U .* conj(U);
Syu = Y .* conj(U);

th   = 1e-15 * max(Suu);     % get rid off extrem small numbers for numerical errors
H2   = Syu ./ Suu;
H2(Suu < th) = NaN;         % in case numbers go under th
H2 = H2(1:floor(N/2));


% create FRD object for bode plot
figure(5)
h = bodeplot(P, Hfrd, opt);
grid on
legend('P(s) (true)','Estimated (FFT)')
op = getoptions(h);
op.Title.String = '';
setoptions(h, op);
sgtitle('Transfer Function S_{yu}/S_{uu} (Welch H_2)');

%% Transfer function according to Welch

Nest    = round(1 / Ts);                 % number of samples per segment
overlap = 0.5;                           % 50% overlap
step    = round(Nest * (1 - overlap));   % segment hop size in samples
starts  = 1:step:(N - Nest + 1);         % start indices of segments
K       = numel(starts);                 % number of segments

w = ones(Nest,1);   % Rectwindow
w1 = hann(Nest,'periodic');               % window of length Nest (or ones(Nest,1))

% Preallocate FFT results
U_all = zeros(Nest, K);                  % Preparation of Vector of Input
Y_all = zeros(Nest, K);                  % Preparation of Vector of Outpur

for k = 1:K
    idx   = starts(k):(starts(k)+Nest-1);
    u_seg = u(idx) .* w;                 % windowed input segment
    y_seg = y(idx) .* w;                 % windowed output segment

    u_seg_win = u(idx) .* w1;                 % windowed input segment
    y_seg_win = y(idx) .* w1;                 % windowed output segment

    % FFT with length = Nest (no nfft parameter)
    U_all(:,k) = fft(u_seg);            % Load of fft from segment
    Y_all(:,k) = fft(y_seg);            % Load of fft form segment
    U_win(:,k) = fft(u_seg_win);            % Load of fft from segment
    Y_win(:,k) = fft(y_seg_win);            % Load of fft form segment
end

% Keep only positive frequencies up to Nyquist
% For real signals the FFT spectrum is symmetric, the second half is just the mirror.
pos_idx = 1:floor(Nest/2)+1;           % indices of unique frequency bins  
U_pos   = U_all(pos_idx,:);            % Load just relevant Data from fft
Y_pos   = Y_all(pos_idx,:);            % Load just relevant Data from fft
U_pos_win   = U_win(pos_idx,:);            % Load just relevant Data from fft
Y_pos_win   = Y_win(pos_idx,:);            % Load just relevant Data from fft


% Frequency axis (0 ... Nyquist, Δf = fs/Nest)
freq   = (0:floor(Nest/2))' * (fs / Nest);   % Hz
w_rad  = 2*pi*freq;                          % rad/s for FRD/Bode

% --- Welch (coherent) estimator: mean of complex quotients Hk = Y/U ---
Hk_coh    = mean(Y_pos,2) ./ (mean(U_pos,2));          % complex FRF per segment
Hk_coh_win = mean(Y_pos_win,2) ./ (mean(U_pos_win,2));          % complex FRF per segment
Hfrd_coh  = frd(Hk_coh, w_rad);
Hfrd_coh_win  = frd(Hk_coh_win, w_rad);

% --- Welch (incoherent) estimator: H2 = Syu / Suu ---
Suu       = mean(abs(U_pos).^2, 2);
Suu_win = mean(abs(U_pos_win).^2, 2);
Syu       = mean(Y_pos .* conj(U_pos), 2);
Syu_win = mean(Y_pos_win.* conj(U_pos_win), 2);
th        = 1e-12 * max(Suu);                % guard to avoid division by ~0
H_welch   = Syu ./ max(Suu, th);
H_welch_win = Syu_win ./ max(Suu_win, th);
Hfrd_welch = frd(H_welch, w_rad);
Hfrd_welch_win = frd(H_welch_win, w_rad);

figure(6)
p6 = bodeplot(P, Hfrd_welch, Hfrd_coh, opt);   % Handle holen
grid on
legend('P(s) (true)', 'Welch H_2 (S_{yu}/S_{uu})', 'Coherent mean(Y/U)', ...
       'Location','SouthWest');

op = getoptions(p6);            
op.Title.String = '';           
setoptions(p6, op);             
sgtitle('Transfer Function after Welch Coherent'); 

figure(7)
p7 = bodeplot(P, Hfrd_welch_win, Hfrd_welch, opt);
grid on
legend('P(s) (true)', 'Welch with Window', 'Welch', ...
       'Location','SouthWest');
op = getoptions(h);
op.Title.String = '';
setoptions(p7, op);
sgtitle('Transfer Function after Welch Incoherent');
